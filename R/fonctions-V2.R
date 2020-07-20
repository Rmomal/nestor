#############

#' FitSparsePCA
#'
#' Fit sparse PCA on a grid of alpha
#' @param Y a data frame
#' @param r the number of missing actors
#' @param minV the minimum number of neighbors for each missing actor
#' @param alphaGrid the grid for alpha
#'
#' @return a list
#' @export
#'
#' @examples #blabla
FitSparsePCA <- function(Y, r=1,minV=1, alphaGrid=10^(seq(-4, 0, by=.1))){
  # estimate of Sigma: empirical variance of the rotated scores  + diagonal risidual variances
  n <- nrow(Y); p <- ncol(Y); alphaNb <- length(alphaGrid); nbDifferentH=-1
  while(nbDifferentH!=r){
    sPCA <- list()

    for(a in 1:alphaNb){
      sPCA[[a]] <- spca(Y, k=r, alpha=alphaGrid[a], verbose=FALSE)
      sPCA[[a]]$Sigma <- cov(sPCA[[a]]$scores%*%t(sPCA[[a]]$transform))

      resVar <- (n-1)*apply(Y - sPCA[[a]]$scores %*% t(sPCA[[a]]$transform), 2, var)/n
      sPCA[[a]]$Sigma <- sPCA[[a]]$Sigma  + diag(resVar)
      sPCA[[a]]$df <- 1 + sum(sPCA[[a]]$loadings!=0)
      sPCA[[a]]$loglik <- sum(dmvnorm((Y), sigma=sPCA[[a]]$Sigma, log=TRUE))
      sPCA[[a]]$bic <- sPCA[[a]]$loglik - log(n)*sPCA[[a]]$df/2
    }
    df<-unlist(lapply(sPCA, function(sPca){sPca$df}))

    good<-do.call(rbind,lapply(sPCA, function(spca){
      vec_col<-apply(spca$loadings, 2,function(col){
        (sum(col!=0)>minV && sum(col!=0)<p)})
      return(sum(vec_col)==r)
    }))
    good2<- do.call(rbind,lapply(sPCA, function(spca){# missing actors should be different
      vec_col<-lapply(1:ncol(spca$loadings), function(col){
        which(spca$loadings[,col]!=0)})
      return(length(unique(vec_col))==r)}))

    # Selects alpha via pseudo-BIC
    loglik <- unlist(lapply(sPCA, function(sPca){sPca$loglik}))
    bic <- unlist(lapply(sPCA, function(sPca){sPca$bic}))
    aOpt <- which(bic==max(bic[good & good2]))
    # Find the cliques
    alphaOpt <- alphaGrid[aOpt]
    sPcaOpt <- sPCA[[aOpt]]

    cliques <- lapply(1:ncol(sPcaOpt$loadings), function(col){
      which(sPcaOpt$loadings[,col]!=0)
    })
    nbDifferentH<-length(unique(cliques))
  }
  return(list(sPcaOpt=sPcaOpt, alphaGrid=alphaGrid, alphaOpt=alphaOpt,
              loglik=loglik, bic=bic, cliques=cliques))
}

#' boot_FitSparsePCA
#'
#' Finds initial cliques using a spca on bosstraps sub-samples
#'
#' @param Y data
#' @param B number of bootstrap samples
#' @param r number of missing actors
#' @param minV minimum number of neighbors of missing actors
#' @param cores number of cores for possible parallel computation
#' @param unique should unique results only be displayed ?
#'
#' @return
#' @export
#'
#' @examples
boot_FitSparsePCA<-function(Y, B,r, minV=1,cores=1, unique=TRUE){
  cliqueList<-mclapply(1:B, function(x){
    n=nrow(Y); v=0.8; n.sample=round(0.8*n, 0)
    ech=sample(1:n,n.sample,replace = FALSE)
    Y.sample=Y[ech,]
    c=FitSparsePCA(Y.sample,r=r,minV=minV)$cliques
    return(c)
  }, mc.cores=cores)

  nb_occ<-tabulate(match(cliqueList,unique(cliqueList)))
  if(unique){
    cliqueList<-unique(cliqueList)
  }
  return(list(cliqueList=cliqueList,nb_occ=nb_occ) )
}

#' init_blockmodels
#'
#'Find initial cliques using blockmodels on the initial network (possibly inferred using EMtree or VEMtree with r=0)
#' @param k number of groups
#' @param counts data
#' @param sigma_obs original covariance matrix estimate
#' @param MO original observed means estimate
#' @param SO original observed marginal variances estimate
#' @param alpha tempering parameter
#'
#' @return a list of possible cliques
#' @export
#'
#' @examples
init_blockmodels<-function(k, counts, sigma_obs, MO, SO, alpha=0.1){
  init=initVEM(counts = counts,initviasigma=NULL, cov2cor(sigma_obs),MO,r = 0)
  Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
  #--- fit VEMtree with 0 missing actor
  resVEM0<- tryCatch(VEMtree(counts,MO,SO,MH=MHinit,omegainit,Winit,Wginit, eps=1e-3, alpha=alpha,
                    maxIter=100, plot=FALSE,print.hist=FALSE, verbatim = FALSE,trackJ=FALSE),
                    error=function(e){e}, finally={})
  if(length(resVEM0)>3){
     sbm.0 <- BM_bernoulli("SBM_sym",1*(resVEM0$Pg>0.5), plotting="", verbosity=0)
  }else{
    p=ncol(counts)
    resEM0 = EMtree(cov2cor(sigma_obs))
    sbm.0 <- BM_bernoulli("SBM_sym",1*(resEM0$edges_prob>2/p), plotting="", verbosity=0)
  }
  sbm.0$estimate()
  paramEstimSBMPoisson <- extractParamBM(sbm.0,k)
  #--- extract k groups from inferred probabilities
  clique=list()
  clique$cliqueList= lapply(1:k, function(z){
    list(which(paramEstimSBMPoisson$Z==z))
  })
  return(clique)
}
#' extractParamBM
#'
#'internal function for initialization with blockmodels
#' @param BMobject
#' @param Q
#'
#' @return
#'
#' @examples
extractParamBM <- function(BMobject,Q){
  model <- BMobject$model_name
  membership_name <-  BMobject$membership_name
  res <- list()
  if (model == 'bernoulli') { res$alpha <- BMobject$model_parameters[Q][[1]]$pi}
  if (model == 'bernoulli_multiplex') { res$alpha <- BMobject$model_parameters[Q][[1]]$pi}
  if (model == 'poisson') {
    res$alpha <- log(BMobject$model_parameters[Q][[1]]$lambda)
    res$lambda <- BMobject$model_parameters[Q][[1]]$lambda
  }
  if (model == 'poisson_covariates') {
    res$lambda <- BMobject$model_parameters[Q][[1]]$lambda
    res$alpha <- log(BMobject$model_parameters[Q][[1]]$lambda)
    res$theta <-  BMobject$model_parameters[Q][[1]]$beta
  }
  if (model == 'bernoulli_covariates') { ### a v??rifier???
    res$alpha <- BMobject$model_parameters[Q][[1]]$pi
    res$theta <-  BMobject$model_parameters[Q][[1]]$beta
  }
  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    res$tau <-  BMobject$memberships[[Q]]$Z
    res$Z <- apply(res$tau, 1, which.max)
    n <- nrow(BMobject$memberships[[Q]]$Z)
    res$pi <-  colSums(BMobject$memberships[[Q]]$Z)/n
    res$Q <- length(res$pi)
  }
  if (membership_name == 'LBM'){
    res$tauRow <-  BMobject$memberships[[Q]]$Z1
    res$tauCol <-  BMobject$memberships[[Q]]$Z2

    res$ZRow <- apply(res$tauRow, 1, which.max)
    res$ZCol <- apply(res$tauCol, 1, which.max)
    nRow <- nrow(BMobject$memberships[[Q]]$Z1)
    nCol <- nrow(BMobject$memberships[[Q]]$Z2)
    res$piRow <-  colSums(BMobject$memberships[[Q]]$Z1)/nRow
    res$piCol <-  colSums(BMobject$memberships[[Q]]$Z2)/nCol
    res$Q <- c(length(res$piRow ),length(res$piCol))
    names(res$Q) <- c('QRow','QCol')
  }

  ########## ordering
  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    o <- switch(model,
                poisson = order(res$lambda %*% matrix(res$pi,ncol = 1),decreasing = TRUE),
                bernoulli  =  order(res$alpha %*% matrix(res$pi,ncol = 1),decreasing = TRUE),
                1:res$Q
    )
    res$pi <- res$pi[o]
    res$alpha <- res$alpha[o,o]
    res$tau <- res$tau[,o]
    res$Z <- apply(res$tau, 1, which.max)
    if (model == 'poisson') {res$lambda <- res$lambda[o,o]}
  }
  if (membership_name == 'LBM'){
    oRow <- switch(model,
                   poisson = order(res$lambda %*% matrix(res$piCol,ncol = 1),decreasing = TRUE),
                   bernoulli  =  order(res$alpha %*% matrix(res$piCol,ncol = 1),decreasing = TRUE),
                   1:res$Q[1])
    oCol <- switch(model,
                   poisson = order(c(matrix(res$piRow,nrow = 1) %*% res$lambda),decreasing = TRUE),
                   bernoulli  =  order(c(matrix(res$piRow,nrow = 1) %*% res$alpha),decreasing = TRUE),
                   1:res$Q[2])
    res$piRow <- res$piRow[oRow]
    res$piCol <- res$piCol[oCol]
    res$alpha <- res$alpha[oRow,oCol]
    res$tauRow <- res$tauRow[,oRow]
    res$tauCol <- res$tauCol[,oCol]

    if(is.vector(res$tauCol)){res$tauCol = matrix(res$tauCol,ncol=1)}
    if(is.vector(res$tauRow)){res$tauRow = matrix(res$tauRow,ncol=1)}

    res$ZRow <- apply(res$tauRow, 1, which.max)
    res$ZCol <- apply(res$tauCol, 1, which.max)

    if (model == 'poisson') {res$lambda <- res$lambda[oRow,oCol]}
  }
  return(res)
}


#' initEM
#'
#'Initialize Sigma and Omega by taking the principal component of cliques as initial value for the hidden variables.
#'  The PCA is done on the correaltion matrix, the variance of hidden variables are set to the empirical variances of the pca.
#'  The corresponding precision term in omega is set to 1 for identifiability reasons.
#' @param Sigma variance-covariance matrix (p x p)
#' @param cst small constant for positive definiteness
#' @param cliqueList list of found cliques, given as vectorss of nodes indices
#'
#' @return  \itemize{
#' \item{Sigma0: }{initial value of complete covariance matrix ((p+r)x(p+r) matrix)}
#' \item{K0: }{initial value of complete precision matrix ((p+r)x(p+r) matrix)}
#' \item{clique: }{vector containing the indices of the nodes in the clique}}
#' @export
#'
#' @examples
initEM <- function(Sigma = NULL, cst=1.1, cliqueList) {
  r=length(cliqueList) ; p=ncol(Sigma);H=(p+1):(p+r)
  Corr <- cov2cor(Sigma); sigma <- sqrt(diag(Sigma))
  coef <- matrix(0, p, r)
  sapply(seq_along(cliqueList), function(c){
    coef[, c] <<- rep(0, p);
    if(length(cliqueList[[c]])>1){ # no pca if only one neighbor
      pca <-eigen(cov2cor(Corr[cliqueList[[c]], cliqueList[[c]]]))
      coef[cliqueList[[c]], c] <<- pca$vectors[, 1]
    }else{
      coef[,c]<<-Corr[cliqueList[[c]],]
    }
  })
  # Recontructing Sigma
  CorrFull <- rbind(cbind(Corr, Corr%*%coef), cbind(t(coef)%*%Corr, cst*t(coef)%*%Corr%*%coef))
  # Initialising Omega
  OmegaFull <- solve(CorrFull)
  if(r>1) OmegaFull[H,H]<-diag(diag(OmegaFull[H,H]))
  return(list( Sigma0 = CorrFull, K0 = OmegaFull, cliquelist = cliqueList))
}



#' initVEM
#'
#'wraper to initialize all parameters for the variational inference
#' @param counts
#' @param cliquelist
#' @param sigma_obs
#' @param MO
#' @param r
#'
#' @return
#' @export
#'
#' @examples
initVEM<-function(counts,cliquelist,sigma_obs,MO,r){
  p=ncol(counts)
  n=nrow(counts)
  # Tree
  Wginit <- matrix(1, p+r, p+r); Wginit =Wginit / sum(Wginit)
  Winit <- matrix(1, p+r, p+r); Winit =Winit / sum(Winit)
  diag(Wginit) = 0;diag(Winit) = 0
  # Z
  if(r!=0){
    initial.param<-initEM(sigma_obs,n=n,cliqueList = (cliquelist),cst=1.05, pca=TRUE) # quick and dirty modif for initEM to take a covariance matrix as input
    omegainit=initial.param$K0
    MHinit<-sapply(cliquelist, function(clique){
      if(length(clique)>1){
         res=matrix(rowMeans(MO[,clique]), n, 1)
      }else{   res=matrix(MO[,clique], n, 1)}
      return(res)
    })

  }else{#init with no missing actors
    omegainit=solve(sigma_obs)
    MHinit=NULL
  }
  return(list(Wginit= Wginit, Winit= Winit, omegainit=omegainit,MHinit=MHinit))
}

## simulate data and parameters
#' missing_from_scratch
#'
#'Function to simulate data and parameters
#' @param n number of samples
#' @param p number of observed species
#' @param r number of missing actors
#' @param type type of conditional dependency structure (graph), either "erdos", "cluster", or "scale-free".
#' @param plot should the simulated graph be plotted ? default to FALSE
#' @param dens density parameter for cluster and erdos graphs. For erdos graphs, this corresponds to the edges probability of connectance
#'
#' @return \itemize{
#' \item{Y:}{}
#' \item{UH:}{}
#' \item{Sigma:}{}
#' \item{Omega:}{}
#' \item{G:}{}
#' \item{TC:}{}
#' \item{H:}{}}
#' @export
#'
#' @examples
missing_from_scratch<-function(n,p,r,type,plot=FALSE, dens=2/p){
  #generate a graph and data Y and U
  data=data_from_scratch(type = type,p = p+r,n = n,signed = FALSE,dens = dens,v = 0)
  omega=data$omega
  G=1*(omega!=0)

  # remove r missing actors
  hidden=which(diag(omega)%in%sort(diag(omega), decreasing = TRUE)[1:r])[r]
  G=G[c(setdiff(1:(p+r), hidden), hidden),c(setdiff(1:(p+r), hidden), hidden)]
  diag(G)=0
  trueClique=lapply(hidden, function(h){ which(omega[h,-h]!=0)})
  group=1*(diag(omega)==diag(omega)[hidden][1])
  if(plot){
    g=draw_network(1*(omega==1),groupes=group,layout="nicely",curv=0,nb=2,
                   pal="black",nodes_label =1:(p+r))$G
    print(g)
  }
  #compute final parameters R and Omega
  if(r!=0){
    Kh  <- omega[hidden,hidden]
    Ko  <- omega[-hidden,-hidden]
    Koh <- omega[-hidden,hidden]
    Km  <- Ko - Koh %*%solve(Kh)%*% t(Koh)
    sigmaO=solve(Km)
    counts=data$Y[,-hidden]
    UH=data$U[,hidden]
  }else{
    group=NULL
    counts=data$Y
    UH=NULL
    sigmaO=solve(omega)
  }
  R=cov2cor(solve(omega))
  omega=solve(R) # inverse of a correlation matrix for the normalized formulation

  return(list(Y=counts, UH=UH, Sigma=sigmaO, Omega=omega,G=G ,TC=trueClique, H=hidden))
}



#################
# For OPTIM

#=====
# VE step

#' computeWg
#'
#'Central function for the update of the variational edges weights.
#' @param Rho Correlation matrix
#' @param Omega Market matrix filled with all the precision values common to all spanning trees
#' @param W Edge weight matrix
#' @param r number of missing actors
#' @param n number of samples
#' @param alpha tempering parameter
#' @param hist should the histogram of the log-values of OO blocs and OH blocs be printed ?
#' @param verbatim controls verbosity
#'
#' @return The variational edge weights matrix
#' @export
#'
#' @examples
computeWg<-function(Rho,Omega,W,r,n, alpha, hist=FALSE, verbatim=FALSE ){
  q=ncol(Rho); p=q-r; O = 1:p ;   binf=exp(-20) ; bsup=exp(30)
  Wg<-matrix(0,q,q)
  logWg<-matrix(0,q,q)
  if(r!=0){   H = (p+1):q   }

  ## update
  null = union(which(W==0),which((1-Rho^2)==0))
  logWg[-null]<-log(W[-null])-alpha*n*(0.5*log((1-Rho^2)[-null])+(Omega*Rho)[-null])
  diag(logWg) = 0
  #--- centrage
  gammaO=logWg[O,O]
  if(r!=0) gammaOH=logWg[O,H]
  if(hist){
    par(mfrow=c(2,1))
    hist(gammaO, breaks=20,main="O")
    hist(gammaOH, breaks=20,main="OH")
  }
  logWg[-null]=logWg[-null]-mean(logWg[-null])

  #--- trimming
  Wg = exp(logWg)
  Wg[null]=0
  if(r!=0) Wg[H,H]=0
  diag(Wg)=0
  vec_null<-apply(Wg,2,function(x){
    vec=(x==0)
    return(sum(vec))
  })

  return(list(Wg=Wg ))
}

#=====
# M step


#' comp_Omega
#'
#'Update of the precision terms
#' @param Pg Edge probability matrix
#' @param Rho Correlation matrix
#' @param p number of observed species
#'
#' @return the market matrix filled with all the precision values common to all spanning trees
#' @export
#'
#' @examples
comp_Omega<-function(Pg,Rho,p){
  q=ncol(Rho) ; hidden=(q!=p)
  # diagonal
  quantity<-Pg*Rho^2/(1-Rho^2)
  diag(quantity)=NA
  vectorSuml=colSums(quantity, na.rm=TRUE)
  OmegaDiag =vectorSuml+1
  #off-diagonal
  Omega = -Rho/(1-Rho^2)
  diag(Omega)=OmegaDiag
  if(hidden){ #model assumption
    H=(p+1):q
    if(length(H)>1) Omega[H,H]<-diag(diag(Omega[H,H]))
  }
  return(Omega)
}


#=====
#Lower bound

#' LowerBound
#'
#'Computes the lower bound.
#' @param Pg
#' @param Omega
#' @param M
#' @param S
#' @param W
#' @param Wg
#' @param p
#' @param logSTW
#' @param logSTWg
#'
#' @return a vector containing the different terms of the lower bound
#' @export
#'
#' @examples
LowerBound<-function(Pg ,Omega, M, S, W, Wg,p, logSTW, logSTWg){
  n=nrow(M) ; q=nrow(Omega) ; O=1:p ; r=q-p
  hidden = (q!=p)
  if(hidden) H=(p+1):q
  Rho=cov2cor((1/n)*(t(M)%*%M+diag(colSums(S))))
  phi=1-Rho^2
  diag(phi)=0
  #Egh lop (Z |T)
  t1<-(-n*0.25)*sum(Pg *log( phi +(phi<1e-16) )) #manip pour 0 sur la diagonale
  t2<-(-n*0.5) *sum((Pg+diag(q))*Omega*Rho)
  t3<- n*0.5* sum(log(diag(Omega)))  - q*n*0.5*log(2*pi)
  T1<-t1+t2+t3

  # Eglog(p) - Eg log(g)
  T2<-0.5*sum(Pg * (log(W+(W==0)) - log(Wg+(Wg==0)) )) - logSTW+ logSTWg

  #Eh log h(Z), reste constant car omegaH fixé à 1 pour identifiabilité
  # T3<- 0.5*sum(apply(S,2,function(x){ log(sum(x))}))+q*n*0.5*(1+log(2*pi))
  T3<- 0.5*sum(log(S))+q*n*0.5*(1+log(2*pi))

  J=T1+T2+T3
  if(is.nan(J)) browser()
  return(c(J=J, T1=T1, T2=T2,T3=T3))
}

True_lowBound<-function(Y, M,S,theta,X, W, Wg, Pg, omega){
  p=ncol(Y)
  logSTW=logSumTree(W)$det
  logSTWg=logSumTree(Wg)$det
  partJ<-LowerBound(Pg = Pg, omega=omega, M=M, S=S,W=W, Wg=Wg,p, logSTW, logSTWg)[1]
  partY<-sum(
    -exp(X%*%t(theta) + M[,1:p]+S[,1:p]/2)+ Y*(X%*%t(theta)+M[,1:p])-
      lgamma(Y+1)
  )
  TrueJ<-as.numeric(partJ)+partY
  return(TrueJ)
}

# Lowerbound correction
part_JPLN<-function(mat_var,EhZZ,n, var=TRUE){
  if(var){
    partJPLN=-n*0.5*(det.fractional(mat_var, log=TRUE)) - 0.5*sum(EhZZ*solve(mat_var))
  }else{# si on donne une matrice de précision
    partJPLN=n*0.5*(det.fractional(mat_var, log=TRUE)) - 0.5*sum(EhZZ*(mat_var))
  }
  return(partJPLN)
}
getJcor<-function(vem,p,eig.tol=1e-14,eps=1e-14){
  Omega=vem$Omega
  q=ncol(Omega)
  O=1:p ; H=(p+1):q ; r=length(H) ; n=nrow(vem$M)
  EhZZ=t(vem$M[,O])%*%vem$M[,O] + diag(colSums(vem$S[,O]))
  sigTilde = (1/n)*EhZZ
  EgO=vem$Pg*vem$Omega+diag(diag(vem$Omega))
  if(q==p){
    r=0
    EgOm = EgO[O,O]
  }else{
    if(r==1){
      EgOm = EgO[O,O] - matrix(EgO[O,H],p,r)%*%matrix(EgO[H,O],r,p)/EgO[H,H]
    }else{
      EgOm = EgO[O,O] - matrix(EgO[O,H],p,r)%*%solve(EgO[H,H])%*%matrix(EgO[H,O],r,p)
    }
  }
  JPLN_SigT = part_JPLN(sigTilde,EhZZ=EhZZ,n)
  JPLN_EgOm = part_JPLN(EgOm,EhZZ=EhZZ,n, var=FALSE)
  diffJPLN = JPLN_SigT-JPLN_EgOm

  det=det.fractional(EgOm, log=TRUE)
  Jcor=tail(vem$lowbound$J,1)+diffJPLN
  delta=norm(sigTilde-EgOm, type="F")
  return(diffJPLN)#c(Jcor=Jcor, diff=diffJPLN,detEg=det, delta=delta)
}

########################################
# edges weight and probability
exactMeila<-function (W,r){ # for edges weight beta
  p = nrow(W) ; index=1
  # cat(paste0(" condL=",signif(cond_lap(W,1),3)))
  L = Laplacian(W)[-index,-index]
  # if(cond_lap(W,1)<1e-16){browser()
  #  message("proj Mei")
  #   L=as.matrix(nearPD(L,eig.tol = 1e-14, posd.tol = 1e-14)$mat)
  # }

  Mei =inverse.gmp(L)
  Mei = rbind(c(0, diag(Mei)),
              cbind(diag(Mei),
                    (diag(Mei) %o% rep(1, p - 1) + rep(1, p - 1) %o% diag(Mei) - 2 * Mei)
              )
  )
  Mei = 0.5 * (Mei + t(Mei))
  if(sum(Mei<0)!=0) stop("unstable Laplacian") # browser()
  return(Mei=Mei)
}

Kirshner <- function(W,r, it1, verbatim=FALSE){# for edges probability from weights W (Kirshner (07) formulas)
 p = nrow(W);   L = Laplacian(W)[-1,-1]
  # if(cond_lap(W, 1)<1e-16 && min(Re(eigen(L)$values))<0){
  #   message("L no PD")
  # }
  K = inverse.gmp(L)
  K =  rbind(c(0, diag(K)),
             cbind(diag(K), (diag(K)%o%rep(1, p-1) + rep(1, p-1)%o%diag(K) - 2*K)))
  K = .5*(K + t(K))
  P = W * K
  P = .5*(P + t(P))
  if(!it1){
    if(sum(P<(-1e-16))!=0){ #browser()
      stop("Instabilities leading to neg. proba")
      #P[P<1e-16]=0
    }
  }
  if(verbatim){
    if(length(P[P>1])!=0) cat(paste0(" / range(P-1>0)= ",min(signif(P[P>1]-1,1))," ; ", max(signif(P[P>1]-1,1))," / "))
  }
 # P[P>1]=1
  P[P<1e-16]=0 # numerical zero
  return(P)
}

logSumTree<-function(W){# exact computation of matrix tree log determinant
  index=1;  max.prec=FALSE
  mat=Laplacian(W)[-index, -index]
  output=det.fractional(mat, log=TRUE)
  if(output==log(.Machine$double.xmax )){
    max.prec=TRUE
  #  message("max.prec!")
  }
  return(list(det=output,max.prec=max.prec))
}

EdgeProba <- function(W, verbatim=FALSE, p.min=1e-16){
  logWcum = logSumTree(W)$det
  if(!isSymmetric(W))  cat('Pb: W non symmetric!')
  p = nrow(W); P = matrix(0, p, p)
  #core of computation
  sapply(1:(p-1),
         function(j){
           sapply((j+1):p,
                  function(k){
                    W_jk = W; W_jk[j, k] = W_jk[k, j] = 0 #kills kj edge in W_kj
                    P[k, j] <<- 1 - exp(logSumTree(W_jk)$det - logWcum )
                    P[j, k] <<- P[k, j]
                  })
         })
  # if(sum(is.na(P))!= 0) browser()
  if(length(which(P<p.min))!=0){

    P[which(P<p.min)]= 0
  }

  diag(P)=0
  return(P)
}

#===========
VE<-function(MO,SO,SH,Omega,W,Wg,MH,Pg,logSTW,logSTWg,eps, alpha,it1, verbatim,trackJ=FALSE, hist=FALSE){
  #--Setting up
  t1=Sys.time()
  n=nrow(MO); q=ncol(Omega) ;  p=ncol(MO);  O=1:ncol(MO); trim=FALSE ;
  hidden=(q!=p)
  if(hidden){
    H=(p+1):ncol(Omega); r=length(H)
    S<-cbind(SO,SH)
    M=cbind(MO,MH)
  }else{
    S=SO
    r=0
    M=MO
  }
  if(trackJ) LB0=LowerBound(Pg = Pg, Omega=Omega, M=M, S=S,W=W, Wg=Wg,p, logSTW=logSTW,logSTWg=logSTWg)[1]

  #-- Updates
  #--- MH
  if(hidden){
    if(r>1){
      MH.new<- (-MO) %*% (Pg[O,H] * Omega[O,H])%*% diag(1/diag(Omega)[H])
    }else{
      MH.new<- (-MO) %*% (Pg[O,H] * Omega[O,H])/diag(Omega)[H]
    }
   # if(diag(t(MH.new)%*%MH.new)<n){#filtre numérique pour la mise à jour de MH
      MH=MH.new
    #}else{ message("no MH update")}
    vec_null<-apply(as.matrix(MH,n,r),2,function(x){
      vec=(x==0)
      return(sum(vec))
    })
    M=cbind(MO,MH)

    #SH <-matrix(rep(pmax((n-diag(t(MH)%*%MH))/n,0),n),n,r, byrow = TRUE)
    SH <-matrix(1/(diag(Omega)[H]),n,r, byrow = TRUE)
    S<-cbind(SO,SH)
    if(trackJ) LB1=c(LowerBound(Pg = Pg, Omega=Omega, M=M, S=S,W=W, Wg=Wg,p, logSTW=logSTW,logSTWg=logSTWg),"MH")
  }else{if(trackJ)LB1=LB0}

  #--- Wg
  Rho=cov2cor((1/n)*(t(M)%*%M+diag(colSums(S))))
  if(max(abs(F_Sym2Vec(Rho)))>1){ browser()
    message("trim Rho")
    Rho[Rho>1]=1
    Rho[Rho<(-1)]=-1}
  compWg= computeWg(Rho, Omega, W, r, n, alpha,  hist=hist, verbatim=verbatim)
  Wg.new = compWg$Wg

  logSTWg.tot=logSumTree(Wg.new)
  logSTWg.new=logSTWg.tot$det
  max.prec=logSTWg.tot$max.prec
  if(max.prec){
    Wg.new=Wg.new/10
    logSTWg.tot=logSumTree(Wg.new)
    logSTWg.new=logSTWg.tot$det
    max.prec=logSTWg.tot$max.prec
    if(max.prec) message("max.prec!")
  }
  Pg.new=Kirshner(Wg.new,r,it1, verbatim=verbatim)
  #Pg.new=EdgeProba(Wg.new,verbatim = verbatim,p.min = 0)
  sumP=signif(sum(Pg.new)-2*(q-1),3)

  if(verbatim) cat(paste0(" sumP=", sumP))
  if(trackJ) LB2=c(LowerBound(Pg = Pg.new, Omega=Omega, M=M, S=S,W=W, Wg=Wg.new,p, logSTW=logSTW,logSTWg=logSTWg.new),"Wg")
 # if(!max.prec){
    Wg=Wg.new
    Pg=Pg.new
    logSTWg=logSTWg.new
  #}else{message("max.prec: no Wg update")}



  #-- end
  if(trackJ){
    if(hidden){ LB= rbind(LB1, LB2) }else{ LB=LB2 }
  }else{ LB=NULL}

  res=list(Pg=Pg,Wg=Wg,M=M,S=S,LB=LB,logSTWg=logSTWg,max.prec=max.prec )
  return(res)
}

#===========
Mstep<-function(M, S, Pg, Omega,W, logSTW, logSTWg, eps, Wg, p, trackJ=FALSE){
  n=nrow(S)  ; O=1:p ; q=ncol(Omega) ; iterM=0 ; diff=1
  hidden=(q!=p)
  if(hidden) H=(p+1):q
  Rho = cov2cor((t(M)%*%M+ diag(colSums(S)) )/ n)

  #--- Omega

  Omega=comp_Omega(Pg, Rho,p)
  if(trackJ) LB1=c(LowerBound(Pg = Pg, Omega=Omega, M=M, S=S,W=W, Wg=Wg,p,logSTW=logSTW,logSTWg=logSTWg),"Omega")

  #--- Beta

  Mei=exactMeila(W,r)
  logW.new = matrix(0,q,q) ;null=which(Pg==0)
  logW.new[-null]= log(Pg[-null]) - log(Mei[-null] )
  logW.new[-null]=logW.new[-null]-mean(logW.new[-null]) #centrage
  W.new = exp(logW.new)
  W.new[null]=0
  W.new[W.new< 1e-16] = 0 # numeric zeros
  if(hidden) W.new[H,H]=0
  diag(W.new)=0
  W=W.new
  #if(sum(is.na(W))!=0) browser()

  logSTW.tot=logSumTree(W)
  logSTW=logSTW.tot$det
  if(trackJ) LB2=c(LowerBound(Pg = Pg, Omega=Omega, M=M, S=S,W=W, Wg=Wg,p,logSTW=logSTW,logSTWg=logSTWg),"W")
  max.prec=logSTW.tot$max.prec
  if(max.prec){
    W=W.new/10
    logSTW.tot=logSumTree(W.new)
    logSTW.new=logSTW.tot$det
    max.prec=logSTW.tot$max.prec
    if(max.prec) message("max.prec!")
  }

  if(trackJ){LB=rbind(LB1,LB2)}else{LB=LowerBound(Pg = Pg, Omega=Omega, M=M, S=S,W=W, Wg=Wg,p,logSTW=logSTW,logSTWg=logSTWg)}
  res=list(W=W, Omega=Omega, LB=LB , logSTW=logSTW,max.prec=max.prec)

  return(res)
}

#===========
VEMtree<-function(counts,MO,SO,MH,upsi_init,W_init,Wg_init, maxIter=20,eps=1e-2, alpha,
                  verbatim=TRUE, plot=FALSE, print.hist=FALSE, trackJ=FALSE){
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO)
  hidden=!is.null(MH)
  if(hidden){
    H=(p+1):ncol(upsi_init); r=length(H)
    SH <-matrix(1/(diag(upsi_init)[H]),n,r, byrow = TRUE)
    #SH <-matrix(rep(pmax((n-diag(t(MH)%*%MH))/n,0),n),n,r, byrow = TRUE)
    M=cbind(MO,MH) ; S=cbind(SO, SH)
  }else{
    r=0
    SH= NULL
    M=MO ; S=SO
  }
  Omega=upsi_init;  W=W_init;  Wg=Wg_init; Pg=matrix(0.5, ncol(W),ncol(W))
  iter=0 ; lowbound=list()
  diffW=c(); diffUpsi=c(); diffWg=c(); diffPg=c(); diffJ=c(); diffMH=c()
  diffWiter=1 ; diffJ=1 ; J=c()
  max.prec=FALSE; projL=FALSE
  t1=Sys.time()
  logSTW=logSumTree(W)$det
  logSTWg=logSumTree(Wg)$det
  #(diffW[iter] > eps) ||
  #&& (abs(diffJ)>eps)
  #
  while(  (diffUpsi[iter] > eps) &&   (iter < maxIter) || iter<2 ){
    iter=iter+1
    if(verbatim) cat(paste0("\n Iter n°", iter))
    #--- VE
    resVE<-VE(MO=MO,SO=SO,SH=SH,Omega=Omega,W=W,Wg=Wg,MH=MH,Pg=Pg,logSTW,logSTWg,eps=1e-3,
              it1=(iter==1),verbatim=verbatim, alpha=alpha, trackJ=trackJ, hist=print.hist)
    M=resVE$M
    S=resVE$S
    Pg.new=resVE$Pg
    Wg.new=resVE$Wg
    if(hidden){
      SH<-matrix(S[,H],n,r)
      MH.new<-matrix(M[,H],n,r)
      diffMH[iter]<-abs(max(MH.new-MH))
      MH=MH.new
    }
    if(resVE$max.prec) max.prec=TRUE
    diffWg[iter]<-abs(max(Wg.new-Wg))
    diffPg[iter]<-abs(max(Pg.new-Pg))
    Wg=Wg.new
    Pg=Pg.new
    logSTWg=resVE$logSTWg
    #--- M

    resM<-Mstep(M=M,S=S,Pg=Pg, Omega=Omega,W=W,logSTW,logSTWg, trackJ=trackJ,eps=1e-3 ,Wg=Wg, p=p)
    W.new=resM$W
    #diffW[iter]=abs(max(log(W.new+(W.new==0))-log(W+(W==0))))
    diffW[iter]=abs(max(W.new-W))
    diffWiter=diffW[iter]
    W=W.new
    Omega.new=resM$Omega
    logSTW=resM$logSTW
    if(resM$max.prec) max.prec=TRUE
    diffUpsi[iter]=abs(max((Omega.new)-(Omega)))
    Omega=Omega.new

    if(trackJ){
      Jiter=as.numeric(tail(resM$LB[,1],1))
    }else{ Jiter=as.numeric(resM$LB[1]) }
    if(trackJ){
      lowbound[[iter]] = rbind( resVE$LB, resM$LB)
    }else{ lowbound[[iter]] =   resM$LB }
    diffJ=Jiter - tail(J,1)
    J<-c(J, Jiter)
  }
  ########################
  resVE<-VE(MO=MO,SO=SO,SH=SH,Omega=Omega,W=W,Wg=Wg,logSTW,logSTWg,MH=MH,Pg=Pg,eps=1e-3,it1=(iter==1),
            verbatim=verbatim, alpha=alpha, trackJ=trackJ, hist=print.hist)
  M=resVE$M
  S=resVE$S
  Pg=resVE$Pg
  Wg=resVE$Wg
  logSTWg=resVE$logSTWg
  if(resVE$max.prec) max.prec=TRUE
  if(hidden){
    SH<-matrix(S[,H],n,r)
    MH<-matrix(M[,H],n,r)
  }
  #--- M
  resM<-Mstep(M=M,S=S,Pg=Pg, Omega=Omega,W=W,logSTW,logSTWg,trackJ=trackJ, eps=1e-3 ,Wg=Wg, p=p )
  if(resM$max.prec) max.prec=TRUE

  W=resM$W
  Omega=resM$Omega

  ##########################

  if(trackJ){lowbound[[iter+1]] = rbind( resVE$LB, resM$LB)
  }else{lowbound[[iter+1]] =resM$LB}
  lowbound=data.frame(do.call(rbind,lowbound))
  lowbound[,-ncol(lowbound)]<-apply(lowbound[,-ncol(lowbound)],2,function(x) as.numeric(as.character(x)))
  if(trackJ){colnames(lowbound)[ncol(lowbound)] = "parameter"}
  else{lowbound$parameter="complete"}
  features<-data.frame(diffPg=diffPg, diffW=diffW, diffUpsi=diffUpsi, diffWg=diffWg)

  t2=Sys.time()
  time=t2-t1
  if(verbatim) cat(paste0("\nVEMtree ran in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal weights difference: ",round(diffW[iter],7)))
  if(plot){
    g1<-features  %>%  rowid_to_column() %>%
      gather(key, values, -rowid) %>%
      ggplot(aes(rowid,values, color=key))+ geom_point()+geom_line() + facet_wrap(~key, scales="free")+
      labs(x="",y="", title="Parameters")+ mytheme.dark("")+guides(color=FALSE)
    if(trackJ){
      g2<- lowbound %>% rowid_to_column() %>%  gather(key,value,-rowid,-parameter) %>%
        ggplot(aes(rowid,value, group=key))+geom_line()+geom_point(aes(color=as.factor(parameter)), size=2, alpha=0.8)+
        facet_wrap(~key, scales="free")+  labs(x="sub-iteration",y="", title="Lower bound and components")+mytheme.dark("")
    }else{ g2<- lowbound %>% rowid_to_column() %>%
      ggplot(aes(rowid,J ))+geom_line()+geom_point(aes(color=as.factor(parameter)),size=2, alpha=0.8)+
      labs(x="iteration",y="", title="Lower bound")+mytheme+ scale_color_manual("",values="#2976d6")+
      guides(color=FALSE)}
    grid.arrange(g1,g2, ncol=1)
  }
  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,Omega=Omega, lowbound=lowbound, features=features,
              finalIter=iter, time=time,max.prec=max.prec))
}




List.VEM<-function(cliquesObj, counts, sigma_obs, MO,SO, r,alpha, cores,maxIter,eps, nobeta,
                   trackJ,save=FALSE){
  p=ncol(counts) ; O=1:p ; n=nrow(counts)

  #--- run all initialisations with parallel computation
  cliques=unique(cliquesObj$cliqueList)
  list<-mclapply(seq_along(cliques), function(num){
    #init
    c=cliquesObj$cliqueList[[num]]
    init=initVEM(counts = counts, initviasigma=c, sigma_obs,MO,r = r)
    Wginit= init$Wginit; Winit= init$Winit; omegainit=init$omegainit ; MHinit=init$MHinit
    #run VEMtree

    VEM<- tryCatch({VEMtree(counts=counts,MO=MO,SO=SO,MH=MHinit,upsi_init=omegainit,W_init=Winit,
                            Wg_init=Wginit, eps=eps, alpha=alpha,maxIter=maxIter,
                            verbatim = FALSE, print.hist=FALSE, trackJ=trackJ)},
                   error=function(e){e},finally={})
    VEM$clique=c
    VEM$nbocc=cliquesObj$nb_occ[num]
    VEM$nbocc_vec=cliquesObj$nb_occ
    if(save){
      saveRDS(VEM, paste0("/Users/raphaellemomal/simulations/Fatala_missing/v2.r",r,"_num",num, ".rds"))
    }
    return(VEM)
  }, mc.cores=cores)
  return(list)
}


