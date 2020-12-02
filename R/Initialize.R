#' Fit sparse PCA on a grid of alpha
#' @param M Gaussian proxy for the original dataset. Possibly obtained using PLNmodels.
#' @param r Number of missing actors.
#' @param min.size Minimal number of neighbors for each missing actor.
#' @param alphaGrid Grid for parsity controlling parameter.
#'
#' @return \itemize{
#' \item{sPcaOpt:}{ optimal spca object.}
#' \item{alphaOpt:}{ best alpha value among the provided grid.}
#' \item{loglik:}{ vector of log likelihood obtained for each value of alpha.}
#' \item{bic:}{ vecotr of BIC values.}
#' \item{cliques:}{ optimal clique of neighbors.}}
#' @export
#' @importFrom sparsepca spca
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats cov var
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' TrueClique=data$TC
#' PLNfit=norm_PLN(data$Y)
#' findclique=FitSparsePCA(PLNfit$MO,r=1)
#' initClique=findclique$cliques
#' TrueClique
#' initClique
FitSparsePCA <- function(M, r=1,min.size=1, alphaGrid=10^(seq(-4, 0, by=.1))){
  # estimate of Sigma: empirical variance of the rotated scores  + diagonal risidual variances
  n <- nrow(M); p <- ncol(M); alphaNb <- length(alphaGrid); nbDifferentH=-1
  while(nbDifferentH!=r){
    sPCA <- list()

    for(a in 1:alphaNb){
      sPCA[[a]] <- spca(M, k=r, alpha=alphaGrid[a], verbose=FALSE)
      sPCA[[a]]$Sigma <- cov(sPCA[[a]]$scores%*%t(sPCA[[a]]$transform))

      resVar <- (n-1)*apply(M - sPCA[[a]]$scores %*% t(sPCA[[a]]$transform), 2, var)/n
      sPCA[[a]]$Sigma <- sPCA[[a]]$Sigma  + diag(resVar)
      sPCA[[a]]$df <- 1 + sum(sPCA[[a]]$loadings!=0)
      sPCA[[a]]$loglik <- sum(mvtnorm::dmvnorm((M), sigma=sPCA[[a]]$Sigma, log=TRUE))
      sPCA[[a]]$bic <- sPCA[[a]]$loglik - log(n)*sPCA[[a]]$df/2
    }
    df<-unlist(lapply(sPCA, function(sPca){sPca$df}))

    good<-do.call(rbind,lapply(sPCA, function(spca){
      vec_col<-apply(spca$loadings, 2,function(col){
        (sum(col!=0)>min.size && sum(col!=0)<p)})
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
  return(list(sPcaOpt=sPcaOpt,  alphaOpt=alphaOpt, loglik=loglik, bic=bic, cliques=cliques))
}

#' Select the k first components and their complement as initial cliques
#'
#'This function aims at efficiently exploring the space of likely cliques when only one missing actor is estimated.
#'It focuses and the first two principal component axes found, and their complements.
#' @param M Gaussian proxy for the original dataset. Possibly obtained using PLNmodels.
#' @param k Number of principal components of sparse PCA to keep.
#'
#' @return A list of 2 x k cliques.
#' @export
#'
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' PLNfit=norm_PLN(data$Y)
#' data$TC
#' complement_spca(PLNfit$MO,k=2)
complement_spca<-function(M,k){
  p=ncol(M)
  cliques_spca<-FitSparsePCA(M, r=k)$cliques
  complement=lapply(cliques_spca, function(clique){setdiff(1:p,clique)})
  four_clique=lapply(c(cliques_spca,complement), function(cl) list(cl))
  return(four_clique)
}

#'Finds initial cliques using a sparse PCA on bootstraps sub-samples
#'
#' @param M Gaussian proxy for the original dataset. Possibly obtained using PLNmodels.
#' @param B Number of bootstrap samples.
#' @param r Number of missing actors.
#' @param min.size Minimum number of neighbors of missing actors.
#' @param cores Number of cores for parallel computation (uses mclapply, not available for Windows).
#' @param unique Boolean for keeping only unique results.
#'
#' @return \itemize{
#' \item{cliqueList:}{ list of all cliques found. Each element is itself a list of size r.}
#' \item{nb_occ:}{ vector of the number of times each cliques has been found by sPCA.}}
#' @export
#' @importFrom parallel mclapply
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' PLNfit=norm_PLN(data$Y)
#' boot_FitSparsePCA(PLNfit$MO, B=100, r=1)
boot_FitSparsePCA<-function(M, B,r, min.size=1,cores=1, unique=TRUE){
  cliqueList<-mclapply(1:B, function(x){
    n=nrow(M); v=0.8; n.sample=round(0.8*n, 0)
    ech=sample(1:n,n.sample,replace = FALSE)
    M.sample=M[ech,]
    c=FitSparsePCA(M.sample,r=r,min.size=min.size)$cliques
    return(c)
  }, mc.cores=cores)

  nb_occ<-tabulate(match(cliqueList,unique(cliqueList)))
  if(unique){
    cliqueList<-unique(cliqueList)
  }
  return(list(cliqueList=cliqueList,nb_occ=nb_occ) )
}

#' Runs PLN function from the PLNmodels package and normalized the outputs
#' @param Y Count dataset (n x p).
#' @param X Matrix of covariates (n x d).
#' @param O Matrix of offsets (n x p).
#' @return \itemize{
#' \item{MO:}{ normalized means (n x p).}
#' \item{SO:}{ normalized marginal variances (n x p).}
#' \item{sigma_O:}{ correlation matrix (p x p) computed from the variance-covariance matrix estimated by PLNmodels, and corresponding to observed variables.}
#' \item{theta:}{ matrix of the covariates regression coefficients (p x (d+1) including the intercept).}}
#' @export
#' @importFrom PLNmodels PLN
#' @importFrom stats cov2cor
#' @examples  n=100 # 100 samples
#' p=10  # 10 species
#' r=1 # 1 missing actor
#' data=generate_missing_data(n=n,p=p,r=r,type="scale-free", plot=TRUE)
#' X=data.frame(X1=rnorm(n), X2=runif(n))
#' normPLNfit<-norm_PLN(data$Y,X)
#' str(normPLNfit)
norm_PLN<-function(Y, X=NULL, O=NULL){
  n=nrow(Y)
  p=ncol(Y)
  if(is.null(O)) O=matrix(1, n, p)
  if(is.null(X)){
    PLNfit<-PLN(Y~1+offset(log(O)), control=list(trace=0))
  }else{
    X=as.matrix(X)
    PLNfit<-PLN(Y~1+X+offset(log(O)), control=list(trace=0))
  }
  MO<-PLNfit$var_par$M
  SO<-PLNfit$var_par$S
  sigma_obs=cov2cor(PLNfit$model_par$Sigma)
  theta=PLNfit$model_par$Theta
  #-- normalize the PLN outputs
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  return(list(MO=MO, SO=SO, sigma_O=sigma_obs, theta=theta))
}
#' Find initial cliques using blockmodels on the initial marginalized network
#'
#' This function aims at clustering the marginal network, inferred from edges probabilities obtained with either `nestorFit()` or `EMtree()`.
#' @param Y Count dataset.
#' @param sigma_O PLNmodels output: covariance matrix estimate.
#' @param MO PLNmodels output: observed means estimate.
#' @param SO PLNmodels output: observed marginal variances estimate.
#' @param k Number of groups to find.
#' @param alpha Tempering parameter.
#' @param cores Number of cores for parallel computation (uses mclapply, not available for Windows).
#'
#' @return A list of cliques
#' @export
#' @importFrom blockmodels BM_bernoulli
#' @importFrom EMtree EMtree
#' @importFrom stats cov2cor
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' data$TC
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' #-- initialize with blockmodels
#' init_blockmodels(data$Y,sigma_O, MO, SO, k=2 )
init_blockmodels<-function(Y, sigma_O, MO, SO, k=3, alpha=0.1, cores=1){
  init=initVEM(Y = Y,cliqueList=NULL, cov2cor(sigma_O),MO,r = 0)
  #--- fit nestorFit with 0 missing actor
  resVEM0<- tryCatch(nestorFit(Y,MO,SO,initList=init, eps=1e-3, alpha=alpha, maxIter=100,verbatim = 0),
                     error=function(e){e}, finally={})
  if(length(resVEM0)>3){
    sbm.0 <- BM_bernoulli("SBM_sym",1*(resVEM0$Pg>0.5), plotting="", verbosity=0,ncores=cores)
  }else{
    p=ncol(Y)
    resEM0 = EMtree(cov2cor(sigma_O))
    sbm.0 <- BM_bernoulli("SBM_sym",1*(resEM0$edges_prob>2/p), plotting="", verbosity=0,ncores=1)
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
#' internal function for initialization with blockmodels
#' @param BMobject object from blickmodels
#' @param k number of desired groups
#' @noRd
extractParamBM <- function(BMobject,k){
  model <- BMobject$model_name
  membership_name <-  BMobject$membership_name
  res <- list()
  if (model == 'bernoulli') { res$alpha <- BMobject$model_parameters[k][[1]]$pi}
  # if (model == 'bernoulli_multiplex') { res$alpha <- BMobject$model_parameters[k][[1]]$pi}
  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    res$tau <-  BMobject$memberships[[k]]$Z
    res$Z <- apply(res$tau, 1, which.max)
    n <- nrow(BMobject$memberships[[k]]$Z)
    res$pi <-  colSums(BMobject$memberships[[k]]$Z)/n
    res$k <- length(res$pi)
  }
  ########## ordering
  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    o <- switch(model,
                #  poisson = order(res$lambda %*% matrix(res$pi,ncol = 1),decreasing = TRUE),
                bernoulli  =  order(res$alpha %*% matrix(res$pi,ncol = 1),decreasing = TRUE),
                1:res$k
    )
    res$pi <- res$pi[o]
    res$alpha <- res$alpha[o,o]
    res$tau <- res$tau[,o]
    res$Z <- apply(res$tau, 1, which.max)
  }
  return(res)
}

#' Find initial cliques using mclust on the estimated correlation matrix
#'
#' @param Sigma A covariance matrix estimate (p x p).
#' @param r Number of missing actors.
#' @param n.noise Quantity of noise for mclust, set to 3 x p by default.
#'
#' @return A list of size r, with initial cliques of neighbors for each missing actor.
#' @importFrom stats prcomp runif
#' @importFrom useful cart2pol pol2cart
#' @importFrom mclust Mclust map mclustBIC
#' @importFrom dplyr mutate select
#' @export
#'
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' EMtree::draw_network(data$G,layout="nicely",curv=0,btw_rank=1,groupes=c(rep(1,10),2),
#'  nodes_size=c(5,5),nodes_label=1:11, pal_nodes= c("#adc9e0","#e7bd42"),
#'  pal_edges = "#31374f")$G
#' PLNfit<-norm_PLN(data$Y)
#' Sigma_hat=PLNfit$sigma_O
#' #-- original true clique
#' data$TC
#' #-- clique found by mclust
#' clique_mclust=init_mclust(Sigma_hat, r=1)
#' clique_mclust

init_mclust<-function(Sigma,r, n.noise=NULL){
  p=ncol(Sigma) ; ok=FALSE
  if(is.null(n.noise)) n.noise=3*p
  # extract PCA axes
  Scomp=stats::prcomp(Sigma,scale. = TRUE)
  data=data.frame(Scomp$rotation[,1:2]%*%diag(Scomp$sdev[1:2]))
  # transform to polar coordinates in half polar circle
  datapolar=useful::cart2pol(x=data[,1],y=data[,2])[,1:2]
  datapolar_half=datapolar %>% dplyr::mutate(theta2=ifelse(.data$theta>pi,.data$theta-pi,.data$theta)) %>%
    dplyr::select(r,theta2)
  colnames(datapolar_half)[1]="radius"
  while(!ok){
    # add noise in polar coords
    radius <- sqrt(stats::runif(n.noise))
    theta2 <- stats::runif(n.noise, 0, pi)
    datapolarall=rbind(datapolar_half,cbind(radius,theta2))
    # back transform to cartesian coordinates
    newdata=useful::pol2cart(datapolarall$radius,datapolarall$theta2)[,1:2]
    noiseInit<-sample(c(T,F), size=ncol(Sigma), replace=T, prob=c(3, 1))
    # run mclust to find r groups among noise
    clust= tryCatch({
      mclust::Mclust(data=newdata, initialization = list(noise=noiseInit),  G=r, verbose = FALSE)
    }, error = function(e) {#if mclust fails, draw new noise data
      message("new noise")
      radius <- sqrt(stats::runif(n.noise))
      theta2 <- stats::runif(n.noise, 0, pi)
      datapolarall=rbind(datapolar_half,cbind(radius,theta2))
      newdata=useful::pol2cart(datapolarall$radius,datapolarall$theta2)[,1:2]
      mclust::Mclust(data=newdata, initialization = list(noise=noiseInit),  G=r, verbose = FALSE)
    }, finally = { })
    # extract probable memberships and build final list of cliques
    groups<-mclust::map(clust$z)[1:p]
    cliques<-lapply(1:r, function(c){
      indices= which(groups==c)
      if(length(indices)>2) ok<<-TRUE
      return(indices)
    })
  }

  return(cliques)
}

#'  Initialize Sigma and Omega using an initial clique of neighbors of the missing actor
#'
#'  Initialize Sigma and Omega by taking the principal component of cliques as initial value for the hidden variables.
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
#' @importFrom stats cov2cor
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' Sigma=data$Sigma
#' initClique=FitSparsePCA(data$Y,r=1)$cliques
#' initOmega(Sigma=Sigma, cliqueList=initClique)
initOmega <- function(Sigma = NULL,  cliqueList,cst=1.1) {
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
  attr(OmegaFull,"dimnames")=NULL
  return(list( Sigma0 = CorrFull, K0 = OmegaFull, cliquelist = cliqueList))
}



#' Initialize all parameters for the variational inference
#'
#' @param Y Count dataset.
#' @param cliqueList List of size r of initial neighbors for each missing actor.
#' @param sigma_O PLNmodels output: estimated observed bloc of the variance-covariance matrix.
#' @param MO PLNmodels output: estimated variational mean values of the latent parameters corresponding to observed species.
#' @param r Number of missing actors.
#'
#' @return Computes the following initial matrices:
#' \itemize{
#' \item{Wginit:}{ variational edges weights matrix.}
#' \item{Winit:}{ edges weights matrix.}
#' \item{omegainit:}{ precision matrix.}
#' \item{MHinit:}{ variational mean values for the hidden latent Gaussian parameters.}}
#' @export
#'
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' sigma_O=PLNfit$sigma_O
#' #-- find initial clique using your favorite method
#' findclique=FitSparsePCA(MO,r=1)
#' initClique=findclique$cliques
#' #-- initialize the VEM
#' initList=initVEM(Y=data$Y,cliqueList=initClique,sigma_O=sigma_O, MO=MO,r=1 )
#' str(initList)
initVEM<-function(Y,cliqueList,sigma_O,MO,r,random.init=FALSE){
  p=ncol(Y)
  n=nrow(Y)
  q=p+r
  # Tree
  # Wginit <- matrix(1, p+r, p+r); Wginit =Wginit / sum(Wginit)
  # Winit <- matrix(1, p+r, p+r); Winit =Winit / sum(Winit)
  # diag(Wginit) = 0;diag(Winit) = 0
  mean.val= (round( q*(q-1)*10^(-300/(q-1)))+1) /(q*(q-1))
  Wginit = matrix(mean.val, q,q);  diag(Wginit)=0

  if(random.init){ # try different starting points for this EM
    Wginit = matrix(runif(n=q*q, min=0.9*mean.val,max=1.1*mean.val ),q,q)
    Wginit=t(Wginit)%*%Wginit/2
    diag(Wginit)=0
   }
  Winit=Wginit
  # Gaussian layer U
  if(r!=0){
    initial.param<-initOmega(sigma_O,cliqueList = cliqueList,cst=1.05 )
    omegainit=initial.param$K0
    MHinit<-sapply(cliqueList, function(clique){
      if(length(clique)>1){
        res=matrix(rowMeans(MO[,clique]), n, 1)
      }else{   res=matrix(MO[,clique], n, 1)}
      return(res)
    })

  }else{#init with no missing actors
    omegainit=solve(sigma_O)
    MHinit=NULL
  }
  return(list(Wginit= Wginit, Winit= Winit, omegainit=omegainit,MHinit=MHinit))
}

#' Heuristic for an upper value of tempering parameter alpha
#'
#' @param q Total number of variables (observed and unobserved).
#' @param n Number of samples.
#' @param sup_val maximal value of the complete estimated correlation matrix (default 0.8).
#' @param delta Maximal value allowed for a determinant (default maximal machine precision).
#'
#' @return An upper value of tempering parameter alpha.
#' @export
#'
#' @examples q=15
#' n=50
#' alphaMax(q,n)
alphaMax<-function(q,n,sup_val=0.8,delta=.Machine$double.xmax){
  x=sup_val
  C= -0.5*log(1-x^2)+x^2/(1-x^2)
  max_value=((1/(q-1))*log(delta)-log(q-1))/(C*n)
  return(max_value)
}
