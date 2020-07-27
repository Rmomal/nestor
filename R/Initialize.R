#' FitSparsePCA
#'
#' Fit sparse PCA on a grid of alpha
#' @param Y a data frame
#' @param r the number of missing actors
#' @param minV the minimum number of neighbors for each missing actor
#' @param alphaGrid the grid for alpha
#'
#' @return \itemize{
#' \item{sPcaOpt}{ the optimal spca object }
#' \item{alphaOpt}{ the best alpha value among the provided grid}
#' \item{loglik}{ vector of log likelihood obtained for each value of alpha}
#' \item{bic}{ vecotr of BIC values}
#' \item{cliques}{ optimal clique of neighbors}}
#' @export
#' @importFrom sparsepca spca
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats cov var
#' @examples data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' TrueClique=data$TC
#' findclique=FitSparsePCA(data$Y,r=1)
#' initClique=findclique$cliques
#' TrueClique
#' initClique
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
      sPCA[[a]]$loglik <- sum(mvtnorm::dmvnorm((Y), sigma=sPCA[[a]]$Sigma, log=TRUE))
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
  return(list(sPcaOpt=sPcaOpt,  alphaOpt=alphaOpt, loglik=loglik, bic=bic, cliques=cliques))
}

#' boot_FitSparsePCA
#'
#' Finds initial cliques using a sparse PCA on bootstraps sub-samples
#'
#' @param Y data
#' @param B number of bootstrap samples
#' @param r number of missing actors
#' @param minV minimum number of neighbors of missing actors
#' @param cores number of cores for possible parallel computation
#' @param unique should unique results only be displayed ?
#'
#' @return \itemize{
#' \item{cliqueList:}{ a list of all possible initial cliques of neighbors. Each element is of size r}
#' \item{nb_occ:}{ vector of the number of times each cliques has been found by sPCA.}}
#' @export
#' @importFrom parallel mclapply
#' @examples data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' boot_FitSparsePCA(data$Y, B=20, r=1)
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

#' norm_PLN
#'
#' Runs PLN function from the PLNmodels package and normalized the outputs
#' @param Y count dataset
#'
#' @return \itemize{
#' \item{MO}{ Normalized means}
#' \item{SO}{ Normalized marginal variances}
#' \item{sigma_obs}{Vairance-covariance matrix Sigma as estimated by PLNmodels}}
#' @export
#' @importFrom PLNmodels PLN
#'
#' @examples  data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' norm_PLN(data$Y)
norm_PLN<-function(Y){
  n=nrow(Y)
  p=ncol(Y)
  PLNfit<-PLN(Y~1)
  MO<-PLNfit$var_par$M
  SO<-PLNfit$var_par$S
  sigma_obs=PLNfit$model_par$Sigma
  #-- normalize the PLN outputs
  D=diag(sigma_obs)
  matsig=(matrix(rep(1/sqrt(D),n),n,p, byrow = TRUE))
  MO=MO*matsig
  SO=SO*matsig^2
  return(list(MO=MO, SO=SO, sigma_obs=sigma_obs))
}
#' init_blockmodels
#'
#'Find initial cliques using blockmodels on the initial network (possibly inferred using EMtree or VEMtree with r=0)
#' @param Y count data
#' @param sigma_obs original covariance matrix estimate
#' @param MO original observed means estimate
#' @param SO original observed marginal variances estimate
#' @param k number of groups
#' @param poisson boolean for the choice of model of blockmodel. If FALSE, runs bernoulli.
#' @param alpha tempering parameter
#'
#' @return a list of possible cliques
#' @export
#' @importFrom blockmodels BM_bernoulli
#' @importFrom EMtree EMtree
#' @importFrom stats cov2cor
#' @examples data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_obs=PLNfit$sigma_obs
#' #-- initialize with blockmodels
#' init_blockmodels(data$Y,sigma_obs, MO, SO, k=2 )
init_blockmodels<-function(Y, sigma_obs, MO, SO, k=3,poisson=FALSE, alpha=0.1){
  init=initVEM(Y = Y,cliqueList=NULL, cov2cor(sigma_obs),MO,r = 0)
  #--- fit VEMtree with 0 missing actor
  resVEM0<- tryCatch(VEMtree(Y,MO,SO,initList=init, eps=1e-3, alpha=alpha,
                             maxIter=100, plot=FALSE,print.hist=FALSE, verbatim = FALSE,trackJ=FALSE),
                     error=function(e){e}, finally={})
  if(length(resVEM0)>3){
    sbm.0 <- BM_bernoulli("SBM_sym",1*(resVEM0$Pg>0.5), plotting="", verbosity=0)
  }else{
    p=ncol(Y)
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
  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    res$tau <-  BMobject$memberships[[Q]]$Z
    res$Z <- apply(res$tau, 1, which.max)
    n <- nrow(BMobject$memberships[[Q]]$Z)
    res$pi <-  colSums(BMobject$memberships[[Q]]$Z)/n
    res$Q <- length(res$pi)
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
  return(res)
}


#' initOmega
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
#' @importFrom stats cov2cor
#' @examples data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=TRUE)
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
  return(list( Sigma0 = CorrFull, K0 = OmegaFull, cliquelist = cliqueList))
}



#' Initialize all parameters for the variational inference
#'
#' @param Y count data matrix
#' @param cliqueList list of initial neighbors for the missing actors
#' @param sigma_obs estimated observed bloc of the variance-covariance matrix
#' @param MO estimated mean values of the latent parameters corresponding to observed species
#' @param r number of missing actors
#'
#' @return initVEM computes the following initial matrices:
#' \itemize{
#' \item{Wginit:}{ Variational edges weights matrix}
#' \item{Winit:}{ Edges weights matrix}
#' \item{omegainit:}{ Precision matrix}
#' \item{MHinit:}{ Mean values for the hidden latent Gaussian parameters}}
#' @export
#'
#' @examples data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=TRUE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' sigma_obs=PLNfit$sigma_obs
#' #-- find initial clique
#' findclique=FitSparsePCA(data$Y,r=1)
#' initClique=findclique$cliques
#' #-- initialize the VEM
#' initVEM(Y=data$Y,cliqueList=initClique,sigma_obs=sigma_obs, MO=MO,r=1 )
initVEM<-function(Y,cliqueList,sigma_obs,MO,r){
  p=ncol(Y)
  n=nrow(Y)
  # Tree
  Wginit <- matrix(1, p+r, p+r); Wginit =Wginit / sum(Wginit)
  Winit <- matrix(1, p+r, p+r); Winit =Winit / sum(Winit)
  diag(Wginit) = 0;diag(Winit) = 0
  # Gaussian layer U
  if(r!=0){
    initial.param<-initOmega(sigma_obs,cliqueList = cliqueList,cst=1.05 )
    omegainit=initial.param$K0
    MHinit<-sapply(cliqueList, function(clique){
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
