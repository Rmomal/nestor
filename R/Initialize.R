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
#' @examples
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
