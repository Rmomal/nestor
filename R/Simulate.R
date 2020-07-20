generator_PLN<-function(Sigma,covariates=NULL, n=50){
  # ajout d'une constante, par rapport à EMtree::generator_PLN
  p<-ncol(Sigma)
  if(!is.null(covariates)){
    n<-nrow(covariates)
    c<-ncol(covariates)
    string<-paste0("~", paste(colnames(covariates), collapse=" + "))
    formula<-as.formula(string)
    m<- model.matrix(formula,covariates)[,-1]
    mc<-ncol(m)
    beta<-matrix(runif(p*mc),mc,p)
    prod=m %*% beta+2
  }else{
    prod=2
  }
  D<-diag(Sigma)
  R<-cov2cor(Sigma)
  U<- rmvnorm(n, rep(0,p), R)
  matsig=(matrix(rep(sqrt(D),n),n,p, byrow = TRUE))
  Y = matrix(rpois(n*p, exp(U*matsig+prod )), n, p)

  return(list(Y=Y, U=U))
}
generator_param<-function(G,signed=FALSE,v=0){
  lambda = 1;  p=ncol(G);  sumlignes=rowSums(matrix(G,p,p));
  if(sum(sumlignes==0)!=0) sumlignes[sumlignes==0]=0.1
  D=diag(sumlignes+v)
  if(signed){
    Gsign = F_Vec2Sym(F_Sym2Vec(G * matrix(2*rbinom(p^2, 1, .3)-1, p, p)))
    omega = lambda*D - Gsign
    while(min(eigen(omega)$values) < 1e-10 & lambda<1e3){
      lambda = 1.1*lambda
      omega = lambda*D - Gsign
    }
  }else{
    omega = lambda*D + G
    while (min(eigen(omega)$values) < 1e-10){
      lambda = 1.1*lambda
      omega =lambda*D + G
    }
  }
  sigma = (solve(omega))
  sim=list(sigma=sigma,omega=omega,cste=lambda)
  return(sim)
}
generator_graph<-function(p = 20, graph = "tree", dens=0.3, r=2){# prob = 0.1,
  theta = matrix(0, p, p)
  if (graph == "cluster") {
    theta<-SimCluster(p,3,dens,r)
  }
  if (graph == "scale-free") {
    theta = huge.generator(d=p,graph="scale-free", verbose=FALSE)$theta

  }
  if(graph=="tree"){
    theta<-SpannTree(p)
  }
  if(graph=="erdos"){
    theta<- erdos(p=p,prob=dens)
  }
  return(theta = Matrix(theta, sparse = TRUE))
}

data_from_scratch<-function(type, p=20,n=50, ratio=10, covariates=NULL,
                            dens=5/p, signed=FALSE,v=0,draw=FALSE){
  graph<- generator_graph(graph=type,p=p,dens=dens,r=ratio)
  param<-generator_param(G=as.matrix(graph),signed=signed,v=v)
  data<-generator_PLN(param$sigma,covariates,n)
  Y=data$Y
  U=data$U
  if(draw){
    g=as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="nicely")+
      geom_edge_link()+
      geom_node_point(size=2, color="blue")
    print(g)
  }
  return(list(Y=Y,U=U,omega= param$omega))
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
