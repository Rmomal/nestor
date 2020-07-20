
#' missing_from_scratch
#'
#'Function to simulate data and parameters including missing actors. Setting r=1 is advised.
#' @param n number of samples
#' @param p number of observed species
#' @param r number of missing actors
#' @param type type of conditional dependency structure (graph), either "erdos", "cluster", or "scale-free".
#' @param plot should the simulated graph be plotted ? default to FALSE
#' @param dens density parameter for cluster and erdos graphs. For erdos graphs, this corresponds to the edges probability of connectance
#'
#' @return \itemize{
#' \item{Y:}{ Observed count data}
#' \item{UH:}{ Hidden latent Gaussian parameters}
#' \item{Sigma:}{ Observed bloc of the variance-covariance matrix}
#' \item{Omega:}{ Observed bloc of the precision matrix}
#' \item{G:}{ The generated conditional dependency structure}
#' \item{TC:}{ The true clique of neighbors of the missing actors (a list if r>1)}
#' \item{H:}{  Indexes corresponding to missing actors in the original data}}
#' @export
#'
#' @examples
missing_from_scratch<-function(n,p,r=1,type,plot=FALSE, dens=2/p){
  #generate a graph and data Y and U
  norm_data=data_from_scratch(type = type,p = p+r,n = n,signed = FALSE,dens = dens,v = 0)
  omega=norm_data$omega
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
    counts=norm_data$data$Y[,-hidden]
    UH=norm_data$data$U[,hidden]
  }else{
    group=NULL
    counts=norm_data$data$Y
    UH=NULL
    sigmaO=solve(omega)
  }
  R=cov2cor(solve(omega))
  omega=solve(R) # inverse of a correlation matrix for the normalized formulation

  return(list(Y=counts, UH=UH, Sigma=sigmaO, Omega=omega,G=G ,TC=trueClique, H=hidden))
}
