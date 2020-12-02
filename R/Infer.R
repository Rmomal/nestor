#=====
# VE step

#' Updates the variational edges weights inside the VEM.
#' @param Rho Correlation matrix.
#' @param Omega Matrix filled with precision values common to all spanning trees.
#' @param W Edge weight matrix.
#' @param r Number of missing actors.
#' @param n Number of samples.
#' @param alpha Tempering parameter.
#' @param hist Prints the histogram of the log-values of OO and OH blocs if TRUE.
#'
#' @return The variational edge weights matrix.
#' @export
#' @importFrom graphics par
computeWg<-function(Rho,Omega,W,r,n, alpha, hist=FALSE, min.val, max.val ){
  q=ncol(Rho); p=q-r; O = 1:p
  Wg<-matrix(0,q,q)
  logWg<-matrix(0,q,q)
  if(r!=0){   H = (p+1):q   }

  ## update
  null = union(which(W==0),which((1-Rho^2)==0))
  logWg[-null]<-log(W[-null])-alpha*n*(0.5*log((1-Rho^2)[-null])+(Omega*Rho)[-null])
  diag(logWg) = 0
#browser()
  #--- centrage
  gammaO=logWg[O,O]
  if(r!=0) gammaOH=logWg[O,H]
  if(hist){
    graphics::par(mfrow=c(2,1))
    hist(gammaO, breaks=20,main="O")
    hist(gammaOH, breaks=20,main="OH")
  }

#   browser()
 # logWg[-null]=logWg[-null]-mean(logWg[-null])+(q-2)*log(q)/(q-1)
  # logWg[-null][logWg[-null]<min.val]=min.val
  # logWg[-null][logWg[-null]>max.val]=max.val
  #--- trimming
  Wg = exp(logWg)
  Wg[null]=0
  if(r!=0) Wg[H,H]=0
  diag(Wg)=0
  vec_null<-apply(Wg,2,function(x){
    vec=(x==0)
    return(sum(vec))
  })

  return(Wg)
}



#' Updates the precision terms
#' @param Pg Edge probability matrix.
#' @param Rho Correlation matrix.
#' @param p Number of observed species.
#'
#' @return The matrix filled with the precision values common to all spanning trees.
#' @export
computeOmega<-function(Pg,Rho,p){
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

# Lowerbound correction
# part_JPLN<-function(mat_var,EhZZ,n, var=TRUE){
#   if(var){
#     partJPLN=-n*0.5*(det.fractional(mat_var, log=TRUE)) - 0.5*sum(EhZZ*solve(mat_var))
#   }else{# si on donne une matrice de precision
#     partJPLN=n*0.5*(det.fractional(mat_var, log=TRUE)) - 0.5*sum(EhZZ*(mat_var))
#   }
#   return(partJPLN)
# }
# getJcor<-function(vem,p,eig.tol=1e-14,eps=1e-14){
#   Omega=vem$Omega
#   q=ncol(Omega)
#   O=1:p ; H=(p+1):q ; r=length(H) ; n=nrow(vem$M)
#   EhZZ=t(vem$M[,O])%*%vem$M[,O] + diag(colSums(vem$S[,O]))
#   sigTilde = (1/n)*EhZZ
#   EgO=vem$Pg*vem$Omega+diag(diag(vem$Omega))
#   if(q==p){
#     r=0
#     EgOm = EgO[O,O]
#   }else{
#     if(r==1){
#       EgOm = EgO[O,O] - matrix(EgO[O,H],p,r)%*%matrix(EgO[H,O],r,p)/EgO[H,H]
#     }else{
#       EgOm = EgO[O,O] - matrix(EgO[O,H],p,r)%*%solve(EgO[H,H])%*%matrix(EgO[H,O],r,p)
#     }
#   }
#   JPLN_SigT = part_JPLN(sigTilde,EhZZ=EhZZ,n)
#   JPLN_EgOm = part_JPLN(EgOm,EhZZ=EhZZ,n, var=FALSE)
#   diffJPLN = JPLN_SigT-JPLN_EgOm
#
#   det=det.fractional(EgOm, log=TRUE)
#   Jcor=tail(vem$lowbound$J,1)+diffJPLN
#   delta=norm(sigTilde-EgOm, type="F")
#   return(diffJPLN)
# }

#' Makes a symmetric matrix from the vector made of its lower tirangular part
#'
#' @param A.vec A vector
#' @noRd
F_Vec2Sym <- function(A.vec){
  n = (1+sqrt(1+8*length(A.vec)))/2
  A.mat = matrix(0, n, n)
  A.mat[lower.tri(A.mat)] = A.vec
  A.mat = A.mat + t(A.mat)
  return(A.mat)
}


#' Makes a vector from the lower triangular par of a symmetric matrix
#'
#' @param A.mat A square matrix
#' @noRd
F_Sym2Vec <- function(A.mat){
  return(A.mat[lower.tri(A.mat)])
}

#' Calculates the Meila matrix using exact computation
#'
#' @param W A weight matrix.
#' @return The Meila matrix.
#' @export
exactMeila<-function (W){ # for edges weight beta
  p = nrow(W) ; index=1
  L = EMtree::Laplacian(W)[-index,-index]
  Mei =inverse.gmp(L)
  Mei = rbind(c(0, diag(Mei)),
              cbind(diag(Mei),
                    (diag(Mei) %o% rep(1, p - 1) + rep(1, p - 1) %o% diag(Mei) - 2 * Mei)
              )
  )
  Mei = 0.5 * (Mei + t(Mei))
  if(sum(Mei<0)!=0) stop("unstable Laplacian")
  return(Mei)
}

#' Computes edges probability from weights W (Kirshner (07) formulas)
#'
#' @param W A weight matrix.
#' @param it1 Checks if nestorFit is at it first iteration.
#' @param verbatim Displays unstable probabilities if set to 2.
#'
#' @return The matrix of edges probabilities.
#' @export
Kirshner <- function(W, it1, verbatim=0){
  p = nrow(W);   L = EMtree::Laplacian(W)[-1,-1]
  #K = inverse.gmp(L)
  K=solve(L)
  K =  rbind(c(0, diag(K)),
             cbind(diag(K), (diag(K)%o%rep(1, p-1) + rep(1, p-1)%o%diag(K) - 2*K)))
  K = .5*(K + t(K))
  P = W * K
  P = .5*(P + t(P))
  if(!it1){ # allow instabilities at first iteration
    if(sum(P<(-1e-10))!=0){ stop("Instabilities leading to neg. proba") }
  }
  if(verbatim==2){
    if(length(P[P>1])!=0) cat(paste0(" / range(P-1>0)= ",min(signif(P[P>1]-1,1))," ; ", max(signif(P[P>1]-1,1))," / "))
  }
  P[P<1e-16]=0 # numerical zero
  return(P)
}

#' Exact computation of matrix tree log determinant
#'
#' @param W A weight matrix.
#'
#' @return \itemize{
#' \item{det:}{ exact log-determinant of W Laplacian matrix.}
#' \item{max.prec:}{ boolean tracking the reach of maximal precision during computation.}}
#' @export
logSumTree<-function(W){
  index=1;  max.prec=FALSE
  mat=EMtree::Laplacian(W)[-index, -index]
# output=det.fractional(mat, log=TRUE)
  output=log(det(mat))
  if(output==log(.Machine$double.xmax )) max.prec=TRUE
  return(list(det=output,max.prec=max.prec))
}


#' Computes the lower bound.
#' @param Pg Edge probability matrix (qxq).
#' @param Omega Matrix with precision values common to all spanning trees (qxq).
#' @param M Matrix of variational means (nxq).
#' @param S Matrix of variational marginal variances (nxq).
#' @param W Edges weight matrix.
#' @param Wg Variational edges weight matrix.
#' @param p Number of observed species.
#' @param logSTW Log of the matrix tree run on the W matrix.
#' @param logSTWg Log of the matrix tree run on the Wg matrix.
#'
#' @return A vector containing the different terms of the lower bound. More precisely :
#' \itemize{
#' \item{J:}{ value of the lower bound.}
#' \item{T1:}{ expectation of log p(Z | T).}
#' \item{T2:}{ expectation of the difference in log scale of the tree density and its variational approximation.}
#' \item{T3:}{ Variational entropy of the latent gaussian parameters Z.} }
#' @export
#' @importFrom stats cov2cor
LowerBound<-function(Pg ,Omega, M, S, W, Wg,p, logSTW, logSTWg){
  n=nrow(M) ; q=nrow(Omega) ; O=1:p ; r=q-p
  hidden = (q!=p)
  if(hidden) H=(p+1):q
  Rho=cov2cor((1/n)*(t(M)%*%M+diag(colSums(S))))
  phi=1-Rho^2
  diag(phi)=0
  #Egh log p(Z |T)
  t1<-(-n*0.25)*sum(Pg *log( phi +(phi<1e-16) ))
  t2<-(-n*0.5) *sum((Pg+diag(q))*Omega*Rho)
  t3<- n*0.5* sum(log(diag(Omega)))  - q*n*0.5*log(2*pi)
  T1<-t1+t2+t3

  # Eglog(p) - Eg log(g)
  T2<-0.5*sum(Pg * (log(W+(W==0)) - log(Wg+(Wg==0)) )) - logSTW+ logSTWg

  #Eh log h(Z)
  T3<- 0.5*sum(log(S))+q*n*0.5*(1+log(2*pi))

  J=T1+T2+T3

  return(c(J=J, T1=T1, T2=T2,T3=T3))
}
#===========
#' Computes the variational expectation step of the algorithm
#'
#' @param MO Matrix of observed means.
#' @param SO Matrix of observed marginal variances.
#' @param SH Matrix of observed hidden variances.
#' @param Omega matrix containing the precision terms of precision matrices faithful ot a tree.
#' @param W Edges weights matrix.
#' @param Wg Variational edges weights matrix.
#' @param MH Matrix of hidden means.
#' @param Pg Edges probabilities matrix.
#' @param logSTW Log of the Matrix Tree quantity of the W matrix.
#' @param logSTWg Log of the Matrix Tree quantity of the Wg matrix.
#' @param alpha Tempering parameter.
#' @param it1 Checks if nestorFit is at its first iteration.
#' @param verbatim Displays verbose if set to 2.
#' @param trackJ Boolean for evaluating the lower bound at each parameter update.
#' @param hist Boolean for printing edges weights histogram at each iteration.
#' @return Quantities required by the Mstep funciton:
#' \itemize{
#' \item{Pg:}{ edges probabilities.}
#' \item{Wg:}{ edges variational weights.}
#' \item{M:}{ variational means.}
#' \item{S:}{ variational marginal variances.}
#' \item{LB:}{ lower bound values.}
#' \item{logSTWg:}{ log of the Matrix Tree quantity of the Wg matrix.}
#' \item{max.prec:}{ boolean tracking the reach of maximal precision during computation.}}
#' @importFrom stats cov2cor
#' @export
VEstep<-function(MO,SO,SH,Omega,W,Wg,MH,Pg,logSTW,logSTWg, alpha,it1, verbatim,trackJ=FALSE, hist=FALSE, min.val, max.val){
  #--Setting up
  n=nrow(MO); q=ncol(Omega) ;  p=ncol(MO);  O=1:ncol(MO); trim=FALSE ;
  hidden=(q!=p)
  if(hidden){
    H=(p+1):ncol(Omega); r=length(H)
    S<-cbind(SO,SH)
    M=cbind(MO,MH)
  }else{
    S=SO;    r=0;    M=MO
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
    MH=MH.new
    vec_null<-apply(as.matrix(MH,n,r),2,function(x){
      vec=(x==0)
      return(sum(vec))
    })
    M=cbind(MO,MH)
    SH <-matrix(1/(diag(Omega)[H]),n,r, byrow = TRUE)
    S<-cbind(SO,SH)
    if(trackJ) LB1=c(LowerBound(Pg = Pg, Omega=Omega, M=M, S=S,W=W, Wg=Wg,p, logSTW=logSTW,logSTWg=logSTWg),"MH")
  }else{if(trackJ) LB1=LB0}

  #--- Wg
  Rho=cov2cor((1/n)*(t(M)%*%M+diag(colSums(S))))
  if(max(abs(F_Sym2Vec(Rho)))>1){
    message("trim Rho")
    Rho[Rho>1]=1
    Rho[Rho<(-1)]=-1}
  Wg.new= computeWg(Rho, Omega, W, r, n, alpha,  hist=hist, min.val=min.val, max.val=max.val)
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
  Pg.new=Kirshner(Wg.new, it1=it1,verbatim=verbatim)
  sumP=signif(sum(Pg.new)-2*(q-1),3)
  if(verbatim==2) cat(paste0(" sumP=", sumP))
  if(trackJ) LB2=c(LowerBound(Pg = Pg.new, Omega=Omega, M=M, S=S,W=W, Wg=Wg.new,p, logSTW=logSTW,logSTWg=logSTWg.new),"Wg")

  Wg=Wg.new ;  Pg=Pg.new;  logSTWg=logSTWg.new
  #-- end
  if(trackJ){
    if(hidden){ LB= rbind(LB1, LB2) }else{ LB=LB2 }
  }else{ LB=NULL}

  res=list(Pg=Pg,Wg=Wg,M=M,S=S,LB=LB,logSTWg=logSTWg,max.prec=max.prec )
  return(res)
}

Meila <- function(W){
  if(!isSymmetric(W)){cat('Pb: W non symmetric!')}
  p = nrow(W) ; index=1
  L = EMtree::Laplacian(W)[-index,-index]
  Mei =solve(L)
  Mei = rbind(c(0, diag(Mei)),
              cbind(diag(Mei),
                    (diag(Mei) %o% rep(1, p - 1) + rep(1, p - 1) %o% diag(Mei) - 2 * Mei)
              )
  )
  Mei = 0.5 * (Mei + t(Mei))

  return(Mei)
}



F_NegLikelihood <- function(beta.vec, P,sum.constraint){
  M = Meila(F_Vec2Sym(beta.vec))
  lambda = SetLambda(P, M,sum.constraint)
  return(- sum(F_Sym2Vec(P)*(log(beta.vec+(beta.vec==0)))) +
           logSumTree(F_Vec2Sym(beta.vec))$det+
           lambda*(sum(beta.vec)-sum.constraint/2))
}

F_NegGradient_Trans <- function(gamma, P,sum.constraint){
  beta=exp(gamma)
  beta[gamma==0]=0
  M = Meila(F_Vec2Sym(beta))
  lambda = SetLambda(P, M,sum.constraint)
  #gradient with log transformation
  return((- F_Sym2Vec(P) + beta*(F_Sym2Vec(M) + lambda)))
}
F_NegLikelihood_Trans <- function(gamma, P,sum.constraint){
  #gamma=gamma-mean(gamma)
  M = Meila(F_Vec2Sym(exp(gamma)))

  lambda = SetLambda(P, M,sum.constraint)

  suppressWarnings(
    res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma)))) )+
      logSumTree(F_Vec2Sym(exp(gamma)))$det+
      lambda*(sum(exp(gamma))-sum.constraint/2))
  # cat("like val is... ",res," !!\n Detail: a=",
  # -sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) ,"/b=",
  # log(SumTree(F_Vec2Sym(exp(gamma)))),"/c=",
  #   lambda*(sum(exp(gamma))-sum.constraint/2),"\n"
  # )
  #cat(paste0("\nSumTree=",SumTree(F_Vec2Sym(exp(gamma)))))


  return( res)
}
sum.constraint.inf<-function(p,min.order=300){
  round( p*(p-1)*10^(-min.order/(p-1)))+1
}
#########################################################################
SetLambda <- function(P, M,sum.constraint=1, eps = 1e-6, start=1){
  # F.x has to be increasing. The target value is 0

  F.x <- function(x){
    if(x!=0){
      sum.constraint - sum(P / (x+M))
    }else{
      sum.constraint - (2*sum(P[upper.tri(P)] / M[upper.tri(M)]))
    }
  }
  i=1
  x.min = ifelse(F.x(0) >0,-sort(F_Sym2Vec(M))[1]+1e-10,max(-1e-4, -min(F_Sym2Vec(M))/2));
  while(F.x(x.min)>0 && i<length(unique(M))){
    i=i+1
    x.min=-sort(F_Sym2Vec(M))[i]+1e-10
  }
  if(F.x(x.min)>0) stop("Could not set lambda.")
  x.max = start
  while(F.x(x.max)<0){x.max = x.max * 10}
  x = (x.max+x.min)/2
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  while(abs(x.max-x.min) > eps){
    if(f > 0) {
      x.max = x
      f.max = f
    } else{
      x.min = x
      f.min = f
    }
    x = (x.max+x.min)/2;
    f = F.x(x)
  }

  return(x)
}
#===========
#' Computes the maximization step of the algorithm
#'
#' @param M Matrix of variational means.
#' @param S Matrix of variational marginal variances.
#' @param Pg Matrix of edges probabilities.
#' @param Omega Matrix gathering precision terms common to all spanning tree structures.
#' @param W Edges weights matrix.
#' @param Wg Variational edge weights matrix.
#' @param p Number of observed species.
#' @param logSTW Log of the Matrix Tree quantity of the W matrix.
#' @param logSTWg Log of the Matrix Tree quantity of the Wg matrix.
#' @param trackJ Boolean for evaluating the lower bound at each parameter update.
#' @return Quantities required by the VEstep function:
#' \itemize{
#' \item{W:}{ edge weights matrix.}
#' \item{Omega:}{ Matrix gathering precision terms common to all spanning tree structures.}
#' \item{LB:}{ lower bound values.}
#' \item{logSTW:}{ log of the Matrix Tree quantity of the W matrix.}
#' \item{max.prec:}{ boolean tracking the reach of maximal precision during computation.}}
#' @importFrom stats cov2cor
#' @export
Mstep<-function(M, S, Pg, Omega,W, Wg, p,logSTW, logSTWg, min.val, max.val, trackJ=FALSE){
  n=nrow(S)  ; O=1:p ; q=ncol(Omega) ; iterM=0 ; diff=1
  hidden=(q!=p)
  if(hidden){ H=(p+1):q
  r=q-p}
  Rho = cov2cor((t(M)%*%M+ diag(colSums(S)) )/ n)

  #--- Omega
  Omega=computeOmega(Pg, Rho,p)
  if(trackJ) LB1=c(LowerBound(Pg = Pg, Omega=Omega, M=M, S=S,W=W, Wg=Wg,p,logSTW=logSTW,logSTWg=logSTWg),"Omega")

  #--- Beta
  # Mei=exactMeila(W,r)
  # logW.new = matrix(0,q,q) ;null=which(Pg==0)
  # logW.new[-null]= log(Pg[-null]) - log(Mei[-null] )
  # logW.new[-null]=logW.new[-null]-mean(logW.new[-null]) #centrage
  # W.new = exp(logW.new)
  # W.new[null]=0
  # W.new[W.new< 1e-16] = 0 # numeric zeros
  init=F_Sym2Vec(W)
  #gradient ascent in log scale
  gamma_init=log(init+(init==0))
  gamma_init[init==0]=0
  sum.weights=sum.constraint.inf(q)

  gamma = stats::optim(gamma_init, F_NegLikelihood_Trans, gr=F_NegGradient_Trans,method='L-BFGS-B',
                        Pg,sum.weights, control=list(trace=0,maxit=500,  pgtol=1e-4, factr=1e+10),
                       lower=rep(min.val, q*(q-1)/2),upper=rep(max.val, q*(q-1)/2))$par

  # cat(round(mean(gamma),2))
  W.new=F_Vec2Sym(exp(gamma))
  if(hidden) W.new[H,H]=0
  diag(W.new)=0
  W=W.new

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


#' Core function of nestorFit
#'
#' @param Y Count dataset.
#' @param MO Estimated means from norm_PLN.
#' @param SO Estimated marginal vairances from norm_PLN.
#' @param initList Result list from initVEM.
#' @param maxIter Maximal number of iterations.
#' @param eps Convergence precision parameter.
#' @param alpha Tempering parameter, default to 0.1.
#' @param verbatim Integer controlling verbosity in three levels, starting at 0.
#' @param print.hist Prints edges weights histograms at each step if TRUE.
#' @param trackJ Computes the lower bound at each parameter update if TRUE. Otherwise, the lower bound is only computed at each new VE step.
#'
#' @return \itemize{
#' \item{M:}{ estimated means.}
#' \item{S:}{ estimated marginal variances.}
#' \item{Pg:}{ edges probabilities.}
#' \item{Wg:}{ edges variational weights.}
#' \item{W:}{ edges weights.}
#' \item{Omega:}{ matrix filled with precision terms common to all spanning trees.}
#' \item{lowbound:}{ table containing the lowerbound trajectory.}
#' \item{features:}{ table containing the parametes trajectory.}
#' \item{finalIter:}{ number of iterations until convergence was reached.}
#' \item{time:}{ running time of the VEM.}
#' \item{max.prec:}{ boolean for reach of maximal precision reached during the VEM fit.}}
#' @export
#' @importFrom magrittr %>%
#' @importFrom tibble rowid_to_column
#' @importFrom tidyr gather
#' @importFrom gridExtra grid.arrange
#' @importFrom utils tail
#' @examples data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' #-- initialize with true clique for example
#' initClique=data$TC
#' #-- initialize the VEM
#' initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=1 )
#' #-- run core function nestorFit
#' fit=nestorFit(data$Y, MO,SO, initList=initList, maxIter=5,verbatim=1)
#' str(fit)

nestorFit<-function(Y,MO,SO,initList, maxIter=20,eps=1e-2, alpha=0.1, verbatim=1,
                   print.hist=FALSE, trackJ=FALSE){
  n=nrow(MO);  p=ncol(MO);  O=1:ncol(MO)
  MH=initList$MHinit;omegainit=initList$omegainit
  Winit=initList$Winit;Wginit=initList$Wginit
  hidden=!is.null(MH)
  if(hidden){
    H=(p+1):ncol(omegainit); r=length(H)
    SH <-matrix(1/(diag(omegainit)[H]),n,r, byrow = TRUE)
    M=cbind(MO,MH) ; S=cbind(SO, SH)
  }else{
    r=0;    SH= NULL;    M=MO ; S=SO
  }
  Omega=omegainit;  W=Winit;  Wg=Wginit; Pg=matrix(0.5, ncol(W),ncol(W))
  iter=0 ; lowbound=list()
  diffW=c(); diffOmega=c(); diffWg=c(); diffPg=c(); diffJ=c(); diffMH=c()
  diffWiter=1 ; diffJ=1 ; J=c()
  max.prec=FALSE; projL=FALSE
  t1=Sys.time()
  logSTW=logSumTree(W)$det
  logSTWg=logSumTree(Wg)$det
  q=p+r
  min.val= (-(q-2)*log(p)+round(log(.Machine$double.xmin))+10)/(q-1)
  max.val= (-(q-2)*log(p)+round(log(.Machine$double.xmax))-10)/(q-1)
  while(  (diffOmega[iter] > eps) &&   (iter < maxIter) || iter<3 ){
    iter=iter+1
    if(verbatim==2) cat(paste0("\n Iter n°", iter))
    #--- VE
    resVE<-VEstep(MO=MO,SO=SO,SH=SH,Omega=Omega,W=W,Wg=Wg,MH=MH,Pg=Pg,logSTW,logSTWg,
              it1=(iter==1),verbatim=verbatim, alpha=alpha, trackJ=trackJ, hist=print.hist,
              min.val, max.val)
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

    resM<-Mstep(M=M,S=S,Pg=Pg, Omega=Omega,W=W,logSTW,logSTWg, trackJ=trackJ,Wg=Wg, p=p,min.val=min.val,
                max.val=max.val)
    W.new=resM$W
    diffW[iter]=abs(max(W.new-W))
    diffWiter=diffW[iter]
    W=W.new
    Omega.new=resM$Omega
    logSTW=resM$logSTW
    if(resM$max.prec) max.prec=TRUE
    diffOmega[iter]=abs(max((Omega.new)-(Omega)))
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
  resVE<-VEstep(MO=MO,SO=SO,SH=SH,Omega=Omega,W=W,Wg=Wg,logSTW,logSTWg,MH=MH,Pg=Pg, it1=(iter==1),
            verbatim=verbatim, alpha=alpha, trackJ=trackJ, hist=print.hist,
            min.val, max.val)
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
  resM<-Mstep(M=M,S=S,Pg=Pg, Omega=Omega,W=W,logSTW,logSTWg,trackJ=trackJ, Wg=Wg, p=p, min.val=min.val,
              max.val=max.val)
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
  features<-data.frame(diffPg=diffPg, diffW=diffW, diffOmega=diffOmega, diffWg=diffWg)

  t2=Sys.time()
  time=t2-t1
  if(verbatim!=0) cat(paste0("\nnestor ran in ",round(time,3), attr(time, "units")," and ", iter," iterations."))

  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,Omega=Omega, lowbound=lowbound, features=features,
              finalIter=iter, time=time,max.prec=max.prec))
}




#' Run function nestorFit on a list of initial cliques using parallel computation (mclapply)
#'
#' @param cliqueList List containing all initial cliques to be tested.
#' @param Y Count dataset.
#' @param sigma_O Result of PLN estimation: variance covariance matrix of observed data.
#' @param MO Result of PLN estimation: means matrix of observed data.
#' @param SO Result of PLN estimation: marginal variances matrix of observed data.
#' @param r Number of hidden variables.
#' @param alpha Tempering parameter.
#' @param cores Number of cores for parallel computation (uses mclapply, not available for Windows).
#' @param maxIter Maximal number of iterations of the algorithm.
#' @param eps Convergence precision parameter.
#' @param trackJ Boolean for the lower bound estimation at each parameter update instead of each step.
#'
#' @return A list containing the fit of nestorFit for every clique contained in cliqueList.
#' @export
#'
#' @examples  data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_O=PLNfit$sigma_O
#' #-- find a list of initial cliques
#' findcliqueList=boot_FitSparsePCA(MO, B=5, r=1)
#' cliqueList=findcliqueList$cliqueList
#' length(cliqueList)
#' #-- run List.nestorFit
#' fitList=List.nestorFit(cliqueList,data$Y, sigma_O, MO,SO, r=1)
#' length(fitList)
#' str(fitList[[1]])
List.nestorFit<-function(cliqueList, Y, sigma_O, MO,SO, r,alpha=0.1, cores=1,maxIter=20,eps=1e-3, trackJ=FALSE){
  p=ncol(Y) ; O=1:p ; n=nrow(Y)
  #--- run all initialisations with parallel computation
  list<-mclapply(cliqueList, function(c){
    #init
    init=initVEM(Y = Y, cliqueList=c, sigma_O,MO,r = r)
   #run nestorFit
    VEM<- tryCatch({nestorFit(Y=Y,MO=MO,SO=SO,initList=init, eps=eps, alpha=alpha,maxIter=maxIter,
                            verbatim = 0, print.hist=FALSE, trackJ=trackJ)},
                   error=function(e){e},finally={})
    VEM$clique=c
    return(VEM)
  }, mc.cores=cores)
  return(list)
}


