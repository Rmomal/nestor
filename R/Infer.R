#=====
# VE step

#' computeWg
#'
#'Central function for the update of the variational edges weights inside the VEM.
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


#' computeOmega
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
#   }else{# si on donne une matrice de précision
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
#' @param A.vec
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
#' @param A.mat
#' @noRd
F_Sym2Vec <- function(A.mat){
  return(A.mat[lower.tri(A.mat)])
}

#' Calculates the Meila matrix using exact computation
#'
#' @param W a weight matrix
#' @param r the number of missing actors
#' @return the Meila matrix
#' @noRd
exactMeila<-function (W,r){ # for edges weight beta
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
  return(Mei=Mei)
}

#' Computes edges probability from weights W (Kirshner (07) formulas)
#'
#' @param W a weight matrix
#' @param r number of missing actors
#' @param it1 check if nestor is at it first iteration
#' @param verbatim boolean controling verbosity
#'
#' @return matrix of edges probabilities
#' @noRd
Kirshner <- function(W,r, it1, verbatim=FALSE){
  p = nrow(W);   L = EMtree::Laplacian(W)[-1,-1]
  K = inverse.gmp(L)
  K =  rbind(c(0, diag(K)),
             cbind(diag(K), (diag(K)%o%rep(1, p-1) + rep(1, p-1)%o%diag(K) - 2*K)))
  K = .5*(K + t(K))
  P = W * K
  P = .5*(P + t(P))
  if(!it1){ # allow instabilities at first iteration
    if(sum(P<(-1e-10))!=0){ stop("Instabilities leading to neg. proba") }
  }
  if(verbatim){
    if(length(P[P>1])!=0) cat(paste0(" / range(P-1>0)= ",min(signif(P[P>1]-1,1))," ; ", max(signif(P[P>1]-1,1))," / "))
  }
  P[P<1e-16]=0 # numerical zero
  return(P)
}

#' exact computation of matrix tree log determinant
#'
#' @param W a weight matrix
#'
#' @return \itemize{
#' \item{det}{ exact log-determinant}
#' \item{max.prec}{ boolean tracking the reach of maximal precision during computation}}
#' @noRd
logSumTree<-function(W){#
  index=1;  max.prec=FALSE
  mat=EMtree::Laplacian(W)[-index, -index]
  output=det.fractional(mat, log=TRUE)
  if(output==log(.Machine$double.xmax )) max.prec=TRUE
  return(list(det=output,max.prec=max.prec))
}


#' LowerBound
#'
#'Computes the lower bound.
#' @param Pg Edge probability matrix (qxq)
#' @param Omega market matrix with precision values (qxq)
#' @param M Matrix of variational means (nxq)
#' @param S Matrix of variational marginal variances (nxq)
#' @param W Edges weight matrix
#' @param Wg Variational edges weight matrix
#' @param p number of observed species
#' @param logSTW log of the matrix tree run on the W matrix
#' @param logSTWg log of the matrix tree run on the Wg matrix
#'
#' @return a vector containing the different terms of the lower bound
#' @export
#' @importFrom stats cov2cor
#' @examples
LowerBound<-function(Pg ,Omega, M, S, W, Wg,p, logSTW, logSTWg){
  n=nrow(M) ; q=nrow(Omega) ; O=1:p ; r=q-p
  hidden = (q!=p)
  if(hidden) H=(p+1):q
  Rho=cov2cor((1/n)*(t(M)%*%M+diag(colSums(S))))
  phi=1-Rho^2
  diag(phi)=0
  #Egh lop (Z |T)
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
#' @param MO matrix of observed means
#' @param SO matrix of observed marginal variances
#' @param SH matrix of observed hidden variances
#' @param Omega matrix containing the precision terms of precision matrices faithful ot a tree
#' @param W edges weights matrix
#' @param Wg variational edges weights matrix
#' @param MH matrix of hidden means
#' @param Pg edges probabilities matrix
#' @param logSTW log-sum tree of the W matrix
#' @param logSTWg log-sum tree of the Wg matrix
#' @param alpha tempering parameter
#' @param it1 check if nestor is at its first iteration
#' @param verbatim boolean controling verbosity
#' @param trackJ boolean for evaluating the lower bound at each parameter update
#' @param hist boolean for printing edges weights histogram at each iteration
#' @importFrom stats cov2cor
#' @noRd
VE<-function(MO,SO,SH,Omega,W,Wg,MH,Pg,logSTW,logSTWg, alpha,it1, verbatim,trackJ=FALSE, hist=FALSE){
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
  Pg.new=Kirshner(Wg.new,r, it1=it1,verbatim=verbatim)
  sumP=signif(sum(Pg.new)-2*(q-1),3)
  if(verbatim) cat(paste0(" sumP=", sumP))
  if(trackJ) LB2=c(LowerBound(Pg = Pg.new, Omega=Omega, M=M, S=S,W=W, Wg=Wg.new,p, logSTW=logSTW,logSTWg=logSTWg.new),"Wg")

  Wg=Wg.new ;  Pg=Pg.new;  logSTWg=logSTWg.new
  #-- end
  if(trackJ){
    if(hidden){ LB= rbind(LB1, LB2) }else{ LB=LB2 }
  }else{ LB=NULL}

  res=list(Pg=Pg,Wg=Wg,M=M,S=S,LB=LB,logSTWg=logSTWg,max.prec=max.prec )
  return(res)
}

#===========
#' Title
#'
#' @param M matrix of means
#' @param S matrix of marginal variances
#' @param Pg matrix of edges probabilities
#' @param Omega matrix gathering the precision terms
#' @param W edges weights matrix
#' @param Wg variational weights matrix
#' @param p number of observed species
#' @param logSTW log-sum tree of the W matrix
#' @param logSTWg log-sum tree of the Wg matrix
#' @param trackJ boolean for evaluating the lower bound at each parameter update
#' @importFrom stats cov2cor
#' @nord
Mstep<-function(M, S, Pg, Omega,W, Wg, p,logSTW, logSTWg,  trackJ=FALSE){
  n=nrow(S)  ; O=1:p ; q=ncol(Omega) ; iterM=0 ; diff=1
  hidden=(q!=p)
  if(hidden){ H=(p+1):q
  r=q-p}
  Rho = cov2cor((t(M)%*%M+ diag(colSums(S)) )/ n)

  #--- Omega
  Omega=computeOmega(Pg, Rho,p)
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
#' Core function of nestor.
#'
#' @param Y count dataset
#' @param MO estimated means from norm_PLN
#' @param SO estimated marginal vairances from norm_PLN
#' @param initList result from initVEM
#' @param maxIter maximun number of iterations
#' @param eps convergence precision parameter
#' @param alpha tempering parameter, default to 0.1
#' @param verbatim boolean for verbosity
#' @param plot plots parameters convergence and lower bound trajectory if TRUE
#' @param print.hist prints edges weights histograms at each step if TRUE
#' @param trackJ computes the lower bound at each parameter update if TRUE. Otherwise, the lower bound is only computed at each new VE step.
#'
#' @return \itemize{
#' \item{M}{ estimated means}
#' \item{S}{ estimated marginal variances}
#' \item{Pg}{ edges probabilities}
#' \item{Wg}{edges variational weights}
#' \item{W}{ edges weights}
#' \item{Omega}{ market matrix filled with precision terms common to all spanning trees}
#' \item{lowbound}{ table containing the lowerbound trajectory}
#' \item{features}{ table containing the parametes trajectory}
#' \item{finalIter}{ number of iterations until convergence was reached}
#' \item{time}{ running time of the VEM}
#' \item{max.prec}{ boolean for reach of maximal precision reached during the VEM fit}}
#' @export
#' @importFrom magrittr %>%
#' @importFrom tibble rowid_to_column
#' @importFrom tidyr gather
#' @importFrom gridExtra grid.arrange
#' @importFrom utils tail
#' @examples data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_obs=PLNfit$sigma_obs
#' #-- initialize with true clique for example
#' initClique=data$TC
#' #-- initialize the VEM
#' initList=initVEM(data$Y,cliqueList=initClique,sigma_obs, MO,r=1 )
#' #-- run core function nestor
#' with_trackJ=nestor(data$Y, MO,SO, initList=initList, maxIter=5)
nestor<-function(Y,MO,SO,initList, maxIter=20,eps=1e-2, alpha=0.1, verbatim=TRUE,
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

  while(  (diffOmega[iter] > eps) &&   (iter < maxIter) || iter<2 ){
    iter=iter+1
    if(verbatim) cat(paste0("\n Iter n°", iter))
    #--- VE
    resVE<-VE(MO=MO,SO=SO,SH=SH,Omega=Omega,W=W,Wg=Wg,MH=MH,Pg=Pg,logSTW,logSTWg,
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

    resM<-Mstep(M=M,S=S,Pg=Pg, Omega=Omega,W=W,logSTW,logSTWg, trackJ=trackJ,Wg=Wg, p=p)
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
  resVE<-VE(MO=MO,SO=SO,SH=SH,Omega=Omega,W=W,Wg=Wg,logSTW,logSTWg,MH=MH,Pg=Pg, it1=(iter==1),
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
  resM<-Mstep(M=M,S=S,Pg=Pg, Omega=Omega,W=W,logSTW,logSTWg,trackJ=trackJ, Wg=Wg, p=p )
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
  if(verbatim) cat(paste0("\nnestor ran in ",round(time,3), attr(time, "units")," and ", iter," iterations.",
                          "\nFinal weights difference: ",round(diffW[iter],7)))

  return(list(M=M,S=S,Pg=Pg,Wg=Wg,W=W,Omega=Omega, lowbound=lowbound, features=features,
              finalIter=iter, time=time,max.prec=max.prec))
}




#' Run function nestor on a list of initial cliques using parallel computation (mclapply)
#'
#' @param cliqueList a list containing all initial cliques to be tested
#' @param Y count dataset
#' @param sigma_obs result of PLN estimation: variance covariance matrix of observed data
#' @param MO result of PLN estimation: means matrix of observed data
#' @param SO result of PLN estimation: marginal variances matrix of observed data
#' @param r number of hidden variables
#' @param alpha tempering parameter
#' @param cores number of cores
#' @param maxIter maximal number of iterations of the algorithm
#' @param eps convergence precision parameter
#' @param trackJ boolean for the estimation of the lower bound at each parameter update
#'
#' @return a list containing the fit of nestor for every cliques contained in cliqueList
#' @export
#'
#' @examples  data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=FALSE)
#' PLNfit<-norm_PLN(data$Y)
#' MO<-PLNfit$MO
#' SO<-PLNfit$SO
#' sigma_obs=PLNfit$sigma_obs
#' #-- find a list of initial cliques
#' findcliqueList=boot_FitSparsePCA(data$Y, B=10, r=1)
#' cliqueList=findcliqueList$cliqueList
#' #-- run List.nestor
#' List.nestor(cliqueList,data$Y, sigma_obs, MO,SO, r=1)
List.nestor<-function(cliqueList, Y, sigma_obs, MO,SO, r,alpha=0.1, cores=1,maxIter=20,eps=1e-3, trackJ=FALSE){
  p=ncol(Y) ; O=1:p ; n=nrow(Y)
  #--- run all initialisations with parallel computation
  list<-mclapply(cliqueList, function(c){
    #init
    init=initVEM(Y = Y, cliqueList=c, sigma_obs,MO,r = r)
   #run nestor
    VEM<- tryCatch({nestor(Y=Y,MO=MO,SO=SO,initList=init, eps=eps, alpha=alpha,maxIter=maxIter,
                            verbatim = FALSE, print.hist=FALSE, trackJ=trackJ)},
                   error=function(e){e},finally={})
    VEM$clique=c
    return(VEM)
  }, mc.cores=cores)
  return(list)
}


