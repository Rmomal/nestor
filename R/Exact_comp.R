#' Exact computation of the determinant using the Gauss pivot method
#'
#' @param A a square positive definite matrix
#' @param trace verbosity
#' @param inverse compute inverse matrix
#'
#' @return \itemize{
#' \item{A}{ the original matrix}
#' \item{det.factors}{ factors of the determinant}
#' \item{A.inv}{ inverse of A if computed}}
#' @noRd
pivot.fractional <- function(A, trace=FALSE, inverse=TRUE) {
  A<-rcdd::d2q(as.matrix(A))
  # trace: afficher ou pas la progression
  if (trace){cat("\nPIVOT DE GAUSS POUR ECHELONNER UNE MATRICE")}
  n <- nrow(A) # nombre de ligne de A
  if (inverse) {
    A<-cbind(A,diag(rep(1,n)))
    m<-2*n    # nombre de colonnes
  } else m<-n
  det.factors<-rep("1",n+1) # La derniere case est pour le signe (suivant
  # le nombre de permutation des lignes

  for (i in 1:(n)) { # on considere le ieme pivot
    if(trace){cat("\n - ITERATION",i)}

    # VERIFICATION DE L'EXISTENCE D'UN PIVOT NON NUL DANS LA COLONE COURANTE
    if (A[i,i] == "0") { # est-ce que le pivot est nul ?
      # Le nouveau pivot candidat est le plus grand element dans le reste de la colonne
      A.coli<-rcdd::qabs(A[(i+1):n,i])
      j <- which(rcdd::qmax(A.coli)==A.coli) + i
      if (A[j,i] != "0") { # si cet element est non nul, on peut s'en servir comme pivot
        if (trace) {cat("\n\t+ Echange des lignes",i,j)}
        A[c(i,j),] <- A[c(j,i),] # echange des ligne i et j
        det.factors[n+1]<-rcdd::qxq(det.factors[n+1],"-1")
      } else { # sinon on n'a pas trouve de pivot non nul: on s'arrete la
        return(A)
      }
    }
    # ECHELONNEMENT DE LA COLONNE COURANTE
    # Normalisation de la ligne du pivot
    det.factors[i]<-A[i,i]
    A[i,]<-rcdd::qdq(A[i,],rep(A[i,i],m))   # C'est la seule division du programme...

    if (inverse) {set <- setdiff(1:n,i)  # Alors on reduit la matrice
    }  else {set <- setdiff(1:n,1:i) }
    if (trace) {cat("\n\t+ Elimination de la variable",i, "dans les lignes", set)}
    for (j in set) {
      A[j, ] <- rcdd::qmq(A[j, ] , rcdd::qxq(rep(A[j,i],m), A[i, ])) # A[i,i]=1
    }
    # AFFICHAGE DE L'ETAT COURANT DU SYSTEME
    if (trace) {
      cat("\n\t+ Etat du systeme:\n")
      print(A)
    }
  }
  if (inverse) {A.inv<-A[,(n+1):(2*n)]
  } else {A.inv=NULL}

  return(list(A=A,
              det.factors=det.factors,
              A.inv=A.inv))
}




#' Get numerator and denominators
#'
#' @param factors
#' @noRd
# logprod<-function(factors){
#   list.of.frac<-stringr::str_split(factors,"/")
#   nd<-matrix(unlist(lapply(list.of.frac,function(x){
#     if (length(x)>1) {c(x[[1]], x[[2]])} else { c(x[[1]], 1)}})),nrow=2)
#   nd<-rcdd::qabs(nd)
#   nd[1,]<-nd[1,order(stringr::str_length(nd[1,]))]
#   nd[2,]<-nd[2,order(stringr::str_length(nd[2,]))]
#   return( sum(log(abs( rcdd::q2d(rcdd::qdq(sort(nd[1,]), sort(nd[2,])) )))))
# }


# inverse.fractional<-function(A){
#   A.inv<-pivot.fractional(A,trace=FALSE,inverse=TRUE)$A.inv
#   return(rcdd::q2d(A.inv))
# }


#' Exact log-determinant computation
#'
#' @param A a square positive definite matrix
#' @param log should result be in log scale
#'
#' @return the exact determinant of A
#' @export
det.fractional<-function(A,log=TRUE){
  factors<-pivot.fractional(A,trace=FALSE,inverse=FALSE)$det.factors
  if(log){
    output<-sum(log(abs(rcdd::q2d(factors))))
    if(output >= log(.Machine$double.xmax)){
      warning("logmax machine precision reached ! \n")
      output <- log(.Machine$double.xmax )
    }

    if(output <= log(.Machine$double.xmin)){
      warning("logmin machine precision reached ! \n")
      output <- log(.Machine$double.xmin)
    }
  }else{
    output <-rcdd::q2d(rcdd::qprod(factors))
    if(output >= .Machine$double.xmax){
      warning("max machine precision reached ! \n")
      output <- .Machine$double.xmax
    }

    if(output <= .Machine$double.xmin){
      warning("min machine precision reached ! \n")
      output <- .Machine$double.xmin
    }
  }
  return(output)
}


#' Computes exact inverses using gmp package
#'
#' @param A squared positive definit matrix
#'
#' @return the exact inverse of A
#' @export
inverse.gmp<-function(A){
  p<-ncol(A)
  A.inv<-matrix(as.double(solve(gmp::as.bigq(A))),p,p)
  return(A.inv)
}
