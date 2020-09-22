
p=4
Rho=Pg=matrix(0.5, p,p)
simpleW = matrix(1,p,p)
resTree=matrix(2/p,p,p) ; diag(resTree)=0
resOmega=matrix(-2/3,p,p)
diag(resOmega)=3/6+1
resOmega[3,4]=resOmega[4,3]=0
test_that("Kirshner", {
  expect_equal(Kirshner(simpleW,r=0, it1=TRUE, verbatim=0), resTree)
})
test_that("Meila", {
  expect_equal( exactMeila(simpleW, r=0), resTree)
})
test_that("logSumTree", {
  expect_equal( exp(logSumTree(simpleW)$det), 16)
})
test_that("computeOmega", {
  expect_equal(computeOmega(Pg, Rho, p=2), resOmega)
})
#----
p=20
Rho=matrix(0.5,p,p)
omega=matrix(0.5,p,p)
W=exp(matrix(runif(25,1,2),p,p)); W=W%*%t(W) ;diag(W)=0
null=which(W==0)
Wg=computeWg(Rho, omega, W, r=1, n=100, alpha=0, hist=TRUE)
test_that("computeWg", {
  expect_equal(sum(abs(log(Wg[-null])-(log(W[-null]) - mean(log(W[-null]))))>1e-15), 0)
})
#----
test_that("Vec-Sym", {
  expect_equal(F_Vec2Sym(F_Sym2Vec(W)), W)
})
#----
library(EMtree)
W=simpleW
logSTW=logSumTree(W)$det
logSTWg=logSumTree(W)$det
data=data_from_scratch("erdos", p = 4, n = 100)
Y=data$data
omega=data$omega
PLNfit=PLN(Y~1)
M=PLNfit$var_par$M
S=PLNfit$var_par$S
bound=LowerBound(Pg, omega, M,S, W, Wg=W, p=4, logSTW, logSTWg)
test_that("bound", {
  expect_equal(as.numeric(bound["J"]), as.numeric(bound["T1"]+bound["T2"]+bound["T3"]))
})
#---- nestor
set.seed(1)
p=10
n=100
#- r0
r=0
data=generate_missing_data(n=n,p=p,r=r,type="scale-free", plot=FALSE)
PLNfit<-norm_PLN(data$Y)
MO<-PLNfit$MO
SO<-PLNfit$SO
sigma_O=PLNfit$sigma_O
#-- initialize
initClique=data$TC
initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=r)
noMA=nestor(data$Y,MO,SO,initList, maxIter=20,eps=1e-2, alpha=0.1,
       verbatim=1, print.hist=FALSE, trackJ=FALSE)
test_that("nestor-r0", {
  expect_equal(dim(noMA$M),dim(noMA$S))
})
#- r1
r=1
data=generate_missing_data(n=n,p=p,r=r,type="scale-free", plot=FALSE)
PLNfit<-norm_PLN(data$Y)
MO<-PLNfit$MO
SO<-PLNfit$SO
sigma_O=PLNfit$sigma_O
#-- initialize
initClique=data$TC
initList=initVEM(data$Y,cliqueList=initClique,sigma_O, MO,r=r)
oneMA=nestor(data$Y,MO,SO,initList, maxIter=20,eps=1e-2, alpha=0.1,
            verbatim=1, print.hist=FALSE, trackJ=FALSE)
test_that("nestor-r1", {
  expect_equal(dim(oneMA$M),dim(oneMA$S))
})
test_that("nestor-r1", {
  expect_equal(dim(oneMA$Pg),dim(oneMA$Omega))
})
oneMAtrack=nestor(data$Y,MO,SO,initList, maxIter=20,eps=1e-2, alpha=0.1,
             verbatim=2, print.hist=FALSE, trackJ=TRUE)
test_that("nestor-r1-track", {
  expect_equal(unique(oneMAtrack$lowbound$parameter),as.factor(c("MH" , "Wg", "Omega", "W")))
})

#---- List.nestor
cliqueList=boot_FitSparsePCA(data$Y, 10, r, min.size = 1, cores = 1, unique = TRUE)$cliqueList
listFit=List.nestor(cliqueList, data$Y, sigma_O, MO,SO, r,alpha=0.1, cores=1,
            maxIter=20,eps=1e-3, trackJ=FALSE)
test_that("listNestor", {
  expect_equal(length(listFit),length(cliqueList))
})


