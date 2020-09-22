set.seed(1)
p=10 ; n=100 ; r=1
data=generate_missing_data(n=n,p=p,r=r,type="scale-free", plot=TRUE)
k=2
compspca=complement_spca(data$Y,k=k)
test_that("complement", {
  expect_length(compspca,2*k)
})

PLNfit<-norm_PLN(data$Y)
MO<-PLNfit$MO
SO<-PLNfit$SO
sigma_O=PLNfit$sigma_O
#-- initialize with blockmodels
blockinit=init_blockmodels(data$Y,sigma_O, MO, SO, k=k )
blockinit$cliqueList
test_that("blockmodels", {
  expect_length(blockinit$cliqueList,k)
})


clique_mclust=init_mclust(PLNfit$sigma_O, r=1)
test_that("mclust", {
  expect_length(clique_mclust,1)
})

#-- alpha
q=15
n=100
alpha=alphaMax(q,n)
test_that("alpha", {
  expect_equal(alpha<0.25,TRUE)
})

cliques=list(list(blockinit$cliqueList[[1]][[1]],blockinit$cliqueList[[2]][[1]]))
initTest=initVEM(data$Y,cliqueList = list(1,2),MO = MO,sigma_O = sigma_O,r=2)
test_that("initVEMr2", {
  expect_equal(ncol(initTest$omegainit),p+2)
})
