set.seed(1)
p=10
n=100
#- r0
r=1
data=generate_missing_data(n=n,p=p,r=r,type="scale-free", plot=TRUE)
test_that("gener", {
  expect_equal(ncol(data$Y),ncol(data$Sigma))
})
test_that("hidden", {
  expect_length(data$H,r)
})

r=0
data=generate_missing_data(n=n,p=p,r=r,type="scale-free", plot=TRUE)
test_that("gener", {
  expect_equal(ncol(data$Y),ncol(data$Sigma))
})
test_that("hidden", {
  expect_length(data$H,r)
})
