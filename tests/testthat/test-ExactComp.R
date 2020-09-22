set.seed(1)
A=matrix(runif(25), 5, 5)
logdet=det.fractional(A,log=TRUE)
nonlogdet=det.fractional(A,log=FALSE)
pivot.inv=pivot.fractional(A, trace=TRUE, inverse=TRUE)

test_that("inv", {
  expect_equal(dim(A), dim(pivot.inv$A.inv))
})
test_that("pivotinv", {
  expect_length(pivot.inv,3)
})
test_that("nonlogdet", {
  expect_length(nonlogdet,1)
})
