
<!-- README.md is generated from README.Rmd. Please edit that file -->
nestor
======

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/Rmomal/nestor.svg?branch=master)](https://travis-ci.org/Rmomal/nestor) [![Codecov test coverage](https://codecov.io/gh/Rmomal/nestor/branch/master/graph/badge.svg)](https://codecov.io/gh/Rmomal/nestor?branch=master) <!-- badges: end -->

The goal of nestor is to infer network of species conditional dependencies from abundances while accounting for possible missing actors in the data.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("Rmomal/nestor")
```

Example
-------

``` r
library(nestor)
data=missing_from_scratch(n=100,p=10,r=1,type="scale-free", plot=TRUE)
PLNfit<-norm_PLN(data$Y)
MO<-PLNfit$MO
SO<-PLNfit$SO
sigma_obs=PLNfit$sigma_obs
data$TC
#-- find initial clique
findclique=init_blockmodels(data$Y,sigma_obs, MO, SO, k=3 )
initClique=findclique$cliqueList
#-- initialize the VEM
initList=initVEM(data$Y,cliquelist=initClique[[2]],sigma_obs, MO,r=1 )
#-- run one clique from cliqueList
fit=nestor(data$Y, MO,SO, initList=initList,alpha=0.2, maxIter=15, plot=TRUE, trackJ=FALSE)
fit$Pg
# run all cliques
List.nestor(initClique, data$Y,sigma_obs, MO,SO,r=1)
```
