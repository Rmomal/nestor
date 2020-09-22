
<!-- README.md is generated from README.Rmd. Please edit that file -->
nestor: Network inference from Species counTs with missing actORs.
==================================================================

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/Rmomal/nestor.svg?branch=master)](https://travis-ci.org/Rmomal/nestor) [![Codecov test coverage](https://codecov.io/gh/Rmomal/nestor/branch/master/graph/badge.svg)](https://codecov.io/gh/Rmomal/nestor?branch=master) <!-- badges: end -->

> `nestor` is an R package for the inference of species interaction networks from their observed abundances, while accounting for possible unobserved missing actors in the data. It is an implementation of the tree-based VEM algorithm described in the preprint <http://arxiv.org/abs/2007.14299>.

Installation
------------

### EMtree dependency

`nestor` uses functions from the R package `EMtree` which development version is available from [GitHub](https://github.com/)

``` r
devtools::install_github("Rmomal/EMtree")
```

### CRAN dependencies

``` r
required_CRAN <- c("utils", "stats", "ROCR","graphics", "mvtnorm", "parallel", 
                   "gridExtra", "reshape2"  ,"ggplot2", "magrittr", "dplyr", 
                   "tidyr", "tibble", "blockmodels", "mclust", "PLNmodels")
not_installed_CRAN <- setdiff(required_CRAN, rownames(installed.packages()))
if (length(not_installed_CRAN) > 0) install.packages(not_installed_CRAN)
```

### Installation of `nestor`

You can install the development version from [GitHub](https://github.com/) with:

``` r
devtools::install_github("Rmomal/nestor")
```
