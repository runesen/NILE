
<!-- README.md is generated from README.Rmd. Please edit that file -->
NILE
====

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![Travis build status](https://travis-ci.org/runesen/NILE.svg?branch=master)](https://travis-ci.org/runesen/NILE) [![Codecov test coverage](https://codecov.io/gh/runesen/NILE/branch/master/graph/badge.svg)](https://codecov.io/gh/runesen/NILE?branch=master) [![R build status](https://github.com/runesen/NILE/workflows/R-CMD-check/badge.svg)](https://github.com/runesen/NILE/actions) <!-- badges: end -->

The goal of `NILE` is to provide an implementation of the statistical estimator proposed by Christiansen et al. (2020). NILE is an estimator applicable in a nonlinear instrumental variables (IV) setting and achieves distribution generalization with linear extensions. For further details refer to the paper Christiansen et al. (2020 <https://arxiv.org/abs/20>??.?????).

Installation
------------

<!-- You can install the released version of NILE from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("NILE") -->
<!-- ``` -->
You can install the the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("runesen/NILE")
```

Example
-------

This is a basic example which shows the idea behind NILE. Let us start by importing the `NILE` library.

``` r
library(NILE)
library(splines)
```

Suppose that the true (uknown) functional relationship between *X* and *Y* is defined by a spline that linearly extrapolates outside the training data.

``` r
# true (unknown) functional relationship X -> Y
# (linearly extrapolating beyond "extrap")
fX <- function(x, extrap, beta){
  bx <- splines::ns(x, knots = seq(from=extrap[1], to=extrap[2],
                                   length.out=(n.splines.true+1))[
                                     -c(1,n.splines.true+1)],
                    Boundary.knots = extrap)
  bx%*%beta
}
```

Let us generate the data with the model ...

``` r
# data generating model
n <- 200
n.splines.true <- 4
set.seed(2)
beta0 <- runif(n.splines.true, -1,1)
alphaA <- alphaEps <- alphaH <- 1/sqrt(3)
A <- runif(n,-1,1)
H <- runif(n,-1,1)
X <- alphaA*A + alphaH*H + alphaEps*runif(n,-1,1)
Y <- fX(x=X,extrap=c(-.7,.7), beta=beta0) + .3*H + .2*runif(n,-1,1)

x.new <- seq(-2,2,length.out=100)
f.new <- fX(x=x.new,extrap=c(-.7,.7), beta=beta0)
plot(X,Y, pch=20)
lines(x.new,f.new,col="#0072B2",lwd=3)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Given the observed dataset, let us now fit NILE.

``` r
fit <- NILE(Y, # response
            X, # predictors (so far, only 1-dim supported)
            A, # anchors (1 or 2-dim, although 2-dim is experimental so far)
            lambda.star = "test", # (0 = OLS, Inf = IV, (0,Inf) =
            # nonlinear anchor regression, "test" = NILE)
            test = "tsls.over.ols",
            intercept = TRUE,
            df = 50, # number of splines used for X -> Y
            p.min = 0.05, # level at which test for lambda is performed
            x.new = x.new, # values at which predictions are required
            plot=TRUE, # diagnostics plots
            f.true = function(x) fX(x,c(-.7,.7), beta0), # if supplied, the
            # true causal function is added to the plot
            par.x = list(lambda=NULL, # positive smoothness penalty for X -> Y,
                         # if NULL, it is chosen by CV to minimize out-of-sample
                         # AR objective
                         breaks=NULL, # if breaks are supplied, exactly these
                         # will be used for splines basis
                         num.breaks=20, # will result in num.breaks+2 splines,
                         # ignored if breaks is supplied.
                         n.order=4 # order of splines
                         ),
            par.a = list(lambda=NULL, # positive smoothness penalty for fit of
                         # residuals onto A. If NULL, we first compute the OLS
                         # fit of Y onto X,
                                      # and then choose lambdaA by CV to
                         # minimize the out-of-sample MSE for predicting
                         # the OLS residuals
                         breaks=NULL, # same as above
                         num.breaks=4, # same as above
                         n.order=4 # same as above
                         ))
#> [1] "lambda.cv.a =  1.42932324587573"
#> [1] "lambda.cv.x =  0.0015933292003293"
#> [1] "lambda.star.p.uncorr =  1.45538091659546"
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" style="display: block; margin: auto;" />

References
----------

Christiansen, R., Pfister, N., Jakobsen, M. E., Gnecco, N., & Peters, J. (2020). *The difficult task of distribution generalization in nonlinear models*.
