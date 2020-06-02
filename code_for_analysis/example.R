## The NILE estimator
# true (unknown) functional relationship X -> Y
# (linearly extrapolating beyond "extrap")
n.splines.true <- 4
fX <- function(x, extrap, beta){
  bx <- splines::ns(x, knots = seq(from=extrap[1], to=extrap[2],
                                   length.out=(n.splines.true+1))[
                                     -c(1,n.splines.true+1)],
                    Boundary.knots = extrap)
  bx%*%beta
}

# data generating model
n <- 200
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
lines(x.new,f.new,col="blue",lwd=3)

## FIT!
fit <- NILE(Y, # response
            X, # predictors (so far, only 1-dim supported)
            A, # anchors (1 or 2-dim, although 2-dim is experimental so far)
            lambda.star = "test", # (0 = OLS, Inf = IV, (0,Inf) =
            # nonlinear anchor regression, "test" = NILE)
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
