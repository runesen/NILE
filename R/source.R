#' NILE estimator
#'
#' Method for estimating a nonlinear causal relationship \eqn{X\to Y}{X -> Y}
#' that is assumed to linearly extrapolate outside of the support
#' of \eqn{X}.
#'
#'
#'
#' @author Rune Christiansen \email{krunechristiansen@@math.ku.dk}
#'
#'
#' @param Y A numeric vector with observations from the target variable
#' @param X A numeric vector with observations from the predictor variable
#' @param A A numeric vector with observations from the exogenous variable
#' @param lambda.star either a positive numeric, \code{Inf} or \code{"test"};
#' weight which determines the relative importance of the OLS loss and the TSLS
#' loss, see details.
#' @param intercept logical; indicating whether an intercept should be included into
#' the regression model
#' @param df positive integer; number of basis splines used
#' to model the nonlinear function \eqn{X\to Y}{X -> Y}
#' @param x.new numeric; x-values at which predictions are desried
#' @param test character; !!!
#' @param p.min numeric between 0 and 1; the significance level at
#' which the test which determines lambda.star should be tested.
#' @param plot logical; diagnostic plots
#' @param f.true real-valued function of one variable; If the groundtruth is known,
#' it can be suplied and will be included in diagnostic plots.
#' @param par.x a list of different parameters determining the \eqn{B}-spline regression
#' of \eqn{Y} onto \eqn{X}
#' \itemize{
#' \item \code{breaks} numeric vector; knots at which B-spline is placed
#' \item \code{num.breaks} positive interger; number of knots used (ignored if breaks is supplied)
#' \item \code{n.order} positive integer; order of B-spline basis (default is 4, which corresponds
#' to cubic splines)
#' \item \code{pen.degree} positive integer; 1 corresponds to a penalty on the first,
#' and 2 (default) corresponds to a penalty on the second order derivative.
#' }
#' @param par.a a list of different parameters determining the B-spline regression
#' of the residuals \eqn{Y - B(X)\beta} onto the residuals \eqn{A- ???} !!!
#' \itemize{
#' \item \code{breaks} same as above
#' \item \code{num.breaks} same as above
#' \item \code{n.order} same as above
#' \item \code{pen.degree} same as above
#' }
#' @param par.cv A list with parameters for the cross-validation procedure.
#' \itemize{
#' \item \code{num.folds} either "leave-one-out" or positive integer; number of CV folds (default = 10)
#' \item \code{optim} one of "optimize" (default) or "grid.search"; optimization method used for cross-validation.
#' \item \code{n.grid} positive integer; number of grid points used in grid search
#' }
#'
#' @details The NILE estimator can be used to learn a nonlinear causal influence
#' of a real-valued predictor \eqn{X} on a real-valued response variable \eqn{Y}.
#' It exploits an instrumental variable setting; it assumes that the
#' variable \eqn{A} is a valid instrument.
#' The estimator uses B-splines to estimate the nonlinear relationship. It further
#' assumes that the causal function extrapolates lienarly outside of the empirical support
#' of \eqn{X}, and can therefore be used to obtain causal predictions (that is, predicted
#' values for \eqn{Y} under confounding-removing interventions on \eqn{X}) even for values of
#' \eqn{X} which lie outside the training support.
#'
#' On a more technical side, the NILE estimator proceeds as follows. First, two B-splines
#' \eqn{B = (B_1,...,B_{df})}{B = (B_1,...,B_df)} and \eqn{C = (C_1, ..., C_{df})}{C = (C_1, ..., C_df)}
#' are constructed, which span the values of \eqn{X} and \eqn{A}, respectively.
#' These give rise to the loss functions OLS(\eqn{\theta}), which
#' corresponds to the MSE for the prediction residuals \eqn{Y - \theta^T B}, and
#' TSLS(\eqn{\theta}),
#' which are the fitted values of the spline-regression of the residuals
#' \eqn{Y - \theta^T B}
#' onto the spline basis \eqn{C}. The NILE estimator the estimates \eqn{\theta}
#' by minimizing the objective function
#' \deqn{OLS(\theta) + \texttt{lambda.star}\times TSLS(\theta) + PEN(\theta),}{
#' OLS(\theta) + lambda.star * TSLS(\theta) + PEN(\theta),}
#' where PEN(\eqn{\theta}) is a quadratic penalty term which enforces smoothness.
#'
#' The parameter \code{lambda.star} is chosen as the largest positive value for which the
#' corresponding solution yields prediction residuals which pass a test for vanishing
#' product moment with all basis functions \eqn{C_1(A), \dots C_{df}(A)}{C_1(A), ... C_df(A)}.
#'
#'
#' @return An object of class AR, containing the following elements
#' \itemize{
#' \item \code{coefficients} estimated splines coefficients for the relationship \eqn{X\to Y}{X -> Y}
#' \item \code{residuals} prediction residuals
#' \item \code{fitted.values} fitted values
#' \item \code{betaA} estimated spline coefficients for the regression of the prediction
#' residuals onto the variable \eqn{A}
#' \item \code{fitA} fitted values for the regression of the prediction
#' residuals onto the variable \eqn{A}
#' \item \code{Y} the response variable
#' \item \code{BX} basis object holding all relevant information about the B-spline
#' basis used for the regression \eqn{X\to Y}{X -> Y}
#' \item \code{BA} basis object holding all relevant information about the B-spline
#' basis used for the regression of the residuals onto A
#' \item \code{ols} OLS loss at the estimated parameters
#' \item \code{iv} TSLS loss at the estimated parameters
#' \item \code{ar} NILE loss at the estimated parameters, that is,
#' \eqn{OLS(\theta) + \texttt{lambda.star}\times TSLS(\theta)}{OLS +
#' lambda.star * TSLS}
#' \item \code{intercept} was an intercept supplied?
#' \item \code{lambdaX} penalty parameter used for the spline regression
#' \eqn{X\to Y}{X -> Y}
#' \item \code{lambdaA} penalty parameter used for the regression of residuals
#' onto \eqn{A}
#' \item \code{lambda.star} the (estimated) value for lambda.star
#' \item \code{pred} predicted values at the supplied vector x.new
#' }
#'
#' @references \insertRef{rune2020}{NILE}
#'
#'
#' @examples
#'    n.splines.true <- 4
#'    fX <- function(x, extrap, beta){
#'      bx <- splines::ns(x, knots = seq(from=extrap[1], to=extrap[2],
#'                                       length.out=(n.splines.true+1))[
#'                                         -c(1,n.splines.true+1)],
#'                        Boundary.knots = extrap)
#'      bx%*%beta
#'    }
#'
#'    # data generating model
#'    n <- 200
#'    set.seed(2)
#'    beta0 <- runif(n.splines.true, -1,1)
#'    alphaA <- alphaEps <- alphaH <- 1/sqrt(3)
#'    A <- runif(n,-1,1)
#'    H <- runif(n,-1,1)
#'    X <- alphaA * A + alphaH * H + alphaEps*runif(n,-1,1)
#'    Y <- fX(x=X,extrap=c(-.7,.7), beta=beta0) + .3 * H + .2 * runif(n,-1,1)
#'
#'    x.new <- seq(-2,2,length.out=100)
#'    f.new <- fX(x=x.new,extrap=c(-.7,.7), beta=beta0)
#'    plot(X,Y, pch=20)
#'    lines(x.new,f.new,col="#0072B2",lwd=3)
#'
#'    ## FIT!
#'    fit <- NILE(Y, # response
#'                X, # predictors (so far, only 1-dim supported)
#'                A, # anchors (1 or 2-dim, although 2-dim is experimental so far)
#'                lambda.star = "test", # (0 = OLS, Inf = IV, (0,Inf) =
#'                # nonlinear anchor regression, "test" = NILE)
#'                test = "penalized",
#'                intercept = TRUE,
#'                df = 50, # number of splines used for X -> Y
#'                p.min = 0.05, # level at which test for lambda is performed
#'                x.new = x.new, # values at which predictions are required
#'                plot=TRUE, # diagnostics plots
#'                f.true = function(x) fX(x,c(-.7,.7), beta0), # if supplied, the
#'                # true causal function is added to the plot
#'                par.x = list(lambda=NULL, # positive smoothness penalty for X -> Y,
#'                             # if NULL, it is chosen by CV to minimize out-of-sample
#'                             # AR objective
#'                             breaks=NULL, # if breaks are supplied, exactly these
#'                             # will be used for splines basis
#'                             num.breaks=20, # will result in num.breaks+2 splines,
#'                             # ignored if breaks is supplied.
#'                             n.order=4 # order of splines
#'                ),
#'                par.a = list(lambda=NULL, # positive smoothness penalty for fit of
#'                             # residuals onto A. If NULL, we first compute the OLS
#'                             # fit of Y onto X,
#'                             # and then choose lambdaA by CV to
#'                             # minimize the out-of-sample MSE for predicting
#'                             # the OLS residuals
#'                             breaks=NULL, # same as above
#'                             num.breaks=4, # same as above
#'                             n.order=4 # same as above
#'                ))
#'
#' @import fda cvTools ggplot2 MASS expm stats
#' @export
NILE <- function(Y, X, A,
                 lambda.star = "test",
                 intercept = TRUE,
                 df = 100,
                 x.new = NULL,
                 test = "penalized",
                 p.min = 0.05,
                 plot = TRUE,
                 f.true = NULL,
                 par.x = list(),
                 par.a = list(),
                 par.cv = list(num.folds = 10, optim = "optimize")
){

  n <- length(Y)
  dA <- ifelse(is.vector(A), 1, ncol(A))

  # set missing parameters to default
  if(!exists("intercept")){
    intercept <- TRUE
  }
  if(!exists("df")){
    df <- 100
  }
  if(!exists("test")){
    test <- "penalized"
  }
  if(!exists("p.min")){
    p.min <- .05
  }
  if(!exists("pen.degree",par.x)){
    par.x$pen.degree <- 2
  }
  if(!exists("n.order",par.x)){
    par.x$n.order <- 4
  }
  if(!exists("num.breaks",par.x)){
    par.x$num.breaks <- df-2
  }
  if(!exists("pen.degree",par.a)){
    par.a$pen.degree <- rep(2,dA)
  }
  if(!exists("n.order",par.a)){
    par.a$n.order <- rep(4,dA)
  }
  if(!exists("num.breaks",par.a)){
    # par.a$num.breaks <- rep(round(20/dA^2),dA)
    par.a$num.breaks <- par.x$num.breaks
  }
  if(!exists("num.folds",par.cv)){
    par.cv$num.folds <- 10
  }
  if(!exists("optim",par.cv)){
    par.cv$optim <- "optimize"
  }
  if(!exists("n.grid",par.cv)){
    par.cv$n.grid <- 20
  }

  if(test == "penalized"){
    test.statistic <- test.statistic.penalized
    quantile.function <- quantile.function.penalized
  }
  if(test == "tsls.over.ols"){
    test.statistic <- test.statistic.tsls.over.ols
    quantile.function <- quantile.function.tsls.over.ols
  }

  # in all other cases, the OLS obective will force small residuals, so intercept in regression on A not needed
  intercept.A <- lambda.star == Inf

  # create basis object for A
  BA <- create.base(x=A, num.breaks=par.a$num.breaks, breaks=par.a$breaks,
                    n.order=par.a$n.order, pen.degree=par.a$pen.degree, lambda=par.a$lambda,
                    intercept = intercept.A, type = "pspline")
  if(is.null(BA$lambda)){
    # choose lambdaA to minimize out-of-sample MSE for OLS residuals
    # R <- lm(Y~X)$residuals
    R <- Y-mean(Y)
    BA$lambda <- pi # arbitrary value
    BA <- cv.lambda(Y=R, BX=BA, BA=BA, lambda.star=0, par.cv) # minimizes out-of-sample AR objective, here with lambda.star=0, i.e., OLS error
    print(paste("lambda.cv.a = ", BA$lambda))
  }

  # create basis object for X
  BX <- create.base(x=X, num.breaks=par.x$num.breaks, breaks=par.x$breaks,
                    n.order=par.x$n.order, pen.degree=par.x$pen.degree, lambda=par.x$lambda,
                    intercept = intercept, type = "pspline")
  if(is.numeric(lambda.star)){
    if(is.null(BX$lambda)){
      # choose lambdaX to minimize out-of-sample AR-objective
      BX <- cv.lambda(Y, BX, BA, lambda.star, par.cv)
      print(paste("lambda.cv.x = ", BX$lambda))
    }
  }

  if(lambda.star == "test"){
    if(is.null(BX$lambda)){
      # first choose lambdaX to minimize MSPE
      BX <- cv.lambda(Y, BX, BA, lambda.star=0, par.cv)
      print(paste("lambda.cv.x = ", BX$lambda))
      lambda.star <- binary.search.lambda(Y, BX, BA, p.min, test.statistic, quantile.function)
      print(paste("lambda.star.p.uncorr = ", lambda.star))
    }
    BX$lambda <- (1+lambda.star)*BX$lambda
  }

  # fit AR model
  out <- ar.fit(Y, BX, BA, lambda.star, p.min)
  obj <- ar.objective(Y, BX, BA, lambda.star, out$coefficients)
  out$ols <- obj$ols
  out$iv <- obj$iv
  out$ar <- obj$ar
  out$intercept <- intercept
  out$lambdaA <- BA$lambda
  out$lambdaX <- BX$lambda
  out$lambda.star <- lambda.star
  out$A <- A
  if(!is.null(x.new)){
    out$pred <- predict(out, x.new)
  }
  if(plot) plot(out, x.new, f.true)
  out
}


create.base <- function(x, num.breaks=3, breaks=NULL, n.order=4,
                        pen.degree=2, lambda=NULL, basis=NULL,
                        intercept = TRUE, type = "pspline"){

  d <- ifelse(is.vector(x),1,ncol(x))

  if(type == "pspline"){
    if(is.null(basis)){ # if a basis is supplied, simply store the new data x
      if(d==1){
        if(is.null(breaks)){
          xmin <- min(x)
          xmax <- max(x)
          # breaks <- quantile(x, probs = seq(0,1,length.out = num.breaks))
          breaks <- seq(xmin, xmax, length.out = num.breaks)
        }
        basis <- create.bspline.basis(norder=n.order, breaks=breaks)
      }
      if(d==2){
        if(is.null(breaks)){
          breaks <- list(quantile(x[,1], probs = seq(0,1,length.out = num.breaks[1])),
                         quantile(x[,2], probs = seq(0,1,length.out = num.breaks[2])))
        }
        basis <- list(create.bspline.basis(norder=n.order[1], breaks=breaks[[1]]),
                      create.bspline.basis(norder=n.order[2], breaks=breaks[[2]]))
      }
    }
  }
  structure(list(x=x, basis=basis, n.order=n.order, pen.degree=pen.degree, lambda=lambda,
                 breaks=breaks, dim.predictor=d, intercept = intercept, type=type),
            class="base")
}


predict.base <- function(object, xnew, ...){
  a <- c(list(x = xnew), object[c("basis", "pen.degree", "lambda",
                                  "intercept", "type")])
  object <- do.call("create.base", a)
  design.mat(object)
}


# penalty matrix lambda*K
penalty.mat <- function(object){
  if(object$type == "pspline"){
    d <- object$dim.predictor
    basis <- object$basis
    lambda <- object$lambda
    pen.degree <- object$pen.degree
    if(d == 1){
      # K <- getbasispenalty(basis,Lfdobj=pen.degree)
      K <- getbasispenalty(basis)
      pm <- lambda*K
    }
    if(d == 2){ # lambda1 and lambda2 penalize the integrated squared partial
                # second derivatives wrt x1 and x2, respectively
      K1 <- getbasispenalty(basis[[1]],Lfdobj=pen.degree[1])
      K2 <- getbasispenalty(basis[[2]],Lfdobj=pen.degree[2])

      K1 <- kronecker(diag(basis[[2]]$nbasis), K1)
      K2 <- kronecker(K2, diag(basis[[1]]$nbasis))

      pm <- lambda[1]*K1 + lambda[2]*K2
    }
  }
  if(object$type=="linear"){
    Z <- design.mat(object)
    pm <- 0*diag(ncol(Z))
  }
  pm
}


# design matrix (basis functions evaluated at data)
design.mat <- function(object){
  if(object$type == "pspline"){
    d <- object$dim.predictor
    basis <- object$basis
    x <- object$x
    if(d==1){
      Z <- getbasismatrix(x, basis, nderiv=0, returnMatrix=FALSE)
    }
    if(d==2){ # Tensor product spline basis
      Z1 <- getbasismatrix(x[,1], basis[[1]], nderiv=0, returnMatrix=FALSE)
      Z2 <- getbasismatrix(x[,2], basis[[2]], nderiv=0, returnMatrix=FALSE)
      Z <- matrix(0, ncol=ncol(Z1)*ncol(Z2), nrow=nrow(Z1))
      i <- 1
      for(j in 1:ncol(Z1)) {
        for(k in 1:ncol(Z2)) {
          Z[,i] <- Z1[,j]*Z2[,k]
          i <- i+1
        }
      }
    }
  }
  if(object$type=="linear"){
    d <- object$dim.predictor
    Z <- matrix(object$x, ncol=d)
    intercept <- object$intercept
    if(intercept){
      Z <- cbind(rep(1, nrow(Z)), Z)
    }
  }
  Z
}


ar.objective <- function(Y, BX, BA, lambda.star, beta){

  ZX <- design.mat(BX)
  lKX <- penalty.mat(BX)
  ZA <- design.mat(BA)
  lKA <- penalty.mat(BA)

  # 'hat matrix' for the fit onto A's spline basis.
  WA <- ZA%*%solve(t(ZA)%*%ZA+lKA, t(ZA))

  ols <- t(Y-ZX%*%beta)%*%(Y-ZX%*%beta)
  pen <- t(beta)%*%lKX%*%beta
  iv <- t(Y-ZX%*%beta)%*%t(WA)%*%WA%*%(Y-ZX%*%beta)

  if(lambda.star < Inf){
    ar <- ols + lambda.star*iv
  } else{
    ar <- iv
  }
  list(ols=ols, ar=ar, iv=iv, pen=pen)
}


ar.fit <- function(Y, BX, BA, lambda.star, p.min){

  n <- length(Y)
  ZX <- design.mat(BX)
  lKX <- penalty.mat(BX)
  ZA <- design.mat(BA)
  lKA <- penalty.mat(BA)

  WA <- ZA%*%solve(t(ZA)%*%ZA + lKA, t(ZA))

  betaX <- solve(t(ZX)%*%ZX + lKX,t(ZX)%*%Y) # OLS
  if(lambda.star!=0){
    if(lambda.star==Inf){ # IV
      betaX <- solve(t(ZX)%*%t(WA)%*%WA%*%ZX + lKX, t(ZX)%*%t(WA)%*%WA%*%Y)
    } else{ # AR
      betaX <- solve(t(ZX)%*%(diag(n)+lambda.star*t(WA)%*%WA)%*%ZX+lKX,
                     t(ZX)%*%(diag(n)+lambda.star*t(WA)%*%WA)%*%Y)
    }
  }
  fitX <- ZX%*%betaX
  resX <- Y - fitX

  betaA <- solve(t(ZA)%*%ZA + lKA, t(ZA)%*%resX)
  fitA <- ZA%*%betaA

  structure(list(coefficients = betaX, residuals = resX, fitted.values = fitX,
                 betaA = betaA, fitA = fitA, Y=Y, BX=BX, BA=BA),
            class = "AR")
}


test.statistic.penalized <- function(Y, BX, BA, beta){

  ZX <- design.mat(BX)
  ZA <- design.mat(BA)
  pA <- ncol(ZA)
  lambdaA <- BA$lambda
  KA <- penalty.mat(BA) / lambdaA
  R <- Y - ZX %*% beta
  Rhat <- ar.fit(R, BA, BA, 0, 0)$fitted.values
  s2hat <- var(R-Rhat)
  Tn <- sum(Rhat^2) # as in Chen et al.
  e <- eigen(t(ZA)%*%ZA)
  ZAtZA.inv.sqrt <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)
  # range(solve(t(ZA)%*%ZA) - ZAtZA.inv.sqrt%*%ZAtZA.inv.sqrt)
  tmp <- ZAtZA.inv.sqrt %*% KA %*% ZAtZA.inv.sqrt
  diag(tmp) <- diag(tmp) + 1e-10
  ee <- eigen(tmp)
  U <- ee$vectors
  S <- diag(ee$values)
  M <- ZA %*% ZAtZA.inv.sqrt %*% U
  H <- M %*% solve(diag(pA) + lambdaA * S) %*% t(M)
  ## check
  # abs(Tn - t(R) %*% (H%^%2) %*% R)
  (Tn - s2hat*sum(diag(H%^%2))) / (s2hat * sqrt(2 * sum(diag(H%^%4))))
}


quantile.function.penalized <- function(BA, p.min){
  qnorm(1-p.min)
}

test.statistic.tsls.over.ols <- function(Y, BX, BA, beta){

  ZX <- design.mat(BX)
  lKX <- penalty.mat(BX)
  ZA <- design.mat(BA)
  lKA <- penalty.mat(BA)
  n <- nrow(ZX)

  # 'hat matrix' for the fit onto A's spline basis.
  WA <- ZA%*%solve(t(ZA)%*%ZA+lKA, t(ZA))
  ols <- t(Y-ZX%*%beta)%*%(Y-ZX%*%beta)
  iv <- t(Y-ZX%*%beta)%*%t(WA)%*%WA%*%(Y-ZX%*%beta)

  n * iv / ols
}

quantile.function.tsls.over.ols <- function(BA, p.min){
  ZA <- design.mat(BA)
  dA <- ncol(ZA)
  qchisq(1-p.min,df=dA, ncp=0,lower.tail = TRUE,log.p = FALSE)
}



binary.search.lambda <- function(Y, BX, BA, p.min, test.statistic, quantile.function, eps = 1e-6){

  q <- quantile.function(BA, p.min)
  lmax <- 2
  lmin <- 0
  betaX <- ar.fit(Y, BX, BA, lmax)$coefficients
  t <- test.statistic(Y, BX, BA, betaX)

  while(t > q){
    lmin <- lmax
    lmax <- lmax^2
    betaX <- ar.fit(Y, BX, BA, lmax)$coefficients
    t <- test.statistic(Y, BX, BA, betaX)
  }

  Delta <- lmax-lmin

  while(Delta > eps || Accept == FALSE){
    l <- (lmin+lmax)/2
    betaX <- ar.fit(Y, BX, BA, l)$coefficients
    t <- test.statistic(Y, BX, BA, betaX)
    if(t > q){
      Accept <- FALSE
      lmin <- l
    }
    else {
      Accept <- TRUE
      lmax <- l
    }
    Delta <- lmax-lmin
  }
  lmax
}


predict.AR <- function(object, xnew, ...){
  BX <- object$BX # fitted base object
  basis <- BX$basis
  x <- BX$x
  # spline basis evalutated at data (i.e., design matrix)
  ZX <- getbasismatrix(x, basis, nderiv=0, returnMatrix=FALSE)
  # derivatives of splines evaluated at data
  ZXprime <- getbasismatrix(x, basis, nderiv=1, returnMatrix=FALSE)
  beta <- object$coefficients # fitted coefficients
  w.min <- which.min(x)
  w.max <- which.max(x)
  xmin <- min(x)
  xmax <- max(x)
  f.xmin <- (ZX%*%beta)[w.min]
  f.xmax <- (ZX%*%beta)[w.max]
  fprime.xmin <- (ZXprime%*%beta)[w.min]
  fprime.xmax <- (ZXprime%*%beta)[w.max]
  w.below <- which(xnew <= xmin)
  w.within <- which(xnew > xmin & xnew < xmax)
  w.above <- which(xnew >= xmax)
  out <- numeric(length(xnew))
  #####
  ZXnew.within <- predict(BX, xnew[w.within], ...) # design matrix for new data
  out[w.within] <- ZXnew.within%*%beta # prediction
  out[w.below] <- f.xmin + fprime.xmin*(xnew[w.below]-xmin)
  out[w.above] <- f.xmax + fprime.xmax*(xnew[w.above]-xmax)
  out
}


plot.AR <- function(object, x.new, f.true=NULL){
  X <- NULL
  x <- NULL
  fhat <- NULL
  ftrue <- NULL
  true <- NULL

  Y <- object$Y
  A <- object$A
  BX <- object$BX
  betaX <- object$coefficients
  BA <- object$BA
  betaA <- object$betaA
  dA <- BA$dim.predictor

  dfX <- data.frame(X=BX$x, Y=Y, A = A, fitted=object$fitted.values,
                    true=ifelse(rep(is.null(f.true), length(Y)), NA,
                                f.true(BX$x)))

  if(1){ # new, experimental
    if(object$BX$type == "pspline"){
      xmin <- min(object$BX$breaks)
      xmax <- max(object$BX$breaks)
    }
    if(object$BX$type == "linear"){
      xmin <- min(object$BX$x)
      xmax <- max(object$BX$x)
    }
    if(!is.null(x.new)){
      xmin <- min(xmin, min(x.new))
      xmax <- max(xmax, max(x.new))
    }
    xseq <- seq(xmin, xmax, length.out = 100)
    dfpred <- data.frame(x = xseq, fhat = predict(object, xseq),
                         ftrue = ifelse(rep(is.null(f.true), length(xseq)), NA,
                                        f.true(xseq)))
    pX <- ggplot2::ggplot(dfX, aes(X, Y)) + geom_point(aes(col=A)) +
      coord_cartesian(xlim = range(xseq)) +
      #geom_line(aes(X,fitted),col="red") +
      geom_line(data=dfpred, aes(x,fhat), col="red") +
      theme_bw() + ggtitle(paste("OLS objective = ", round(object$ols,2)))
    if(!is.null(f.true)) pX <- pX +
      geom_line(data=dfpred, aes(x,ftrue), col="darkgreen")
  }

  if(0){ # old, works
    pX <- ggplot2::ggplot(dfX, aes(X, Y)) + geom_point() +
      geom_line(aes(X,fitted),col="red") +
      theme_bw() + ggtitle(paste("OLS objective = ", round(object$ols,2)))
    if(!is.null(f.true)) pX <- pX + geom_line(aes(X,true),col="darkgreen")
  }

  if(dA == 1){
    dfA <- data.frame(A = BA$x,
                      residuals = object$residuals,
                      fitted=object$fitA)
    pA <- ggplot2::ggplot(dfA, aes(A, residuals)) + geom_point() +
      geom_line(aes(A,fitted),col="red") + theme_bw() +
      ggtitle(paste("IV objective = ", round(object$iv,2)))
  }
  if(dA == 2){
    A1 <- BA$x[,1]
    A2 <- BA$x[,2]
    dfA <- expand.grid(A1=seq(min(A1),max(A1),length.out = 100),
                       A2=seq(min(A2),max(A2),length.out = 100))
    ZA <- predict(BA, as.matrix(dfA))
    dfA$fitted <- ZA%*%betaA

    dA1 <- (sort(unique(dfA$A1))[2]-sort(unique(dfA$A1))[1])/2
    dA2 <- (sort(unique(dfA$A2))[2]-sort(unique(dfA$A2))[1])/2

    pA <- ggplot(dfA, aes(xmin=A1-dA1,xmax=A1+dA1,ymin=A2-dA2,ymax=A2+dA2,
                          fill=fitted), col="transparent") +
      geom_rect() +
      scale_fill_gradient2(low="darkblue", high="gold", mid="darkgray") +
      theme_bw() +
      ggtitle(paste("IV objective = ", round(object$iv,2)))
  }

  gridExtra::grid.arrange(pX, pA, widths = c(1, .8), ncol=2)
}


cv.lambda <- function(Y, BX, BA, lambda.star, par.cv){

  n <- length(Y)
  num.folds <- par.cv$num.folds
  optim <- par.cv$optim
  n.grid <- par.cv$n.grid
  X <- BX$x
  dX <- BX$dim.predictor
  A <- BA$x
  dA <- BA$dim.predictor

  ## Construct folds for CV
  if(num.folds=="leave-one-out"){
    folds <- vector("list", n)
    for(i in 1:n){
      folds[[i]] <- i
    }
    num.folds <- length(folds)
  } else if(is.numeric(num.folds)){
    if(num.folds>1){
      folds_tmp <- cvFolds(n, K=num.folds, type="interleaved")
      folds <- vector("list", num.folds)
      for(i in 1:num.folds){
        folds[[i]] <- folds_tmp$subsets[folds_tmp$which == i]
      }
    } else{
      stop("num.folds should be at least 2")
    }
  } else{
    stop("num.folds was specified incorrectly")
  }


  BX.train <- vector("list", num.folds)
  BX.val <- vector("list", num.folds)
  BA.train <- vector("list", num.folds)
  BA.val <- vector("list", num.folds)
  Y.train <- vector("list", num.folds)
  Y.val <- vector("list", num.folds)
  for(i in 1:num.folds){
    train <- (1:n)[-folds[[i]]]
    val <- (1:n)[folds[[i]]]
    # base objects, each holding data corresponding to their respective CV fold
    BX.train[[i]] <- create.base(x=matrix(X,ncol=dX)[train,], basis=BX$basis,
                                 intercept=BX$intercept, type=BX$type)
    BX.val[[i]] <- create.base(x=matrix(X,ncol=dX)[val,], basis=BX$basis,
                               intercept=BX$intercept, type=BX$type)
    BA.train[[i]] <- create.base(x=matrix(A,ncol=dA)[train,], basis=BA$basis,
                                 lambda=BA$lambda, intercept=BA$intercept,
                                 type=BA$type)
    BA.val[[i]] <- create.base(x=matrix(A,ncol=dA)[val,], basis=BA$basis,
                               lambda = BA$lambda, intercept=BA$intercept,
                               type=BA$type)
    Y.train[[i]] <- Y[train]
    Y.val[[i]] <- Y[val]
  }

  # Initialize upper and lower bound for penalty (see 5.4.1 in Functional
  # Data Analysis)
  ZX <- design.mat(BX)
  BX$lambda <- rep(1,dX)
  KX <- penalty.mat(BX)
  ZXnorm <- sum(diag(t(ZX)%*%ZX))
  KXnorm <- sum(diag(KX))

  r <- ZXnorm/KXnorm
  lower.spar <- 0
  upper.spar <- 2
  lower.lambda <- r*256^(3*lower.spar-1)
  upper.lambda <- r*256^(3*upper.spar-1)

  cost.function <- function(spar){
    # here spar is 1-dim (i.e., same penalty for all dimensions).
    # Currently, implementation is unstable if spar higher-dimensional
    lambda <- rep(r*256^(3*spar-1),dX)
    cost <- 0
    for(i in 1:num.folds){
      # impose current penalty lambda
      BX.train[[i]]$lambda <- lambda
      BX.val[[i]]$lambda <- lambda
      # fit on training da
      betahat <- ar.fit(Y=Y.train[[i]], BX=BX.train[[i]], BA=BA.train[[i]],
                        lambda.star)$coefficients
      # AR-objective evaluated on test data
      cost.ar <- ar.objective(Y.val[[i]], BX.val[[i]], BA.val[[i]],
                              lambda.star, betahat)$ar
      cost <- cost + cost.ar
    }
    cost/num.folds
  }

  if(optim=="optimize"){
    solutions.optim <- optimize(cost.function, c(lower.spar, upper.spar))
    spar <- as.numeric(solutions.optim$minimum)
    lambda <- rep(r*256^(3*spar-1),dX)
  }
  if(optim=="grid.search"){
    spar.vec <- seq(lower.spar, upper.spar, length.out=n.grid)
    cost.vec <- rep(NA, length(spar.vec))
    for(k in 1:length(spar.vec)){
      cost.vec[k] <- cost.function(spar.vec[k])
    }
    index.best.spar <- length(cost.vec)
    lowest.cost <- cost.vec[index.best.spar]
    for(j in (index.best.spar-1):1){
      if(cost.vec[j] < lowest.cost){
        index.best.spar <- j
        lowest.cost <- cost.vec[index.best.spar]
      }
    }
    spar <- spar.vec[index.best.spar]
    lambda <- r*256^(3*spar-1)
  }

  # Check whether lambda is attained at boundaries
  if(min(abs(lambda - lower.lambda)) < 10^(-16)){
    warning(paste("There was at least one case in which CV yields",
                  "a lambda at the lower boundary."))
  }
  if(min(abs(lambda - upper.lambda)) < 10^(-16)){
    warning(paste("There was at least one case in which CV yields",
                  "a lambda at the upper boundary."))
  }

  BX$lambda <- lambda
  BX
}
