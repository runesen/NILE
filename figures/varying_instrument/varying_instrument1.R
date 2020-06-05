## Non-linear IV estimation. 
SIM <- 1:10
onServer <- TRUE
if(!onServer){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  library(ggplot2)
}
library(NILE)
library(R.utils)
library(splines)
library(np)
# X has fixed variance 1/3
n <- 200
fact <- 1 # three times the desired variance of X
df.X <- 50
## missing settings
alphaA.seq <- sqrt(seq(0,2/3,length.out=5))
alphaH.seq <- rep(sqrt(fact/3),5)
alphaEps.seq <- sqrt(fact - alphaA.seq^2 - alphaH.seq^2)
## check
alphaA.seq^2 + alphaH.seq^2 + alphaEps.seq^2
## quantiles for different settings
set.seed(1)
n.exp <- length(alphaA.seq)
qX.mat <- sapply(1:n.exp, function(i){
  N <- 100000
  alphaA <- alphaA.seq[i]
  alphaH <- alphaH.seq[i]
  alphaEps <- alphaEps.seq[i]
  A <- runif(N,-1,1)
  H <- runif(N,-1,1)
  X <- alphaA*A + alphaH*H + alphaEps*runif(N,-1,1) 
  q <- quantile(X, probs = c(.05,.95))
  c(-1,1)*(q[2]-q[1])/2
})
suppX.mat <- sapply(1:n.exp, function(i){
  alphaA <- alphaA.seq[i]
  alphaH <- alphaH.seq[i]
  alphaEps <- alphaEps.seq[i]
  c(-1,1)*(alphaA + alphaH + alphaEps)
})

n.splines.true <- 4
fX <- function(x, qx, beta){
  bx <- ns(x, knots = seq(from=qx[1], to=qx[2], length.out=(n.splines.true+1))[-c(1,n.splines.true+1)], Boundary.knots = qx)
  bx%*%beta
}

pred.frame <- NULL
meta.frame <- NULL
x.new <- seq(-2,2,length.out = 100)
k <- max(SIM)*n.exp
for(i in SIM){
  for(j in 1:n.exp){
    k <- k+1
    set.seed(k)
    print(paste("sim = ", i, ", alphaseq = ", j))
    alphaA <- alphaA.seq[j]
    alphaEps <- alphaEps.seq[j]
    alphaH <- alphaH.seq[j]
    qX <- qX.mat[,j]
    suppX <- suppX.mat[,j]
    
    beta <- runif(n.splines.true, -1,1)
    A <- runif(n,-1,1)
    H <- runif(n,-1,1)
    X <- alphaA*A + alphaH*H + alphaEps*runif(n,-1,1) 
    Y <- fX(X,qX,beta) + .3*H + .2*runif(n,-1,1)
    f.new <- fX(x.new,qX,beta)
    
    ols <- tryCatch({
      t1 <- Sys.time()
      res <- withTimeout({NILE(Y=Y,X=X,A=A,lambda.star=0,df=df.X,x.new=x.new,plot=FALSE,par.a=list(lambda=0.1))}, timeout = 180, onTimeout = "silent")
      t2 <- Sys.time()
      if(!is.null(res)){
        out <- list(time = as.numeric(t2-t1), 
                    fit = res$pred, 
                    timeout = FALSE, 
                    error = FALSE)
      }
      if(is.null(res)){
        out <- list(time = 180, 
                    fit = rep(NA, length(x.new)), 
                    timeout = TRUE, 
                    error = FALSE)
      }
      out
    }, 
    error=function(e){
      print(paste("ERROR OLS: ", e))
      list(time = NA, fit = rep(NA, length(x.new)), timeout = NA, error = TRUE)
    })
    
    
    
    
    nile <- tryCatch({
      t1 <- Sys.time()
      res <- withTimeout({NILE(Y=Y,X=X,A=A,lambda.star="test",df=df.X,x.new=x.new,plot=FALSE)}, timeout = 180, onTimeout = "silent")
      t2 <- Sys.time()
      if(!is.null(res)){
        out <- list(time = as.numeric(t2-t1), 
                    lambda = res$lambda.star, 
                    fit = res$pred, 
                    timeout = FALSE, 
                    error = FALSE)
      }
      if(is.null(res)){
        out <- list(time = 180, 
                    lambda = NA, 
                    fit = rep(NA, length(x.new)), 
                    timeout = TRUE, 
                    error = FALSE)
      }
      out
    }, 
    error=function(e){
      print(paste("ERROR NILE: ", e))
      list(time = NA, lambda = NA, fit = rep(NA, length(x.new)), timeout = NA, error = TRUE)
    })
    
    
    
    npregiv <- tryCatch({
      t1 <- Sys.time()
      res <- withTimeout({npregiv(y=Y,z=X,w=A,zeval=x.new)}, timeout = 180, onTimeout = "silent")
      t2 <- Sys.time()
      if(!is.null(res)){
        out <- list(time = as.numeric(t2-t1), 
                    fit = res$phi.eval, 
                    timeout = FALSE, 
                    error = FALSE)
      }
      if(is.null(res)){
        out <- list(time = 180, 
                    fit = rep(NA, length(x.new)), 
                    timeout = TRUE, 
                    error = FALSE)
      }
      out
    }, 
    error=function(e){
      print(paste("ERROR OLS: ", e))
      list(time = NA, fit = rep(NA, length(x.new)), timeout = NA, error = TRUE)
    })
    
    meta.frame.loop <- data.frame(method = c("OLS", "NILE", "NPREGIV"), 
                                  time = c(ols$time, nile$time, npregiv$time), 
                                  timeout = c(ols$timeout, nile$timeout, npregiv$timeout),
                                  error = c(ols$error, nile$error, npregiv$error),
                                  lambda.star = nile$lambda,
                                  alphaA = alphaA, 
                                  alphaH = alphaH, 
                                  alphaEps = alphaEps,
                                  qXmin = min(qX),
                                  qXmax = max(qX),
                                  suppXmin = min(suppX),
                                  suppXmax = max(suppX),
                                  experiment = j, 
                                  sim = i)
    meta.frame <- rbind(meta.frame, meta.frame.loop)
    
    pred.frame.loop <- data.frame(x = x.new, 
                                  fhat = c(ols$fit, nile$fit, npregiv$fit), 
                                  ftrue = rep(f.new, 3),
                                  method = rep(c("OLS", "NILE", "NPREGIV"), each = length(x.new)),
                                  lambda.star = nile$lambda,
                                  alphaA = alphaA, 
                                  alphaH = alphaH, 
                                  alphaEps = alphaEps,
                                  qXmin = min(qX),
                                  qXmax = max(qX),
                                  suppXmin = min(suppX),
                                  suppXmax = max(suppX),
                                  experiment = j, 
                                  sim = i)
    pred.frame <- rbind(pred.frame, pred.frame.loop)
    
    write.table(meta.frame, paste0("varying_instrument_meta", min(SIM), ".txt"), quote = FALSE)
    write.table(pred.frame, paste0("varying_instrument_fit", min(SIM), ".txt"), quote = FALSE)
  }
}