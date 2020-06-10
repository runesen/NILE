## Non-linear IV estimation. 
SIM <- 21:30
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
alphaA <- alphaH <- alphaEps <- sqrt(fact/3)
## quantiles for different settings
set.seed(1)
N <- 100000
A <- runif(N,-1,1)
H <- runif(N,-1,1)
X <- alphaA*A + alphaH*H + alphaEps*runif(N,-1,1) 
q <- quantile(X, probs = c(.05,.95))
qX <- c(-1,1)*(q[2]-q[1])/2
suppX <- c(-1,1)*(alphaA + alphaH + alphaEps)

n.exp <- 5
curv.max.seq <- seq(0, 4, length.out = n.exp)

n.splines.true <- 4
fX <- function(x, qx, beta, curv){
  wlo <- which(x < qx[1])
  wup <- which(x > qx[2])
  curv.term <- numeric(length(x))
  curv.term[wlo] <- .5*curv[1]*(x[wlo]-qx[1])^2
  curv.term[wup] <- .5*curv[2]*(x[wup]-qx[2])^2
  bx <- ns(x, knots = seq(from=qx[1], to=qx[2], length.out=(n.splines.true+1))[-c(1,n.splines.true+1)], Boundary.knots = qx)
  lin.term <- bx%*%beta
  lin.term + curv.term
}

# curv.max <- .5
# curv <- runif(2,min=-curv.max,max=curv.max)
# beta <- runif(n.splines.true, -1,1)
# x.new <- seq(-2,2,length.out = 100)
# f.new <- fX(x.new,qX,beta,curv)
# plot(x.new, f.new)
# abline(v=qx)

pred.frame <- NULL
x.new <- seq(-2,2,length.out = 100)
k <- max(SIM)*n.exp
for(i in SIM){
  for(j in 1:n.exp){
    k <- k+1
    set.seed(k)
    print(paste("sim = ", i, ", curv.max.seq = ", j))
    curv.max <- curv.max.seq[j]
    curv <- runif(2,min=-curv.max,max=curv.max)
    beta <- runif(n.splines.true, -1,1)
    A <- runif(n,-1,1)
    H <- runif(n,-1,1)
    X <- alphaA*A + alphaH*H + alphaEps*runif(n,-1,1) 
    Y <- fX(X,qX,beta,curv) + .3*H + .2*runif(n,-1,1)
    f.new <- fX(x.new,qX,beta,curv)
    
    
    # fit.nile <- tryCatch({withTimeout({NILE(Y=Y,X=X,A=A,lambda.star="PULSE",df=df.X,x.new=x.new,plot=FALSE)$pred}, timeout = 180, onTimeout = "error")}, 
    #                      error=function(e){
    #                        print(paste("ERROR NILE: ", e))
    #                        rep(NA, length(x.new))
    #                      })
    nile <- tryCatch({withTimeout({NILE(Y=Y,X=X,A=A,lambda.star="test",df=df.X,x.new=x.new,plot=FALSE)}, timeout = 180, onTimeout = "error")}, 
                     error=function(e){
                       print(paste("ERROR NILE: ", e))
                       list(pred = rep(NA, length(x.new)), lambda.star = NA)
                     })
    fit.nile <- nile$pred
    lambda.star <- nile$lambda.star
    
    fit.ols <- tryCatch({withTimeout({NILE(Y=Y,X=X,A=A,lambda.star=0,df=df.X,x.new=x.new,plot=FALSE,par.a=list(lambda=0.1))$pred}, timeout = 180, onTimeout = "error")}, 
                        error=function(e){
                          print(paste("ERROR OLS: ", e))
                          rep(NA, length(x.new))
                        })
    fit.iv <- tryCatch({withTimeout({NILE(Y=Y,X=X,A=A,lambda.star=Inf,df=df.X,x.new=x.new,plot=FALSE)$pred}, timeout = 180, onTimeout = "error")}, 
                       error=function(e){
                         print(paste("ERROR IV: ", e))
                         rep(NA, length(x.new))
                       })
    fit.npregiv <- tryCatch({withTimeout({npregiv(y=Y,z=X,w=A,zeval=x.new)$phi.eval}, timeout = 180, onTimeout = "error")}, 
                            error=function(e){
                              print(paste("ERROR NPREGIV: ", e))
                              rep(NA, length(x.new))
                            })
    
    pred.frame.loop <- data.frame(x = x.new, 
                                  fhat = c(fit.ols, fit.nile, fit.iv, fit.npregiv), 
                                  ftrue = rep(f.new, 4),
                                  method = rep(c("OLS", "NILE", "IV", "NPREGIV"), each = length(x.new)),
                                  lambda.star = lambda.star,
                                  alphaA = alphaA, 
                                  alphaH = alphaH, 
                                  alphaEps = alphaEps,
                                  qXmin = min(qX),
                                  qXmax = max(qX),
                                  suppXmin = min(suppX),
                                  suppXmax = max(suppX),
                                  curv.max = curv.max,
                                  sim = i)
    pred.frame <- rbind(pred.frame, pred.frame.loop)
    write.table(pred.frame, paste0("nz_curvature_strong", min(SIM), ".txt"), quote = FALSE)
  }
}