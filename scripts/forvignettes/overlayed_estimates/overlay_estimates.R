## Non-linear IV estimation. 
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
alphaA <- sqrt(fact/3)
alphaH <- sqrt(2*fact/3)
alphaEps <- 0
## quantiles for different settings
set.seed(1)
N <- 100000
A <- runif(N,-1,1)
H <- runif(N,-1,1)
X <- alphaA*A + alphaH*H + alphaEps*runif(N,-1,1) 
q <- quantile(X, probs = c(.05,.95))
qX <- c(-1,1)*(q[2]-q[1])/2
suppX <- c(-1,1)*(alphaA + alphaH + alphaEps)

n.splines.true <- 4
set.seed(4)
beta <- runif(n.splines.true, -1,1)
fX <- function(x=x.new, qx=qX, beta=beta){
  bx <- ns(x, knots = seq(from=qx[1], to=qx[2], length.out=(n.splines.true+1))[-c(1,n.splines.true+1)], Boundary.knots = qx)
  bx%*%beta
}
A <- runif(n,-1,1)
H <- runif(n,-1,1)
X <- alphaA*A + alphaH*H + alphaEps*runif(n,-1,1) 
Y <- fX(X,qX,beta) + .3*H + .2*runif(n,-1,1)
x.new <- seq(-2,2,length.out = 100)
f.new <- fX(x.new, qX, beta)
plot(X,Y)
lines(x.new, f.new)

pred.frame <- NULL
k <- 0 
Nsim <- 20
for(i in 1:Nsim){
  k <- k+1
  set.seed(k)
  print(paste("sim = ", i))
  
  A <- runif(n,-1,1)
  H <- runif(n,-1,1)
  X <- alphaA*A + alphaH*H + alphaEps*runif(n,-1,1) 
  Y <- fX(X,qX,beta) + .3*H + .2*runif(n,-1,1)
  
  fit.ols <- tryCatch({withTimeout({NILE(Y=Y,X=X,A=A,lambda.star=0,df=df.X,x.new=x.new,plot=FALSE,par.a=list(lambda=0.1))$pred}, timeout = 180, onTimeout = "error")}, 
                      error=function(e){
                        print(paste("ERROR OLS: ", e))
                        rep(NA, length(x.new))
                      })
  
  nile <- tryCatch({withTimeout({NILE(Y=Y,X=X,A=A,lambda.star="test",df=df.X,x.new=x.new,plot=FALSE)}, timeout = 180, onTimeout = "error")}, 
                   error=function(e){
                     print(paste("ERROR NILE: ", e))
                     list(pred = rep(NA, length(x.new)), lambda.star = NA)
                   })
  fit.nile <- nile$pred
  lambda.star <- nile$lambda.star
  
  fit.npregiv <- tryCatch({withTimeout({npregiv(y=Y,z=X,w=A,zeval=x.new)$phi.eval}, timeout = 180, onTimeout = "error")}, 
                          error=function(e){
                            print(paste("ERROR NPREGIV: ", e))
                            rep(NA, length(x.new))
                          })
  
  pred.frame.loop <- data.frame(x = x.new, 
                                fhat = c(fit.ols, fit.nile, fit.npregiv), 
                                ftrue = rep(f.new, 3),
                                method = rep(c("OLS", "NILE", "NPREGIV"), each = length(x.new)),
                                lambda.star = lambda.star,
                                alphaA = alphaA, 
                                alphaH = alphaH, 
                                alphaEps = alphaEps,
                                qXmin = min(qX),
                                qXmax = max(qX),
                                suppXmin = min(suppX),
                                suppXmax = max(suppX),
                                sim = i)
  pred.frame <- rbind(pred.frame, pred.frame.loop)
  write.table(pred.frame, "overlay_estimates_strongconf.txt", quote = FALSE)
}