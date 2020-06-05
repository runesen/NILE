setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(R.utils)
library(splines)
library(ggplot2)
theme_set(theme_bw())


alphaA <- alphaH <- alphaEps <- 1/sqrt(3)
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
fX <- function(x=x.new, qx=qX, beta=beta){
  bx <- ns(x, knots = seq(from=qx[1], to=qx[2], length.out=(n.splines.true+1))[-c(1,n.splines.true+1)], Boundary.knots = qx)
  bx%*%beta
}

x.new <- seq(-2,2,length.out = 100)
n <- 200
Nsim <- 18
df.data <- df.ftrue <- NULL
k <- 0
for(i in 1:Nsim){
  k <- k+1
  set.seed(k)
  beta <- runif(n.splines.true, -1,1)
  A <- runif(n,-1,1)
  H <- runif(n,-1,1)
  X <- alphaA*A + alphaH*H + alphaEps*runif(n,-1,1) 
  Y <- fX(X,qX,beta) + .3*H + .2*runif(n,-1,1)
  
  f.new <- fX(x.new,qX,beta)
  df.data <- rbind(df.data, data.frame(x=X, y=Y, sim = i))
  df.ftrue <- rbind(df.ftrue, data.frame(x=x.new, f=f.new, sim = i))
}

p <- ggplot() + 
  geom_point(data=df.data, aes(x,y), alpha = .7) + 
  geom_line(data=df.ftrue, aes(x, f), col = "#009E73", size = 1) + 
  geom_vline(xintercept = qX, lty=2) + 
  facet_wrap(.~sim, nrow=3)
p

pdf("f_samples.pdf", width = 10, height = 5)
p
dev.off()

