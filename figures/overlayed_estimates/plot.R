setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(R.utils)
library(splines)
library(np)
theme_set(theme_bw())
pred.frame <- read.table("overlay_estimates_strongconf.txt", header = TRUE)

fact <- 1 # three times the desired variance of X
## missing settings
alphaA <- sqrt(1/3)
alphaH <- sqrt(2/3)
alphaEps <- 0
## quantiles for different settings
qX <- c(pred.frame$qXmin[1], pred.frame$qXmax[1])
suppX <- c(pred.frame$suppXmin[1], pred.frame$suppXmax[1])

n.splines.true <- 4
fX <- function(x=x.new, qx=qX, beta=beta){
  bx <- ns(x, knots = seq(from=qx[1], to=qx[2], length.out=(n.splines.true+1))[-c(1,n.splines.true+1)], Boundary.knots = qx)
  bx%*%beta
}

n.plt <- 200
set.seed(4)
beta <- runif(n.splines.true, -1,1)
A1.plt <- runif(n.plt,-1,1)
H1.plt <- runif(n.plt,-1,1)
X1.plt <- alphaA*A1.plt + alphaH*H1.plt + alphaEps*runif(n.plt,-1,1) 
Y1.plt <- fX(X1.plt,qX,beta) + .3*H1.plt + .2*runif(n.plt,-1,1)

plt.frame <- data.frame(x=X1.plt, y=Y1.plt)
# ftrue.frame <- data.frame(x=x.pred, ftrue=f.true)

pred.frame$method <- factor(pred.frame$method, levels = c("OLS", "NILE", "NPREGIV"))

p <- ggplot(pred.frame) + 
  geom_point(data=plt.frame, aes(x,y), alpha=.5) +
  geom_vline(xintercept = qX, lty = 2) + 
  geom_line(aes(x, fhat, col = method, group = sim), size = .2) + 
  geom_line(aes(x, ftrue), col = "#009E73", size=1, lty = "dashed") + 
  scale_color_manual(values = c("#D55E00","#0072B2", "#CC79A7")) + 
  # ggtitle("")
  facet_grid(.~method) + 
  coord_cartesian(xlim=c(-2,2), ylim = c(-1, 1)) + 
  theme(legend.position = "none") 
p


pdf("overlay_estimates.pdf", width = 7, height = 2.5)
p
dev.off()


pdf("~/Google_Drive/phd/nonlinearIV/figures/overlay_estimates.pdf", width = 7, height = 2.5)
p
dev.off()
