## impossibility result for X = g(A) + H + eps
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## ---- plot-impossibility ----
library(ggplot2)
library(gridExtra)
library(grid)
theme_set(theme_bw())


join.smoothly <- function(x,x1,x2,f1,f2,f1prime,f2prime){
  y <- c(f1,f2,f1prime,f2prime)
  A <- matrix(c(1,x1,x1^2,x1^3,
                1,x2,x2^2,x2^3,
                0,1,2*x1,3*x1^2,
                0,1,2*x2,3*x2^2), 4,4, byrow=TRUE)
  b <- solve(A,y)
  sapply(x, function(xx) sum(b*(xx^(0:3))))
}

x <- seq(1,2,.01)
# plot(x,join.smoothly(x,0,1,-1,1))

g.global <- function(a) a+1
g.global.prime <- function(a) 1

g <- function(a, supp = c(-1,1), eps = .5, C=5){
  out <- numeric(length(a))
  w.below <- which(a <= min(supp)-eps)
  w.midlow <- which(a > min(supp)-eps & a < min(supp))
  w.within <- which(a >= min(supp) & a <= max(supp))
  w.midup <- which(a > max(supp) & a < max(supp)+eps)
  w.above <- which(a >= max(supp)+eps)
  if(length(w.below)>0) out[w.below] <- out[w.above] <- C
  if(length(w.midlow)>0) out[w.midlow] <- join.smoothly(a[w.midlow],supp[1]-eps,supp[1],C,g.global(supp[1]),0,g.global.prime(supp[2]))
  if(length(w.within)>0) out[w.within] <- g.global(a[w.within])
  if(length(w.midup)>0) out[w.midup] <- join.smoothly(a[w.midup],supp[2], supp[2]+eps,g.global(supp[2]),C, g.global.prime(supp[2]),0)
  if(length(w.above)>0) out[w.above] <- C
  out
}
a <- seq(-2,2,length.out = 1000)
# plot(a,g(a,eps=.5),type="l",lwd=5)

n <- 200
beta <- 1
set.seed(1)
# observational
A <- runif(n,-1,1)
H <- rnorm(n)
X <- g(A, eps=.5, C=5) + H + .5*rnorm(n)
Y <- beta*X + H + .5*rnorm(n)

bOLS <- lm(Y~X)$coefficients

# interventional
Ai <- rep(1.5,n)
Hi <- rnorm(n)
Xi <- g(Ai, eps=.5, C=5) + Hi + .5*rnorm(n)
Yi <- beta*Xi + Hi + .5*rnorm(n)

df <- data.frame(X=c(X,Xi), A=c(A,Ai), Y=c(Y,Yi),
                 setting = factor(rep(c("observational", "interventional"), each=n),
                                  levels = c("observational", "interventional")))
x.seq <- seq(min(df$X), max(df$X), length.out = 100)
df.line <- data.frame(x = rep(x.seq, 2),
                      y = c(bOLS[1] + bOLS[2]*x.seq,
                            x.seq),
                      model = factor(rep(c("candidate", "causal"), each = length(x.seq))))

pXY <- ggplot() +
  geom_point(data=df, aes(X,Y, col=setting, shape = setting), alpha=.7, size=3) +
  scale_color_manual(values = c("black", "red")) +
  geom_line(data=subset(df.line, model == "candidate"), aes(x,y),size=1, col = "blue") +
  geom_line(data=subset(df.line, model == "causal"), aes(x,y),size=1, col = "#009E73", lty = "dashed") +
  annotate(geom = "text", x = -1, y = 5, label = "training data") +
  annotate(geom = "text", x = 5, y = -1, label = "test data") +
  annotate(geom = "text", x = 3, y = 10, label = "candidate model") +
  # annotate(geom = "text", x = 6, y = 1, label = "causal model") +
  annotate(geom = "text", x = 7.5, y = 6.5, label = "f", col = "#009E73") +
  geom_segment(aes(x = -.5, y = 4, xend = 0, yend =2),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 5, y = 0, xend = 5, yend = 2),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 5, y = 10, xend = 6, yend = 10),
               arrow = arrow(length = unit(0.3, "cm"))) +
  # geom_segment(aes(x = 6, y = 2, xend = 6.6, yend = 6.3),
  #              arrow = arrow(length = unit(0.3, "cm"))) +
  xlab("x") + ylab("y") +
  # geom_abline(intercept = bOLS[1], slope = bOLS[2]) +
  theme(legend.position = "none",
        text = element_text(size=15),
        plot.title = element_text(size = 14, hjust=.5))

a.seq <- seq(-2,2,length.out = 1000)
df.A <- data.frame(a = a.seq, g = g(a.seq))
df.A$seq <- factor(sapply(1:nrow(df.A), function(i) ifelse(df.A$a[i] < -1, 1, ifelse(df.A$a[i] <= 1, 2, 3))))

grid.a <- seq(-.9, .9, .2)
hist.a <- sapply(grid.a, function(g) sum(abs(A-g) < .1))
m <- max(hist.a)
hist.a <- 4*hist.a/m
df.hist <- data.frame(a=c(grid.a, 1.55), h = c(hist.a, 4), data = c(rep("1",length(grid.a)), "2"))

pA <- ggplot(df.A, aes(a,g)) +
  geom_segment(aes(x = 1, y = .8, xend = 1.4, yend = .1),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_line(aes(lty = seq), size=1, col = "black") +
  # geom_line(size=1, col = "black") +
  xlab("a") + ylab("") +
  geom_bar(data=df.hist, aes(x=a,y=h,fill=data), alpha=.2, stat = "identity") +
  annotate(geom = "text", x = -1.6, y = 4.5, label = "g") +
  geom_segment(aes(x = -.8, y = 3.5, xend = -1, yend = 3.5),
               arrow = arrow(length = unit(0.3, "cm"), angle = 90)) +
  geom_segment(aes(x = .8, y = 3.5, xend = 1, yend = 3.5),
               arrow = arrow(length = unit(0.3, "cm"), angle = 90)) +
  annotate(geom = "text", x = 0, y = 3.5, label = "training support") +
  annotate(geom = "text", x = .9, y = 1, label = "intervention") +
  scale_fill_manual(values= c("black", "red")) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed")) +
  theme(legend.position = "none",
        text = element_text(size=15),
        plot.title = element_text(size = 14, hjust=.5))

all.A <- arrangeGrob(pA, pXY, ncol = 2,
                     top = textGrob("support-extending interventions on A", gp=gpar(fontsize=15)))

# # grid.arrange(pA, pXY, ncol=2)
#
# pdf("../figures/impossibility_AtoX_nonlinear.pdf", width = 8, height = 3.5)
# grid.arrange(pA, pXY, ncol=2, top = "support-extending interventions on A")
# dev.off()



#######################
## interventions on X
#######################
f.global <- function(x) -x
f.global.prime <- function(x) -1

f <- function(a, supp = c(-1,1), eps = 1, C=c(-5,5)){
  out <- numeric(length(a))
  w.below <- which(a <= min(supp)-eps)
  w.midlow <- which(a > min(supp)-eps & a < min(supp))
  w.within <- which(a >= min(supp) & a <= max(supp))
  w.midup <- which(a > max(supp) & a < max(supp)+eps)
  w.above <- which(a >= max(supp)+eps)
  if(length(w.below)>0) out[w.below] <- C[1]
  if(length(w.midlow)>0) out[w.midlow] <- join.smoothly(a[w.midlow],supp[1]-eps,supp[1],C[1],f.global(supp[1]),20,f.global.prime(supp[2]))
  if(length(w.within)>0) out[w.within] <- f.global(a[w.within])
  if(length(w.midup)>0) out[w.midup] <- join.smoothly(a[w.midup],supp[2],supp[2]+eps,f.global(supp[2]),C[2], f.global.prime(supp[2]),20)
  if(length(w.above)>0) out[w.above] <- C[2]
  out
}
x <- seq(-2,2,length.out = 1000)
# plot(x,f(x,eps=1),type="l",lwd=5)

n <- 200
beta <- 1
set.seed(1)
# observational
A <- runif(n,-1,1)
H <- runif(n,-1,1)
X <- .4*H + .4*A + .2*runif(n,-1,1)
Y <- f(X) + .5*H + .5*rnorm(n)

Ai <- runif(n,-1,1)
Hi <- runif(n,-1,1)
Xi <- runif(n,1.5,2)
Yi <- f(Xi) + Hi + .5*rnorm(n)


df1 <- data.frame(X=c(X,Xi), A=c(A,Ai), Y=c(Y,Yi),
                 setting = factor(rep(c("observational", "interventional"), each=n),
                                  levels = c("observational", "interventional")))
x.seq <- seq(-2,2,length.out = 1000)
df1.line <- data.frame(x = rep(x.seq, 2),
                      y = c(-x.seq, f(x.seq)),
                      model = factor(rep(c("candidate", "causal"), each = length(x.seq))))

pXY1 <- ggplot() +
  geom_point(data=df1, aes(X,Y, col=setting, shape = setting), alpha=.7, size=3) +
  scale_color_manual(values = c("black", "red")) +
  geom_line(data=subset(df1.line, model == "candidate"), aes(x,y),size=1, col = "blue") +
  geom_line(data=subset(df1.line, model == "causal"), aes(x,y),size=1, col = "#009E73", lty = "dashed") +
  annotate(geom = "text", x = 0, y = 3.5, label = "training data") +
  annotate(geom = "text", x = 1, y = 5.5, label = "test data") +
  annotate(geom = "text", x = -1.2, y = 5.5, label = "candidate model") +
  annotate(geom = "text", x = -1.5, y = -1, label = "f", col = "#009E73") +
  ggtitle("support-extending interventions on X") +
  xlab("x") + ylab("y") +
  # annotate(geom = "text", x = 6, y = 1, label = "causal model") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 1.5),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 1.2, y = 5, xend = 1.6, yend = 4),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = -1.2, y = 5, xend = -1.6, yend = 2.5),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = -.8, y = -5, xend = -1, yend = -5),
               arrow = arrow(length = unit(0.3, "cm"), angle = 90)) +
  geom_segment(aes(x = .8, y = -5, xend = 1, yend = -5),
               arrow = arrow(length = unit(0.3, "cm"), angle = 90)) +
  geom_segment(aes(x = 1.5, y = -5, xend = 2, yend = -5),
               arrow = arrow(length = unit(0.3, "cm"), angle = 90)) +
  geom_segment(aes(x = 2, y = -5, xend = 1.5, yend = -5),
               arrow = arrow(length = unit(0.3, "cm"), angle = 90)) +
  annotate(geom = "text", x = 0, y = -5, label = "training support") +
  coord_cartesian(ylim=c(-5.5,6)) +
  annotate(geom = "text", x = 1.5, y = -3.5, label = "test support") +
  # geom_segment(aes(x = 1, y = -3.5, xend = 1.5, yend = -4),
  #              arrow = arrow(length = unit(0.3, "cm"))) +
  # xlab("x") + ylab("y") +
  # geom_abline(intercept = bOLS[1], slope = bOLS[2]) +
  theme(legend.position = "none",
        text = element_text(size=15),
        plot.title = element_text(size = 14, hjust=.5))

grid.arrange(pXY1, all.A,ncol=2, widths = c(1,2))

## ---- private ----
pdf("../figures/impossibility_all.pdf", width = 12, height = 3.5)
dev.off()
