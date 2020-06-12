setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## ---- plot-varying-inst ----
library(ggplot2)
library(reshape2)
library(gridExtra)
theme_set(theme_bw())
pred.frame <- rbind(read.table("varying_instrument_fit1.txt", header = TRUE),
                    read.table("varying_instrument_fit11.txt", header = TRUE),
                    read.table("varying_instrument_fit21.txt", header = TRUE),
                    read.table("varying_instrument_fit31.txt", header = TRUE),
                    read.table("varying_instrument_fit41.txt", header = TRUE),
                    read.table("varying_instrument_fit51.txt", header = TRUE),
                    read.table("varying_instrument_fit61.txt", header = TRUE),
                    read.table("varying_instrument_fit71.txt", header = TRUE),
                    read.table("varying_instrument_fit81.txt", header = TRUE),
                    read.table("varying_instrument_fit91.txt", header = TRUE))

pred.frame <- subset(pred.frame, experiment <= 3)

n <- nrow(pred.frame)
n.x <- length(unique(pred.frame$x))
n.m <- length(unique(pred.frame$method))
n.sim <- length(unique(pred.frame$sim))
n.exp <- n/(n.x*n.m*n.sim)

var_xi_Y <- (1/3)*(.3^2+.2^2)

# pred.frame <- subset(pred.frame, method != "IV")
pred.frame$sq.err <- (pred.frame$fhat-pred.frame$ftrue)^2
pred.frame$x.abs <- round(abs(pred.frame$x),2)
pred.frame$alphaA <- round(pred.frame$alphaA,3)
pred.frame$alphaH <- round(pred.frame$alphaH,3)
pred.frame$alphaEps <- round(pred.frame$alphaEps,3)
pred.frame$qXmax <- round(pred.frame$qXmax,3)
pred.frame$lambda.star <- round(pred.frame$lambda.star,3)
error.frame.x.abs <- aggregate(sq.err ~ method + alphaA + qXmax + x.abs + sim, pred.frame, mean)
error.frame.x.abs <- error.frame.x.abs[order(error.frame.x.abs$method, error.frame.x.abs$alphaA, error.frame.x.abs$sim),]
n.x.abs <- length(unique(error.frame.x.abs$x.abs))
error.frame.x.abs$worst.case.mse <- sapply(1:nrow(error.frame.x.abs), function(i) max(error.frame.x.abs$sq.err[(1+floor((i-1)/n.x.abs)*n.x.abs):i])) + var_xi_Y

error.frame.x.abs$alphaA <- factor(error.frame.x.abs$alphaA,
                                   levels = sort(unique(error.frame.x.abs$alphaA)),
                                   labels = c(expression(paste("(", alpha[A], ",", alpha[H], ",", alpha[epsilon], ")  =  (", 0,",", sqrt(1/3),",", sqrt(2/3), ")")),
                                              expression(paste("(", alpha[A],",", alpha[H],",", alpha[epsilon], ")  =  (", sqrt(1/6),",", sqrt(1/3),",", sqrt(1/2), ")")),
                                              expression(paste("(", alpha[A], ",",alpha[H],",", alpha[epsilon], ")  =  (", sqrt(1/3),",", sqrt(1/3),",", sqrt(1/3), ")"))))

error.frame.x.abs.avg <- aggregate(worst.case.mse ~ method + alphaA + qXmax + x.abs, error.frame.x.abs, mean)
tmp <- subset(error.frame.x.abs.avg, method == "OLS")
tmp$worst.case.mse <- var_xi_Y
tmp$method <- "ORACLE"
error.frame.x.abs.avg <- rbind(error.frame.x.abs.avg, tmp)

error.frame.x.abs.avg$sim <- max(error.frame.x.abs$sim+1)
error.frame.x.abs.avg$sq.err <- NA
error.frame.x.abs.avg$avg <- TRUE
error.frame.x.abs$avg <- FALSE
error.frame.x.abs <- rbind(error.frame.x.abs, error.frame.x.abs.avg)

error.frame.x.abs$method <- factor(error.frame.x.abs$method,
                                   levels = c("OLS", "NILE", "NPREGIV", "ORACLE"))


qX.frame <- subset(error.frame.x.abs.avg, method == "OLS" & x.abs == 2)
rects1 <- qX.frame
rects1$xstart <- 0
rects1$xend <- rects1$qXmax
rects2 <- qX.frame
rects2$xstart <- rects2$qXmax
rects2$xend <- 2
rects <- rbind(rects1, rects2)
rects$col <- factor(rep(1:2, each=n.exp))


p <- ggplot() +
  # geom_density(data=df.plt, aes(x=x.plt),alpha=.5) +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=col), alpha=.1) +
  geom_line(data=error.frame.x.abs, aes(x.abs, worst.case.mse, col=method, lty = method, group = interaction(method,sim,avg), alpha=avg, size = avg)) +
  #geom_line(data=error.frame.x.abs.avg, aes(x.abs, worst.case.mse, col=method, lty = method), size=1) +
  # facet_wrap(. ~ alphaA + alphaH + alphaEps, ncol=5) +
  # geom_vline(aes(xintercept = qXmax),lty=2) +
  # geom_hline(yintercept = var_xi_Y,lty=2, col = "black") +
  facet_grid(. ~ alphaA, labeller = label_parsed) +
  # ggtitle(expression(paste("constant instrument strength alpha_A = ", 1/sqrt(3)))) +
  # ggtitle(expression(paste("constant instrument strength ", alpha[A], " = ", sqrt(1/3)))) +
  # gtitle("generalization performance for different confounding strength") +
  scale_size_manual(values = c(.2, 1)) +
  scale_alpha_manual(values = c(.2, 1)) +
  scale_color_manual(values = c("#D55E00",
                                #"#0072B2",
                                "#0072B2", "#CC79A7", "#009E73")) +
  scale_linetype_manual(values = c(rep("solid",3), "dashed")) +
  scale_fill_manual(values = c("black", "transparent")) +
  coord_cartesian(xlim = c(0,2), ylim = c(0,.3)) +
  xlab("generalization to interventions up to strength") + ylab("worst-case MSE")  +
  theme(plot.title = element_text(hjust=0.5),
        #legend.position = "none",
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA)) +
  guides(fill="none", lty= "none", size = "none", alpha = "none")
p



pdf("varying_instrument_with_variability.pdf", width = 8, height = 2.5)
p
dev.off()

## ---- private ----
pdf("~/Google_Drive/phd/nonlinearIV/figures/varying_instrument_with_variability.pdf", width = 8, height = 2.5)
p
dev.off()

