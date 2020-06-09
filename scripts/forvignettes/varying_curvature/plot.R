setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
theme_set(theme_bw())
pred.frame <- rbind(read.table("nz_curvature_strong1.txt", header = TRUE), 
                    read.table("nz_curvature_strong11.txt", header = TRUE), 
                    read.table("nz_curvature_strong21.txt", header = TRUE), 
                    read.table("nz_curvature_strong31.txt", header = TRUE), 
                    read.table("nz_curvature_strong41.txt", header = TRUE), 
                    read.table("nz_curvature_strong51.txt", header = TRUE), 
                    read.table("nz_curvature_strong61.txt", header = TRUE), 
                    read.table("nz_curvature_strong71.txt", header = TRUE), 
                    read.table("nz_curvature_strong81.txt", header = TRUE), 
                    read.table("nz_curvature_strong91.txt", header = TRUE))


n <- nrow(pred.frame)
n.x <- length(unique(pred.frame$x))
n.m <- length(unique(pred.frame$method))
n.sim <- length(unique(pred.frame$sim))
n.exp <- n/(n.x*n.m*n.sim)

var_xi_Y <- (1/3)*(.3^2+.2^2)

# pred.frame <- subset(pred.frame, method != "IV")
pred.frame$sq.err <- (pred.frame$fhat-pred.frame$ftrue)^2
pred.frame$x.abs <- round(abs(pred.frame$x),2)
pred.frame$curv.max <- round(pred.frame$curv.max,1)
error.frame.x.abs <- aggregate(sq.err ~ method + curv.max + qXmax + x.abs + sim, pred.frame, mean)
error.frame.x.abs <- error.frame.x.abs[order(error.frame.x.abs$method, error.frame.x.abs$curv.max, error.frame.x.abs$sim, error.frame.x.abs$x.abs),]
n.x.abs <- length(unique(error.frame.x.abs$x.abs))
error.frame.x.abs$worst.case.mse <- sapply(1:nrow(error.frame.x.abs), function(i) max(error.frame.x.abs$sq.err[(1+floor((i-1)/n.x.abs)*n.x.abs):i])) + var_xi_Y

error.frame.x.abs$curv.max <- factor(error.frame.x.abs$curv.max, 
                                         levels = sort(unique(error.frame.x.abs$curv.max)), 
                                         labels = sapply(sort(unique(error.frame.x.abs$curv.max)), function(c) paste0("max. curvature = ", c)))

error.frame.x.abs.avg <- aggregate(worst.case.mse ~ method + curv.max + qXmax + x.abs, error.frame.x.abs, mean)
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
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=col), alpha=.1) +
  geom_line(data=error.frame.x.abs, aes(x.abs, worst.case.mse, col=method, lty = method, group = interaction(method,sim,avg), alpha=avg, size = avg)) + 
  facet_grid(. ~ curv.max) + 
  scale_size_manual(values = c(.15, .85)) + 
  scale_alpha_manual(values = c(.15, .85)) + 
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


pdf("varying_extrapolation_curvature_with_variability.pdf", width = 10, height = 2.2)
p
dev.off()

pdf("~/Google_Drive/phd/nonlinearIV/figures/varying_extrapolation_curvature_with_variability.pdf", width = 10, height = 2.2)
p
dev.off()



