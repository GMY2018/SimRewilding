## ---------------------------------------
## Graph and tables for simulation study
## ---------------------------------------

library(ggplot2)
library(RColorBrewer)
# library(gridExtra)


## (1) Heat maps for offline changepoint detection (AMOC) result

## Trend change
# Load the saved simulation result
load("sim_result_trend_change_AMOC.RData")
# load("sim_result_trend_change_AMOC_AR1.RData")  
# Use the same method for plotting the AR1 scenario

npre = 24
npost.vec = seq(24, 120, by=12) 
level = 10
trend.increase = round(level * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 
simN = 1000

# The test statistics
Tmax.mat = sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.v1[[i]][[j]][[1]][1, ]
  })
})

# The idenfied changepoint
Cpt.mat = sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.v1[[i]][[j]][[1]][2, ]
  })
})

# The delay in effect (simulated)
delay.mat = sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult.v1[[i]][[j]][[2]]
  })
})

# The critical values of the test statistics 
nl = round(length(npost.vec) / 2)
criticalVal = rep(getCritical.v1[1:2], times=c(nl, length(npost.vec)-nl))
criticalMat = rep(criticalVal, times=length(trend.increase))

# True positive rate (TPR): how many time the cpt is correctly detected
# H0: there is no cpt at tau, H1: there is a cpt at tau
Cpt.mat2 = Cpt.mat
Cpt.mat2[Tmax.mat < criticalMat] = -100
# give a large negative value to mean no cpt is detected
Delay.mat = apply(delay.mat, MARGIN=2, FUN=rep, each=length(npost.vec))
gap.mat = abs(Cpt.mat2 - (Delay.mat + npre))
Correct.count = apply(gap.mat <= 12, MARGIN=1, FUN=sum)
TPR.mat = matrix(Correct.count / simN, nrow=length(npost.vec),
                 ncol=length(trend.increase))

# False discovery rate (FDR): among all positive claims, how many get the cpt wrongly
# H0: there is no cpt at tau, H1: there is a cpt at tau
Positive.count = apply(Tmax.mat > criticalMat, MARGIN=1, FUN=sum)
Wrong.count = Positive.count - Correct.count
FDR.mat = matrix(Wrong.count / Positive.count, nrow=length(npost.vec),
                 ncol=length(trend.increase))

## Create the heat maps
trend.percent = c(0.05, seq(0.1, 1, by=0.1)) * 100   # for the coordinates of the maps

plot.df = expand.grid(x=as.factor(npost.vec),
                      y=as.factor(trend.percent[1:7]))
# If we use numeric values in grid expansion, there will be uneven intervals,
# and this will result in space in the heat maps. So we use factors
plot.df$TPR = as.vector(TPR.mat[, 1:7])
plot.df$FDR = as.vector(FDR.mat[, 1:7])

p1 = ggplot(data=plot.df, aes(x=x, y=y, fill=TPR)) + 
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd") + 
  labs(x = "Post intervention length", y = "Size of change (%)") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey95"),
        # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "white", fill=NA),
        legend.position = 'right',
        legend.title = element_text(size=10),
        axis.title = element_text(size=10))

p2 = ggplot(data=plot.df, aes(x=x, y=y, fill=FDR)) + 
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd") + 
  labs(x = "Post intervention length", y = "Size of change (%)",) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey95"),
        # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "white", fill=NA),
        legend.position = 'right',
        legend.title = element_text(size=10),
        axis.title = element_text(size=10))

ggsave(filename="Sim_trend_AMOC_large_1.png", plot=p1, 
       device = "png", width=4, height=3, units="in")
ggsave(filename="Sim_trend_AMOC_large_2.png", plot=p2, 
       device = "png", width=4, height=3, units="in")


## Distribution change
load("sim_result_distribution_change_AMOC.RData")

npre = 24
npost.vec = seq(24, 120, by=12) 
mu = 10
trend.increase.mu = round(mu * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 
simN = 1000

# The test statistics
Tmax.mat = sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase.mu), FUN=function(j) {
    simResult.mu1[[i]][[j]][[1]][1, ]
  })
})

# The identified changepoint
Cpt.mat = sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase.mu), FUN=function(j) {
    simResult.mu1[[i]][[j]][[1]][2, ]
  })
})

# The delay in effect
delay.mat = sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase.mu), FUN=function(j) {
    simResult.mu1[[i]][[j]][[2]]
  })
})

# The critical values of the test statistics 
nl = round(length(npost.vec) / 2)
criticalVal = rep(getCritical.mu1[1:2], times=c(nl, length(npost.vec)-nl))
criticalMat = rep(criticalVal, times=length(trend.increase.mu))

# True positive rate (TPR): how many time the cpt is correctly detected
# H0: there is no cpt at tau, H1: there is a cpt at tau
Cpt.mat2 = Cpt.mat
Cpt.mat2[Tmax.mat < criticalMat] = -100
# give a large negative value to mean no cpt is detected
Delay.mat = apply(delay.mat, MARGIN=2, FUN=rep, each=length(npost.vec))
gap.mat = abs(Cpt.mat2 - (Delay.mat + npre))
Correct.count = apply(gap.mat <= 12, MARGIN=1, FUN=sum)
TPR.mat = matrix(Correct.count / simN, nrow=length(npost.vec),
                 ncol=length(trend.increase.mu))

# False discovery rate (FDR): among all positive claims, how many get the cpt wrongly
# H0: there is no cpt at tau, H1: there is a cpt at tau
Positive.count = apply(Tmax.mat > criticalMat, MARGIN=1, FUN=sum)
Wrong.count = Positive.count - Correct.count
FDR.mat = matrix(Wrong.count / Positive.count, nrow=length(npost.vec),
                 ncol=length(trend.increase.mu))

## Create the heat maps
trend.percent.mu = c(0.05, seq(0.1, 1, by=0.1)) * 100

plot.df = expand.grid(x=as.factor(npost.vec), 
                      y=as.factor(trend.percent.mu[1:7]))
plot.df$TPR = as.vector(TPR.mat[, 1:7])
plot.df$FDR = as.vector(FDR.mat[, 1:7])

p1 = ggplot(data=plot.df, aes(x=x, y=y, fill=TPR)) + 
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd") + 
  labs(x = "Post intervention length", y = "Size of change (%)") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey95"),
        # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "white", fill=NA),
        legend.position = 'right',
        legend.title = element_text(size=10),
        axis.title = element_text(size=10))

p2 = ggplot(data=plot.df, aes(x=x, y=y, fill=FDR)) + 
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd") + 
  labs(x = "Post intervention length", y = "Size of change (%)",) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey95"),
        # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "white", fill=NA),
        legend.position = 'right',
        legend.title = element_text(size=10),
        axis.title = element_text(size=10))

ggsave(filename="Sim_distribution_AMOC_location_1.png", plot=p1, 
       device = "png", width=4, height=3, units="in")
ggsave(filename="Sim_distribution_AMOC_location_2.png", plot=p2, 
       device = "png", width=4, height=3, units="in")



## (2) Tables for online changepoint detection result

## Detecting change with the forecast error method
load("sim_result_trend_change_Forecast.RData")
# load("sim_result_distribution_change_Forecast.RData")
# Use the same method for the distribution change scenario.

npre = 24
npost.vec = seq(24, 120, by=12) 
level = 10
trend.increase = round(level * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 
simN = 1000

# The estimated changepoint
Cpt.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult[[i]][1, j]
  })
})

# The detection delay
Time.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult[[i]][2, j]   
  })
})

# The delay in effect
Delay.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult[[i]][3, j]  
  })
})

# True positive rate (TPR)
Cpt.true = Delay.mat + npre
correctMat = abs(Cpt.mat - Cpt.true) <= 12
Detect.mat = apply(correctMat, MARGIN=1, FUN=sum, na.rm=TRUE) / simN

# Mean and median delay, standard deviation of delay
delayMat = Time.mat + npre - Cpt.true
Delay.median = apply(delayMat, MARGIN=1, FUN=median, na.rm=TRUE)
Delay.mean = apply(delayMat, MARGIN=1, FUN=function(x) {
  x.fin = x[!is.infinite(x)]   # remove the infinite values (i.e., no detection)
  mean(x.fin, na.rm=TRUE)
})
Delay.std = apply(delayMat, MARGIN=1, FUN=function(x) {
  x.fin = x[!is.infinite(x)]
  sd(x.fin, na.rm=TRUE)
})

# Count of no detection
NoDetection = apply(delayMat, MARGIN=1, FUN=function(x) {
  sum(is.infinite(x))
})

# The summary table
summary.mat = rbind(Detect.mat, Delay.median, Delay.mean, Delay.std, NoDetection)
colnames(summary.mat) = paste0("trend_", trend.increase)
rownames(summary.mat) = c("True_pos", "median_delay", "mean_delay", "std_delay",
                          "no_detection")
summary.mat = round(summary.mat, 2)


## Detecting change with the BOCPD method
load("sim_result_trend_change_BOCPD.RData")
# load("sim_result_distribution_change_BOCPD.RData")
# Use the same method for the distribution change scenario.

npre = 24
npost.vec = seq(24, 120, by=12) 
level = 10
trend.increase = round(level * c(0.05, seq(0.1, 1, by=0.1)) / 120, 4) 
simN = 1000

# The estimated changepoint
Cpt.mat = sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult[[i]][1, j]
  })
})

# The detection time
Time.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult[[i]][2, j]
  })
})

# The delay in effect
Delay.mat =  sapply(1:simN, FUN=function(i) {
  sapply(1:length(trend.increase), FUN=function(j) {
    simResult[[i]][3, j]  
  })
})

# True positive rate (TPR)
Cpt.true = Delay.mat + npre
correctMat = abs(Cpt.mat - Cpt.true) <= 12
Detect.mat = apply(correctMat, MARGIN=1, FUN=sum, na.rm=TRUE) / simN

# Mean and median delay, standard deviation of delay
Cpt.mat2 = Cpt.mat
Cpt.mat2[!is.na(Cpt.mat)] = 1
delayMat = Time.mat * Cpt.mat2 - (Delay.mat + npre)
Delay.median = apply(delayMat, MARGIN=1, FUN=median, na.rm=TRUE)
Delay.mean = apply(delayMat, MARGIN=1, FUN=mean, na.rm=TRUE)
Delay.std = apply(delayMat, MARGIN=1, FUN=sd, na.rm=TRUE)

# Count of no detection
NoDetection = apply(Cpt.mat, MARGIN=1, FUN=function(x) {
  sum(is.na(x))
})

summary.mat = rbind(Detect.mat, Delay.median, Delay.mean, Delay.std, NoDetection)
colnames(summary.mat) = paste0("trend_", trend.increase)
rownames(summary.mat) = c("True_pos", "median_delay", "mean_delay", "std_delay",
                          "no_detection")
summary.mat = round(summary.mat, 2)

