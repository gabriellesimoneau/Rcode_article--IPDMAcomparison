library(lme4)
library(MASS)
library(copula)
setwd("/Users/gabriellesimoneau/Dropbox/article - Master/Biometrical/presentableRcode/RcodeGIT_IPDMAcomparison/R")
source("simulate_example.R")
source("poisson_functions.R")

## As the PHQ9 dataset cannot be publish to this day, we simulated a similar dataset.
## We will use this simulated dataset to demonstrate how to reproduce:
##      - Figure 1: Individual ROC curves
##      - Table 3: Parameters estimates with bivariate approach
##      - Table 5: Parameters estimates with Poisson approach
##      - Table 2: Estimates of pooled sensitivity and specificity from all models
## Results from the simulations. How to reproduce:
##      - Figure 2: results of simulation study
##      - Figure 3: threshold effect simulation study
##      - Table 6: results simulation correlation parameters




#### --------------- Simulate dataset ---------------####
## similar simulation mechanism used in simulation study

example_data <- read.csv("/Users/gabriellesimoneau/Dropbox/article - Master/Biometrical/presentableRcode/RcodeGIT_IPDMAcomparison/Data/simulated_dataset.csv")
Bdata <- data.for.bivariate.example(example_data)
Pdata <- data.for.poisson.example(example_data)




#### --------------- Figure 1 ---------------####
## individual ROC curves from the simulated data ##

D1 <- Bdata[which(Bdata$study==1),]
D1.y <- c(D1$prop[which(D1$d1==1)], 0)
D1.x <- c(D1$prop[which(D1$d1==0)], 0)

D <- rbind(D1.x,D1.y)
for(i in 2:13)
{
  Di <- Bdata[which(Bdata$study==i),]
  D.y <- c(Di$prop[which(Di$d1==1)], 0)
  D.x <- c(Di$prop[which(Di$d1==0)], 0)
  D <- rbind(D, D.x, D.y)
}

par(mar=c(4.1, 4.1, 2.1, 2.1), xpd = TRUE)
design <- c(2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1)
color <- c("darkgrey", "black", "darkgrey", "black", "darkgrey", "black", "darkgrey", "black", "darkgrey", "black", "darkgrey", "black")
plot(D1.x, D1.y, type = "l", lty = 1, xlab = "1-Specificity", ylab = "Sensitivity", cex.axis = 1, cex.lab = 1, lwd = 1.2)
points(D1.x[7], D1.y[7], pch = 16, col = "white")
points(D1.x[14], D1.y[14], pch = 16, col = "white")
points(D1.x[21], D1.y[21], pch = 16, col = "white")
text(D1.x[7], D1.y[7], label = "1", cex = 0.6)
text(D1.x[14], D1.y[14], label = "1", cex = 0.6)
text(D1.x[21], D1.y[21], label = "1", cex = 0.6)
points(0, 0, pch = 18, col = "black", cex = 1)
points(1, 1, pch = 18, col = "black", cex = 1)

lab <- seq(2, 13, by = 1)
count <- 1
for(i in seq(3, 25, by = 2))
{
  lines(D[i,], D[(i+1),], type = "l", lty = design[(i-1)/2], col = color[(i-1)/2])
  points(D[i,7], D[(i+1),7], pch = 16, col = "white")
  points(D[i,14], D[(i+1),14], pch = 16, col = "white")
  points(D[i,21], D[(i+1),21], pch = 16, col = "white")
  text(x=D[i,7], D[(i+1),7], label = lab[count], cex = 0.6)
  text(x=D[i,14], D[(i+1),14], label = lab[count], cex = 0.6)
  text(x=D[i,21], D[(i+1),21], label = lab[count], cex = 0.6)
  count <- count+1
}




#### --------------- Table 3 ---------------####
## parameter estimates with the Bivariate approach ##

# run one model per threshold
cut7 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==7, weights = n)
cut7s <- summary(cut7)

cut8 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==8, weights = n)
cut8s <- summary(cut8)

cut9 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==9, weights = n)
cut9s <- summary(cut9)

cut10 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==10, weights = n)
cut10s <- summary(cut10)

cut11 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==11, weights = n)
cut11s <- summary(cut11)

cut12 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==12, weights = n)
cut12s <- summary(cut12)

cut13 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==13, weights = n)
cut13s <- summary(cut13)

cut14 <- glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family = binomial, nAGQ = 1, data = Bdata, subset = cutoff==14, weights = n)
cut14s <- summary(cut14)

# logit sensitivity and 1-specificity, standarho_dis error logit sensitivity, standard error logit specificity for thresholds 7 to 14
logit_sens_bivariate <- c(fixef(cut7)[2], fixef(cut8)[2], fixef(cut9)[2], fixef(cut10)[2], fixef(cut11)[2], fixef(cut12)[2], fixef(cut13)[2], fixef(cut14)[2])
logit_1_spec_bivariate <- c(fixef(cut7)[1], fixef(cut8)[1], fixef(cut9)[1], fixef(cut10)[1], fixef(cut11)[1], fixef(cut12)[1], fixef(cut13)[1], fixef(cut14)[1])
sd_logit_sens_bivariate <- c(cut7s$coef[2,2], cut8s$coef[2,2], cut9s$coef[2,2], cut10s$coef[2,2], cut11s$coef[2,2], cut12s$coef[2,2], cut13s$coef[2,2], cut14s$coef[2,2])
sd_logit_1_spec_bivariate  <- c(cut7s$coef[1,2], cut8s$coef[1,2], cut9s$coef[1,2], cut10s$coef[1,2], cut11s$coef[1,2], cut12s$coef[1,2], cut13s$coef[1,2], cut14s$coef[1,2])
corr <- c(cut7s$varcor$study[3]/(sqrt(cut7s$varcor$study[1])*sqrt(cut7s$varcor$study[4])), cut8s$varcor$study[3]/(sqrt(cut8s$varcor$study[1])*sqrt(cut8s$varcor$study[4])), cut9s$varcor$study[3]/(sqrt(cut9s$varcor$study[1])*sqrt(cut9s$varcor$study[4])), cut10s$varcor$study[3]/(sqrt(cut10s$varcor$study[1])*sqrt(cut10s$varcor$study[4])), cut11s$varcor$study[3]/(sqrt(cut11s$varcor$study[1])*sqrt(cut11s$varcor$study[4])), cut12s$varcor$study[3]/(sqrt(cut12s$varcor$study[1])*sqrt(cut12s$varcor$study[4])), cut13s$varcor$study[3]/(sqrt(cut13s$varcor$study[1])*sqrt(cut13s$varcor$study[4])), cut14s$varcor$study[3]/(sqrt(cut14s$varcor$study[1])*sqrt(cut14s$varcor$study[4])))




#### --------------- Table 5 ---------------####

Pois_pestimates <- poissonModel(Pdata)

# estimates hazard, frailty variance 
hazards_diseased <- Pois_pestimates$lambda1
hazards_healthy <- Pois_pestimates$lambda0
frailty_variance_1 <- Pois_pestimates$frailty1
frailty_variance_0 <- Pois_pestimates$frailty0

# set up the data to estimate rho_thres and rho_dis: keep only information threshold 7 to 14
Pdata$lam <- c(rep(Pois_pestimates$lambda1,13), rep(Pois_pestimates$lambda0,13))
Pdata$mu <- Pdata$r*Pdata$lam
Pdata2 <- Pdata[which(Pdata$cutoff %in% 7:14),]

# stabilize estimation of rho_thres (but less valid...)
Pdata2$y <- ceiling(Pdata2$y)
Pdata2$r <- ceiling(Pdata2$r)

# estimate the correlation parameter rho_thres
rho_thres <- optimize(rho_thres_estimation, interval = c(0.01,0.99), th1 = 1/Pois_pestimates$frailty1, th0 = 1/Pois_pestimates$frailty0, data = Pdata2, maximum = TRUE)$maximum

# estimate the correlation parameter rho_dis
ab <- optimize(rho_dis_estimation, interval = c(0, sqrt(1/Pois_pestimates$frailty1*1/Pois_pestimates$frailty0)), th1 = 1/Pois_pestimates$frailty1, th0 = 1/Pois_pestimates$frailty0, data = Pdata2, maximum = TRUE)$maximum
rho_dis <- ab*sqrt(Pois_pestimates$frailty1*Pois_pestimates$frailty0)

# estimate the standard errors via parametric bootstrap
# define correlation matrix + margins specification for parametric bootstrap
R <- c(rho_thres^(1:13), rho_dis*rho_thres^(0:13), rho_thres^(1:12), rho_dis*rho_thres, rho_dis*rho_thres^(0:12), rho_thres^(1:11), rho_dis*rho_thres^(2:1), rho_dis*rho_thres^(0:11), rho_thres^(1:10), rho_dis*rho_thres^(3:1), rho_dis*rho_thres^(0:10), rho_thres^(1:9), rho_dis*rho_thres^(4:1), rho_dis*rho_thres^(0:9), rho_thres^(1:8), rho_dis*rho_thres^(5:1), rho_dis*rho_thres^(0:8), rho_thres^(1:7), rho_dis*rho_thres^(6:1), rho_dis*rho_thres^(0:7), rho_thres^(1:6), rho_dis*rho_thres^(7:1), rho_dis*rho_thres^(0:6), rho_thres^(1:5), rho_dis*rho_thres^(8:1), rho_dis*rho_thres^(0:5), rho_thres^(1:4), rho_dis*rho_thres^(9:1), rho_dis*rho_thres^(0:4), rho_thres^(1:3), rho_dis*rho_thres^(10:1), rho_dis*rho_thres^(0:3), rho_thres^(1:2), rho_dis*rho_thres^(11:1), rho_dis*rho_thres^(0:2), rho_thres, rho_dis*rho_thres^(12:1), rho_dis*rho_thres^(0:1), rho_dis*rho_thres^(13:0), rho_thres^(1:13), rho_thres^(1:12), rho_thres^(1:11), rho_thres^(1:10), rho_thres^(1:9), rho_thres^(1:8), rho_thres^(1:7), rho_thres^(1:6), rho_thres^(1:5), rho_thres^(1:4), rho_thres^(1:3), rho_thres^(1:2), rho_thres)
margins <- c(list(list(shape = 1/frailty_variance_1, rate = 1/frailty_variance_1))[rep(1,14)],list(list(shape = 1/frailty_variance_0, rate = 1/frailty_variance_0))[rep(1,14)])
copula_input <- mvdc(normalCopula(R, dim = 28, dispstr="un"),  margins = rep("gamma", 28), paramMargins = margins)

results_SE <- as.data.frame(matrix(NA, nrow = 500, ncol = 36))
colnames(results_SE) <- c("lam1_7", "lam1_8", "lam1_9", "lam1_10", "lam1_11", "lam1_12", "lam1_13", "lam1_14", "lam0_7", "lam0_8", "lam0_9", "lam0_10", "lam0_11", "lam0_12", "lam0_13", "lam0_14", "frailty1", "frailty0", "rho_thres", "rho_dis", "sens7", "sens8", "sens9", "sens10", "sens11", "sens12", "sens13", "sens14", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12", "spec13", "spec14")
set.seed(2072016)
for(i in 1:500)
{
  results_SE[i,] <- parametric_boot(data = Pdata, copula_input = copula_input)
}

SE_hazard_diseased <- apply(results_SE[,1:8], 2, sd)
SE_hazard_healthy <- apply(results_SE[,9:16], 2, sd)
SE_variance_frailty_1 <- sd(results_SE[,17])
SE_variance_frailty_0 <- sd(results_SE[,18])
SE_rho_thres <- sd(results_SE[,19])
SE_rho_dis <- sd(results_SE[,20])




#### --------------- Table 2 ---------------####

## first two columns: results from the Bivariate approach
# extract pooled sensitivity, pooled specificity for thresholds 7 to 14
pooled_specificity_bivariate <- c(expit(-fixef(cut7)[1]), expit(-fixef(cut8)[1]), expit(-fixef(cut9)[1]), expit(-fixef(cut10)[1]), expit(-fixef(cut11)[1]), expit(-fixef(cut12)[1]), expit(-fixef(cut13)[1]), expit(-fixef(cut14)[1]))
pooled_sensitivity_bivariate <- c(expit(fixef(cut7)[2]), expit(fixef(cut8)[2]), expit(fixef(cut9)[2]), expit(fixef(cut10)[2]), expit(fixef(cut11)[2]), expit(fixef(cut12)[2]), expit(fixef(cut13)[2]), expit(fixef(cut14)[2]))

# confidence intervals
CI.sens.L <- c(expit(logit_sens_bivariate - 1.96*sd_logit_sens_bivariate))
CI.sens.U <- c(expit(logit_sens_bivariate + 1.96*sd_logit_sens_bivariate))
CI.spec.L <- 1-c(expit(logit_1_spec_bivariate - 1.96*sd_logit_1_spec_bivariate))
CI.spec.U <- 1-c(expit(logit_1_spec_bivariate + 1.96*sd_logit_1_spec_bivariate))

## last two columns: results from the Poisson approach
# extract pooled sensitivity, pooled specificity for thresholds 7 to 14
pooled_specificity_poisson <- Pois_pestimates$pooled_specificity
pooled_sensitivity_poisson <- Pois_pestimates$pooled_sensitivity
SE_pooled_spec_poisson <- apply(results_SE[,29:36], 2, sd)
SE_pooled_sens_poisson <- apply(results_SE[,21:28], 2, sd)

# confidence intervals
CI.sens_poisson <- pooled_sensitivity_poisson + c(-1,1)*1.96*SE_pooled_sens_poisson 
CI.spec_poisson <- pooled_specificity_poisson + c(-1,1)*1.96*SE_pooled_spec_poisson 




#### --------------- Figure 2 ---------------####

## import csv with results
biv <- read.csv("bivariate_all.csv", header = TRUE)
pois <- read.csv("poisson_all.csv", header = TRUE)
simulation_scenarios <- read.csv("simulation_scenarios.R")
#biv <- merge(biv, simulation_scenarios, by = "scenario")
#pois <- merge(pois, simulation_scenarios, by = "scenario")

## focus only on bias of sensitivity and specificity
true_sens <-  round(seq(0.90, 0.20, length=8), 2)
true_spec <- round(seq(0.20, 0.99, length=8), 2)

## the dataframe figure2 contains all info to construct figure 2
figure2 <- as.data.frame(matrix(NA, nrow = 2*8*2*16))
colnames(figure2) <- c("scenario", "model", "sens", "threshold", "bias", "abs_bias", "mse")
figure2$scenario <- rep(seq(1:16), each = 2*16)
figure2$sens <- rep(c(rep(1,8), rep(0,8)), 2*16)
figure2$threshold <- rep(seq(1:8), 2*2*16)
figure2$model <- c(rep("biv", 8*2*16), rep("pois", 8*2*16))

count <- 1
for(i in 1:16)
{
  sce <- biv[which(biv$scenario == i),]
  bias <- apply(sce[,1:16], 1, function(x) x-c(true_sens, true_spec))
  figure2$bias[count:count+15] <- apply(bias, 2, mean)
  abs_bias <- apply(sce[,1:16], 1, function(x) abs(x-c(true_sens, true_spec)))
  figure2$abs_bias[count:count+15] <- apply(abs_bias, 2, mean)
  mse_1 <- apply(sce[,1:16], 1, function(x) (x-c(true_sens, true_spec))^2)
  figure2$mse[count:count+15] <- apply(mse_1, 2, mean)
  count <- count + 16
}

for(i in 1:16)
{
  sce <- pois[which(pois$scenario == i),]
  bias <- apply(sce[,1:16], 1, function(x) x-c(true_sens, true_spec))
  figure2$bias[count:count+15] <- apply(bias, 2, mean)
  abs_bias <- apply(sce[,1:16], 1, function(x) abs(x-c(true_sens, true_spec)))
  figure2$abs_bias[count:count+15] <- apply(abs_bias, 2, mean)
  mse_1 <- apply(sce[,1:16], 1, function(x) (x-c(true_sens, true_spec))^2)
  figure2$mse[count:count+15] <- apply(mse_1, 2, mean)
  count <- count + 16
}

figure2 <- merge(figure2, simulation_scenarios, by = "scenario")

## define vectors of what will be plotted
figure2_sens <- figure2[which(figure2$sens == 1),]
figure2_spec <- figure2[which(figure2$sens == 1),]

##  Bias Sensitivity
# threshold 1 to 8
bias_sens_thresBiv <- c(mean(figure2_sens$bias[which(figure2_sens$threshold == 1 & figure2_sens$model == "biv")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 2 & figure2_sens$model == "biv")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 3 & figure2_sens$model == "biv")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 4 & figure2_sens$model == "biv")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 5 & figure2_sens$model == "biv")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 6 & figure2_sens$model == "biv")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 7 & figure2_sens$model == "biv")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 8 & figure2_sens$model == "biv")]))
bias_sens_thresPois <- c(mean(figure2_sens$bias[which(figure2_sens$threshold == 1 & figure2_sens$model == "pois")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 2 & figure2_sens$model == "pois")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 3 & figure2_sens$model == "pois")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 4 & figure2_sens$model == "pois")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 5 & figure2_sens$model == "pois")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 6 & figure2_sens$model == "pois")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 7 & figure2_sens$model == "pois")]), mean(figure2_sens$bias[which(figure2_sens$threshold == 8 & figure2_sens$model == "pois")]))
# variance RE low and high
bias_sens_varianceBiv <- c(mean(figure2_sens$bias[which(figure2_sens$var_sens==0.1 & figure2_sens$model == "biv")]),mean(figure2_sens$bias[which(figure2_sens$var_sens==0.25 & figure2_sens$model == "biv")]))
bias_sens_variancePois <- c(mean(figure2_sens$bias[which(figure2_sens$var_sens==0.1 & figure2_sens$model == "pois")]),mean(figure2_sens$bias[which(figure2_sens$var_sens==0.25 & figure2_sens$model == "pois")]))
# correlation rho dis low and high
bias_sens_disBiv <- c(mean(figure2_sens$bias[which(figure2_sens$rho_dis==0.25 & figure2_sens$model == "biv")]),mean(figure2_sens$bias[which(figure2_sens$rho_dis==0.75 & figure2_sens$model == "biv")]))
bias_sens_disPois <- c(mean(figure2_sens$bias[which(figure2_sens$rho_dis==0.25 & figure2_sens$model == "pois")]),mean(figure2_sens$bias[which(figure2_sens$rho_dis==0.75 & figure2_sens$model == "pois")]))
# correlation rho thres low and high
bias_sens_rthBiv <- c(mean(figure2_sens$bias[which(figure2_sens$rho_thres==0.25 & figure2_sens$model == "biv")]),mean(figure2_sens$bias[which(figure2_sens$rho_thres==0.75 & figure2_sens$model == "biv")]))
bias_sens_rthPois <- c(mean(figure2_sens$bias[which(figure2_sens$rho_thres==0.25 & figure2_sens$model == "pois")]),mean(figure2_sens$bias[which(figure2_sens$rho_thres==0.75 & figure2_sens$model == "pois")]))
# sample size diseased small large
bias_sens_ndBiv <- c(mean(figure2_sens$bias[which(figure2_sens$ND==50 & figure2_sens$model == "biv")]),mean(figure2_sens$bias[which(figure2_sens$ND==100 & figure2_sens$model == "biv")]))
bias_sens_ndPois <- c(mean(figure2_sens$bias[which(figure2_sens$ND==50 & figure2_sens$model == "pois")]),mean(figure2_sens$bias[which(figure2_sens$ND==100 & figure2_sens$model == "pois")]))

## Bias Specificity
bias_spec_thresBiv <- c(mean(figure2_spec$bias[which(figure2_spec$threshold == 1 & figure2_spec$model == "biv")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 2 & figure2_spec$model == "biv")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 3 & figure2_spec$model == "biv")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 4 & figure2_spec$model == "biv")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 5 & figure2_spec$model == "biv")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 6 & figure2_spec$model == "biv")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 7 & figure2_spec$model == "biv")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 8 & figure2_spec$model == "biv")]))
bias_spec_thresPois <- c(mean(figure2_spec$bias[which(figure2_spec$threshold == 1 & figure2_spec$model == "pois")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 2 & figure2_spec$model == "pois")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 3 & figure2_spec$model == "pois")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 4 & figure2_spec$model == "pois")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 5 & figure2_spec$model == "pois")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 6 & figure2_spec$model == "pois")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 7 & figure2_spec$model == "pois")]), mean(figure2_spec$bias[which(figure2_spec$threshold == 8 & figure2_spec$model == "pois")]))

bias_spec_varianceBiv <- c(mean(figure2_spec$bias[which(figure2_spec$var_spec==0.1 & figure2_spec$model == "biv")]),mean(figure2_spec$bias[which(figure2_spec$var_spec==0.25 & figure2_spec$model == "biv")]))
bias_spec_variancePois <- c(mean(figure2_spec$bias[which(figure2_spec$var_spec==0.1 & figure2_spec$model == "pois")]),mean(figure2_spec$bias[which(figure2_spec$var_spec==0.25 & figure2_spec$model == "pois")]))

bias_spec_disBiv <- c(mean(figure2_spec$bias[which(figure2_spec$rho_dis==0.25 & figure2_spec$model == "biv")]),mean(figure2_spec$bias[which(figure2_spec$rho_dis==0.75 & figure2_spec$model == "biv")]))
bias_spec_disPois <- c(mean(figure2_spec$bias[which(figure2_spec$rho_dis==0.25 & figure2_spec$model == "pois")]),mean(figure2_spec$bias[which(figure2_spec$rho_dis==0.75 & figure2_spec$model == "pois")]))

bias_spec_rthBiv <- c(mean(figure2_spec$bias[which(figure2_spec$rho_thres==0.25 & figure2_spec$model == "biv")]),mean(figure2_spec$bias[which(figure2_spec$rho_thres==0.75 & figure2_spec$model == "biv")]))
bias_spec_rthPois <- c(mean(figure2_spec$bias[which(figure2_spec$rho_thres==0.25 & figure2_spec$model == "pois")]),mean(figure2_spec$bias[which(figure2_spec$rho_thres==0.75 & figure2_spec$model == "pois")]))

bias_spec_ndBiv <- c(mean(figure2_spec$bias[which(figure2_spec$ND==50 & figure2_spec$model == "biv")]),mean(figure2_spec$bias[which(figure2_spec$ND==100 & figure2_spec$model == "biv")]))
bias_spec_ndPois <- c(mean(figure2_spec$bias[which(figure2_spec$ND==50 & figure2_spec$model == "pois")]),mean(figure2_spec$bias[which(figure2_spec$ND==100 & figure2_spec$model == "pois")]))

## absolute bias sensitivity
abs_bias_sens_thresBiv <- c(mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 1 & figure2_sens$model == "biv")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 2 & figure2_sens$model == "biv")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 3 & figure2_sens$model == "biv")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 4 & figure2_sens$model == "biv")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 5 & figure2_sens$model == "biv")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 6 & figure2_sens$model == "biv")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 7 & figure2_sens$model == "biv")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 8 & figure2_sens$model == "biv")]))
abs_bias_sens_thresPois <- c(mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 1 & figure2_sens$model == "pois")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 2 & figure2_sens$model == "pois")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 3 & figure2_sens$model == "pois")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 4 & figure2_sens$model == "pois")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 5 & figure2_sens$model == "pois")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 6 & figure2_sens$model == "pois")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 7 & figure2_sens$model == "pois")]), mean(figure2_sens$abs_bias[which(figure2_sens$threshold == 8 & figure2_sens$model == "pois")]))

abs_bias_sens_varianceBiv <- c(mean(figure2_sens$abs_bias[which(figure2_sens$var_sens==0.1 & figure2_sens$model == "biv")]),mean(figure2_sens$abs_bias[which(figure2_sens$var_sens==0.25 & figure2_sens$model == "biv")]))
abs_bias_sens_variancePois <- c(mean(figure2_sens$abs_bias[which(figure2_sens$var_sens==0.1 & figure2_sens$model == "pois")]),mean(figure2_sens$abs_bias[which(figure2_sens$var_sens==0.25 & figure2_sens$model == "pois")]))

abs_bias_sens_disBiv <- c(mean(figure2_sens$abs_bias[which(figure2_sens$rho_dis==0.25 & figure2_sens$model == "biv")]),mean(figure2_sens$abs_bias[which(figure2_sens$rho_dis==0.75 & figure2_sens$model == "biv")]))
abs_bias_sens_disPois <- c(mean(figure2_sens$abs_bias[which(figure2_sens$rho_dis==0.25 & figure2_sens$model == "pois")]),mean(figure2_sens$abs_bias[which(figure2_sens$rho_dis==0.75 & figure2_sens$model == "pois")]))

abs_bias_sens_rthBiv <- c(mean(figure2_sens$abs_bias[which(figure2_sens$rho_thres==0.25 & figure2_sens$model == "biv")]),mean(figure2_sens$abs_bias[which(figure2_sens$rho_thres==0.75 & figure2_sens$model == "biv")]))
abs_bias_sens_rthPois <- c(mean(figure2_sens$abs_bias[which(figure2_sens$rho_thres==0.25 & figure2_sens$model == "pois")]),mean(figure2_sens$abs_bias[which(figure2_sens$rho_thres==0.75 & figure2_sens$model == "pois")]))

abs_bias_sens_ndBiv <- c(mean(figure2_sens$abs_bias[which(figure2_sens$ND==50 & figure2_sens$model == "biv")]),mean(figure2_sens$abs_bias[which(figure2_sens$ND==100 & figure2_sens$model == "biv")]))
abs_bias_sens_ndPois <- c(mean(figure2_sens$abs_bias[which(figure2_sens$ND==50 & figure2_sens$model == "pois")]),mean(figure2_sens$abs_bias[which(figure2_sens$ND==100 & figure2_sens$model == "pois")]))

## absolute bias specificity
abs_bias_spec_thresBiv <- c(mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 1 & figure2_spec$model == "biv")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 2 & figure2_spec$model == "biv")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 3 & figure2_spec$model == "biv")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 4 & figure2_spec$model == "biv")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 5 & figure2_spec$model == "biv")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 6 & figure2_spec$model == "biv")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 7 & figure2_spec$model == "biv")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 8 & figure2_spec$model == "biv")]))
abs_bias_spec_thresPois <- c(mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 1 & figure2_spec$model == "pois")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 2 & figure2_spec$model == "pois")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 3 & figure2_spec$model == "pois")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 4 & figure2_spec$model == "pois")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 5 & figure2_spec$model == "pois")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 6 & figure2_spec$model == "pois")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 7 & figure2_spec$model == "pois")]), mean(figure2_spec$abs_bias[which(figure2_spec$threshold == 8 & figure2_spec$model == "pois")]))

abs_bias_spec_varianceBiv <- c(mean(figure2_spec$abs_bias[which(figure2_spec$var_spec==0.1 & figure2_spec$model == "biv")]),mean(figure2_spec$abs_bias[which(figure2_spec$var_spec==0.25 & figure2_spec$model == "biv")]))
abs_bias_spec_variancePois <- c(mean(figure2_spec$abs_bias[which(figure2_spec$var_spec==0.1 & figure2_spec$model == "pois")]),mean(figure2_spec$abs_bias[which(figure2_spec$var_spec==0.25 & figure2_spec$model == "pois")]))

abs_bias_spec_disBiv <- c(mean(figure2_spec$abs_bias[which(figure2_spec$rho_dis==0.25 & figure2_spec$model == "biv")]),mean(figure2_spec$abs_bias[which(figure2_spec$rho_dis==0.75 & figure2_spec$model == "biv")]))
abs_bias_spec_disPois <- c(mean(figure2_spec$abs_bias[which(figure2_spec$rho_dis==0.25 & figure2_spec$model == "pois")]),mean(figure2_spec$abs_bias[which(figure2_spec$rho_dis==0.75 & figure2_spec$model == "pois")]))

abs_bias_spec_rthBiv <- c(mean(figure2_spec$abs_bias[which(figure2_spec$rho_thres==0.25 & figure2_spec$model == "biv")]),mean(figure2_spec$abs_bias[which(figure2_spec$rho_thres==0.75 & figure2_spec$model == "biv")]))
abs_bias_spec_rthPois <- c(mean(figure2_spec$abs_bias[which(figure2_spec$rho_thres==0.25 & figure2_spec$model == "pois")]),mean(figure2_spec$abs_bias[which(figure2_spec$rho_thres==0.75 & figure2_spec$model == "pois")]))

abs_bias_spec_ndBiv <- c(mean(figure2_spec$abs_bias[which(figure2_spec$ND==50 & figure2_spec$model == "biv")]),mean(figure2_spec$abs_bias[which(figure2_spec$ND==100 & figure2_spec$model == "biv")]))
abs_bias_spec_ndPois <- c(mean(figure2_spec$abs_bias[which(figure2_spec$ND==50 & figure2_spec$model == "pois")]),mean(figure2_spec$abs_bias[which(figure2_spec$ND==100 & figure2_spec$model == "pois")]))

## MSE sensitivity
mse_sens_thresBiv <- c(mean(figure2_sens$mse[which(figure2_sens$threshold == 1 & figure2_sens$model == "biv")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 2 & figure2_sens$model == "biv")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 3 & figure2_sens$model == "biv")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 4 & figure2_sens$model == "biv")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 5 & figure2_sens$model == "biv")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 6 & figure2_sens$model == "biv")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 7 & figure2_sens$model == "biv")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 8 & figure2_sens$model == "biv")]))
mse_sens_thresPois <- c(mean(figure2_sens$mse[which(figure2_sens$threshold == 1 & figure2_sens$model == "pois")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 2 & figure2_sens$model == "pois")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 3 & figure2_sens$model == "pois")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 4 & figure2_sens$model == "pois")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 5 & figure2_sens$model == "pois")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 6 & figure2_sens$model == "pois")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 7 & figure2_sens$model == "pois")]), mean(figure2_sens$mse[which(figure2_sens$threshold == 8 & figure2_sens$model == "pois")]))

mse_sens_varianceBiv <- c(mean(figure2_sens$mse[which(figure2_sens$var_sens==0.1 & figure2_sens$model == "biv")]),mean(figure2_sens$mse[which(figure2_sens$var_sens==0.25 & figure2_sens$model == "biv")]))
mse_sens_variancePois <- c(mean(figure2_sens$mse[which(figure2_sens$var_sens==0.1 & figure2_sens$model == "pois")]),mean(figure2_sens$mse[which(figure2_sens$var_sens==0.25 & figure2_sens$model == "pois")]))

mse_sens_disBiv <- c(mean(figure2_sens$mse[which(figure2_sens$rho_dis==0.25 & figure2_sens$model == "biv")]),mean(figure2_sens$mse[which(figure2_sens$rho_dis==0.75 & figure2_sens$model == "biv")]))
mse_sens_disPois <- c(mean(figure2_sens$mse[which(figure2_sens$rho_dis==0.25 & figure2_sens$model == "pois")]),mean(figure2_sens$mse[which(figure2_sens$rho_dis==0.75 & figure2_sens$model == "pois")]))

mse_sens_rthBiv <- c(mean(figure2_sens$mse[which(figure2_sens$rho_thres==0.25 & figure2_sens$model == "biv")]),mean(figure2_sens$mse[which(figure2_sens$rho_thres==0.75 & figure2_sens$model == "biv")]))
mse_sens_rthPois <- c(mean(figure2_sens$mse[which(figure2_sens$rho_thres==0.25 & figure2_sens$model == "pois")]),mean(figure2_sens$mse[which(figure2_sens$rho_thres==0.75 & figure2_sens$model == "pois")]))

mse_sens_ndBiv <- c(mean(figure2_sens$mse[which(figure2_sens$ND==50 & figure2_sens$model == "biv")]),mean(figure2_sens$mse[which(figure2_sens$ND==100 & figure2_sens$model == "biv")]))
mse_sens_ndPois <- c(mean(figure2_sens$mse[which(figure2_sens$ND==50 & figure2_sens$model == "pois")]),mean(figure2_sens$mse[which(figure2_sens$ND==100 & figure2_sens$model == "pois")]))

## MSE specificity
mse_spec_thresBiv <- c(mean(figure2_spec$mse[which(figure2_spec$threshold == 1 & figure2_spec$model == "biv")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 2 & figure2_spec$model == "biv")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 3 & figure2_spec$model == "biv")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 4 & figure2_spec$model == "biv")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 5 & figure2_spec$model == "biv")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 6 & figure2_spec$model == "biv")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 7 & figure2_spec$model == "biv")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 8 & figure2_spec$model == "biv")]))
mse_spec_thresPois <- c(mean(figure2_spec$mse[which(figure2_spec$threshold == 1 & figure2_spec$model == "pois")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 2 & figure2_spec$model == "pois")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 3 & figure2_spec$model == "pois")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 4 & figure2_spec$model == "pois")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 5 & figure2_spec$model == "pois")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 6 & figure2_spec$model == "pois")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 7 & figure2_spec$model == "pois")]), mean(figure2_spec$mse[which(figure2_spec$threshold == 8 & figure2_spec$model == "pois")]))

mse_spec_varianceBiv <- c(mean(figure2_spec$mse[which(figure2_spec$var_spec==0.1 & figure2_spec$model == "biv")]),mean(figure2_spec$mse[which(figure2_spec$var_spec==0.25 & figure2_spec$model == "biv")]))
mse_spec_variancePois <- c(mean(figure2_spec$mse[which(figure2_spec$var_spec==0.1 & figure2_spec$model == "pois")]),mean(figure2_spec$mse[which(figure2_spec$var_spec==0.25 & figure2_spec$model == "pois")]))

mse_spec_disBiv <- c(mean(figure2_spec$mse[which(figure2_spec$rho_dis==0.25 & figure2_spec$model == "biv")]),mean(figure2_spec$mse[which(figure2_spec$rho_dis==0.75 & figure2_spec$model == "biv")]))
mse_spec_disPois <- c(mean(figure2_spec$mse[which(figure2_spec$rho_dis==0.25 & figure2_spec$model == "pois")]),mean(figure2_spec$mse[which(figure2_spec$rho_dis==0.75 & figure2_spec$model == "pois")]))

mse_spec_rthBiv <- c(mean(figure2_spec$mse[which(figure2_spec$rho_thres==0.25 & figure2_spec$model == "biv")]),mean(figure2_spec$mse[which(figure2_spec$rho_thres==0.75 & figure2_spec$model == "biv")]))
mse_spec_rthPois <- c(mean(figure2_spec$mse[which(figure2_spec$rho_thres==0.25 & figure2_spec$model == "pois")]),mean(figure2_spec$mse[which(figure2_spec$rho_thres==0.75 & figure2_spec$model == "pois")]))

mse_spec_ndBiv <- c(mean(figure2_spec$mse[which(figure2_spec$ND==50 & figure2_spec$model == "biv")]),mean(figure2_spec$mse[which(figure2_spec$ND==100 & figure2_spec$model == "biv")]))
mse_spec_ndPois <- c(mean(figure2_spec$mse[which(figure2_spec$ND==50 & figure2_spec$model == "pois")]),mean(figure2_spec$mse[which(figure2_spec$ND==100 & figure2_spec$model == "pois")]))

#### figure 2 - start here
okgrey=c("gray40")
par(mfrow=c(3,2),mar=c(4.1,4.1,1.1,1.1))

# bias sens
plot(x = rep(1,2), y = c(min(bias_sens_thresBiv), max(bias_sens_thresBiv)), xlim = c(0,15), ylim = c(-0.015,0.02), type = "l", xlab = "Factors", ylab = "Mean bias of sensitivity", xaxt = "n", col = "black", cex.lab=1, cex.axis=0.8)
abline(0,0)
axis(1, at = c(1.5,4.5,7.5,10.5,13.5), labels = c("Threshold","Heterogeneity","rho(dis)","rho(thres)","n1"), cex.axis = 1)
lines(x = rep(2,2), y = c(min(bias_sens_thresPois), max(bias_sens_thresPois)), col = okgrey)
points(x = rep(1,8), y = bias_sens_thresBiv, pch="-", col="black")
text(x = rep(0.75,8), y = bias_sens_thresBiv, labels = c("1","2","3","4","5","6","7","8"), cex = 1, col = "black")
points(x = rep(2,8), y = bias_sens_thresPois, pch = "-", col = okgrey)
text(x = rep(1.75,8), y = bias_sens_thresPois, labels = c("1","2","3","4","5","6","7","8"), cex = 1, col = okgrey)

lines(x = rep(4,2), y = c(min(bias_sens_varianceBiv), max(bias_sens_varianceBiv)), col = "black")
lines(x = rep(5,2), y = c(min(bias_sens_variancePois), max(bias_sens_variancePois)), col = okgrey)
points(x = rep(4,2), y = bias_sens_varianceBiv, pch = "-", col = "black")
text(x = rep(4,2), y = c(bias_sens_varianceBiv[1]-0.00117, bias_sens_varianceBiv[2]+0.0013), labels = c("low","high"), cex = 1, col = "black")
points(x = rep(5,2), y = bias_sens_variancePois, pch = "-", col = okgrey)
text(x = rep(5,2), y = c(bias_sens_variancePois[1]-0.00117, bias_sens_variancePois[2]+0.0013), labels = c("low","high"), cex = 1, col = okgrey)

lines(x = rep(7,2), y = c(min(bias_sens_disBiv), max(bias_sens_disBiv)), col = "black")
lines(x = rep(8,2), y = c(min(bias_sens_disPois), max(bias_sens_disPois)), col = okgrey)
points(x = rep(7,2), y = bias_sens_disBiv, pch = "-", col = "black")
text(x = rep(7,2), y = c(bias_sens_disBiv[1]+0.0013, bias_sens_disBiv[2]-0.00117), labels = c("0.25","0.75"), cex = 1, col = "black")
points(x = rep(8,2), y = bias_sens_disPois, pch = "-", col = okgrey)
text(x = rep(8,2), y = c(bias_sens_disPois[1]+0.0013, bias_sens_disPois[2]-0.00117), labels = c("0.25","0.75"), cex = 1, col = okgrey)

lines(x = rep(10,2), y = bias_sens_rthBiv, col="black")
lines(x = rep(11,2), y = bias_sens_rthPois, col = okgrey)
points(x = rep(10,2), y = bias_sens_rthBiv, pch = "-", col = "black")
text(x = rep(10,2), y = c(bias_sens_rthBiv[1]+0.0013, bias_sens_rthBiv[2]-0.00117), labels = c("0.25","0.75"), cex = 1, col = "black")
points(x = rep(11,2), y = bias_sens_rthPois, pch = "-", col = okgrey)
text(x = rep(11,2), y = c(bias_sens_rthPois[1]-0.00117, bias_sens_rthPois[2]+0.0013), labels = c("0.25","0.75"), cex = 1, col = okgrey)

lines(x = rep(13,2), y = bias_sens_ndBiv, col = "black")
lines(x = rep(14,2), y = bias_sens_ndPois, col = okgrey)
points(x = rep(13,2), y = bias_sens_ndBiv, pch = "-", col = "black")
text(x = rep(13,2), y = c(bias_sens_ndBiv[1]+0.0013, bias_sens_ndBiv[2]-0.00117), labels = c("50","100"), cex = 1, col = "black")
points(x = rep(14,2), y = bias_sens_ndPois, pch = "-", col = okgrey)
text(x = rep(14,2), y = c(bias_sens_ndPois[1]-0.00117, bias_sens_ndP[2]+0.0013), labels = c("50","100"), cex = 1, col = okgrey)



#### figure 2 - end here









