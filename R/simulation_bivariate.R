## simulation bivariate
source("functions_general.R")
library(lme4)
set.seed(21)

## scenario parameters
sc <- read.csv("simulation_scenarios.csv")

## allocate space to save scenario results
results_bivariate <- vector(mode = "list", length = 16)

for(i in 1:16) # loop over scenarios
{
  results_bivariate[[i]] <- matrix(NA, nrow = 1000, ncol = 40)
  for(j in 1:1000) # number of simulated scenarios per scenario
  {
    ## simulate data + set up data for bivariate model
    simulD <- simulation(rhodis = sc$rho_dis[i], rhothres = sc$rho_thres[i], sigma_sens = sc$var_sens[i], sigma_spec = sc$var_spec[i], nd = sc$ND[i], nnd = 400, sdd = sc$SD.ND.[i])
    Bdata <- data.for.bivariate.simulation(x0 = simulD$x0, x1 = simulD$x1, nthres = 8)
    
    ## apply Poisson to all cutoff and extract sensitivty specificity var frailty rhothres rhodis
    sensitivity <- rep(NA, 8)
    specificity <- rep(NA, 8)
    var.sens <- rep(NA, 8)
    var.spec <- rep(NA, 8)
    corr <- rep(NA, 8)
    
    for(k in 1:8)
    {
      mod <- try(glmer(prop ~ -1 + d0 + d1 + (-1 + d0 + d1|study), family=binomial, nAGQ=1, data=Bdata, subset = cutoff == k, weight=ni))
      sensitivity[k] <- expit(fixef(mod)[2])
      specificity[k] <- expit(-fixef(mod)[1])
      var.sens[k] <- summary(mod)$varcor$study[4]
      var.spec[k] <- summary(mod)$varcor$study[1]
      corr[k] <- summary(mod)$varcor$study[3]/(sqrt(var.sens[k])*sqrt(var.spec[k]))
    }
    
    ## save estimates
    results_bivariate[[i]][j,] <- c(sensitivity, specificity, var.sens, var.spec, corr)
    print(j)
  }
  results_bivariate[[i]] <- as.data.frame(results_bivariate[[i]])
  colnames(results_bivariate[[i]]) <- c("sens1", "sens2", "sens3", "sens4", "sens5", "sens6", "sens7", "sens8", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "varsens1", "varsens2", "varsens3", "varsens4", "varsens5", "varsens6", "varsens7", "varsens8", "varspec1", "varspec2", "varspec3", "varspec4", "varspec5", "varspec6", "varspec7", "varspec8", "corr1", "corr2", "corr3", "corr4", "corr5", "corr6", "corr7", "corr8")
  name <- paste("bivariate_scenario", paste(i), ".csv", sep = "")
  write.csv(results_bivariate[[i]], file = name)
}

## write csv with all results
all <- as.data.frame(rbind(results_bivariate[[1]], results_bivariate[[2]], results_bivariate[[3]], results_bivariate[[4]], results_bivariate[[5]], results_bivariate[[6]], results_bivariate[[7]], results_bivariate[[8]], results_bivariate[[9]], results_bivariate[[10]], results_bivariate[[11]], results_bivariate[[12]], results_bivariate[[13]], results_bivariate[[14]], results_bivariate[[15]], results_bivariate[[16]]))
colnames(all) <- c("sens1", "sens2", "sens3", "sens4", "sens5", "sens6", "sens7", "sens8", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "varsens1", "varsens2", "varsens3", "varsens4", "varsens5", "varsens6", "varsens7", "varsens8", "varspec1", "varspec2", "varspec3", "varspec4", "varspec5", "varspec6", "varspec7", "varspec8", "corr1", "corr2", "corr3", "corr4", "corr5", "corr6", "corr7", "corr8")
all$scenario <- rep(seq(1:16), each = 1000)
write.csv(all, file = "simulation_results_bivariate.csv")









