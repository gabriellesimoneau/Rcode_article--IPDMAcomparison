## simulation Poisson
source("functions_general.R")
source("functions_poisson.R")
library(MASS)


## scenario parameters
sc <- read.csv("simulation_scenarios.csv")

## allocate space to save scenario results
results_poisson <- vector(mode = "list", length = 16)

## Run each for loop separately to speed to computational time
## seed for each scenario:
##    1,2,3,4 - set.seed(345)
##    5,6,7,8 - set.seed(5714)
##    9,10,11,12 - set.seed(514)
##    13,14,15,16 - set.seed(170)

set.seed(345)
for(i in 1:4) # loop over scenarios
{
  results_poisson[[i]] <- matrix(NA, nrow = 1000, ncol = 20)
  for(j in 1:1000) # number of simulated scenarios per scenario
  {
    ## simulate data + set up data for bivariate model
    simulD <- simulation(rhodis = sc$rho_dis[i], rhothres = sc$rho_thres[i], sigma_sens = sc$var_sens[i], sigma_spec = sc$var_spec[i], nd = sc$ND[i], nnd = 400, sdd = sc$SD.ND.[i])
    Pdata <- data.for.poisson.simulation(x0 = simulD$x0, x1 = simulD$x1, nthres = 8)
    
    ## apply Poisson model
    cutD <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==1, init.theta =10))
    cutND <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==0, init.theta = 10))
    if (class(cutD) != "try-error" & class(cutND) != "try-error") 
    {
      lam1 <- as.numeric(c(exp(cutD$coeff[1]), exp(cutD$coeff[2]), exp(cutD$coeff[3]), exp(cutD$coeff[4]), exp(cutD$coeff[5]), exp(cutD$coeff[6]), exp(cutD$coeff[7]), exp(cutD$coeff[8])))
      lam0 <- as.numeric(c(exp(cutND$coeff[1]), exp(cutND$coeff[2]), exp(cutND$coeff[3]), exp(cutND$coeff[4]), exp(cutND$coeff[5]), exp(cutND$coeff[6]), exp(cutND$coeff[7]), exp(cutND$coeff[8])))
      ksi <- c(1/cutD$theta, 1/cutND$theta) 
      theta <- c(cutD$theta, cutND$theta)
    
      sensitivity <- cumprod(1-lam1)
      specificity <- 1-cumprod(1-lam0)
    
      Pdata$lam <- c(rep(lam1,13), rep(lam0,13))
      Pdata$mu <- Pdata$r*Pdata$lam
    
      rho_thres <- optimize(rho_thres_estimation, interval = c(0.001,0.99), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      ab <- optimize(rho_dis_estimation, interval = c(0, sqrt(theta[1]*theta[2])), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      rho_dis <- ab*sqrt(ksi[1]*ksi[2])
    }
    ## save estimates
    results_poisson[[i]][j,] <- c(sensitivity, specificity, ksi, rho_thres, rho_dis)
    print(j)
  }
  results_i <- as.data.frame(results_poisson[[i]])
  colnames(results_i) <- c("sens1", "sens2", "sens3", "sens4", "sens5", "sens6", "sens7", "sens8", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "frail1", "frail0", "rhothres", "rhodis")
  name <- paste("poisson_scenario", paste(i), ".csv", sep = "")
  write.csv(results_i, file = name)
}

set.seed(5714)
for(i in 5:8) # loop over scenarios
{
  results_poisson[[i-4]] <- matrix(NA, nrow = 1000, ncol = 20)
  for(j in 1:1000) # number of simulated scenarios per scenario
  {
    ## simulate data + set up data for bivariate model
    simulD <- simulation(rhodis = sc$rho_dis[i], rhothres = sc$rho_thres[i], sigma_sens = sc$var_sens[i], sigma_spec = sc$var_spec[i], nd = sc$ND[i], nnd = 400, sdd = sc$SD.ND.[i])
    Pdata <- data.for.poisson.simulation(x0 = simulD$x0, x1 = simulD$x1, nthres = 8)
    
    ## apply Poisson model
    cutD <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==1, init.theta =10))
    cutND <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==0, init.theta = 10))
    if (class(cutD) != "try-error" & class(cutND) != "try-error") 
    {
      lam1 <- as.numeric(c(exp(cutD$coeff[1]), exp(cutD$coeff[2]), exp(cutD$coeff[3]), exp(cutD$coeff[4]), exp(cutD$coeff[5]), exp(cutD$coeff[6]), exp(cutD$coeff[7]), exp(cutD$coeff[8])))
      lam0 <- as.numeric(c(exp(cutND$coeff[1]), exp(cutND$coeff[2]), exp(cutND$coeff[3]), exp(cutND$coeff[4]), exp(cutND$coeff[5]), exp(cutND$coeff[6]), exp(cutND$coeff[7]), exp(cutND$coeff[8])))
      ksi <- c(1/cutD$theta, 1/cutND$theta) 
      theta <- c(cutD$theta, cutND$theta)
      
      sensitivity <- cumprod(1-lam1)
      specificity <- 1-cumprod(1-lam0)
      
      Pdata$lam <- c(rep(lam1,13), rep(lam0,13))
      Pdata$mu <- Pdata$r*Pdata$lam
      
      #rho_thres <- optimize(rho_thres_estimation, interval = c(0.001,0.99), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      ab <- optimize(rho_dis_estimation, interval = c(0, sqrt(theta[1]*theta[2])), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      rho_dis <- ab*sqrt(ksi[1]*ksi[2])
    }
    ## save estimates
    results_poisson[[i-4]][j,] <- c(sensitivity, specificity, ksi, 0, rho_dis)
    print(j)
  }
  results_i <- as.data.frame(results_poisson[[i-4]])
  colnames(results_i) <- c("sens1", "sens2", "sens3", "sens4", "sens5", "sens6", "sens7", "sens8", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "frail1", "frail0", "rhothres", "rhodis")
  name <- paste("poisson_scenario", paste(i-4), ".csv", sep = "")
  write.csv(results_i, file = name)
}

set.seed(514)
for(i in 9:12) # loop over scenarios
{
  results_poisson[[i-8]] <- matrix(NA, nrow = 1000, ncol = 20)
  for(j in 1:1000) # number of simulated scenarios per scenario
  {
    ## simulate data + set up data for bivariate model
    simulD <- simulation(rhodis = sc$rho_dis[i], rhothres = sc$rho_thres[i], sigma_sens = sc$var_sens[i], sigma_spec = sc$var_spec[i], nd = sc$ND[i], nnd = 400, sdd = sc$SD.ND.[i])
    Pdata <- data.for.poisson.simulation(x0 = simulD$x0, x1 = simulD$x1, nthres = 8)
    
    ## apply Poisson model
    cutD <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==1, init.theta =10))
    cutND <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==0, init.theta = 10))
    if (class(cutD) != "try-error" & class(cutND) != "try-error") 
    {
      lam1 <- as.numeric(c(exp(cutD$coeff[1]), exp(cutD$coeff[2]), exp(cutD$coeff[3]), exp(cutD$coeff[4]), exp(cutD$coeff[5]), exp(cutD$coeff[6]), exp(cutD$coeff[7]), exp(cutD$coeff[8])))
      lam0 <- as.numeric(c(exp(cutND$coeff[1]), exp(cutND$coeff[2]), exp(cutND$coeff[3]), exp(cutND$coeff[4]), exp(cutND$coeff[5]), exp(cutND$coeff[6]), exp(cutND$coeff[7]), exp(cutND$coeff[8])))
      ksi <- c(1/cutD$theta, 1/cutND$theta) 
      theta <- c(cutD$theta, cutND$theta)
      
      sensitivity <- cumprod(1-lam1)
      specificity <- 1-cumprod(1-lam0)
      
      Pdata$lam <- c(rep(lam1,13), rep(lam0,13))
      Pdata$mu <- Pdata$r*Pdata$lam
      
      #rho_thres <- optimize(rho_thres_estimation, interval = c(0.001,0.99), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      ab <- optimize(rho_dis_estimation, interval = c(0, sqrt(theta[1]*theta[2])), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      rho_dis <- ab*sqrt(ksi[1]*ksi[2])
    }
    ## save estimates
    results_poisson[[i-8]][j,] <- c(sensitivity, specificity, ksi, 0, rho_dis)
    print(j)
  }
  results_i <- as.data.frame(results_poisson[[i-8]])
  colnames(results_i) <- c("sens1", "sens2", "sens3", "sens4", "sens5", "sens6", "sens7", "sens8", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "frail1", "frail0", "rhothres", "rhodis")
  name <- paste("poisson_scenario", paste(i-8), ".csv", sep = "")
  write.csv(results_i, file = name)
}

set.seed(170)
for(i in 13:16) # loop over scenarios
{
  results_poisson[[i-12]] <- matrix(NA, nrow = 1000, ncol = 20)
  for(j in 1:1000) # number of simulated scenarios per scenario
  {
    ## simulate data + set up data for bivariate model
    simulD <- simulation(rhodis = sc$rho_dis[i], rhothres = sc$rho_thres[i], sigma_sens = sc$var_sens[i], sigma_spec = sc$var_spec[i], nd = sc$ND[i], nnd = 400, sdd = sc$SD.ND.[i])
    Pdata <- data.for.poisson.simulation(x0 = simulD$x0, x1 = simulD$x1, nthres = 8)
    
    ## apply Poisson model
    cutD <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==1, init.theta =10))
    cutND <- try(glm.nb(y ~ -1 + offset(log(r)) + as.factor(cutoff), weights = n, data = Pdata, subset = dis==0, init.theta = 10))
    if (class(cutD) != "try-error" & class(cutND) != "try-error") 
    {
      lam1 <- as.numeric(c(exp(cutD$coeff[1]), exp(cutD$coeff[2]), exp(cutD$coeff[3]), exp(cutD$coeff[4]), exp(cutD$coeff[5]), exp(cutD$coeff[6]), exp(cutD$coeff[7]), exp(cutD$coeff[8])))
      lam0 <- as.numeric(c(exp(cutND$coeff[1]), exp(cutND$coeff[2]), exp(cutND$coeff[3]), exp(cutND$coeff[4]), exp(cutND$coeff[5]), exp(cutND$coeff[6]), exp(cutND$coeff[7]), exp(cutND$coeff[8])))
      ksi <- c(1/cutD$theta, 1/cutND$theta) 
      theta <- c(cutD$theta, cutND$theta)
      
      sensitivity <- cumprod(1-lam1)
      specificity <- 1-cumprod(1-lam0)
      
      Pdata$lam <- c(rep(lam1,13), rep(lam0,13))
      Pdata$mu <- Pdata$r*Pdata$lam
      
      #rho_thres <- optimize(rho_thres_estimation, interval = c(0.001,0.99), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      ab <- optimize(rho_dis_estimation, interval = c(0, sqrt(theta[1]*theta[2])), th1 = theta[1], th0 = theta[2], data = Pdata, maximum = TRUE)$maximum
      rho_dis <- ab*sqrt(ksi[1]*ksi[2])
    }
    ## save estimates
    results_poisson[[i-12]][j,] <- c(sensitivity, specificity, ksi, 0, rho_dis)
    print(j)
  }
  results_i <- as.data.frame(results_poisson[[i-12]])
  colnames(results_i) <- c("sens1", "sens2", "sens3", "sens4", "sens5", "sens6", "sens7", "sens8", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "frail1", "frail0", "rhothres", "rhodis")
  name <- paste("poisson_scenario", paste(i-12), ".csv", sep = "")
  write.csv(results_i, file = name)
}

## write csv with all results
files <-  list.files(pattern = "*.csv")
all <- do.call(rbind, lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE)))
all$scenario <- rep(seq(1:16), each = 1000)
write.csv(all, file = "poisson_all.csv")
