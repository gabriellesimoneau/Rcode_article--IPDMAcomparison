### functions to estimate poisson model from dataset

## fit poisson model and get pooled estimates
poissonModel <- function(data)
{
  # data needs to be set up in a certain way.
  # Each line gives information on sensitivity OR specificity for one threshold and one study
  # Dimension of data = number of studies X number of thresholds X 2
  # Column names
  #       - study: study id
  #       - cutoff: threshold considered
  #       - dis: indicator for true disease status, 1 = truly diseased, 0 = healthy
  #       - y: number of subjects (diseased subject if dis=1, healthy otherwise) with test results = cutoff - 1 
  #       - r: number of subjects (diseased subject if dis=1, healthy otherwise) with test results >= cutoff - 1 
  #       - n: total number of diseased subject if dis=1 or healthy subjects if dis=0 in that study
  
  # negative binomial models
  poisD <- data[which(data$dis==1),]
  poisND <- data[which(data$dis==0),] 
  cutD <- glm.nb(poisD$y ~ -1 + offset(log(poisD$r)) + as.factor(poisD$cutoff), weights = poisD$n)
  cutND <- glm.nb(poisND$y ~ -1 + offset(log(poisND$r)) + as.factor(poisND$cutoff), weights = poisND$n)
  
  # save estimates
  lam1 <- as.numeric(c(exp(cutD$coeff[1]), exp(cutD$coeff[2]), exp(cutD$coeff[3]), exp(cutD$coeff[4]), exp(cutD$coeff[5]), exp(cutD$coeff[6]), exp(cutD$coeff[7]), exp(cutD$coeff[8]), exp(cutD$coeff[9]), exp(cutD$coeff[10]), exp(cutD$coeff[11]), exp(cutD$coeff[12]), exp(cutD$coeff[13]), exp(cutD$coeff[14])))
  lam0 <- as.numeric(c(exp(cutND$coeff[1]), exp(cutND$coeff[2]), exp(cutND$coeff[3]), exp(cutND$coeff[4]), exp(cutND$coeff[5]), exp(cutND$coeff[6]), exp(cutND$coeff[7]), exp(cutND$coeff[8]), exp(cutND$coeff[9]), exp(cutND$coeff[10]), exp(cutND$coeff[11]), exp(cutND$coeff[12]), exp(cutND$coeff[13]), exp(cutND$coeff[14])))
  lambda <- c(lam1,lam0)
  ksi <- c(1/cutD$theta, 1/cutND$theta) # estimates fraily
  
  # corresponding sens/spec for th 7-14
  sens7 <- 1*(1-lam1[1])*(1-lam1[2])*(1-lam1[3])*(1-lam1[4])*(1-lam1[5])*(1-lam1[6])*(1-lam1[7])
  sens8 <- sens7*(1-lam1[8])
  sens9 <- sens8*(1-lam1[9])
  sens10 <- sens9*(1-lam1[10])
  sens11 <- sens10*(1-lam1[11])
  sens12 <- sens11*(1-lam1[12])
  sens13 <- sens12*(1-lam1[13])
  sens14 <- sens13*(1-lam1[14])
  
  spec7 <- 1-(1*(1-lam0[1])*(1-lam0[2])*(1-lam0[3])*(1-lam0[4])*(1-lam0[5])*(1-lam0[6])*(1-lam0[7]))
  spec8 <- 1-((1-spec7)*(1-lam0[8]))
  spec9 <- 1-((1-spec8)*(1-lam0[9]))
  spec10 <- 1-((1-spec9)*(1-lam0[10]))
  spec11 <- 1-((1-spec10)*(1-lam0[11]))
  spec12 <- 1-((1-spec11)*(1-lam0[12]))
  spec13 <- 1-((1-spec12)*(1-lam0[13]))
  spec14 <- 1-((1-spec13)*(1-lam0[14]))

  # extract lam1, lam0, frailty1, frailty0, sensitivity 7to14, specificity 7to14
  return(list(lambda1 = lam1, lambda0 = lam0, frailty1 = ksi[1], frailty0 = ksi[2], pooled_sensitivity = round(c(sens7, sens8, sens9, sens10, sens11, sens12, sens13, sens14), 2), pooled_specificity = round(c(spec7, spec8, spec9, spec10, spec11, spec12, spec13, spec14), 2)))
}




## functions to estimate correlation parameter rho_thres

# function to repeat vector v n times (following M.Fiocco)
rowrep <- function(v,n)
{
  # Input:
  # 	v: vector
  #	n: number of rows in the matrix
  # Output:
  # matrix of repeated rows, dimension: length(v)X n
  return(t(colrep(v,n)))
}

# function to repeat vector v n times (following M.Fiocco)
colrep <- function(v,n)
{
  # Input:
  # 	v: vector
  #	n: number of columns in the matrix
  # Output:
  # matrix of repeated columns, dimension: length(v)X n
  return(matrix(rep(v,n),length(v),n))
}

# function to calculate the log-likelihood contribution of one pair of observation (following M.Fiocco)
loglike.jk <- function(y1, y2, mu1, mu2, rho, theta, j, k)
{
  # joint distribution of (yijd,yikd)
  mu12 <- mu1+mu2
  rho.jk <- rho^(abs(j-k))
  Pist1 <- dnbinom(y1:0, size = theta*(1-rho.jk), mu = mu1*(1-rho.jk))
  Pist2 <- dnbinom(y2:0, size = theta*(1-rho.jk), mu = mu2*(1-rho.jk))
  P1 <- colrep(Pist1, y2+1)
  P2 <- rowrep(Pist2, y1+1)
  outerm <-outer(0:y1, 0:y2, "+") 
  outerv <-as.vector(outerm)
  helpv <- as.vector(colrep(0:y1, y2+1)) 
  P3 <- matrix(dnbinom(outerv, size = theta*rho.jk, mu = mu12*rho.jk), y1+1, y2+1)
  P4 <- matrix(dbinom(helpv, outerv, mu1/mu12), y1+1, y2+1)
  P <- P1*P2*P3*P4  	           
  return(logP = log(sum(P)))
}

rt <- function(rho, data, th1, th0)
{
  ## Arguments
  # rho: correlation parameter, to maximize
  # data: dataset as returned by the functions data.for.poisson.example or data.for.poisson.simulation
  # th1: inverse of estimates of frailty1 when running the poisson model with the function poissonModel
  # th0: inverse of estimates of frailty0 when running the poisson model with the function poissonModel
  
  loglike <- 0
  for(i in seq(1:13)) 
  {
    datai <- data[which(data$study==i), ]
    dataD <- datai[which(datai$dis==1), ]
    dataND <- datai[which(datai$dis==0), ]
    for(j in 7:13)
    {
      for(k in (j+1):14)
      {
        yij1 <- dataD$y[j-6]
        yik1 <- dataD$y[k-6]
        muij1 <- dataD$mu[j-6]
        muik1 <- dataD$mu[k-6]
        yij0 <- dataND$y[j-6]
        yik0 <- dataND$y[k-6]
        muij0 <- dataND$mu[j-6]
        muik0 <- dataND$mu[k-6]
        #print(i)
        loglike <- loglike + loglike.jk(yij1, yik1, muij1, muik1, rho, th1, j, k) + loglike.jk(yij0, yik0, muij0, muik0, rho, th0, j, k)
      }
    }
  }
  return(-loglike)
}

rt_simulation <- function(rho, data, th1, th0)
{
  loglike <- 0
  for(i in seq(1:13))
  {
    datai <- data[which(data$study==i),]
    dataD <- datai[which(datai$dis==1),]
    dataND <- datai[which(datai$dis==0),]
    for(j in 1:7)
    {
      for(k in (j+1):8)
      {
        yij1 <- dataD$y[j]
        yik1 <- dataD$y[k]
        muij1 <- dataD$mu[j]
        muik1 <- dataD$mu[k]
        yij0 <- dataND$y[j]
        yik0 <- dataND$y[k]
        muij0 <- dataND$mu[j]
        muik0 <- dataND$mu[k]
        #print(i)
        loglike <- loglike + loglike.jk(yij1, yik1, muij1, muik1, rho, th1, j, k) + loglike.jk(yij0, yik0, muij0, muik0, rho, th0, j, k)
      }
    }
  }
  return(-loglike)
}




## estimate correlation parameter rho_dis
rd <- function(ab, data, th1, th0) ## similar function to estimate rho_dis, same input
{
  like <- 0
  for(i in 1:13) # loop over 13 studies
  {
    bloc <- data[which(data$study==i),]
    blocD <- bloc[which(bloc$dis==1),]
    blocND <- bloc[which(bloc$dis==0),]
    for(j in 1:8) # 8 thresholds
    {
      p_j <- 0
      yij1 <- blocD$y[j]
      yij0 <- blocND$y[j]
      muij1 <- blocD$mu[j]/th1
      muij0 <- blocND$mu[j]/th0
      alphaa.b <- th1-ab
      alphab.a <- th0-ab
      for(l in 0:yij1)
      {
        for(m in 0:yij0)
        {
          p_j <- p_j + dnbinom(l, mu = muij1*alphaa.b, size = alphaa.b) * dnbinom(m, mu = muij0*alphab.a, size = alphab.a) * dnbinom(yij1+yij0-l-m, mu = (muij1+muij0)*ab, size = ab) * dbinom(yij1-l, size = yij1+yij0-l-m, p = muij1/(muij1+muij0))
        }
      }
      like <- like + log(p_j)
    }
  }  
  return(-like)
}




## functions to estimate the standard error of the poisson estimate via parametric bootstrap

parametric_boot <- function(R, copula_input, data)
{
  # Given theta and rho, generate N independent copies z_i*=c(z_i1*,...,z_i8*) from the multivariate gamma
  # N=13 studies
  generate.gamma <- rMvdc(13, copula_input)
  
  # Given lambda and r, derive mu_ij=lambda_ij*r_ij = column mu 
  # generate y_i*=c(y_i1*,...,y_i8*) with y_ij*~Pois(mu_ij*z_ij*) independent.
  data$zi <- as.vector(c(t(generate.gamma)[1:14,],t(generate.gamma)[15:28,]))
  for(i in 1:nrow(data)) data$y[i] <- rpois(1, lambda = data$mu[i]*data$zi[i])
  
  # apply poisson model
  esti <- poissonModel(data)
  data2 <- data[which(data$cutoff %in% 1:14),]
  data2$lam <- c(rep(esti$lambda1,13),rep(esti$lambda0,13))
  data2$mu <- data2$r*data2$lam
  theta1 <- 1/esti$frailty1
  theta0 <- 1/esti$frailty0
  
  # estimate rho_thres and rho_dis
  rho_thres <- optimize(rt, interval = c(0.01,0.99), th1 = theta1, th0 = theta0, data = data2)
  ab <- optimize(rd, interval = c(0, sqrt(theta1*theta0)), th1 = theta1, th0 = theta0, data = data2)$minimum
  rho_dis <- ab*sqrt(theta1*theta0)
  
  # return all parameters for which we need SD
  return(c(lam1[7:14], lam0[7:14], esti$frailty1, esti$frailty0, rho_thres, rho_dis, esti$pooled_sensitivity, esti$pooled_specificity))
}
