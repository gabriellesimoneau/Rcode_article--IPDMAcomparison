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




## the following 3 functions are required for the estimation of rho_dis

# function to repeat vector v n times (written by M.Fiocco)
rowrep <- function(v,n)
{
  # Arguments
  # v: vector
  #	n: number of rows in the matrix
  # Output:
  # matrix of repeated rows, dimension: length(v)X n
  return(t(colrep(v, n)))  
}

# function to repeat vector v n times (written by M.Fiocco)
colrep <- function(v, n)  
{
  # Arguments
  # v: vector
  #	n: number of columns in the matrix
  # Output:
  # matrix of repeated columns, dimension: length(v)X n
  return(matrix(rep(v,n), length(v), n))
}

# function to calculate the log-likelihood contribution of one pair of observation (written M.Fiocco): useful for the estimation of rho_dis
loglike.dis <- function(y1, y2, mu1, mu2, theta1, theta2, ab) 
{
  # joint distribution of (yijd,yikd)
  mu1t <- mu1/theta1
  mu2t <- mu2/theta2
  alphaa.b <- theta1-ab
  alphab.a <- theta2-ab
  mu12t <- mu1t+mu2t
  if(alphaa.b>0) 
  {
    Pist1 <- dnbinom(y1:0, size = alphaa.b, mu = mu1t*alphaa.b)
  }
  if(alphaa.b<=0) 
  {
    alphaa.b <- 0.000001
    Pist1 <- dnbinom(y1:0, size = alphaa.b, mu = mu1t*alphaa.b)
  }
  Pist2 <- dnbinom(y2:0, size = alphab.a, mu = mu2t*alphab.a)
  P1 <- colrep(Pist1, y2+1)
  P2 <- rowrep(Pist2, y1+1)
  outerm <-outer(0:y1, 0:y2, "+")#dim=y1+1Xy2+1
  outerv <-as.vector(outerm)
  helpv <- as.vector(colrep(0:y1, y2+1)) 
  P3 <- matrix(dnbinom(outerv, size = ab, mu = mu12t*ab), y1+1, y2+1)
  P4 <- matrix(dbinom(helpv, outerv, mu1t/mu12t), y1+1, y2+1)
  P <- P1*P2*P3*P4  	           
  return(logP = log(sum(P)))
}




## function to estimate the correlation parameter rho_thres GOOD

rho_thres_estimation <- function(rho, data, th0, th1)
{
  ## Arguments
  # rho: correlation parameter, to maximize
  # data: dataset as returned by the functions data.for.poisson.example or data.for.poisson.simulation
  #       with 2 additional columns: lam and mu
  # th1: inverse of estimates of frailty1 when running the poisson model with the function poissonModel
  # th0: inverse of estimates of frailty0 when running the poisson model with the function poissonModel
  
  like <- 0
  for(i in 1:13) # Loop over the 13 studies in the example dataset (and simulation)
  {
    bloc <- data[which(data$study==i),]
    blocD <- bloc[which(bloc$dis==1),]
    blocND <- bloc[which(bloc$dis==0),]
    for(j in 1:7) 
    {
      if(dim(blocD)[1]>j)
      {
        for(k in (j+1):dim(blocD)[1])
        {
          p_sens <- 0
          yij1 <- blocD$y[j]
          yik1 <- blocD$y[k]
          rhojk <- rho^abs(j-k)
          muij1 <- blocD$mu[j]
          muik1 <- blocD$mu[k]
          theta1 <- th1
          for(l in 0:yij1)
          {
            for(m in 0:yik1)
            {
              p_sens <- p_sens + dnbinom(l, mu = muij1*(1-rhojk), size = theta1*(1-rhojk))*dnbinom(m, mu = muik1*(1-rhojk), size = theta1*(1-rhojk))*dnbinom(yij1+yik1-l-m, mu = (muij1+muik1)*rhojk, size = theta1*rhojk)*dbinom(yij1-l, size = yij1+yik1-l-m, p = muij1/(muij1+muik1))
            }
          }
        }
        like <- like + log(p_sens)
      }
      if(dim(blocND)[1]>j)
      {
        for(k in (j+1):dim(blocND)[1])
        {
          p_spec <- 0
          yij0 <- blocND$y[j]
          yik0 <- blocND$y[k]
          rhojk <- rho^abs(j-k)
          muij0 <- blocND$mu[j]
          muik0 <- blocND$mu[k]
          theta0 <- th0
          for(l in 0:yij0)
          {
            for(m in 0:yik0)
            {
              p_spec <- p_spec + dnbinom(l, mu = muij0*(1-rhojk), size = theta0*(1-rhojk))*dnbinom(m, mu = muik0*(1-rhojk), size = theta0*(1-rhojk))*dnbinom(yij0+yik0-l-m, mu = (muij0+muik0)*rhojk, size = theta0*rhojk)*dbinom(yij0-l, size = yij0+yik0-l-m, p = muij0/(muij0+muik0))  
            }
          }
          like <- like + log(p_spec)
        }
      }
    } 
  } 
  return(like)
}




## function to estimate correlation parameter rho_dis GOOD

rho_dis_estimation <- function(ab, data, th1, th0) 
{
  ## Arguments
  # rho: correlation parameter, to maximize
  # data: dataset as returned by the functions data.for.poisson.example or data.for.poisson.simulation
  #       with 2 additional columns: lam and mu
  # th1: inverse of estimates of frailty1 when running the poisson model with the function poissonModel
  # th0: inverse of estimates of frailty0 when running the poisson model with the function poissonModel
  
  loglike <- 0
  for(i in 1:13) 
  {
    bloc <- data[which(data$study==i), ]
    blocD <- bloc[which(bloc$dis==1), ]
    blocND <- bloc[which(bloc$dis==0), ]
    sj <- min(dim(blocD)[1], dim(blocND)[1])
    for(j in 1:sj) 
    {              
      p_j <- 0
      yij1 <- blocD$y[j]
      yij0 <- blocND$y[j]
      muij1 <- blocD$mu[j]
      muij0 <- blocND$mu[j]
      theta1 <- th1
      theta2 <- th0
      loglike <- loglike + loglike.dis(yij1, yij0, muij1, muij0, theta1, theta2, ab)
    }
  }  
  return(loglike)
}




## functions to estimate the standard error of the poisson estimate via parametric bootstrap

parametric_boot <- function(copula_input, data)
{
  ## Arguments
  #
  
  # Given theta and rho, generate N independent copies z_i*=c(z_i1*,...,z_i8*) from the multivariate gamma
  # N=13 studies
  generate.gamma <- rMvdc(13, copula_input)
  
  # Given lambda and r, derive mu_ij=lambda_ij*r_ij = column mu 
  # generate y_i*=c(y_i1*,...,y_i8*) with y_ij*~Pois(mu_ij*z_ij*) independent.
  data$zi <- as.vector(c(t(generate.gamma)[1:14,], t(generate.gamma)[15:28,]))
  for(i in 1:nrow(data)) data$y[i] <- rpois(1, lambda = data$mu[i]*data$zi[i])
  
  # apply poisson model
  esti <- poissonModel(data)
  data2 <- data[which(data$cutoff %in% 7:14),]
  data2$lam <- c(rep(esti$lambda1[7:14], 13), rep(esti$lambda0[7:14], 13))
  data2$mu <- data2$r*data2$lam
  theta1 <- 1/esti$frailty1
  theta0 <- 1/esti$frailty0
  
  # estimate rho_thres and rho_dis
  rho_thres <- optimize(rho_thres_estimation, interval = c(0.01,0.99), th1 = theta1, th0 = theta0, data = data2, maximum = TRUE)$maximum
  ab <- optimize(rho_dis_estimation, interval = c(0, sqrt(theta1*theta0)), th1 = theta1, th0 = theta0, data = data2, maximum = TRUE)$maximum
  rho_dis <- ab*sqrt(esti$frailty1*esti$frailty0)
  
  # return all parameters for which we need SD
  return(c(esti$lambda1[7:14], esti$lambda0[7:14], esti$frailty1, esti$frailty0, rho_thres, rho_dis, esti$pooled_sensitivity, esti$pooled_specificity))
}
