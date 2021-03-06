library(rje)
library(mvtnorm)


## function to set up the data for application of the bivariate method in example data analysis

data.for.bivariate.application <- function(original_data)
{
  ## Argument
  # original_data: dataset containing each patient diagnostic score. It should include 4 columns "study", "score", "status" and "weight"
  #                column "study": a study id
  #                column "score": the diagnostic score of a patient
  #                column "status": indicates the true disease status (truly healthy or truly diseased) as given by some reference / gold standard test
  #                column "weight": if weighting was used in a study. If no weight, set all entries of "weight" to 1
  original_data$score <- as.integer(original_data$score)
  phqBiv <- as.data.frame(matrix(NA, nrow = 728, ncol = 6))
  colnames(phqBiv) <- c("cutoff", "study", "POS", "n", "d1", "d0")
  phqBiv$cutoff <- rep(seq(0,27), each = 26)
  phqBiv$study <- rep(rep(seq(1,13), each = 2), 28)
  phqBiv$d1 <- rep(c(1,0), 364)
  phqBiv$d0 <- rep(c(0,1), 364)
  
  # compute n = c(total # diseased, total # healthy) per study
  size <- c()
  for(i in 1:13)
  {
    size <- c(size, sum(original_data$weight[which(original_data$study==i & original_data$status==1)]), sum(original_data$weight[which(original_data$study==i & original_data$status==0)]))
  }
  phqBiv$n <- rep(size, 28)
  
  # compute POS+prop
  aug <- 0
  phqBiv$POS <- rep(0, 728)
  ss <- c()
  for(i in 1:13)
  {
    stu <- original_data[which(original_data$study==i), ]
    stud1 <- stu[which(stu$status==1), ]
    stud0 <- stu[which(stu$status==0), ]
    ss <- c(ss, nrow(stud1), nrow(stud0))
    k <- -1
    for(j in seq(1, 703, 26))
    {
      phqBiv$POS[j+aug] <- sum(stud1$weight[which(stud1$score>k)])
      phqBiv$POS[j+1+aug] <- sum(stud0$weight[which(stud0$score>k)])
      k <- k+1
    }
    aug <- aug+2
  }
  phqBiv$prop <- phqBiv$POS/phqBiv$n
  phqBiv$ss <- rep(ss, 28) # weighted sample size
  return(phqBiv)
}




## function to set up the data for application of the poisson method in example data analysis

data.for.poisson.application <- function(original_data)
{
  ## Argument
  # original_data: dataset containing each patient diagnostic score. It should include 4 columns "study", "score", "status" and "weight"
  #                column "study": a study id
  #                column "score": the diagnostic score of a patient
  #                column "status": indicates the true disease status (truly healthy or truly diseased) as given by some reference / gold standard test
  #                column "weight": if weighting was used in a study
  original_data$score <- as.integer(original_data$score)
  ydis <- c()
  yndis <- c()
  rdis <- c()
  rndis <- c()
  studyD <- c()
  studyND <- c()
  nd <- c()
  nnd <- c()
  
  for(i in 1:13) # loop over the studies
  {
    phqi <- original_data[which(original_data$study==i), ]
    phqD <- phqi[which(phqi$status==1), ] # diseased
    phqND <- phqi[which(phqi$status==0), ] # healthy
    
    for(j in 0:13)
    {
      ydis <- c(ydis, sum(phqD$weight[which(phqD$score==j)]))
      yndis <- c(yndis, sum(phqND$weight[which(phqND$score==j)]))
      
      rdis <- c(rdis, sum(phqD$weight[which(phqD$score>=j)]))
      rndis <- c(rndis, sum(phqND$weight[which(phqND$score>=j)]))
    }
    
    nd <- c(nd, rep(dim(phqD)[1], 14)) 
    nnd <- c(nnd, rep(dim(phqND)[1], 14))
    
    studyD <- c(studyD, rep(i, 14))
    studyND <- c(studyND, rep(i, 14))
  }
  
  phqPois <- as.data.frame(cbind(c(ydis,yndis), c(rdis,rndis), c(studyD,studyND), c(nd,nnd)))
  colnames(phqPois) <- c("y", "r", "study", "n")
  phqPois$cutoff <- as.factor(rep(seq(1,14), dim(phqPois)[1]/14))
  phqPois$dis <- c(rep(1, dim(phqPois)[1]/2), rep(0, dim(phqPois)[1]/2))
  return(phqPois)
}




## function to set up the data for application of the bivariate method in simulations

data.for.bivariate.simulation <- function(x0, x1, nthres)
{
  ## Arguments
  # x0: a dataframe with # rows = # of studies, and # columns = # of category in the diagnostic test
  #     each element is the number of truly healthy patients in the study (row) with results in that category (column)
  # x0: a dataframe with # rows = # of studies, and # columns = # of category in the diagnostic test
  #     each element is the number of truly diseased patients in the study (row) with results in that category (column)
  # nthres: total number of thresholds in the diagnostic test
  
  num.cat <- as.vector(rbind(apply(x1, 1, function(x) rev(cumsum(rev(x)))), apply(x0, 1, function(x) rev(cumsum(rev(x)))))) # number of TP and FP for each threshold -> to estimate sens./1-spec.
  cutoff <- rep(seq(0, nthres), 13*2) # cutoff 0 to nthres, 0 not interesting
  study <- sort(rep(1:13, 2*(nthres+1))) # identify 13 studies
  bdata <- as.data.frame(cbind(num.cat, study, cutoff)) 
  bdata$ni <- rep(bdata$num.cat[which(bdata$cutoff==0)], each = nthres+1) # total number of truly diseased / healthy by study
  bdata$prop <- bdata$num.cat/bdata$ni # estimates of sensitivity / 1-specificity for each threshold/each study
  bdata$d1 <- rep(c(rep(1, nthres+1), rep(0, nthres+1)), 13) # d1=1 for diseased, =0 healthy
  bdata$d0 <- rep(c(rep(0, nthres+1), rep(1, nthres+1)), 13) # d0=1 for healthy, =0 diseased
  return(bdata)
}




## function to set up the data for application of the poisson method in simulations

data.for.poisson.simulation <- function(x0, x1, nthres)
{
  ## Arguments
  # x0: a dataframe with # rows = # of studies, and # columns = # of category in the diagnostic test
  #     each element is the number of truly healthy patients in the study (row) with results in that category (column)
  # x0: a dataframe with # rows = # of studies, and # columns = # of category in the diagnostic test
  #     each element is the number of truly diseased patients in the study (row) with results in that category (column)
  # nthres: total number of thresholds in the diagnostic test
  
  y <- c(as.vector(t(x1[,1:nthres])), as.vector(t(x0[,1:nthres])))
  r <- c(as.vector(apply(x1, 1, function(x) rev(cumsum(rev(x)))[-length(x)])), as.vector(apply(x0, 1, function(x) rev(cumsum(rev(x)))[-length(x)])))
  n <- c(apply(x1, 1, function(x)rep(sum(x),nthres)), apply(x0, 1, function(x)rep(sum(x),nthres)))
  cutoff <- rep(seq(1, nthres), 13*2) # cutoff 0 to nthres, 0 not interesting
  study <- rep(sort(rep(1:13, nthres)), 2) 
  dis <- c(rep(1,13*nthres), rep(0,13*nthres)) # dis=1 for diseased, dis=0 healthy 
  pdata <- as.data.frame(cbind(study, cutoff, y, r, n, dis))
  return(pdata)
}




### code to simulate diagnostic accuracy data for simulations with fixed number of threshold to 8

simulation <- function(rhodis, rhothres, sigma_sens, sigma_spec, nd, nnd, sdd) 
{
  ## rhodis: the correlation between sensitivity and specificity at each threshold
  ## rhothres: the correlation between sensitivities and specificities across thresholds 
  ## sigma_sens: variance of the random effects associated to sensitivity
  ## sigma_spec: variance of the random effects associated to specificity
  ## nd: mean # truly diseased patients per study 100 or 50
  ## nnd: mean # truly healthy patients per study
  ## sdd: standard deviation for # truly diseased patients 80 or 40
  
  ## specifity true sensitivity and specificity
  sens <- c(0.94, 0.91, 0.88, 0.84, 0.79, 0.74, 0.67, 0.57)
  spec <- c(0.74, 0.79, 0.83, 0.87, 0.89, 0.91, 0.93, 0.95)
  ss.overall <- logit(c(sens, 1-spec))
  
  ## correlation structure for random noise
  mat <- matrix(c(sigma_sens, rhodis*sqrt(sigma_sens*sigma_spec), rhodis*sqrt(sigma_sens*sigma_spec), sigma_spec),2,2)
  matti <- matrix(1, 8, 8)
  matti <- rhothres^(abs(row(matti)-col(matti)))
  Sig <- mat %x% matti
  
  ## generate diagnostic data
  n <- 13 # number of studies
  n1i <- round(rnorm(n, mean = nd, sd = sdd), 0) # number of truly diseased patients for each study from a Normal (nd, sdd) -> mimicks PHQ9
  n0i <- round(rnorm(n, mean = nnd, sd = 150), 0) # number of truly healthy patients for each study from a Normal (nnd, sd=150) -> mimicks PHQ9
  while (any(c(n1i,n0i) <= 0)) # make sure we have n > 0
  {
    n1i <- round(rnorm(n, mean = nd, sd = sdd), 0) 
    n0i <- round(rnorm(n, mean = nnd, sd = 150), 0)
  }
  ss.specific <- rmvnorm(n, ss.overall, Sig) # study-specific logit sens / 1-spec
  x0 <- x1 <- matrix(NA, n, 9) # matrices for number of diseased (x1)/healthy (x0) falling in each 9 category of 
  
  # the diagnostic data for each study
  for(i in 1:n) # for each study
  {
    sensi <- expit(ss.specific[i,1:8]) # sens in study i for all 8 thresholds
    speci <- 1 - expit(ss.specific[i,9:16]) # spec in study i for all 8 thresholds
    # p0 probability of falling in each category of the test for healthy group
    p0 <- c(speci, 1) - c(0, speci)
    # p1 probability of falling in each category of the test for diseased group
    p1 <- c(1,sensi) - c(sensi, 0)
    
    # make sure sens/spec are monotone
    while (any(c(p0,p1) <= 0)) 
    {
      ss.specific[i,] <- rmvnorm(1, ss.overall, Sig)
      sensi <- expit(ss.specific[i,1:8]) 
      speci <- 1 - expit(ss.specific[i,9:16]) 
      p0 <- c(speci, 1) - c(0, speci)
      p1 <- c(1,sensi) - c(sensi, 0)
    }
    # number of healthy patients per category
    x0i <- sample(0:8, n0i[i], replace = TRUE, prob = p0)
    x0[i,] <- table(factor(x0i, levels = 0:8))
    # number of diseased patients per category
    x1i <- sample(0:8, n1i[i], replace = TRUE, prob = p1)
    x1[i,] <- table(factor(x1i, levels = 0:8))
  }
  return(list(x0 = x0, x1 = x1))
}











