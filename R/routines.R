
###################################################################################################
#' Piecewise Exponential hazard function
###################################################################################################
#' @param x,q :	  vector of quantiles.
#' @param p	  :   vector of probabilities.
#' @param n	  :   number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param rate :	vector of rates.
#' @param t	   : vector of the same length as rate, giving the times at which the rate changes. The first element of t should be 0, and t should be in increasing order.
#' @param log :	logical; if TRUE, log function is returned.
#' @return Piecewise Exponential hazard function
#' @export
hpwexp  <- function(x, rate = 1, t = 0, log = FALSE){
  logpdf <- dpexp(x, rate = rate, t = t, log = TRUE)
  logsurv <- ppexp(x, rate = rate, t = t, lower.tail = FALSE, log.p = TRUE)
  if(log) return( logpdf - logsurv )
  else return( exp( logpdf - logsurv ) )
}

#------------------------------------------------------------------
#' Piecewise Exponential cumulative hazard function
#------------------------------------------------------------------
#' @param x,q :	  vector of quantiles.
#' @param p	  :   vector of probabilities.
#' @param n	  :   number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param rate :	vector of rates.
#' @param t	   : vector of the same length as rate, giving the times at which the rate changes. The first element of t should be 0, and t should be in increasing order.
#' @param log :	logical; if TRUE, log function is returned.
#' @return Piecewise Exponential cumulative hazard function
#' @export
chpwexp  <- function(x, rate = 1, t = 0){
  logsurv <- ppexp(x, rate = rate, t = t, lower.tail = FALSE, log.p = TRUE)
  return( - logsurv )
}

#--------------------------------------------------------------------------
#' Simulation from a Piecewise Exponential (piece-wise constant rate)
#--------------------------------------------------------------------------
#' @param seed : seed for simulation
#' @param x,q :	  vector of quantiles.
#' @param p	  :   vector of probabilities.
#' @param n	  :   number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param rate :	vector of rates.
#' @param t	   : vector of the same length as rate, giving the times at which the rate changes. The first element of t should be 0, and t should be in increasing order.
#' @param log :	logical; if TRUE, log function is returned.
#' @return generates random numbers from a Piecewise Exponential hazard function
#' @export
rsim.pwexp <- function(seed, n, rate = 1, t = 0){
  rate <- as.vector(rate);
  t <- as.vector(t);
  set.seed(seed)
  u <- runif(n)
  vals <- -log(1-u)
  # maximum follow-up time
  maxfup <- max(t)
  # cumulative hazard at maxfup
  chmax <- chpwexp(maxfup, rate = rate, t = t)
  # censoring indicators
  ind.cens <- (vals > chmax)
  # simulated times
  sim.times <- rep(0,n)
  # censored times
  sim.times[ind.cens] <- maxfup
  # observed times
  ind.obs <- (1:n)[!ind.cens]
  # Probability Integral Transform based on the cumulative hazard
  for(i in ind.obs){
    chaz <- Vectorize(function(z) chpwexp(z, rate = rate, t = t) - vals[i])
    sim.times[i] <- uniroot(f = chaz, interval = c( 0, maxfup ), maxiter = 10000 )$root
  }
  out <- list(sim.pop = sim.times, status.pop = (1-ind.cens))
  return(out)
}

#' Which element in a vector is one
#' @param vec: a vector
#' @return which values are 1
#' @export
which1 <- function(vec) which(vec == 1)

###################################################################################################################
#' Function to simulate a design matrix. Based on real data on colon and lung (male and female) patients in England.
###################################################################################################################
#------------------------------------------------------------------------------------------------------------------
#' @param n : Sample size
#' @param seed : seed for simulation ("DD-MM-YYYY")
#' @param  years.diag : years of diagnosis in the simulation (scalar or vector)
#' @param  Returns: ID, sex (0 = male, 1 = female), dod (date of diagnosis), IMD (Index of Multiple Deprivation),
#' @param  admin.cens (year of administrative censoring),
#' @param  sex    :  sex ("female", "male")
#' @param  site     : cancer site ("lung" or "colon")
#' @return a design matrix containing the ID, date.diag, dep, gor, age, and fup.time.
#' @export
#------------------------------------------------------------------------------------------------------------------
simDesMatrix <- function(seed = 123, n = 100, admin.cens = "31-12-2015", scale.age = FALSE,
                         site = "colon", sex = "female"){
  #-----------------------------
  # patient ID
  #-----------------------------
  ID <- as.vector(1:n)
  #-----------------------------
  # Seed
  #-----------------------------
  set.seed(seed)
  # dep
  #tempmat <- t(rmultinom(n = n,size = 1,rep(0.2,5)))
  if(site == "colon" & sex == "female") tempmat <- t(rmultinom(n = n,size = 1, prob = c(0.21,0.22,0.22,0.20,0.15)))
  if(site == "colon" & sex == "male")   tempmat <- t(rmultinom(n = n,size = 1, prob = c(0.22,0.22,0.21,0.19,0.16)))
  if(site == "lung" & sex == "female")  tempmat <- t(rmultinom(n = n,size = 1, prob = c(0.14,0.17,0.20,0.24,0.25)))
  if(site == "lung" & sex == "male")   tempmat <- t(rmultinom(n = n,size = 1, prob = c(0.14,0.18,0.20,0.23,0.25)))
  dep <- apply(tempmat,1,which1)
  #-----------------------------
  # gor
  #-----------------------------
  #tempmat <- t(rmultinom(n = n,size = 1,rep(1/9,9)))
  # gor <- vector()
  tempmat <- matrix(0, ncol = 9, nrow = n)
  if(site == "colon" & sex == "female"){
    for(i in 1:n){
      if(dep[i]==1) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.37,12.54,8.15,8.34,9.32,13.22,6.10,27.41,11.55)/100))
      if(dep[i]==2) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.32,11.61,8.92,8.97,10.31,15.78,7.94,17.98,15.17)/100))
      if(dep[i]==3) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.51,11.47,9.00,9.22,9.76,15.07,8.46,17.00,16.51)/100))
      if(dep[i]==4) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(7.20,12.71,8.64,8.54,12.16,11.42,12.41,13.36,13.56)/100))
      if(dep[i]==5) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(10.93,22.33,14.30,5.67,11.40,5.13,17.95,6.34,5.95)/100))
    }}

  if(site == "colon" & sex == "male"){
    for(i in 1:n){
      if(dep[i]==1) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.18,11.90,8.48,9.05,9.66,14.10,5.58,27.47,10.58)/100))
      if(dep[i]==2) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.31,11.18,10.86,10.38,10.82,12.74,6.63,18.97,15.11)/100))
      if(dep[i]==3) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(4.41,10.64,9.39,9.82,10.39,14.19,8.22,16.96,15.98)/100))
      if(dep[i]==4) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(8.20,13.09,10.72,9.48,9.58,12.52,11.24,12.94,12.23)/100))
      if(dep[i]==5) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(11.40,23.03,12.50,7.12,14.99,4.22,15.16,6.37,5.21)/100))
    }}

  if(site == "lung" & sex == "female"){
    for(i in 1:n){
      if(dep[i]==1) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.87,13.63,9.13,8.89,7.88,13.44,6.19,26.15,10.82)/100))
      if(dep[i]==2) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(4.07,13.12,11.12,9.75,9.64,13.87,7.64,17.94,12.85)/100))
      if(dep[i]==3) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(5.46,13.05,12.27,9.64,9.47,12.58,8.80,16.25,12.48)/100))
      if(dep[i]==4) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(9.52,16.18,12.26,9.57,10.05,10.13,11.06,11.22,10.01)/100))
      if(dep[i]==5) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(15.11,25.88,15.62,5.56,10.85,3.56,14.02,4.84,4.56)/100))
    }}

  if(site == "lung" & sex == "male"){
    for(i in 1:n){
      if(dep[i]==1) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.04,13.87,8.44,9.36,8.13,13.91,6.51,26.13,10.61)/100))
      if(dep[i]==2) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(3.59,11.94,10.69,8.98,10.75,14.88,6.52,18.95,13.70)/100))
      if(dep[i]==3) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(4.31,12.63,11.61,9.66,10.22,12.97,8.27,15.76,14.57)/100))
      if(dep[i]==4) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(8.10,14.26,11.20,9.57,10.43,10.48,12.86,12.67,10.43)/100))
      if(dep[i]==5) tempmat[i,] <- t(rmultinom(n = 1,size = 1, prob = c(12.16,24.48,14.91,6.77,12.90,4.02,15.26,5.09,4.41)/100))
    }}
  gor <- apply(tempmat,1,which1)

  #-----------------------------
  # Date of diagnosis
  #-----------------------------
  date.diag0 <- sample(x = 0:364, size = n, replace = TRUE)
  date.diag <- as.Date(date.diag0, origin = "2010-01-01")
  format(date.diag, "%d-%m-%Y")
  #format(date.diag, "%d%b%Y")
  #-------------------------------------------------------------------
  # Age at diagnosis (based on summary statistics from real data)
  #-------------------------------------------------------------------
  if(site == "colon" & sex == "female") points0 <- c(18,32.00,49.61,56.69,65.55,74.86,82.60,87.76,90.34,94.94)
  if(site == "colon" & sex == "male")   points0 <- c(18,37.73,50.92,57.17,64.40,72.72,79.89,85.54,88.12,92.46)
  if(site == "lung" & sex == "female")  points0 <- c(18,43.21,52.49,57.44,64.65,73.18,80.94,86.31,88.98,93.68)
  if(site == "lung" & sex == "male")    points0 <- c(18,44.34,53.77,58.09,65.23,72.90,80.06,85.43,88.24,92.44)
  probs0 <- diff(c(0,0.01,0.05,0.1,0.25,0.5,0.75,0.90,0.95,0.99))
  k <- sample(x = 1:length(probs0), size = n, prob = probs0, replace = TRUE)
  age <-  runif(n, points0[k], points0[k+1])
  if(scale.age) age <- as.vector(scale(age))
  #-----------------------------
  # Maximum follow-up time
  #-----------------------------
  admin.cens <- as.Date(admin.cens, "%d-%m-%Y")
  format(admin.cens, "%d-%m-%Y")
  fup.time <- as.vector(admin.cens - date.diag)/365.25
  #-----------------------------
  # data frame
  #-----------------------------
  data <- data.frame(ID = ID, date.diag = date.diag, dep = dep, gor = gor, age = age, fup.time = fup.time)
  #-----------------------------
  # output
  #-----------------------------
  return(data)
}


###################################################################################################################
#' Function to simulate survival times from the population hazard
###################################################################################################################
#' @param lst : list with time points at which the piecewise hazard jumps (times) and the values of the piecewise hazard (hrates)
#' @param seed: seed for simulation
#' @return a list containing the simulated survival times and the ID
#' @export

sim_pophaz <- function(seed,lst){
  sim.pop <- vector()
  ind.pop <- vector()
  n <- length(lst$times)
  for(i in 1:n){
    times.temp <- cumsum(lst$times[[i]])/365.25
    rates.temp <- lst$hrates[[i]]
    temp.sim <- rsim.pwexp(seed = (seed+i), n = 1, rate = rates.temp, t = times.temp)
    sim.pop[i] <- temp.sim$sim.pop
    ind.pop[i] <- temp.sim$status.pop
  }
  out <- list(sim.pop = sim.pop, status = ind.pop)
}


