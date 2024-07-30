## ------------------------------------------------------------------------------------------------------------------
rm(list=ls())

# Required packages
#library(devtools)
#install_github("FJRubio67/SimLT")
library(SimLT)
library(msm)
library(lubridate)

#------------------------------------------------------------------------
# Life tables
#------------------------------------------------------------------------
LT <- as.data.frame(read.table("England_LT_2010_2015_dep_gor.txt",header = T))
head(LT)


## ------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Design matrix simulation 
#------------------------------------------------------------------------
# Sample size
n <- 1000
# Design matrix
Xcf <- simDesMatrix(seed = 123, n = n, admin.cens = "31-12-2015", scale.age = FALSE,
                    site = "colon", sex = "female")
head(Xcf)


## ------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Follow-up period
#------------------------------------------------------------------------
# fixed points (new year and last day of follow-up)
year.init <- 2010
last.day <- as.Date("2015-12-31")
fixed.times <- c(seq(as.Date("2010-01-01"), as.Date("2015-01-01"), by = "year"),last.day)


## ------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Time points and hazard rates
#------------------------------------------------------------------------

# function to calculate ages at time points
ages <- Vectorize(function(dates){
  out <- time_length(difftime(as.Date(dates), as.Date(dob)), "years")
  return(out)
})

# Variables to be used from the life table
MLT <- LT[,c("X_year","sex","dep","age","gor")]

# Initialising times, dates, and hrates
times <- dates <- hrates <- list()

# sex
sex0 <- 2

# Extracting the hazard rates from the life tables for each patient at "age + t"

#pb = txtProgressBar(min = 0, max = n, initial = 1) # uncomment to track progress
for(i in 1:n){
  # Date of birth
  dob <- as.POSIXlt(as.Date(Xcf[i,2]) - Xcf[i,5]*365.25)
  if(day(dob)==29 & month(dob)==2) day(dob) <- 1; month(dob) <- 3
  # Birthday on 2010
  bday1 <- dob; year(bday1) <- year.init; bday1 <-as.Date(bday1);
  # Date of diagnosis
  dod <- as.Date(Xcf[i,2])
  # Variable times (Birthdays + date of diagnosis)
  var.times <- c(seq(as.Date(bday1),last.day, by = "year"),dod)
  # Unique sorted time points
  time.points <- sort(unique(c(fixed.times,var.times)))
  # Age at each time point
  age.tp <- ages(time.points)
  #*********************************************************
  # Removing dates before the date of diagnosis
  #*********************************************************
  ind <- length(which(time.points==dod))
  if(ind == 1 ){
    # Length of times in days between time points
    times[[i]] <- c(0,as.numeric(diff(time.points)))
    # Length of times in days between time points
    dates[[i]] <- time.points
    # Age at each time point (integer part)
    age.tp <- floor(age.tp)
  }
  if(ind > 1 ){
    # Length of times in days between time points
    times[[i]] <- c(0,as.numeric(diff(time.points))[-c(1:(ind-1))])
    # Length of times in days between time points
    dates[[i]] <- time.points[-c(1:(ind-1))]
    # Age at each time point (integer part)
    age.tp <- floor(age.tp[-c(1:(ind-1))])
  }
  MAT.ID <- cbind(1:length(age.tp),year(dates[[i]]),rep(sex0,length(age.tp)), rep(Xcf[i,"dep"],length(age.tp)),
                  age.tp,rep(Xcf[i,"gor"],length(age.tp)))
  MAT.ID <- as.data.frame(MAT.ID)
  colnames(MAT.ID) <- c("index","X_year", "sex", "dep", "age", "gor")
  hrates[[i]] <- merge(x = MAT.ID, y = LT, by = colnames(MAT.ID)[-1], all.x = TRUE, sort = FALSE)
  hrates[[i]] <- hrates[[i]][order(hrates[[i]]$index),]
  hrates[[i]] <- hrates[[i]]$rate
#  setTxtProgressBar(pb,i) # uncomment to track progress
}

# List containing the simulated survival times, hazard rates, and last date of follow-up
Xcf_rates <- list(times = times, hrates = hrates, dates = date)



## ------------------------------------------------------------------------------------------------------------------
##################################################################################################
# Population survival times
##################################################################################################
simC <- sim_pophaz(seed = 123,lst = Xcf_rates)
sim.pop <- simC$sim.pop
status.pop <- simC$status

# Censoring
table(status.pop)


library(survival)

# Kaplan-Meier estimator for the simulated times
kmsim <- survfit(Surv(sim.pop, status.pop) ~ 1)

# Comparison
plot(kmsim$time, kmsim$surv, type = "l", col = "black", lwd = 2, lty = 1, ylim = c(0,1),
     xlab = "Time", ylab = "Survival", cex.axis = 1.5, cex.lab = 1.5, main = "")


