## --------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

# Required packages
#library(devtools)
#install_github("FJRubio67/SimLT")
library(SimLT)
library(msm)
library(lubridate)
library(data.table)
source("routinesSPLT.R")

#------------------------------------------------------------------------
# Life tables
#------------------------------------------------------------------------
LT <- data.table(read.table("England_LT_2010_2015_dep_gor.txt",header = T))
head(LT)


## --------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Design matrix simulation 
#------------------------------------------------------------------------
# Sample size
n <- 1000
# Design matrix
Xcf <- simDesMatrix(seed = 123, n = n, admin.cens = "31-12-2015", scale.age = FALSE,
                    site = "colon", sex = "female")
head(Xcf)

Xcf$sex <- rep(2, nrow(Xcf))



## --------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Population survival probabilities at specific time points
#------------------------------------------------------------------------

# Time points of interest
tps <- seq(0.1,5, by = 0.5)

# Calculation of the survival probabilities
sps <- spLT(tps = tps, LT = LT, des = Xcf, sdvars = c("dep","gor"), 
     dodn = "date.diag", progress = FALSE, dy = 365.25)

head(sps$psurv)


# Some plots for specific individuals
index <- c(1,200,400,600,800,1000)

# Individual characteristics
Xcf[index,]

# Survival probabilities at times of interest
par(mfrow=c(2,3))
for(i in index){
  plot(sps$tps, sps$psurv[i,], type = "p", ylim = c(0,1), lwd = 2, main = paste("Individual ", i),
       xlab = "time", ylab = "Survival", cex.axis = 1.5, cex.lab = 1.5)
}


