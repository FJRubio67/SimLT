## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

# Required packages
#pak::pak("FJRubio67/SimLT")
library(SimLT)
library(msm)
library(lubridate)

#------------------------------------------------------------------------
# Life tables
#------------------------------------------------------------------------
LT <- as.data.frame(read.table("England_LT_2010_2015_dep_gor.txt",header = T))
head(LT)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Design matrix simulation 
#------------------------------------------------------------------------
# Sample size
n <- 1000
# Design matrix
Xcf <- simDesMatrix(seed = 123, n = n, admin.cens = "31-12-2015", scale.age = FALSE,
                    site = "colon", sex = "female")
head(Xcf)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Follow-up period
#------------------------------------------------------------------------
# fixed points (new year and last day of follow-up)
year.init <- 2010
last.day <- as.Date("2015-12-31")
fixed.times <- c(seq(as.Date("2010-01-01"), as.Date("2015-01-01"), by = "year"),last.day)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------
# Time points and hazard rates
#------------------------------------------------------------------------

# Date of diagnosis and date of birth (29 February is moved to 1 March)
dod <- as.Date(Xcf$date.diag)
dob <- as.POSIXlt(dod - Xcf$age*365.25)
leap <- dob$mday == 29 & dob$mon == 1
dob$mday[leap] <- 1
dob$mon[leap] <- 2
dob <- as.Date(format(dob, "%Y-%m-%d"))

# sex
sex0 <- 2

# Candidate time points for all patients at once:
# date of diagnosis, new years, birthdays, and last day of follow-up
yrs <- year.init:year(last.day)
cand <- data.frame(
  id   = c(1:n, rep(1:n, each = length(yrs)), rep(1:n, each = length(yrs)), 1:n),
  date = c(dod,
           rep(as.Date(paste0(yrs, "-01-01")), n),
           as.Date(paste0(rep(yrs, n), rep(format(dob, "-%m-%d"), each = length(yrs)))),
           rep(last.day, n)))
# Keep the points between diagnosis and the last day of follow-up, sorted by patient and date
cand <- cand[cand$date >= dod[cand$id] & cand$date <= last.day, ]
cand <- unique(cand)
cand <- cand[order(cand$id, cand$date), ]

# Age (completed years) at each time point
cand$age <- year(cand$date) - year(dob[cand$id]) -
  (format(cand$date, "%m%d") < format(dob[cand$id], "%m%d"))

# Extracting the hazard rates from the life table for each patient at "age + t":
# a single hashed lookup on (year, sex, dep, age, gor) for all patients
LT.key <- paste(LT$X_year, LT$sex, LT$dep, LT$age, LT$gor)
key <- paste(year(cand$date), sex0, Xcf$dep[cand$id], cand$age, Xcf$gor[cand$id])
cand$rate <- LT$rate[match(key, LT.key)]

# Dates, lengths of times in days between time points, and hazard rates for each patient
dates  <- unname(split(cand$date, cand$id))
times  <- lapply(dates, function(d) c(0, as.numeric(diff(d))))
hrates <- unname(split(cand$rate, cand$id))

# List containing the time points, hazard rates, and dates
Xcf_rates <- list(times = times, hrates = hrates, dates = dates)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##################################################################################################
# Population survival times
##################################################################################################
simC <- sim_pophaz(seed = 123,lst = Xcf_rates, method = "direct")
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


