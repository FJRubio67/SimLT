# NEEDS TO BE SPEED UP

#######################################################################
# function to calculate (approximate) ages at time points
#######################################################################
ages <- Vectorize(function(dates, dob){
  out <- time_length(difftime(as.Date(dates), as.Date(dob)), "years")
  return(out)
})


#######################################################################
# Function to calculate the population survival at specific time points
# It uses all variables in the life table
#######################################################################
# tps: time points of interest to evaluate the population survival.
#        These values represent the time length after the initial
#        point (date of diagnosis) in.
# LT: data frame containing the life table.
#     It must contain the variables "year", "sex", "age", "rate"
# des: design matrix containing the patient's characteristics.
#      It may contain additional variables not used in the matching step.
#      It must include "age" in years and date of diagnosis.
#      It must include the variable "sex" coded as in the life table.
#      Date of diagnosis must be entered in "YYYY-MM-DD".
# sdvars: socio-demographic variable names additional to "age" and "sex"
#         to be matched to the life tables (e.g. c("dep", "gor") )
# dodn: date of diagnosis variable name in the life table
# progress: progress bar
# dy: days in a year 
# Returns a list containing the survival probabilities as a matrix (psurv)

spLT <- function(tps, LT, des, sdvars = NULL, dodn, progress = FALSE, dy = 365.25){
  des <- as.data.frame(des); 
  times <- as.vector(tps)
  # Number of individuals
  n <- nrow(des)
  # Initialising times, dates, and hrates
  times <- dates <- hrates <- list()
  # Initialising population survival probabilities
  PSPs <- matrix(0, ncol = length(tps), nrow = nrow(des))
  maxt <-  max(tps) # maximum follow-up

  # Date of birth (approximate if no exact age provided)
  dob <- as.POSIXlt(as.Date(des[,dodn]) - des[,"age"]*dy)
  # Transforming date of birth if equal to 29th February  
 #  for(i in 1:n){
  #   if(day(dob[i])==29 & month(dob[i])==2) day(dob[i]) <- 1; month(dob[i]) <- 3
  # }
  temp <- (month(dob)==2 & day(dob)==29)
  dob[temp] <- dob[temp] + days(1)
  # initial year
  year_init <- year(des[,dodn])
  # Last date
  last_day <- as.Date(as.Date(des[,dodn]) + maxt*dy)
  # Birthdays at the initial year
  bday1 <- dob; year(bday1) <- year_init; bday1 <-as.Date(bday1);
  # New years
  ny1 <- as.Date(paste(year_init,"-01-01",sep="")); year(ny1) <- year_init; ny1 <-as.Date(ny1);
  # Date of diagnosis
  dod <- as.Date(des[,dodn])
  
    
  # Extracting the hazard rates from the life tables for each 
  # patient at "age + t"
  if(progress) pb = txtProgressBar(min = 0, max = n, initial = 1) 
  
for(i in 1:n){
  # Birthdays during period of interest
  bdays <- seq(as.Date(bday1[i]),last_day[i], by = "year")
  # New years
  nys <- seq(as.Date(ny1[i]),last_day[i], by = "year") 
  # Time points as dates
  tpsd <- as.Date(as.Date(des[i,dodn]) + tps*dy)
    
  # Birthdays & new years & time points & last day & date of diagnosis
  # Unique sorted time points
  time_points <- sort(unique(c(bdays, nys, tpsd, last_day[i], dod[i])))
  # Age at each time point
  age_tp <- ages(time_points, dob[i])
  #*********************************************************
  # Removing dates before the date of diagnosis
  #*********************************************************
  ind <- which(time_points==dod[i])
    if(ind == 1 ){
      # Length of times in days between time points
      times <- c(0,as.numeric(diff(time_points)))
      # dates starting at date of diagnosis
      dates <- time_points
      # Age at each time point (integer part)
      age_tp <- floor(age_tp)
    }
    if(ind > 1 ){
      # Length of times in days between time points
      times <- c(0,as.numeric(diff(time_points))[-c(1:(ind-1))])
      # dates starting at date of diagnosis
      dates <- time_points[-c(1:(ind-1))]
      # Age at each time point (integer part)
      age_tp <- floor(age_tp[-c(1:(ind-1))])
    }
    
  ltp <- length(age_tp)

  # Creating the ID matrix to merge with life table
  
  # Without additional socio-demographic variables
  if(is.null(sdvars)){
      MAT_ID <- cbind(1:length(age_tp),
                      year(dates),
                      age_tp,
                      rep(des[i,"sex"],ltp))
      MAT_ID <- data.table(MAT_ID)
      colnames(MAT_ID) <- c("index", "year", "age", "sex") 
    }
    
  # With additional socio-demographic variables
  if(!is.null(sdvars)){
    MAT_ID <- cbind(1:ltp,
                    year(dates),
                    age_tp,
                    rep(des[i,"sex"],ltp),
                    do.call("rbind", replicate(ltp, des[i,sdvars], simplify = FALSE)))
    MAT_ID <- data.table(MAT_ID)
    colnames(MAT_ID) <- c("index", "year", "age", "sex", sdvars)
  }
  # Merging life-tables and extracting hazard rates
  hrates <- merge(x = MAT_ID, y = LT, by = colnames(MAT_ID)[-1], 
                  all.x = TRUE, sort = FALSE)
  #hrates <- left_join(x = MAT_ID, y = LT, by = colnames(MAT_ID)[-1])
  hrates <- hrates[order(hrates$index),]
  hrates <- head(hrates$rate,-1)
  times <- times[-1]
  
  # Population survival at all time points 
  allps <- exp( - cumsum( hrates*times/dy  )  )
  
  # Extracting time points of interest 
  ts <- round(c(cumsum(times)/dy), digits = 12)
  MT <- cbind(1:length(ts), ts)
  colnames(MT) <- c("index", "tps")
  mtps <- as.matrix(round(tps, digits = 12))
  colnames(mtps) <- c("tps")
  
  indtps <- merge(x = mtps, y = MT, by.x = "tps", 
                  all.x = TRUE, sort = FALSE)[,2]
  
  # Population survival at time points of interest
  PSPs[i,] <- allps[indtps]
  
    # Progress bar
    if(progress)  setTxtProgressBar(pb,i) 
  }

  out <- list( psurv = as.matrix(PSPs), tps = tps)
  
  return(out)
      
}