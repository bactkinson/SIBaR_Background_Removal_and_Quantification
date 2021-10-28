rm(list=ls())
start_time <- Sys.time()
require(lubridate)
require(tidyverse)
require(mgcv)
require(depmixS4)
require(zoo)
## FUNCTION DEFINITIONS ##
daysFromYearStartConversion <- function(inputVec){
  outputVec <- numeric(length(inputVec))
  se <- as.numeric(difftime(inputVec[year(inputVec)==2017],as.POSIXct("2017-01-01 00:00:00 CDT"),units="days"))
  et <- as.numeric(difftime(inputVec[year(inputVec)==2018],as.POSIXct("2018-01-01 00:00:00 CDT"),units="days"))
  outputVec[year(inputVec)==2017] <- se
  outputVec[year(inputVec)==2018] <- et
  return(outputVec)
}
secondsFromDayStartConversion <- function(inputVec){
  outputVec <- as.numeric(difftime(inputVec,trunc(inputVec,"days"),units="secs"))
  return(outputVec)
}
# Function which applies depmixS4 model to data with time as response covariate
applyTimeDepmix <- function(poll,time){
  fm.temp <- tryCatch(
    {
      trst <- runif(4)
      model <- depmixS4::depmix(poll~time, nstates=2, data=data.frame(poll,time),family=gaussian(link=identity), trstart=trst)
      fm.temp <- depmixS4::fit(model)
    },
    error=function(cond) {
      return(NA)
    }
  )
  return(fm.temp)
}

findMode <- function(vec){
  vals <- unique(vec)
  mode <- vals[1]
  mode.idx <- length(which(vec==mode))
  for(a in 2:length(vals)){
    temp.idx <- length(which(vec==vals[a]))
    if(temp.idx>mode.idx){
      mode <- vals[a]
      mode.idx <- length(which(vec==mode))
    }
  }
  return(mode)
}

smoothData <- function(poll,time,t.interval){
  if(t.interval==0){
    return(list(poll,time))
  } else{
    poll.zoo <- zoo::zoo(poll, order.by = time)
    total.time <- seq(time[1],time[length(time)],1)
    na.zoo <- zoo::zoo(,order.by = total.time)
    combined.zoo <- merge.zoo(poll.zoo,na.zoo,all=TRUE)
    smoothed.zoo <- zoo::rollapply(combined.zoo,t.interval,mean,na.rm=TRUE,fill=NA,partial=TRUE)
    Poll.smooth <- coredata(smoothed.zoo)
    Poll.times.final <- index(smoothed.zoo)
    return(list(Poll.smooth,Poll.times.final))
  }
}
# Function which runs hidden markov models (HMMs) on separate days in time series
# Time must be POSIXct
# Bootstrap iterations describes number of time EM step of HMM fit runs to fit the data
# More bootstrap iterations will result in better fit but will come at computational cost
# Plot bool = 1 if returning plots of state 1 points, state 2 points by is desired
# Partitions points into two groups
# Classifies lower partition as background
# Upper partition as source
partitionPoints <- function(poll,time,lat,long,index,bootstrapIters,transformString){
  poll.days <- lubridate::mday(time)
  ## Extract and store indices of instances in which day changes in breaks
  day.breaks <- which(diff(poll.days)!=0)
  day.diffs <- diff(day.breaks)
  while(length(which(day.diffs<600))!=0){
    idx2takeout <- which(day.diffs<600)[1]
    poll <- poll[-((day.breaks[idx2takeout]+1):day.breaks[idx2takeout+1])]
    time <- time[-((day.breaks[idx2takeout]+1):day.breaks[idx2takeout+1])]
    long <- long[-((day.breaks[idx2takeout]+1):day.breaks[idx2takeout+1])]
    lat <- lat[-((day.breaks[idx2takeout]+1):day.breaks[idx2takeout+1])]
    index <- index[-((day.breaks[idx2takeout]+1):day.breaks[idx2takeout+1])]
    poll.days <- poll.days[-((day.breaks[idx2takeout]+1):day.breaks[idx2takeout+1])]
    day.breaks <- which(diff(poll.days)!=0)
    day.diffs <- diff(day.breaks)  
  }
  ## Preallocation
  poll.s1 <- vector(,)
  poll.s2 <- vector(,)
  s1.list <- list()
  state.vec <- vector(,)
  time.s1 <- list()
  time.s2 <- list()
  total.idx <- length(day.breaks)+1
  # Note: Seed value may inhibit convergence of certain segments of time series due to predetermined values. Consider changing
  # value for convergence issues.
  # Analyze densities of data and time distributions of data
  mod.list <- list()
  for (i in 1:total.idx){
    print(i)
    if (i==1){
      temp.poll <- poll[1:day.breaks[i]]
      temp.times <- time[1:day.breaks[i]]
    } else if (i==(length(day.breaks)+1)){
      temp.poll <- poll[(day.breaks[length(day.breaks)]+1):length(poll)]
      temp.times <- time[(day.breaks[length(day.breaks)]+1):length(time)]
    } else {
      temp.poll <- poll[(day.breaks[i-1]+1):day.breaks[i]]
      temp.times <- time[(day.breaks[i-1]+1):day.breaks[i]]
    }
  ## Fit the HMM to the day i of data
    if(transformString=="log"){
      temp.poll[temp.poll<0] <- 0
      temp.poll <- log(temp.poll+1)
    }
    temp.list <- list(temp.poll,temp.times)
    boot.list <- list()
    for(q in 1:bootstrapIters){boot.list[[q]] <- temp.list}
    print("Begin model fitting")
    mod.list <- lapply(boot.list, function(x) applyTimeDepmix(unlist(x[1]),unlist(x[2])))
    print("End model fitting")
    ll.storage <- numeric(bootstrapIters)
    for(q in 1:bootstrapIters){ll.storage[q] <- tryCatch(
      {
        ll.storage[q] <- logLik(mod.list[[q]])
      },
      error=function(cond){
        return(NA)
      }
    )
    }
    ll.storage <- round(ll.storage,2)
    selected.ll <- findMode(ll.storage)
    ## Extract data corresponding to states 1 and 2
    ## Choose model which returns highest log likelihood
    fm <- mod.list[[which(ll.storage==selected.ll)[1]]]
    resulting.states <- viterbi(fm)[,1]
    state.one <- resulting.states==1
    state.two <- resulting.states==2
    mirror.states <- resulting.states
    ## Let Poll.s1 represent the background and let Poll.s2 represent the signal
    ## Background is defined as whichever set of state points has the lower median
    if (length(temp.poll[state.one])==0){
      poll.s1 <- c(poll.s1,temp.poll[state.two])
      resulting.states[mirror.states==2] <- 1
      resulting.states[mirror.states==1] <- 2
      poll.s2 <- c(poll.s2,temp.poll[state.one])
      time.s1 <- append(time.s1,temp.times[state.two])
      time.s2 <- append(time.s2,temp.times[state.one])
    } else if (length(temp.poll[state.two])==0){
      poll.s1 <- c(poll.s1,temp.poll[state.one])
      poll.s2 <- c(poll.s2,temp.poll[state.two])
      time.s1 <- append(time.s1,temp.times[state.one])
      time.s2 <- append(time.s2,temp.times[state.two])
    } else {
      if (median(temp.poll[state.one],na.rm = TRUE)>median(temp.poll[state.two],na.rm = TRUE)){
        poll.s1 <- c(poll.s1,temp.poll[state.two])
        resulting.states[mirror.states==2] <- 1
        resulting.states[mirror.states==1] <- 2
        poll.s2 <- c(poll.s2,temp.poll[state.one])
        time.s1 <- append(time.s1,temp.times[state.two])
        time.s2 <- append(time.s2,temp.times[state.one])
      } else {
        poll.s1 <- c(poll.s1,temp.poll[state.one])
        poll.s2 <- c(poll.s2,temp.poll[state.two])
        time.s1 <- append(time.s1,temp.times[state.one])
        time.s2 <- append(time.s2,temp.times[state.two])
      }
    }
    state.vec <- c(state.vec,resulting.states)
    s1.list[[i]] <- temp.poll[resulting.states==1]
  }
  final.results <- list(poll,time,lat,long,index,state.vec,poll.s1,time.s1,poll.s2,time.s2,s1.list)
  return(final.results)
}
all.signals <- list()
s1.runs <- list()
n.runs <- 5
## END DEFINITIONS##
# Number of iterations to find global maximum in expectation maximization of depmix likelihood
bootstrap.iters <- 100
# Read in data
## CHANGE FILENAME DEPENDING ON POLLUTANT ##
Poll.times.c1 <- parse_date_time(as.character(unlist(read.csv("C:/Path/LST_Car1NOxFinal.csv", header = FALSE),use.names=FALSE)), orders=c("mdy HMS"), tz=Sys.timezone())
Long.c1 <- as.vector(unlist(read.csv("C:/Path/NOx_Long_Raw_Car1.csv", header = FALSE),use.names=FALSE)) 
Lat.c1 <- as.vector(unlist(read.csv("C:/Path/NOx_Lat_Raw_Car1.csv", header = FALSE),use.names=FALSE))
Index.c1 <- rep(1,length(Long.c1))
Poll.measurements.c1 <- as.vector(unlist(read.csv("C:/Path/NOx_Raw_Car1.csv", header = FALSE),use.names=FALSE)) 
## Make corrections to data ##
## Perform initial NaN purge
idxNA <- (is.na(Poll.times.c1) | is.na(Long.c1) | is.nan(Lat.c1) | is.na(Index.c1) | is.na(Poll.measurements.c1))
Poll.times.c1 <- Poll.times.c1[!idxNA]
Poll.measurements.c1 <- Poll.measurements.c1[!idxNA]
Long.c1 <- Long.c1[!idxNA]
Lat.c1 <- Lat.c1[!idxNA]
Index.c1 <- Index.c1[!idxNA]
## Smooth the data
smoothed.list.c1 <- smoothData(Poll.measurements.c1,Poll.times.c1,30)
Poll.smooth.c1 <- smoothed.list.c1[[1]]
Poll.times.smooth.c1 <- smoothed.list.c1[[2]]
## Take out NaNs and inconsistent times.
idxNA <- (is.na(Poll.times.smooth.c1) | is.na(Poll.smooth.c1))
Poll.times.smooth.c1 <- Poll.times.smooth.c1[!idxNA]
Poll.smooth.c1 <- Poll.smooth.c1[!idxNA]
idxMatch <- Poll.times.smooth.c1 %in% Poll.times.c1
Poll.times.smooth.c1 <- Poll.times.smooth.c1[idxMatch]
Poll.smooth.c1 <- Poll.smooth.c1[idxMatch]
## Define IdxHour if want to look at particular subset of data
idxHour <- hour(Poll.times.c1)>7 & hour(Poll.times.c1)<19
Poll.times.smooth.c1 <- Poll.times.smooth.c1[idxHour]
Poll.smooth.c1 <- Poll.smooth.c1[idxHour]
Long.c1 <- Long.c1[idxHour]
Lat.c1 <- Lat.c1[idxHour]
Index.c1 <- Index.c1[idxHour]
Poll.times.c2 <- parse_date_time(as.character(unlist(read.csv("C:/Path/LST_Car2NOxFinal.csv", header = FALSE),use.names=FALSE)), orders=c("mdy HMS"), tz=Sys.timezone())
Long.c2 <- as.vector(unlist(read.csv("C:/Path/NOx_Long_Raw_Car2.csv", header = FALSE),use.names=FALSE)) 
Lat.c2 <- as.vector(unlist(read.csv("C:/Path/NOx_Lat_Raw_Car2.csv", header = FALSE),use.names=FALSE)) 
Index.c2 <- rep(2,length(Long.c2))
Poll.measurements.c2 <- as.vector(unlist(read.csv("C:/Path/NOx_Raw_Car2.csv", header = FALSE),use.names=FALSE)) 
## Make Corrections to Data ##
## Perform initial nan purge
idxNA <- (is.na(Poll.times.c2) | is.na(Long.c2) | is.nan(Lat.c2) | is.na(Index.c2) | is.na(Poll.measurements.c2))
Poll.times.c2 <- Poll.times.c2[!idxNA]
Poll.measurements.c2 <- Poll.measurements.c2[!idxNA]
Long.c2 <- Long.c2[!idxNA]
Lat.c2 <- Lat.c2[!idxNA]
Index.c2 <- Index.c2[!idxNA]
## Smooth the data
smoothed.list.c2 <- smoothData(Poll.measurements.c2,Poll.times.c2,30)
Poll.smooth.c2 <- smoothed.list.c2[[1]]
Poll.times.smooth.c2 <- smoothed.list.c2[[2]]
## Take out NaNs and inconsistent times.
idxNA <- (is.na(Poll.times.smooth.c2) | is.na(Poll.smooth.c2))
Poll.times.smooth.c2 <- Poll.times.smooth.c2[!idxNA]
Poll.smooth.c2 <- Poll.smooth.c2[!idxNA]
idxMatch <- Poll.times.smooth.c2 %in% Poll.times.c2
Poll.times.smooth.c2 <- Poll.times.smooth.c2[idxMatch]
Poll.smooth.c2 <- Poll.smooth.c2[idxMatch]
## Define idxHour
idxHour <- hour(Poll.times.c2)>4
Poll.times.smooth.c2 <- Poll.times.smooth.c2[idxHour]
Poll.smooth.c2 <- Poll.smooth.c2[idxHour]
Long.c2 <- Long.c2[idxHour]
Lat.c2 <- Lat.c2[idxHour]
Index.c2 <- Index.c2[idxHour]

for (b in 1:n.runs){
c1.res <- partitionPoints(Poll.smooth.c1,Poll.times.smooth.c1,Lat.c1,Long.c1,Index.c1,bootstrap.iters,"log")
c2.res <- partitionPoints(Poll.smooth.c2,Poll.times.smooth.c2,Lat.c2,Long.c2,Index.c2,bootstrap.iters,"log")
## Save state designations
total.state <- c(c1.res[[6]],c2.res[[6]])
total.time <- c(c1.res[[2]],c2.res[[2]])
total.poll <- c(c1.res[[1]],c2.res[[1]])
total.lat <- c(c1.res[[3]],c2.res[[3]])
total.long <- c(c1.res[[4]],c2.res[[4]])
total.index <- c(c1.res[[5]],c2.res[[5]])
total.s1.time <- c(c1.res[[8]],c2.res[[8]])
total.s1.poll <- c(c1.res[[7]],c2.res[[7]])
# total.s2.time <- c(c1.res[[10]],c2.res[[10]])
# total.s2.poll <- c(c1.res[[9]],c2.res[[9]])
s1.Dayseconds <- secondsFromDayStartConversion(total.s1.time) 
s1.yeardays <- daysFromYearStartConversion(total.s1.time)
# s2.Dayseconds <- secondsFromDayStartConversion(total.s2.time)
# s2.yeardays <- daysFromYearStartConversion(total.s2.time)
HMM.dataframe <- data.frame("Poll"= total.s1.poll,"Dayseconds"= s1.Dayseconds,"Yeardays"= s1.yeardays)
rowsToSample <- round(1*nrow(HMM.dataframe))
HMMPtsToFit <- HMM.dataframe[sample(nrow(HMM.dataframe),size = rowsToSample,replace = FALSE),]  
## Create overall data frame
overall.Dayseconds <- secondsFromDayStartConversion(total.time)
overall.yeardays <- daysFromYearStartConversion(total.time)
overall.matrix <- matrix(c(overall.Dayseconds,overall.yeardays),nrow=length(overall.Dayseconds),ncol=2)
## Fit spline to HMM data
print("Begin Smoothing Process")

HMM.smooth <- bam(Poll~s(Dayseconds,Yeardays,k=900,bs='tp'),method="fREML",data = HMMPtsToFit,family=gaussian(link="identity"),discrete=TRUE)
HMM.signal <- predict.bam(HMM.smooth,newdata = data.frame(Dayseconds=overall.matrix[,1],Yeardays=overall.matrix[,2]))
c1.HMM.signal <- HMM.signal[1:length(Poll.measurements.c1)]
c2.HMM.signal <- HMM.signal[(length(Poll.measurements.c1)+1):length(HMM.signal)]
## Save background signal file
# write.csv(HMM.signal,"C:/Path/Total_HMMBackground_NOx_Log_Test2.csv")
# end_time <- Sys.time()
# total_time <- end_time - start_time
# total_time
all.signals[[b]] <- HMM.signal
s1.runs[[b]] <- c1.res[[11]]
}
## Calculate RMSE values between all list entries.
rmse.res <- matrix(,nrow=n.runs,ncol=n.runs)
for (p in 1:n.runs){
  for (q in 1:n.runs){
    rmse.res[p,q] <- Metrics::rmse(all.signals[[p]],all.signals[[q]])
  }
}
diffs.vec <- vector(,)
for (k in 1:(length(s1.runs)-1)){
  for (j in (k+1):length(s1.runs)){
    for (i in 1:length(s1.runs[[j]])){
      if(!identical(s1.runs[[k]][[i]],s1.runs[[j]][[i]])){
        diffs.vec <- c(diffs.vec,i)
      }
    }
  }
}
