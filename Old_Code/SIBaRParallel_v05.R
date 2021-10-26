rm(list=ls())
require(lubridate)
require(tidyverse)
require(mgcv)
require(depmixS4)
require(zoo)
require(parallel)
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
recursiveCorrections <- function(misclassed.split,tol){
  old.series <- list(list(misclassed.split$pollutant,misclassed.split$timestamps,misclassed.split$state))
  classified.res <- F
  length.tolerance <- tol*length(misclassed.split$pollutant)
  while(any(!classified.res)){
    # Continue to perform partition corrections
    new.series <- list()
    idx <- 1
    for(j in 1:length(old.series)){
      if(classified.res[j]==F & (length(old.series[[j]][[1]]) > length.tolerance)){
        output <- partitionCorrection(old.series[[j]][[1]],old.series[[j]][[2]],50)
        new.series[[idx]] <- list(output[[1]],output[[2]],output[[3]])
        new.series[[idx+1]] <- list(output[[4]],output[[5]],output[[6]])
        idx <- idx+2
      } else{
        new.series[[idx]] <- list(old.series[[j]][[1]],old.series[[j]][[2]],old.series[[j]][[3]])
        idx <- idx+1
      }
    }
    classified.res <- logical()
    final.state <- vector(,)
    final.poll <- vector(,)
    final.time <- vector(,)
    for(k in 1:length(new.series)){
      idx <- rep(1,length(new.series[[k]][[1]]))
      if(length(new.series[[k]][[1]]) > length.tolerance){ 
        temp <- fittedLineClassifier(new.series[[k]][[3]],
                                     new.series[[k]][[1]],
                                     new.series[[k]][[2]],
                                     idx,
                                     50,
                                     dir="Blah",
                                     F)
        classified.res[k] <- temp[[1]]
      } else {
        classified.res[k] <- T
      }
    }
    old.series <- new.series
  }
  finalized.states <- vector(,)
  for(i in 1:length(old.series)){
    finalized.states <- c(finalized.states,old.series[[i]][[3]])
  }
  return(finalized.states)
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
# Function which smooths irregularly spaced time series
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
# Function which sets observations taken on days in which total observations < counts to NA.
validDayIdxReturn <- function(time,counts){
  idxs <- vector(,)
  ms <- month(time)
  ds <- day(time)
  df <- data.frame("month"=ms,"day"=ds)
  unique.entries <- unique(df)
  for (i in 1:nrow(unique.entries)){
    matching.idxs <- which((ms==unique.entries[i,1]) & (ds==unique.entries[i,2]))
    if (length(matching.idxs)>=counts){
      idxs <- c(idxs,matching.idxs)
    }
  }
  return(idxs)
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
  idxValid <- validDayIdxReturn(time,600)
  poll <- poll[idxValid]
  time <- time[idxValid]
  lat <- lat[idxValid]
  long <- long[idxValid]
  index <- index[idxValid]
  poll.days <- mday(time)
  day.breaks <- which(diff(poll.days)!=0)
  ## Preallocation
  poll.s1 <- vector(,)
  poll.s2 <- vector(,)
  state.vec <- vector(,)
  time.s1 <- vector(,)
  time.s2 <- vector(,)
  total.idx <- length(day.breaks)+1
  # value for convergence issues.
  # Analyze densities of data and time distributions of data
  mod.list <- list()
  ll.storage <- numeric(bootstrapIters)
  no.cores <- detectCores()-1
  cl <- makeCluster(no.cores)
  clusterExport(cl, varlist = c("applyTimeDepmix"), envir = .GlobalEnv)
  clusterExport(cl, list("bootstrapIters"), envir = environment())
  if(transformString=="log"){
    poll[poll<0] <- 0
    poll <- log(poll+1)
  }
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
    temp.list <- list(temp.poll,temp.times)
    boot.list <- list()
    for(q in 1:bootstrapIters){boot.list[[q]] <- temp.list}
    print("Begin model fitting")
    clusterExport(cl, list("boot.list"), envir = environment())
    mod.list <- parLapply(cl, boot.list, function(x) applyTimeDepmix(unlist(x[1]),unlist(x[2])))
    print("End model fitting")
    ll.storage <- numeric(bootstrapIters)
    for(q in 1:bootstrapIters){ll.storage[q] <- try({logLik(mod.list[[q]])})}
    ## Extract data corresponding to states 1 and 2
    ## Choose model which returns highest log likelihood
    fm <- mod.list[[which.max(ll.storage)[1]]]
    resulting.states <- viterbi(fm)[,1]
    state.one <- resulting.states==1
    state.two <- resulting.states==2
    mirror.states <- resulting.states
    ## Let Poll.s1 represent the background and let Poll.s2 represent the signal
    ## Background is defined as whichever set of state points has the lower median
    if(length(temp.list[[1]][state.one])==0){
      resulting.states[mirror.states==2] <- 1
      resulting.states[mirror.states==1] <- 2
    } 
    if(median(temp.list[[1]][state.one],na.rm = TRUE)>median(temp.list[[1]][state.two],na.rm = TRUE)){
      resulting.states[mirror.states==2] <- 1
      resulting.states[mirror.states==1] <- 2
    } 
    poll.s1 <- c(poll.s1,temp.poll[resulting.states==1])
    poll.s2 <- c(poll.s2,temp.poll[resulting.states==2])
    time.s1 <- c(time.s1,temp.times[resulting.states==1])
    time.s2 <- c(time.s2,temp.times[resulting.states==2])
    state.vec <- c(state.vec,resulting.states)
  }
  stopCluster(cl)
  time.s1 <- as.POSIXct(time.s1,tz="US/Central",origin="1970-01-01")
  time.s2 <- as.POSIXct(time.s2,tz="US/Central",origin="1970-01-01")
  final.results <- list(poll,time,lat,long,index,state.vec,poll.s1,time.s1,poll.s2,time.s2)
  return(final.results)
}

partitionCorrection <- function(poll,time,bootstrapIters){
  halfway.point <- round(length(poll)/2,0)
  poll.a <- poll[1:halfway.point]
  poll.b <- poll[(halfway.point+1):length(poll)]
  time.a <- time[1:halfway.point]
  time.b <- time[(halfway.point+1):length(poll)]
  mod.list.a <- list()
  mod.list.b <- list()
  boot.list.a <- list()
  boot.list.b <- list()
  temp.list.a <- list(poll.a,time.a)
  temp.list.b <- list(poll.b,time.b)
  for(q in 1:bootstrapIters){
    boot.list.a[[q]] <- temp.list.a
    boot.list.b[[q]] <- temp.list.b
  }
  print("Begin model refitting")
  mod.list.a <- lapply(boot.list.a, function(x) applyTimeDepmix(unlist(x[1]),unlist(x[2])))
  mod.list.b <- lapply(boot.list.b, function(x) applyTimeDepmix(unlist(x[1]),unlist(x[2])))
  print("End model refitting")
  ll.storage.a <- numeric(bootstrapIters)
  ll.storage.b <- numeric(bootstrapIters)
  for(q in 1:bootstrapIters){
    ll.storage.a[q] <- try({logLik(mod.list.a[[q]])})
    ll.storage.b[q] <- try({logLik(mod.list.b[[q]])})
  }
  ## Extract data corresponding to states 1 and 2
  ## Choose model which returns highest log likelihood
  fm.a <- mod.list.a[[which.max(ll.storage.a)[1]]]
  fm.b <- mod.list.b[[which.max(ll.storage.b)[1]]]
  states.a <- viterbi(fm.a)[,1]
  state.one.a <- states.a==1
  state.two.a <- states.a==2
  mirror.states.a <- states.a
  states.b <- viterbi(fm.b)[,1]
  state.one.b <- states.b==1
  state.two.b <- states.b==2
  mirror.states.b <- states.b
  
  if(length(temp.list.a[[1]][state.one.a])==0){
    states.a[mirror.states.a==2] <- 1
    states.a[mirror.states.a==1] <- 2
  } 
  if(median(temp.list.a[[1]][state.one.a],na.rm = TRUE)>median(temp.list.a[[1]][state.two.a],na.rm = TRUE)){
    states.a[mirror.states.a==2] <- 1
    states.a[mirror.states.a==1] <- 2
  }
  if(length(temp.list.b[[1]][state.one.b])==0){
    states.b[mirror.states.b==2] <- 1
    states.b[mirror.states.b==1] <- 2
  } 
  if(median(temp.list.b[[1]][state.one.b],na.rm = TRUE)>median(temp.list.b[[1]][state.two.b],na.rm = TRUE)){
    states.b[mirror.states.b==2] <- 1
    states.b[mirror.states.b==1] <- 2
  }
  total.new.states <- c(states.a,states.b)
  return(list(poll.a,time.a,states.a,poll.b,time.b,states.b,total.new.states))
}
## END DEFINITIONS##
# Number of iterations to find global maximum in expectation maximization of depmix likelihood
bootstrap.iters <- 3
# Read in data
## CHANGE FILENAME DEPENDING ON POLLUTANT ##
Poll.times.c1 <- parse_date_time(as.character(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/LST_Car1_NOx_Revised_Final.csv", header = FALSE),use.names=FALSE)), orders=c("mdy HMS"), tz=Sys.timezone())
Long.c1 <- as.vector(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/NOx_Long_Raw_Car1_Revised_Final.csv", header = FALSE),use.names=FALSE))
Lat.c1 <- as.vector(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/NOx_Lat_Raw_Car1_Revised_Final.csv", header = FALSE),use.names=FALSE))
Index.c1 <- rep(1,length(Long.c1))
Poll.measurements.c1 <- as.vector(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/NOx_Raw_Car1_Revised_Final.csv", header = FALSE),use.names=FALSE))
## Make corrections to data ##
## Perform initial NaN purge
idxNA <- (is.na(Poll.times.c1) | is.na(Long.c1) | is.nan(Lat.c1) | is.na(Index.c1) | is.na(Poll.measurements.c1))
Poll.times.c1 <- Poll.times.c1[!idxNA]
Poll.measurements.c1 <- Poll.measurements.c1[!idxNA]
Long.c1 <- Long.c1[!idxNA]
Lat.c1 <- Lat.c1[!idxNA]
Index.c1 <- Index.c1[!idxNA]
## Smooth the data
smoothed.list.c1 <- smoothData(Poll.measurements.c1,Poll.times.c1,0)
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
idxHour <- hour(Poll.times.c1)>4
Poll.times.smooth.c1 <- Poll.times.smooth.c1[idxHour]
Poll.smooth.c1 <- Poll.smooth.c1[idxHour]
Long.c1 <- Long.c1[idxHour]
Lat.c1 <- Lat.c1[idxHour]
Index.c1 <- Index.c1[idxHour]
## Do the same process for Car 2 ##
## CHANGE FILENAME DEPENDING ON POLLUTANT ##
Poll.times.c2 <- parse_date_time(as.character(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/LST_Car2_NOx_Revised_Final.csv", header = FALSE),use.names=FALSE)), orders=c("mdy HMS"), tz=Sys.timezone())
Long.c2 <- as.vector(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/NOx_Long_Raw_Car2_Revised_Final.csv", header = FALSE),use.names=FALSE))
Lat.c2 <- as.vector(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/NOx_Lat_Raw_Car2_Revised_Final.csv", header = FALSE),use.names=FALSE))
Index.c2 <- rep(2,length(Long.c2))
Poll.measurements.c2 <- as.vector(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/NOx_Raw_Car2_Revised_Final.csv", header = FALSE),use.names=FALSE))
# ## Make Corrections to Data ##
# ## Perform initial nan purge
idxNA <- (is.na(Poll.times.c2) | is.na(Long.c2) | is.nan(Lat.c2) | is.na(Index.c2) | is.na(Poll.measurements.c2))
Poll.times.c2 <- Poll.times.c2[!idxNA]
Poll.measurements.c2 <- Poll.measurements.c2[!idxNA]
Long.c2 <- Long.c2[!idxNA]
Lat.c2 <- Lat.c2[!idxNA]
Index.c2 <- Index.c2[!idxNA]
## Smooth the data
smoothed.list.c2 <- smoothData(Poll.measurements.c2,Poll.times.c2,0)
Poll.smooth.c2 <- smoothed.list.c2[[1]]
Poll.times.smooth.c2 <- smoothed.list.c2[[2]]
# ## Take out NaNs and inconsistent times.
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
## Run the partitioning
c1.res <- partitionPoints(Poll.smooth.c1,Poll.times.smooth.c1,Lat.c1,Long.c1,Index.c1,bootstrap.iters,"log")
c2.res <- partitionPoints(Poll.smooth.c2,Poll.times.smooth.c2,Lat.c2,Long.c2,Index.c2,bootstrap.iters,"log")
# ## Save parameter designations
total.state <- c(c1.res[[6]],c2.res[[6]])
total.time <- c(c1.res[[2]],c2.res[[2]])
total.poll <- c(c1.res[[1]],c2.res[[1]])
total.lat <- c(c1.res[[3]],c2.res[[3]])
total.long <- c(c1.res[[4]],c2.res[[4]])
total.index <- c(c1.res[[5]],c2.res[[5]])
## Evaluate misclassification
source('C:/Users/bwa2/Documents/Academic Work/Research/SIBaR Background Removal/SIBaR Scripts/SIBaR Classification Error.R')
misclass.eval <- fittedLineClassifier(total.state,total.poll,total.time,total.index,50,dir="Blah",F)
classified.res <- !misclass.eval[[1]]
# classified.res <- c(TRUE,FALSE,FALSE,TRUE,FALSE)
if(any(classified.res)){
  test.splits <- tibble("state"=total.state,"pollutant"=total.poll,"timestamps"=total.time,"Index"=total.index,"Vector_Loc"=seq(1,length(total.state),1)) %>%
    cbind("Month"=month(total.time)) %>%
    cbind("Day"=mday(total.time)) %>%
    group_split(Month,Day,Index,.keep=FALSE)
  
  misclassed.splits <- tibble("state"=total.state,"pollutant"=total.poll,"timestamps"=total.time,"Index"=total.index,"Vector_Loc"=seq(1,length(total.state),1)) %>%
    cbind("Month"=month(total.time)) %>%
    cbind("Day"=mday(total.time)) %>%
    group_split(Month,Day,Index,.keep=FALSE) %>%
    .[classified.res]
  
  for(n in 1:length(misclassed.splits)){
    print("")
    print(n)
    print("")
    corrected.states <- recursiveCorrections(misclassed.splits[[n]],0.05)
    total.state[misclassed.splits[[n]]$Vector_Loc] <- corrected.states
  }
}
write.csv(total.state,"C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_State_Test.csv")
write.csv(total.poll,"C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Poll_Test.csv")
write.csv(total.time,"C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Time_Test.csv")
write.csv(total.lat,"C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Lat_Test.csv")
write.csv(total.long,"C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Long_Test.csv")
write.csv(total.index,"C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Index_test.csv")
#
