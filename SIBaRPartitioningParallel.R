## Function that performs the partitioning step in SIBaR.
require(parallel)
require(depmixS4)
require(lubridate)
require(zoo)
require(dplyr)
require(tidyverse)


## Wrapper function that performs the SIBaR partitioning routine
## Function already assumes that user wants to correct misclassified series
## subject to a threshold of 50%. 

partitionRoutine <- function(poll,time,bootstrap_iterations,transform_string="none",
                             minTimePts = 600, cores=parallel::detectCores()-2,
                             index=rep(1,length(poll)), threshold = 50,
                             length_tolerance = 0.05){

  initialPartition <- partitionPoints(poll,time,bootstrap_iterations)
  
  ## Initiate the correction routine
  initial_evaluation <- fittedLineClassifier(initialPartition$State,initialPartition$Poll,
                                             initialPartition$Timestamps,initialPartition$Index,
                                             threshold = 50)
  
  ## If there are splits deemed misclassified by the given threshold, recursively
  ## correct using partitionCorrection
  
  
}

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

partitionPoints <- function(poll,time,bootstrapIters,transformString="none",minTimePts=600,
                            cores=parallel::detectCores()-1,
                            index=rep(1,length(poll))){

  idxValid <- validDayIdxReturn(time,minTimePts)
  poll <- poll[idxValid]
  time <- time[idxValid]
  index <- index[idxValid]
  
  ## Preallocation
  poll.s1 <- vector(,)
  poll.s2 <- vector(,)
  state.vec <- vector(,)
  time.s1 <- vector(,)
  time.s2 <- vector(,)
  mod.list <- list()
  ll.storage <- numeric(bootstrapIters)
  
  cl <- makeCluster(cores)
  clusterExport(cl, varlist = c("applyTimeDepmix"), envir = .GlobalEnv)
  clusterExport(cl, list("bootstrapIters"), envir = environment())
  
  if(transformString=="log"){
    poll[poll<0] <- 0
    poll <- log1p(poll)
  }
  
  data_splits <- tibble("Poll"=poll,"Time"=time,"Index"=index,"Month"=month(time),"Day"=mday(time))%>%
    group_split(Month,Day,Index,.keep = "unused")
  
  for (i in 1:length(data_splits)){
    print(i)
    temp.poll <- data_splits[[i]]$Poll
    temp.times <- data_splits[[i]]$Time
  
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
  final.results <- tibble("Timestamps"=time,"Poll"=poll,"States"=state.vec,"Index"=index)
  
  # final.results <- list(poll,time,state.vec,poll.s1,time.s1,poll.s2,time.s2)
  return(final.results)
}

stateEvaluation <- function(meas,state)
{
  if(length(which(is.na(meas)))>(length(meas)/2)){
    return(NA)
  }
  if((length(meas[state==1])==0) | (length(meas[state==2])==0)){
    return(TRUE)
  }
  if(mean(meas[state==1],na.rm=T)<mean(meas[state==2],na.rm=T)){
    return(TRUE)
  } else{ 
    return(FALSE)
  }
}
## 10/25 Rewrite index to be optional factor passed in. Determine suitable 
## default behavior
fittedLineClassifier <- function(state,poll,timestamps,index,threshold,dir,saveGraphsBool){
  data.splits <- tibble("state"=state,"pollutant"=poll,"timestamps"=timestamps,"Index"=index) %>%
    cbind("Month"=month(timestamps)) %>%
    cbind("Day"=mday(timestamps)) %>%
    group_split(Month,Day,Index,.keep=FALSE)
  
  classifier <- logical(length(data.splits))
  for(i in 1:length(data.splits)){
    temp.data <- data.splits[[i]]
    temp.poll <- temp.data$pollutant
    temp.state <- temp.data$state
    temp.times <- temp.data$timestamps
    lagged.state <- lag(temp.state)
    transitions <- temp.state != lagged.state
    trans.times <- temp.times[transitions]
    trans.indices <- which(transitions)
    trans.poll <- numeric(length(trans.times))
    for(j in 1:length(trans.indices)){
      current.index <- trans.indices[j]
      trans.poll[j] <- mean(c(temp.poll[current.index-1],temp.poll[current.index],temp.poll[current.index+1]),na.rm=T)
    }
    temp.df <- data.frame("Time"=trans.times,"Measurement"=trans.poll)
    trans.line.fit <- lm(Measurement~Time,data=temp.df)
    trans.line.preds <- predict(trans.line.fit,newdata=data.frame("Time"=temp.times))
    states.below <- temp.state[temp.poll <= trans.line.preds]
    states.above <- temp.state[temp.poll > trans.line.preds]
    pct.above.misclass <- length(which(states.above==1))/length(states.above)*100
    pct.below.misclass <- length(which(states.below==2))/length(states.below)*100
    if(pct.above.misclass >= threshold | pct.below.misclass >= threshold){
      titlestring <- "Misclassified"
      classifier[i] <- F
    } else{
      titlestring <- "Classified correctly"
      classifier[i] <- T
    }
    if(saveGraphsBool)
    {
      png(filename = paste0(dir,'Day ',i, '.png'),
          width = 480, height = 480, units = "px", pointsize = 12,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,
          type = c("windows"))
      plot(temp.poll~temp.times,col=temp.state,ylab="ln(NOx+1)",xlab="Time",main=titlestring)
      abline(trans.line.fit,col="green",lty=1,lwd=2)
      legend("topright", inset=c(0,0), y.intersp = 1, legend = c("Background", "Source"),  lty = 1, bty = "n", col = c(1,2), cex = 1.2)
      dev.off()  
    }
    
  }
  pct.correct <- length(which(classifier))/length(classifier)*100
  print(pct.correct)
  return(list(classifier,pct.correct))
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

