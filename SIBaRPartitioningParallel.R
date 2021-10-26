## Function that performs the partitioning step in SIBaR.
require(parallel)
require(depmixS4)

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