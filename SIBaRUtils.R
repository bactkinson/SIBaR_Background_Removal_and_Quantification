require(zoo)
require(lubridate)
## Utility functions to make processing data through SIBaR easier. 
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

## Rewrite to remove NAs, match up timestamps inside function call  (10/25)
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
