require(zoo)
require(lubridate)
## Utility functions to make processing data through SIBaR easier. 
## These functions convert POSIXct timestamps to numerical days from year start
## and seconds from start of day.
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

## Function which smooths irregularly spaced time series given a prescribed
## time interval.
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
    idxNA <- (is.na(Poll.smooth) | is.na(Poll.times.final))
    Poll.times.final <- Poll.times.final[!idxNA]
    Poll.smooth <- Poll.smooth[!idxNA]
    idxMatch <- Poll.times.final %in% time
    Poll.times.final <- Poll.times.final[idxMatch]
    Poll.smooth <- Poll.smooth[idxMatch]
    return(list(Poll.smooth,Poll.times.final))
  }
}
