# rm(list=ls())
require(lubridate)
require(zoo)
require(dplyr)
require(tidyverse)

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

rollingWindowClassifier <- function(state,poll,timestamps,index,tw,threshold,dir,saveGraphsBool){
  ## Partition the data
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
    start.time <- temp.times[1]
    end.time <- temp.times[length(temp.times)]
    time.seq <- seq.POSIXt(start.time,end.time,by="sec")
    na.zoo <- zoo(NA,order.by=time.seq)
    temp.zoo <- zoo(matrix(c(temp.poll,temp.state),ncol=2),order.by=temp.times)
    working.zoo <- merge.zoo(na.zoo,temp.zoo)
    working.zoo <- working.zoo[,-1]
    time.window <- tw
    output.bool <- logical(length(time.seq)-time.window+1)
    ## Now to apply stateEvaluation function
    for(j in 1:(length(time.seq)-time.window+1)){
      window <- j:(j+time.window-1)
      windowed.poll <- working.zoo[window,1]
      windowed.state <- working.zoo[window,2]
      output.bool[j] <- stateEvaluation(windowed.poll,windowed.state)
    }
    output.bool <- output.bool[!is.na(output.bool)]
    if((length(which(!output.bool))/length(output.bool)*100) > threshold){
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
      legend("topright", inset=c(0,0), y.intersp = 1, legend = c("Background", "Source"),  lty = 1, bty = "n", col = c(1,2), cex = 1.2)
      dev.off()
    }
  }
  pct.correct <- length(which(classifier))/length(classifier)*100
  print(pct.correct)
  return(pct.correct)
}

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

# state <- read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_State_Revised_Final.csv")[,2]
# 
# poll <- read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Poll_Revised_Final.csv")[,2]
# 
# timestamps <- parse_date_time(as.character(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Time_Revised_Final.csv")[,2]),orders = c("ymd HMS"),tz="US/Central")
# 
# index <- read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Index_Revised_Final.csv")[,2]

# dir <- 'C:/Users/bwa2/Documents/Academic Work/Research/SIBaR Background Removal/SIBaR State Time Series/Daily NOx HMM Time Series Rolling Window Classifer 25/'
# 
# rollingWindowClassifier(state,poll,timestamps,index,900,25,dir)

# dir <- 'C:/Users/bwa2/Documents/Academic Work/Research/SIBaR Background Removal/SIBaR State Time Series/Daily NOx HMM Time Series Fitted Line Classifier 50/'

# output <- fittedLineClassifier(state,poll,timestamps,index,50,dir,F)

# sens.vec <- seq(5,100,5)
# pct.correct <- numeric(length(sens.vec))
# for(i in 1:length(sens.vec)){
#   pct.correct[i] <- fittedLineClassifier(state,poll,timestamps,index,sens.vec[i],dir,F)
# }
# df <- tibble("Threshold_Percent"=sens.vec,"Correct_Class_Percent"=pct.correct)
# ggplot(data=df,aes(Threshold_Percent,Correct_Class_Percent)) + 
#   geom_line(size=1) + 
#   geom_point(size=3,color="red") +
#   labs(title="Sensitivity Analysis for the Fitted Line Classifier",
#        x = "Threshold Percentage (%)",
#        y = "Percentage Classified as Correct (%)") +
#   theme_classic() +
#   expand_limits(y=0)
