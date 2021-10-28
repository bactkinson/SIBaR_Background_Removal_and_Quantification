rm(list=ls())
require(lubridate)
require(mgcv)
require(dplyr)
require(magrittr)
require(splines)
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

state <- read.csv("C:/Path/Total_NOx_State_Corrected.csv")[,2]
poll <- read.csv("C:/Path/Total_NOx_Poll_Corrected.csv")[,2]
timestamps <- parse_date_time(as.character(read.csv("C:/Path/Total_NOx_Time_Corrected.csv")[,2]),orders = c("ymd HMS"),tz="US/Central")
lat <- read.csv("C:/Path/Total_NOx_Lat_Corrected.csv")[,2]
long <- read.csv("C:/Path/Total_NOx_Long_Corrected.csv")[,2]
index <- read.csv("C:/Path/Total_NOx_Index_Corrected.csv")[,2]

total.s1.time <- timestamps[state==1]
total.s1.poll <- poll[state==1]
total.s1.index <- index[state==1]
## Convert LST seconds and days to numerical seconds and days
s1.Dayseconds <- secondsFromDayStartConversion(total.s1.time)
s1.yeardays <- daysFromYearStartConversion(total.s1.time)
s1.dataframe <- tibble("Poll"=total.s1.poll,"Dayseconds"=s1.Dayseconds,"Yeardays"=s1.yeardays,"Month"=month(total.s1.time),"Day"=mday(total.s1.time),"Timestamps"=total.s1.time,"Index"=total.s1.index) %>%
  group_by(Month,Day,Index)
## Create overall data frame for prediction in predict.gam step
overall.dayseconds <- secondsFromDayStartConversion(timestamps)
overall.yeardays <- daysFromYearStartConversion(timestamps)
overall.dataframe <- tibble("Dayseconds"=overall.dayseconds,"Yeardays"=overall.yeardays,"Month"=month(timestamps),"Day"=mday(timestamps),"Index"=index) %>%
  group_by(Month,Day,Index)
  
## Fit spline to HMM data
background.signal <- vector(,)
unique.month.days <- overall.dataframe %>%
  group_indices() %>%
  unique()

for (i in 1:length(unique.month.days)){
  temp.s1.dat <- s1.dataframe %>% filter(cur_group_id()==unique.month.days[i])
  temp.overall.dat <- overall.dataframe %>% filter(cur_group_id()==unique.month.days[i])
  temp.smooth <- lm(Poll~ns(Dayseconds,df=length(unique(hour(Timestamps)))),data=temp.s1.dat)
  temp.signal <- predict(temp.smooth,newdata=temp.overall.dat)
  background.signal <- c(background.signal,as.numeric(temp.signal))
#   # print("Begin Fitting")
#   # temp.smooth <- gam(Poll~te(Dayseconds,Yeardays,k=15,bs='tp'),method="REML",discrete=FALSE,data = temp.s1.dat,family=gaussian(link="identity"))
#   # print("End Fitting")
#   # temp.signal <- predict.gam(temp.smooth,newdata=temp.overall.dat)
#   # background.signal <- c(background.signal,as.numeric(temp.signal))
}

## Save background signal file
write.csv(background.signal,"C:/Path/Total_HMMBackground_NOx_NaturalSpline_ByDay_Corrected.csv")
overall.smooth <- gam(Poll~te(Dayseconds,Yeardays,k=5,bs='tp'),method="REML",data = s1.dataframe,family=gaussian())
overall.signal <- predict.gam(overall.smooth,newdata=overall.dataframe)
write.csv(overall.signal,"C:/Path/Total_HMMBackground_NOx_Corrected.csv")

