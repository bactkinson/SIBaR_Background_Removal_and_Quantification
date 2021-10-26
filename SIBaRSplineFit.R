require(lubridate)
require(mgcv)
require(dplyr)
require(magrittr)
require(splines)

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