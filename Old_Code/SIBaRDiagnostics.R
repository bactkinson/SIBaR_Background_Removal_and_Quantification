rm(list=ls())
require(zoo)
require(lubridate)
smoothData <- function(poll,time,t.interval){
  if(t.interval==0){
    return(list(poll,time))
  } else{
    poll.zoo <- zoo(poll, order.by = time)
    total.time <- seq(time[1],time[length(time)],1)
    na.zoo <- zoo(,order.by = total.time)
    combined.zoo <- merge.zoo(poll.zoo,na.zoo,all=TRUE)
    smoothed.zoo <- rollapply(combined.zoo,t.interval,mean,align="left",na.rm=FALSE,partial=FALSE)
    Poll.smooth <- coredata(smoothed.zoo)
    Poll.times.final <- index(smoothed.zoo)
    return(list(Poll.smooth,Poll.times.final))
  }
}
# Import Moody Tower Data
moody.data <- rbind(
  read.csv('C:/Users/bwa2/Documents/Academic Work/Research/TCEQ Monitor Data/MoodyTower_Jan1_Jan15_2018.csv',colClasses = "character"),
  read.csv('C:/Users/bwa2/Documents/Academic Work/Research/TCEQ Monitor Data/MoodyTower_Jan15_April30_2018.csv',colClasses = "character")
)
# Convert time and date to datetime
moody.data.datetime <- lubridate::parse_date_time(paste(moody.data$dateGMT, moody.data$timeGMT), order="ymd HMS")
# Clean up the data frame. Insert datetime column, take out the date and time columns, and rename the columns to be easier to work with
moody.data$dateGMT <- moody.data.datetime
moody.data.clean <- subset(moody.data,select= -c(timeGMT))
colnames(moody.data.clean)[colnames(moody.data.clean) == "NOx_NO2conc_flag"] <- "NO2_flag"
colnames(moody.data.clean)[colnames(moody.data.clean) == "NOx_NO2conc_value"] <- "NO2_value"
colnames(moody.data.clean)[colnames(moody.data.clean) == "NOx_NOconc_flag"] <- "NO_flag"
colnames(moody.data.clean)[colnames(moody.data.clean) == "NOx_NOconc_value"] <- "NO_value"
rmIdx <- moody.data.clean$NO2_flag=="K" & moody.data.clean$NO_flag=="K"
moody.data.clean <- moody.data.clean[rmIdx,]
attr(moody.data.clean$dateGMT,"tzone") <- "US/Central"
moody.data.clean <- moody.data.clean[-2888,]

## Evaluating NRMSE between SIBaR background, Apte background, Moody tower
# Import background data
require(zoo)
sibar.background <- unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_HMMBackground_NOx_Revised_Final.csv", header = TRUE, row.names = 1),use.names=FALSE)
# Transform the data
sibar.background <- exp(sibar.background)-1
NOx.idx <- unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Index_Revised_Final.csv", header = TRUE, row.names = 1),use.names=FALSE)
NOx.time <- parse_date_time(as.character(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_NOx_Time_Revised_Final.csv", header = TRUE, row.names=1),use.names=FALSE)), orders=c("ymd HMS"), tz=Sys.timezone())
# Filtering data to be between 9 AM and 4 PM. Comment if undesired.
time.filter <- ((hour(NOx.time) > 9) & (hour(NOx.time) < 16) & (month(NOx.time) >= month(moody.data.clean$dateGMT[1])) & (month(NOx.time) < month(moody.data.clean$dateGMT[nrow(moody.data.clean)]))) 
sibar.background <- sibar.background[time.filter]
NOx.idx <- NOx.idx[time.filter]
NOx.time <- NOx.time[time.filter]

sibar.car.1 <- sibar.background[NOx.idx==1]
NOx.time.1 <- NOx.time[NOx.idx==1]
sibar.car.2 <- sibar.background[NOx.idx==2]
NOx.time.2 <- NOx.time[NOx.idx==2]
start.time <- Sys.time()
sibar.mean.1 <- smoothData(sibar.car.1,NOx.time.1,300)

# Take out nans
sibar.mean.1.NOx <- sibar.mean.1[[1]]
sibar.time.1 <- sibar.mean.1[[2]]
na.idx <- is.na(sibar.mean.1.NOx)
sibar.time.1 <- sibar.time.1[!na.idx]
sibar.mean.1.NOx <- sibar.mean.1.NOx[!na.idx]
# Car 2
sibar.mean.2 <- smoothData(sibar.car.2,NOx.time.2,300)
# Take out nans
sibar.mean.2.NOx <- sibar.mean.2[[1]]
sibar.time.2 <- sibar.mean.2[[2]]
na.idx <- is.na(sibar.mean.2.NOx)
sibar.mean.2.NOx <- sibar.mean.2.NOx[!na.idx]
sibar.time.2 <- sibar.time.2[!na.idx]
# Match up dates with one another to compare.
NOx.compare.c1 <- numeric(length(moody.data.clean$dateGMT))
NOx.compare.c2 <- numeric(length(moody.data.clean$dateGMT))

for (i in 1:length(NOx.compare.c1)){
  matchingIdx <- which(moody.data.clean$dateGMT[i]==sibar.time.1)
  if(length(matchingIdx)==0){
      NOx.compare.c1[i] <- NA
  } else {
      NOx.compare.c1[i] <- sibar.mean.1.NOx[matchingIdx]
  }
}

for (i in 1:length(NOx.compare.c2)){
  matchingIdx <- which(moody.data.clean$dateGMT[i]==sibar.time.2)
  if(length(matchingIdx)==0){
      NOx.compare.c2[i] <- NA
  } else {
      NOx.compare.c2[i] <- sibar.mean.2.NOx[matchingIdx]
  }
}

moody.NOx <- as.numeric(moody.data.clean$NO_value)+as.numeric(moody.data.clean$NO2_value)
zoo.one <- zoo::zoo(NOx.compare.c1,order.by=moody.data.clean$dateGMT)
zoo.two <- zoo::zoo(NOx.compare.c2,order.by=moody.data.clean$dateGMT)
zoo.moody <- zoo::zoo(moody.NOx,order.by=moody.data.clean$dateGMT)
zoo.overall <- merge.zoo(zoo.one,zoo.two,zoo.moody)
## Trim double na columns.
na.double <- is.na(zoo.overall[,1]) & is.na(zoo.overall[,2])
zoo.overall <- zoo.overall[!na.double,]
sibar.values <- numeric(nrow(zoo.overall))
for (i in 1:length(sibar.values)){sibar.values[i] <- mean(zoo.overall[i,1:2],na.rm = TRUE)}
sibar.overall <- merge(zoo.overall,sibar.values)

# Do the same process for brantley, apte background signals.
# Import background data
brant.background <- unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_BrantBackground_NOx_Log.csv", header = TRUE, row.names = 1),use.names=FALSE)
# Transform the data
brant.background <- brant.background[time.filter]
brant.Car.1 <- brant.background[NOx.idx==1]
brant.Car.2 <- brant.background[NOx.idx==2]
brant.mean.1 <- smoothData(brant.Car.1,NOx.time.1,300)
# Take out nans
brant.mean.1.NOx <- brant.mean.1[[1]]
brant.time.1 <- brant.mean.1[[2]]
na.idx <- is.na(brant.mean.1.NOx)
brant.time.1 <- brant.time.1[!na.idx]
brant.mean.1.NOx <- brant.mean.1.NOx[!na.idx]
# Car 2
brant.mean.2 <- smoothData(brant.Car.2,NOx.time.2,300)
# Take out nans
brant.mean.2.NOx <- brant.mean.2[[1]]
brant.time.2 <- brant.mean.2[[2]]
na.idx <- is.na(brant.mean.2.NOx)
brant.mean.2.NOx <- brant.mean.2.NOx[!na.idx]
brant.time.2 <- brant.time.2[!na.idx]

# brant.mean.1.NOx <- read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/brant_mean_1_NOx.csv")[,2]
# brant.mean.2.NOx <- read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/brant_mean_2_NOx.csv")[,2]
# brant.time.1 <- parse_date_time(as.character(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/brant_time_1_NOx.csv")[,2]), orders=c("ymd HMS"), tz=Sys.timezone())
# brant.time.2 <- parse_date_time(as.character(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/brant_time_2_NOx.csv")[,2]), orders=c("ymd HMS"), tz=Sys.timezone())

brant.compare.c1 <- numeric(length(moody.data.clean$dateGMT))
brant.compare.c2 <- numeric(length(moody.data.clean$dateGMT))
for (i in 1:length(brant.compare.c1)){
  matchingIdx <- which(moody.data.clean$dateGMT[i]==brant.time.1)
  if(length(matchingIdx)==0){
    brant.compare.c1[i] <- NA
  } else {
    brant.compare.c1[i] <- brant.mean.1.NOx[matchingIdx]
  }
}

for (i in 1:length(brant.compare.c2)){
  matchingIdx <- which(moody.data.clean$dateGMT[i]==brant.time.2)
  if(length(matchingIdx)==0){
    brant.compare.c2[i] <- NA
  } else {
    brant.compare.c2[i] <- brant.mean.2.NOx[matchingIdx]
  }
}

brant.one <- zoo::zoo(brant.compare.c1,order.by=moody.data.clean$dateGMT)
brant.two <- zoo::zoo(brant.compare.c2,order.by=moody.data.clean$dateGMT)

brant.overall <- merge.zoo(brant.one,brant.two,zoo.moody)
## Trim double na columns.
na.double <- is.na(brant.overall[,1]) & is.na(brant.overall[,2])
brant.overall <- brant.overall[!na.double,]
brant.values <- numeric(nrow(brant.overall))
for (i in 1:length(brant.values)){brant.values[i] <- mean(brant.overall[i,1:2],na.rm = TRUE)}
brant.overall <- merge(brant.overall,brant.values)

# Import background data
apte.background <- unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/State Files/Total_ApteBackground_NOx_Log.csv", header = TRUE, row.names = 1),use.names=FALSE)
# Transform the data
apte.background <- apte.background[time.filter]
apte.car.1 <- apte.background[NOx.idx==1]
apte.car.2 <- apte.background[NOx.idx==2]
apte.mean.1 <- smoothData(apte.car.1,NOx.time.1,300)
# Take out nans
apte.mean.1.NOx <- apte.mean.1[[1]]
apte.time.1 <- apte.mean.1[[2]]
na.idx <- is.na(apte.mean.1.NOx)
apte.time.1 <- apte.time.1[!na.idx]
apte.mean.1.NOx <- apte.mean.1.NOx[!na.idx]
# Car 2
apte.mean.2 <- smoothData(apte.car.2,NOx.time.2,300)
# Take out nans
apte.mean.2.NOx <- apte.mean.2[[1]]
apte.time.2 <- apte.mean.2[[2]]
na.idx <- is.na(apte.mean.2.NOx)
apte.mean.2.NOx <- apte.mean.2.NOx[!na.idx]
apte.time.2 <- apte.time.2[!na.idx]

# apte.mean.1.NOx <- read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/apte_mean_1_NOx.csv")[,2]
# apte.mean.2.NOx <- read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/apte_mean_2_NOx.csv")[,2]
# apte.time.1 <- parse_date_time(as.character(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/apte_time_1_NOx.csv")[,2]), orders=c("ymd HMS"), tz=Sys.timezone())
# apte.time.2 <- parse_date_time(as.character(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/Processed Data Files/apte_time_2_NOx.csv")[,2]), orders=c("ymd HMS"), tz=Sys.timezone())

apte.compare.c1 <- numeric(length(moody.data.clean$dateGMT))
apte.compare.c2 <- numeric(length(moody.data.clean$dateGMT))
for (i in 1:length(apte.compare.c1)){
  matchingIdx <- which(moody.data.clean$dateGMT[i]==apte.time.1)
  if(length(matchingIdx)==0){
    apte.compare.c1[i] <- NA
  } else {
    apte.compare.c1[i] <- apte.mean.1.NOx[matchingIdx]
  }
}

for (i in 1:length(apte.compare.c2)){
  matchingIdx <- which(moody.data.clean$dateGMT[i]==apte.time.2)
  if(length(matchingIdx)==0){
    apte.compare.c2[i] <- NA
  } else {
    apte.compare.c2[i] <- apte.mean.2.NOx[matchingIdx]
  }
}

apte.one <- zoo::zoo(apte.compare.c1,order.by=moody.data.clean$dateGMT)
apte.two <- zoo::zoo(apte.compare.c2,order.by=moody.data.clean$dateGMT)
apte.overall <- merge.zoo(apte.one,apte.two,zoo.moody)

## Trim double na columns.
na.double <- is.na(apte.overall[,1]) & is.na(apte.overall[,2])
apte.overall <- apte.overall[!na.double,]
apte.values <- numeric(nrow(apte.overall))
for (i in 1:length(apte.values)){apte.values[i] <- mean(apte.overall[i,1:2],na.rm = TRUE)}
apte.overall <- merge(apte.overall,apte.values)

## Create plots
## Want to create plots comparing background signals on day-by-day basis
apte.zoo.to.plot <- apte.overall[,3:4]
brant.zoo.to.plot <- brant.overall[,3:4]
sibar.zoo.to.plot <- sibar.overall[,3:4]

poll.days <- lubridate::mday(index(apte.zoo.to.plot))
# Extract and store indices of instances in which day changes in breaks
day.breaks <- which(diff(poll.days)!=0)
day.diffs <- diff(day.breaks)
total.idx <- length(day.breaks)+1
apte.rmse.vec <- numeric(total.idx)
apte.mae.vec <- numeric(total.idx)
# Create the plots
for (i in 1:total.idx){
  print(i)
  if (i==1){
    temp.zoo <- apte.zoo.to.plot[1:day.breaks[i],]
  } else if (i==(length(day.breaks)+1)){
    temp.zoo <- apte.zoo.to.plot[(day.breaks[length(day.breaks)]+1):nrow(apte.zoo.to.plot),]
  } else {
    temp.zoo <- apte.zoo.to.plot[(day.breaks[i-1]+1):day.breaks[i],]
  }
  plot.zoo(temp.zoo,plot.type="single",col=1:2,type="p",pch=19,main=paste("Apte Date:",index(temp.zoo)[1]),ylab = "NOx (ppb)")
  legend("topright", inset=c(0,0), y.intersp = 1, legend = c("Moody", "Apte"),  lty = 1, bty = "n", col = c(1,2), cex = .5)
  ## CHANGE FILE PATH DEPENDING ON WHERE YOU WANT PLOTS SAVED
  apte.rmse.vec[i] <- Metrics::rmse(temp.zoo$zoo.moody,temp.zoo$apte.values)
  apte.mae.vec[i] <- Metrics::mae(temp.zoo$zoo.moody,temp.zoo$apte.values)
}

poll.days <- mday(index(sibar.zoo.to.plot))
# Extract and store indices of instances in which day changes in breaks
day.breaks <- which(diff(poll.days)!=0)
day.diffs <- diff(day.breaks)
total.idx <- length(day.breaks)+1
sibar.rmse.vec <- numeric(total.idx)
sibar.mae.vec <- numeric(total.idx)
# Create the plots and calculate daily RMSE, MAE
for (i in 1:total.idx){
  print(i)
  if (i==1){
    temp.zoo <- sibar.zoo.to.plot[1:day.breaks[i],]
  } else if (i==(length(day.breaks)+1)){
    temp.zoo <- sibar.zoo.to.plot[(day.breaks[length(day.breaks)]+1):nrow(sibar.zoo.to.plot),]
  } else {
    temp.zoo <- sibar.zoo.to.plot[(day.breaks[i-1]+1):day.breaks[i],]
  }
  plot.zoo(temp.zoo,plot.type="single",col=1:2,type="p",pch=19,main=paste("SIBaR Date:",index(temp.zoo)[1]),ylab = "NOx (ppb)")
  legend("topright", inset=c(0,0), y.intersp = 1, legend = c("Moody", "SIBaR"),  lty = 1, bty = "n", col = c(1,2), cex = .5)
  ## CHANGE FILE PATH DEPENDING ON WHERE YOU WANT PLOTS SAVED
  sibar.rmse.vec[i] <- Metrics::rmse(temp.zoo$zoo.moody,temp.zoo$sibar.values)
  sibar.mae.vec[i] <- Metrics::mae(temp.zoo$zoo.moody,temp.zoo$sibar.values)
}

poll.days <- mday(index(brant.zoo.to.plot))
# Extract and store indices of instances in which day changes in breaks
day.breaks <- which(diff(poll.days)!=0)
day.diffs <- diff(day.breaks)
total.idx <- length(day.breaks)+1
brant.rmse.vec <- numeric(total.idx)
brant.mae.vec <- numeric(total.idx)
# Create the plots
for (i in 1:total.idx){
  print(i)
  if (i==1){
    temp.zoo <- brant.zoo.to.plot[1:day.breaks[i],]
  } else if (i==(length(day.breaks)+1)){
    temp.zoo <- brant.zoo.to.plot[(day.breaks[length(day.breaks)]+1):nrow(brant.zoo.to.plot),]
  } else {
    temp.zoo <- brant.zoo.to.plot[(day.breaks[i-1]+1):day.breaks[i],]
  }
  plot.zoo(temp.zoo,plot.type="single",col=1:2,type="p",pch=19,main=paste("Brant Date:",index(temp.zoo)[1]),ylab = "NOx (ppb)")
  legend("topright", inset=c(0,0), y.intersp = 1, legend = c("Moody", "Brantley"),  lty = 1, bty = "n", col = c(1,2), cex = .5)
  ## CHANGE FILE PATH DEPENDING ON WHERE YOU WANT PLOTS SAVED
  brant.rmse.vec[i] <- Metrics::rmse(temp.zoo$zoo.moody,temp.zoo$brant.values)
  brant.mae.vec[i] <- Metrics::mae(temp.zoo$zoo.moody,temp.zoo$brant.values)
}

## Numerical metrics evaluation
apte.rmse <- Metrics::rmse(apte.overall$zoo.moody,apte.overall$apte.values)
apte.rmse
lowest.apte <- length(which(apte.rmse.vec<brant.rmse.vec & apte.rmse.vec<sibar.rmse.vec))
lowest.apte

brant.rmse <- Metrics::rmse(brant.overall$brant.values,brant.overall$zoo.moody)
brant.rmse
lowest.brant <- length(which(brant.rmse.vec<apte.rmse.vec & brant.rmse.vec<sibar.rmse.vec))
lowest.brant

sibar.rmse <- Metrics::rmse(sibar.overall$sibar.values,sibar.overall$zoo.moody)
sibar.rmse
lowest.sibar <- length(which(sibar.rmse.vec<apte.rmse.vec & sibar.rmse.vec<brant.rmse.vec))
lowest.sibar

apte.mae <- Metrics::mae(apte.overall$zoo.moody,apte.overall$apte.values)
apte.mae
lowest.apte.mae <- length(which(apte.mae.vec<brant.mae.vec & apte.mae.vec<sibar.mae.vec))
lowest.apte.mae

brant.mae <- Metrics::mae(brant.overall$zoo.moody,brant.overall$brant.values)
brant.mae
lowest.brant.mae <- length(which(brant.mae.vec<apte.mae.vec & brant.mae.vec<sibar.mae.vec))
lowest.brant.mae

sibar.mae <- Metrics::mae(sibar.overall$zoo.moody,sibar.overall$sibar.values)
sibar.mae
lowest.sibar.mae <- length(which(sibar.mae.vec<apte.mae.vec & sibar.mae.vec<brant.mae.vec))
lowest.sibar.mae

length(apte.zoo.to.plot)
length(sibar.zoo.to.plot)
length(brant.zoo.to.plot)

qz <- merge.zoo(apte.zoo.to.plot,sibar.zoo.to.plot,brant.zoo.to.plot)
