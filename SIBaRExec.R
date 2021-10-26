require(tidyverse)
## Source in files containing SIBaR functions
source(paste0(getwd(),'/SIBaRUtils.R'))
source(paste0(getwd(),'/SIBaRPartitioningParallel.R'))
source(paste0(getwd(),'/SIBaRSplineFit.R'))

## To get started, copy and paste the directory that holds demo data next to dir
dir <- getwd()

## Read in the data. If using your own data, be sure to have some process
## which accounts for NA removal as functions won't work with NAs inserted.
# demo_data <- read.csv(paste0(dir,"/DemoData.csv")) %>%
#   mutate("Timestamps"=as.POSIXct(Time,format="%Y-%m-%d %H:%M:%S",tz="US/Central"),.keep="unused") %>%
#   mutate("NOx"=Original,.keep="unused") %>%
#   drop_na()

## Take the timestamps and NOx measurements out of the dataframe, store as 
## individual vectors
# times <- demo_data$Timestamps
# NOx <- demo_data$NOx

Poll.times.c1 <- parse_date_time(as.character(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/LST_Car1_NOx_Revised_Final.csv", header = FALSE),use.names=FALSE)), orders=c("mdy HMS"), tz=Sys.timezone())
Index.c1 <- rep(1,length(Poll.times.c1))
Poll.measurements.c1 <- as.vector(unlist(read.csv("C:/Users/bwa2/Documents/Academic Work/Research/EDF Matlab Files/NOx_Raw_Car1_Revised_Final.csv", header = FALSE),use.names=FALSE))
temp_tibble <- tibble("Times"=Poll.times.c1,"Meas"=Poll.measurements.c1,
                      "Index"=Index.c1) %>%
  drop_na()

## First, we create an initial point partition. We do this by calling the 
## partitionPoints function from the SIBaRPartitioningParallel script file
partitionOutput <- partitionRoutine(temp_tibble$Meas,temp_tibble$Times,2,transform_string = "log",
                                    length_tolerance = 0.2)

## To visualize our point paritions, we can make a simple call to plot.
## The first list entry of partition output are the measurements input into 
## the function, less any measurements corresponding to time periods shorter
## than minTimePts.
## The second list entry are the timestamps corresponding to the measurements
## in the first list entry
## The third list entry are the states returned by the partitioning step.
## 1 = background, 2 = non-background

plot(Poll~Timestamps,data=partitionOutput,
     col=States,
     xlab = "Time",
     ylab = "ln(NOx+1)",
     main = "Visualizing the partitioning step")

## Next, we fit a spline to our partitioned data. 
background <- sibarSplineFit(partitionOutput$Poll,partitionOutput$Timestamps,
                             partitionOutput$States)

## Visualize the background signal
plot(Poll~Timestamps,data=partitionOutput,
     col=States,
     xlab = "Time",
     ylab = "ln(NOx+1)",
     main = "Visualizing the partitioning step")
lines(partitionOutput$Timestamps,background,col="blue",lwd=2)


