require(tidyverse)
## Source in files containing SIBaR functions
source(paste0(getwd(),'/SIBaRUtils.R'))
source(paste0(getwd(),'/SIBaRPartitioningParallel.R'))
source(paste0(getwd(),'/SIBaRSplineFit.R'))

## To get started, copy and paste the directory that holds demo data next to dir
dir <- getwd()

## Read in the data. If using your own data, be sure to have some process
## which accounts for NA removal as functions won't work with NAs inserted.
demo_data <- read.csv(paste0(dir,"/DemoData.csv")) %>%
  mutate("Timestamps"=as.POSIXct(Time,format="%Y-%m-%d %H:%M:%S",tz="US/Central"),.keep="unused") %>%
  mutate("NOx"=Original,.keep="unused") %>%
  drop_na()

## Take the timestamps and NOx measurements out of the dataframe, store as 
## individual vectors
times <- demo_data$Timestamps
NOx <- demo_data$NOx

## First, we create an initial point partition. We do this by calling the 
## partitionPoints function from the SIBaRPartitioningParallel script file
partitionOutput <- partitionRoutine(NOx,times,25,transform_string = "log",
                                    length_tolerance = 0.05)
                                    
## To visualize our point partitions, we can make a simple call to plot.
## The first list entry of partition output are the measurements input into 
## the function, less any measurements corresponding to time periods shorter
## than minTimePts.
## The second list entry are the timestamps corresponding to the measurements
## in the first list entry
## The third list entry are the states returned by the partitioning step.
## 1 = background, 2 = non-background

separate_days <- partitionOutput %>%
  mutate(Month = month(Timestamps)) %>%
  mutate(Day = mday(Timestamps)) %>%
  group_split(Month,Day,.keep = FALSE)

for(i in 1:length(separate_days)){
  plot(Poll~Timestamps,data=separate_days[[i]],
       col=States,
       xlab = "Time",
       ylab = "ln(NOx+1)",
       main = paste0("Day ", i))
}
  
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

