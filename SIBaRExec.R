require(tidyverse)
## Source in files containing SIBaR functions
source(paste0(getwd(),'/SIBaRUtils.R'))
source(paste0(getwd(),'/SIBaRPartitioningParallel.R'))
source(paste0(getwd(),'/SIBaRSplineFit.R'))
source(paste0(getwd(),'/SIBaRRecursiveCorrections.R'))

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
partitionOutput <- partitionPoints(NOx,times,25,transformString = "log")

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

## Run diagnostics on the fit to ensure it's of high quality



## First, we need to log transform the NOx data. We'll create a new tibble
## and work from there

transformed_data <- demo_data %>%
  select(Original,Timestamps) %>%
  mutate(NOx=log(Original+1),.keep="unused")

## Now we come to fitting the HMM model. We utilize the depmix and fit functions 
## from depmixS4 to do this.
set.seed(10)
## We use depmix to create the model.
model <- depmix(NOx~Timestamps, nstates=2, 
                    data=transformed_data,
                    family=gaussian())

## We use fit to fit the model.
model_fit <- fit(model)

## We use the viterbi algorithm, accessed by the viterbi call, to extract 
## decoded states
fitted_states <- viterbi(model_fit)[,1]

## And we plot the results
retransformed_data <- transformed_data %>%
  cbind(fitted_states) %>%
  mutate(NOx=exp(NOx)-1,.keep="unused")

plot(Original~Timestamps,col=State,data=demo_data,
     main = "These are the state designations that appear in the manuscript",
     ylab = "NOx (ppb)",
     ylim = c(-5,100),
     pch=19)
plot(NOx~Timestamps,col=fitted_states,data=retransformed_data,
     main = "These are the state designations we just fit",
     ylab = "NOx (ppb)",
     ylim = c(-5,100),
     pch=19)

## Changing the starting values for the depmix model can change the outcome
## To change the starting values for the transition probabilities, use 
## trstart
trst <- runif(4)
new_model <- depmix(NOx~Timestamps, nstates=2, 
                    data=transformed_data,
                    family=gaussian(), 
                    trstart=trst)
new_fit <- fit(new_model)