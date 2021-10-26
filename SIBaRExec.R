require(depmixS4)
require(tidyverse)

## To get started, copy and paste the directory that holds demo data next to dir
dir <- getwd()

demo_data <- read.csv(paste0(dir,"/DemoData.csv")) %>%
  mutate("Timestamps"=as.POSIXct(Time,format="%Y-%m-%d %H:%M:%S",tz="US/Central"),.keep="unused")

## Our goal is to reproduce Figure 1 in the manuscript
## The results are contained within demo data
plot(Original~Timestamps,col=State,data=demo_data,
     main = "These are the state designations that appear in the manuscript",
     ylab = "NOx (ppb)",
     ylim = c(-5,100),
     pch=19)

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