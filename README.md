# SIBaR_Background_Removal_and_Quantification

This repository contains functions necessary for implmentation of State Informed Background Removal (SIBaR). The method is described in detail in SIBaR: a new method for background quantification and removal, published in Atmospheric Measurement Techniques (https://amt.copernicus.org/articles/14/5809/2021/amt-14-5809-2021.html).

### Overview
Rather than implement a flat or rolling percentile to estimate background in air pollution time series, SIBaR casts the background problem as a time series regime change problem and utilizes time series regime change modeling techniques to solve it. We visualize the car going
through a series of regime changes that informs when it's in a clean state from when it's in a dirty state. We then take all the points designated as the clean state and fit a spline to them.

![Transition Image](/Misc/Transition_Figure.jpg)
*Credits for this figure belong to Michael H. Actkinson.*

### How to Use
The files SIBaRUtils.R, SIBaRSplineFit.R, and SIBaRPartitioningParallel.R contain the functions called to run the SIBaR routine in SIBaRExec.R. Feel free to use SIBaRExec.R as your own template or write your own R script to incorporate the variety of SIBaR
function calls.

Below contains a brief description of the function calls used in SIBaRExec.R.

'partitionRoutine': Wrapper function which calls the remaining functions in SIBaRPartitioningParallel to perform the partitioning step on the data. Automatically implements recursive corrections routine if misclassified data are detected. partitionRoutine requires the following parameters.

'poll': Pollution measurements to be partitioned.

'time': Timestamps corresponding to pollution measurements. Format must be POSIXct. Timestamps can span multiple days.

'bootstrap_iterations': Number of Expectation Maximization algorithm executions performed in the HMM step.

'transform_string': String which indicates whether to transform the data. Default is "none." "log" is available.

'minTimePts': The minimum number of time points needed for data to be partitioned. Data taken on days with less than minTimePts are
removed during the partitioning step. Default is 600.

'cores': Number of cores to be used in parallelization. Utilizes parallel library. Default is parallel::detectCores()-2

'index': Optional factor index to be included (e.g. 1 designates set of measurements from Car 1, 2 designates set of measurements from Car 2).

'threshold': In evaluating misclassification instances, determines the percentage of fit background points lying above line of best fit or percentage of non-background points lying below line of best fit required for misclassification.

'length_tolerance': Parameter (as decimal, between 0 and 1) which determines stopping point for recursive corrections step. Recursive corrections won't happen for time series below this parameter multiplied by the original time series length.

Once we obtain the output from partitionRoutine, we can fit background splines to it using the SIBaRSplineFit function available in the SIBaRSplineFit.R file.






This work was created by Blake Actkinson with input from Robert Griffin and Katherine Ensor from Rice University. 


