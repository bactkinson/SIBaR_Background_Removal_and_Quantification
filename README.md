# SIBaR_Background_Removal_and_Quantification

This repository contains code written in R for implmentation of State Informed Background Removal (SIBaR). The method is described in detail in SIBaR: a new method for background quantification and removal from mobile air pollution measurements, published in Atmospheric Measurement Techniques (https://amt.copernicus.org/articles/14/5809/2021/amt-14-5809-2021.html).

### Overview
Rather than implement a flat or rolling percentile to estimate background in air pollution time series, SIBaR casts the background estimation problem as a time series regime change problem and utilizes Hidden Markov models to solve it. We visualize the car going
through a series of changes between a clean regime and a dirty regime. We then take all the clean regime points and fit a spline to them. The method can be broken into two steps: the partitioning step, or when we identify clean and dirty regimes within our time series, and the spline fitting step, or when we fit a spline across all points taken in the clean regimes determined in the partitioning step.

![Transition Image](/Misc/Transition_Figure.jpg)
*Credits for this figure belong to Michael H. Actkinson.*

### How to Use
The files SIBaRUtils.R, SIBaRSplineFit.R, and SIBaRPartitioningParallel.R all contain functions needed to run the partitioning step and spline fitting step as described in the overview. `partitionRoutine` is a wrapper function which calls functions from SIBaRUtils.R and SIBaRPartitioningParallel.R to perform the partitioning step. 'sibarSplineFit` is a function which takes output from the partitioning step and fits the set of background splines which determine the background signal. An example of these two routines in action is given in SIBaRExec.R. Feel free to use SIBaRExec.R as a template for your own SIBaR applications.

It's worth noting that if the data contain timestamps from separate days, SIBaR will recognize and perform the partitioning and fitting steps on each day separately. For example, if the data contain 100 points taken on 1/1/2000 and 100 points on 1/2/2000, the partitioning step will be performed for the first set of 100 points, then the second set of 100 points. Background splines will then be fit to points taken in the clean regimes of each partitioning step separately.

### `partitionRoutine`

`partitionRoutine` is a wrapper function which calls additional functions in SIBaRUtils.R and SIBaRPartitioningParallel to perform the partitioning step on the data. Automatically implements a recursive corrections routine if data are initially labeled as misclassified. `partitionRoutine` requires the following parameters:

`poll`: Pollution measurements to be partitioned.

`time`: Timestamps corresponding to pollution measurements. Format must be POSIXct. Timestamps can span multiple days.

`bootstrap_iterations`: Number of HMM fits to be executed given randomly chosen starting transition probabilities. Out of the total number of HMM fits given by `bootstrap_iterations`, chooses the model fit with highest likelihood. 

`transform_string`: String which indicates whether to transform the data. Available transformations are "none" and "log". Default is "none."

`minTimePts`: The minimum number of time points needed for data to be partitioned. Data taken on days with less than minTimePts are
removed during the partitioning step. Default is 600.

`cores`: Number of cores to be used in parallelization. Parallelization is implemeted using R parallel library. Default is parallel::detectCores()-2.

`index`: Optional factor index to be included (e.g. 1 designates set of measurements from Car 1, 2 designates set of measurements from Car 2).

`threshold`: In evaluating misclassification instances, determines the percentage of fit background points lying above line of best fit or percentage of non-background points lying below line of best fit required for misclassification. Default is 50%.

`length_tolerance`: Parameter (as decimal, between 0 and 1) which determines stopping point for recursive corrections step. Recursive corrections won't continue if the lengths of the time series block fall below this tolerance. Default is 0.05, or 5% of the original time series length.

`save_misclassifications`: Boolean which determines whether SIBaR saves graphs (as a .png file) associated with initial classification decisions. Time series graphs are titled with either 'classified correctly' if the series categorizations are below the predetermined correction threshold and 'misclassified' if the series are above the predetermined correction threshold specified by `threshold`. Misclassified time series are automatically input into the SIBaR correction routine. If `save_misclassification` is `TRUE`, must be accompanied by `directory`

`directory`: The directory to which SIBaR saves graphs of initial classification decisions.

### `sibarSplineFit`

Takes output from the partitioning step performed in `partitionRoutine` and fits background splines to data taken during clean regimes. Data from clean regimes on multiple days are fit to separate background splines. Calls functions from SIBaRSplineFit.R Parameters are:

`poll`: Pollutant observations from the partitioning step.

`timestamps`: Timestamps from the partitioning step.

`state`: State designations determined in the partitioning step.

`index`: Optional factor index taken from the partitioning step.

If you have any issues in implementing this method or have suggestions for improvement, I'd love to hear about it. Please raise it as an issue ticket on the repository. Alternatively, you can email me at blake.w.actkinson@rice.edu.

This work was created by Blake Actkinson with input from Robert Griffin and Katherine Ensor from Rice University. 


