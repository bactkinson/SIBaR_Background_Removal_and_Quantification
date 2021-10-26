# SIBaR_Background_Removal_and_Quantification

This repository contains functions necessary for implmentation of State Informed Background Removal (SIBaR). The method is described in detail in SIBaR: a new method for background quantification and removal, published in Atmospheric Measurement Techniques (https://amt.copernicus.org/articles/14/5809/2021/amt-14-5809-2021.html).

Overview:
Rather than implement a flat or rolling percentile to estimate background in air pollution time series, SIBaR casts the background problem as a regime change problem. We visualize the car going
through a series of regime changes that informs when it's in a clean state from when it's in a dirty state. We then take all the points designated as the clean state and fit a spline to them.

Insert picture here.

How to use:
The files A.R, B.R, and C.R contain the functions called to run SIBaR in SIBaRDemo.R. Feel free to use SIBaRDemo.R as your own template or write your own R script to incorporate the variety of SIBaR
function calls.

This work was created by Blake Actkinson, Katherine Ensor, and Robert Griffin of Rice University. 


