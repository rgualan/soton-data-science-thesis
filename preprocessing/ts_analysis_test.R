## Load required libraries
library(sp)
library(gstat)
library(spacetime)
library(reshape2)
library(lattice)
library(scales)
source("util/my_helper.R")

## Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
epa <- addDoyField(epa)
epa <- addDowField(epa)
sites <- getSites(epa)


## Simple TS analysis ###################################################
aggregate(Ozone~Station.Code,epa,length)

epa$Ozone100 <- rescale(epa$Ozone, to=c(0, 100)) 
x = ts(epa$Ozone100[epa$Station.Code=="113-1003"], frequency=7, start=c(2016,1))
fit = stl(x, s.window="periodic")
plot(fit)

## Notes:
## No weekly seasonal trend: its contribution is neglectable
## The trend is not appropiate
## Remainder is large