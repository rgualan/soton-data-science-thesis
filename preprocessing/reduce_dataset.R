## Clean the environment
rm(list=ls())

## Load libraries
#library(randomForest)
source("util/my_helper.R")

## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
sites <- getSites(epa)

## Reduce by Latitude ##############################################################
plotStations(F,sites$Station.Code)

sites2 <- sites[sites$Latitude<36,]
plotStations(F,sites2$Station.Code)
dim(sites2)

epa <- epa[epa$Station.Code %in% sites2$Station.Code,]
saveRDS(epa, "data/epa/epa_daily/2016/california_ozone_plus_rcov_4.RDS")
