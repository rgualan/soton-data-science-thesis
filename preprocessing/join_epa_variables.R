## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(mice)
library(Amelia)
library(VIM)
library(gstat)
library(raster)
library(spacetime)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()


## Read data #########################################################################
ozone <- readRDS("data/epa/epa_daily/2016/california_ozone_2.RDS")
temp <- readRDS("data/epa/epa_daily/2016/california_temperature.RDS")
wind <- readRDS("data/epa/epa_daily/2016/california_wind.RDS")
rh <- readRDS("data/epa/epa_daily/2016/california_rh.RDS")
# ozone.sites <- getSites(ozone)
# temp.sites <- getSites(temp)
# wind.sites <- getSites(wind)
# rh.sites <- getSites(rh)
## Read sites
#sites <- readRDS("data/epa/sites/aqs_sites.RDS")
#sites <- sites[sites$State.Name=="California",]


## Merge with the other variables ########################################################
ozone <- fillMissingDates(ozone)
## Implicit filling of missing dates in the other data frames:
cali <- merge(ozone, temp, all.x=T)
cali <- merge(cali, wind, all.x=T)
cali <- merge(cali, rh, all.x=T)
#View(cali)
summary(cali)

sites <- getSites(cali)
cali <- merge(cali,sites[,c("Station.Code","Longitude","Latitude","UTM.X","UTM.Y","Elevation","Location.Setting")])
saveRDS(cali,file="data/epa/epa_daily/2016/california_ozone_plus_cov.RDS")

