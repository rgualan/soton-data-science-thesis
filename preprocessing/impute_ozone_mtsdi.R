## Clean environment
rm(list=ls())

## Load required libraries
library(reshape2)
library(mtsdi)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Sites
sites <- getSites(epa)
#str(epa)

## Simple imputation of Ozone ################################################
covByStation <- aggregate(cbind(Ozone,Temperature)~Station.Code,epa,length)
covByStation <- covByStation[order(covByStation$Ozone),]
head(covByStation)
testStation <- "065-2002" # The one with more missing data (35 NAs)

kn <- getKneighbours(testStation, sites, 5, include_own=T)
epa_sub <- epa[epa$Station.Code %in% kn$Station.Code,] ## Zone of interest (near the target variable)
epa_sub2 <- dcast(epa_sub, Date ~ Station.Code, mean, value.var="Ozone") # average effect of time
epa_sub2 <- epa_sub2[,-1]
epa_sub2 <- edaprep(epa_sub2)
mstats(epa_sub2)

names(epa_sub2) <- paste0("c",1:6)

#f <- ~Temperature+Elevation
f <- ~c1+c2+c3+c4+c5+c6
i <- mnimput(f,epa_sub2,eps=1e-3,ts=TRUE,method="spline")
## DOES NOT WORK!!!

summary(i)



## Example
data(miss)
f <- ~c31+c32+c33+c34+c35
## one-window covariance
i <- mnimput(f,miss,eps=1e-3,ts=TRUE, method="spline",sp.control=list(df=c(7,7,7,7,7)))
summary(i)
## two-window covariances
b<-c(rep("year1",12),rep("year2",12))
ii <- mnimput(f,miss,by=b,eps=1e-3,ts=TRUE, method="spline",sp.control=list(df=c(7,7,7,7,7)))
summary(ii)
plot(i)
