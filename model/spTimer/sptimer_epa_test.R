## demo("nyExample")

## Clean the environment
rm(list=ls())

## Load libraries
library(spTimer)
library(mgcv)
library(spBayes)
library(maps)
library(parallel)
source("util/my_helper.R")

## Global variables ###############################################################
paper <- setupPaper()


## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]
sites <- getSites(epa)
## Features derived from the Date field
epa <- addDoyField(epa)
epa <- addDowField(epa, sinTx = T)
## Standardize variable 
epa$sOzone <- scale(epa$Ozone)[,1]

## Reduce the size of the data for testing purposes
epa <- epa[epa$Station.Code %in% sort(unique(epa$Station.Code))[1:20],]
testStations <- unique(epa$Station.Code)[c(14,4)]
plotStations(F, unique(epa$Station.Code), redIds=testStations)

## Feature selection - through a linear model ############################################
head(epa)
lm.a <- lm(sOzone ~ Temperature+RH+sqrtWind+Elevation+Location.Setting+Doy+Dow.name+Dow.number, epa)
summary(lm.a)
stop("debug")

## Feature selection #####################################################################
#fm <- sOzone ~ 1
#fm <- sOzone ~ Temperature+Elevation+Dow.number+Doy
#fm <- sOzone ~ Temperature + Elevation
#fm <- sOzone ~ Temperature+RH+sqrtWind+Elevation+Location.Setting+Doy+Dow.name+Dow.number+isWeekday

epa.train <- epa[!epa$Station.Code %in% testStations,]
epa.test <- epa[epa$Station.Code %in% testStations,]
epa.test$sOzoneH <- NA  ## For the prediction
    
## Fit model
time.data<-spT.time(t.series=366,segments=1)
ticToc({
  post.gp <- spT.Gibbs(formula = fm, data = epa.train, model = "GP", 
                       coords = ~UTM.X+UTM.Y, #time.data=time.data,
                       newdata=epa.test, newcoords=~UTM.X+UTM.Y,
                       spatial.decay = spT.decay(distribution = Gamm(2,1), tuning = 0.1),
                       report=2)
})
#print(post.gp)
summary(post.gp)

## Prediction
epa.test$sOzoneH <- post.gp$prediction$Mean
plot(epa.test$sOzone, type="l")
lines(epa.test$sOzoneH, col=2)

## Metrics
evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
#a)   0.376   0.613   0.523 -92.008   0.317  -3.133   0.261   0.836  (sOzone ~ Temperature+Elevation+Dow.number+Doy)
#b)   0.253   0.503   0.418 -32.880  -0.033   0.325   0.188   0.886  (sOzone ~ 1)
#c)   0.196   0.443   0.364 -49.568   0.023  -0.226   0.146   0.899  (sOzone ~ Temperature)
#d)   0.299   0.547   0.462 -69.977   0.213  -2.110   0.216   0.867  (sOzone ~ Temperature + RH)
#e)   0.376   0.613   0.522 -93.775   0.316  -3.130   0.261   0.837

## Parameters
#post.gp$parameters
