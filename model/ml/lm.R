## Clean environment
rm(list=ls())

## Libraries
library("openair")
source("util/my_helper.R")

## Read data ########################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- addDateDerivedFeatures(epa)
epa <- transformFeatures(epa)
epa$sOzone <- scale(epa$Ozone)[,1]
sites <- getSites(epa)

epa <- addNeighboursAverage(epa,5)
## Problem:
## Station 023-1005
## It is isolated in the northwest R:-0.18!!!
## Luckily this station is the only problem

## Split data #######################################################################
folds <- getFolds()
epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=1],] 
epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==1],]

## Linear model A ###################################################################
## Withouth including a neighbor combination
fmA <- sOzone ~ Temperature+RH+logRain+sqrtWind+Dew.Point+Water.Evap+
  Geop.Height+Geop.Height.Tropo+Lat.Heat.Flux+ #Heat.Flux
  Tropo.Press+Press.MSL+Vegetation+Longitude+Latitude+Elevation+
  Location.Setting+Doy+Dow.name+Month+Day
lmModelA <- lm(fmA, data=epa.train)
summary(lmModelA)
#plot(lmModel)

epa.test$sOzoneH <- predict(lmModelA, epa.test)
(ma <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH))
plot(epa.test$sOzone,type="l"); lines(epa.test$sOzoneH,col=2)

epa.test$sOzoneRes <- epa.test$sOzone - epa.test$sOzoneH
plot(epa.test$sOzoneRes,type="l")
hist(epa.test$sOzoneRes)


## Linear model B ###################################################################
## Including a neighbor combination
# fm <- sOzone ~ Temperature+RH+logRain+sqrtWind+Dew.Point+Water.Evaporation+
#   Geopotential.Height+Geopotential.Height.Tropo+Latent.Heat.FLux+ #Heat.Flux
#   Tropopause.Press+Press.MSL+Vegetation+Longitude+Latitude+Elevation+
#   Location.Setting+Doy+Dow.name+Month+Day+sqrtWind+logRain+Neighbor
fmB <- sOzone ~ Temperature+Dew.Point+Water.Evap+
  Geop.Height+Geop.Height.Tropo+#Heat.Flux
  Tropo.Press+Press.MSL+Longitude+Latitude+Elevation+
  Location.Setting+Doy+Neighbor
lmModelB <- lm(fmB, data=epa.train)
summary(lmModelB)
#plot(lmModel)

epa.test$sOzoneH <- predict(lmModelB, epa.test)
(mb <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH))
plot(epa.test$sOzone,type="l"); lines(epa.test$sOzoneH,col=2)

## Notes:
## Adding the neighbor increases the R2 from
## And several covariates stop being important
rbind(ma,mb)



