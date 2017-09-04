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
fm <- sOzone ~ Temperature+RH+logRain+sqrtWind+Dew.Point+Water.Evaporation+
  Geopotential.Height+Geopotential.Height.Tropo+Latent.Heat.FLux+ #Heat.Flux
  Tropopause.Press+Press.MSL+Vegetation+Longitude+Latitude+Elevation+
  Location.Setting+Doy+Dow.name+Month+Day
lmModelA <- lm(fm, data=epa.train)
summary(lmModelA)
#plot(lmModel)

epa.test$sOzoneH <- predict(lmModelA, epa.test)
(ma <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH))
plot(epa.test$sOzone,type="l"); lines(epa.test$sOzoneH,col=2)

## Linear model B ###################################################################
## Including a neighbor combination
# fm <- sOzone ~ Temperature+RH+logRain+sqrtWind+Dew.Point+Water.Evaporation+
#   Geopotential.Height+Geopotential.Height.Tropo+Latent.Heat.FLux+ #Heat.Flux
#   Tropopause.Press+Press.MSL+Vegetation+Longitude+Latitude+Elevation+
#   Location.Setting+Doy+Dow.name+Month+Day+sqrtWind+logRain+Neighbor
fm <- sOzone ~ Temperature+Dew.Point+Water.Evaporation+
  Geopotential.Height+Geopotential.Height.Tropo+#Heat.Flux
  Tropopause.Press+Press.MSL+Longitude+Latitude+Elevation+
  Doy+Neighbor #Location.Setting+
lmModelB <- lm(fm, data=epa.train)
summary(lmModelB)
#plot(lmModel)

epa.test$sOzoneH <- predict(lmModelB, epa.test)
(mb <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH))
plot(epa.test$sOzone,type="l"); lines(epa.test$sOzoneH,col=2)

## Notes:
## Adding the neighbor increases the R2 from
## And several covariates stop being important
rbind(ma,mb)



x <- 1:366
# Coldest day?
# dailyTempAgg <- aggregate(Temperature~Doy,epa,mean)
# dailyTempAgg[which.min(dailyTempAgg$Temperature),]
# plot(Temperature~Doy, dailyTempAgg, type="l")
cdayt <- sin((x-1)*(2*pi)/(366))
plot(cdayt, type="l")

# fit.gam <- gam(o8hrmax ~ s(cMAXTMP) + s(WDSP) + s(RH) + s(Longitude, Latitude, k = 10), 
#   data = DataFit)
# pred.gam <- predict(fit.gam, DataValPred, interval = "prediction")