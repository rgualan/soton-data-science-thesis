## Clean environment
rm(list=ls())

## Libraries
source("util/my_helper.R")

## Read data ########################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- addDateDerivedFeatures(epa)
epa <- addNeighboursAverage(epa,5)
epa <- transformFeatures(epa)
epa <- scaleTargetVariable(epa)
sites <- getSites(epa)

## Add combination of closest neighbours as a predictor
## Problem: Station 023-1005. It is isolated in the northwest R:-0.18!!!
## Luckily this station is the only problem
## Notes:
## k was chosen based on the best R2

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


## 10-fold cross validation ############################################################
metrics = c()
epa$sOzoneH <- NA
ticToc({
  for(k in 1:10){
    print(paste("Fold", k,paste(rep("=",50),collapse = "")))
    epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],] 
    epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
    
    epa.test$sOzoneH <- NA
    lmModel <- lm(fmB, data=epa.train)
    epa.test$sOzoneH <- predict(lmModel,epa.test)
    #m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
    m <- getMetricsByStationFromDF(epa.test, "sOzone", "sOzoneH")
    # plot(epa.test$sOzone,type="l")
    # lines(epa.test$OzoneTps,col=2)
    
    ## Output
    metrics <- round(rbind(metrics,m),3)
    for(s in unique(epa.test$Station.Code)){
      epa$sOzoneH[epa$Station.Code==s] <- epa.test$sOzoneH[epa.test$Station.Code==s]      
    }
  }
})
print(metrics[1:10,])
apply(metrics,2,function(x){round(mean(x),3)})
saveRDS(metrics, "output/lm.metrics.RDS")
saveRDS(epa, "output/lm.out.RDS")

# plot(epa$sOzone, type="l")
# lines(epa$sOzoneH, col=2)
# cor(epa$sOzone,epa$sOzoneH)
