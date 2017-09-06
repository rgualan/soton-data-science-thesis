## Clean environment
rm(list=ls())

## Libraries
library(e1071)
source("util/my_helper.R")

## Read data ########################################################################
epa <- readEpaDataset()
epa <- addDateDerivedFeatures(epa)
epa <- addNeighboursAverage(epa,5)
epa <- transformFeatures(epa)
epa <- scaleTargetVariable(epa)
sites <- getSites(epa)

## Formula #######################################################################
fm <- sOzone ~ Temperature+Dew.Point+Water.Evap+
  Geop.Height+Geop.Height.Tropo+
  Tropo.Press+Press.MSL+Longitude+Latitude+Elevation+
  Location.Setting+Doy+Neighbor

## Grid search #################################################################################
## Ref: https://www.svm-tutorial.com/2014/10/support-vector-regression-r/

## Reduced dataset
epa <- epa[epa$Station.Code %in% sample(unique(epa$Station.Code),50), ]

## perform a grid search
gridSearch <- tune(svm, fm,  data = epa,
                   ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9))
)
saveRDS(gridSearch,"output/svr/gridSearch.RDS")
print(tuneResult)
plot(tuneResult)
