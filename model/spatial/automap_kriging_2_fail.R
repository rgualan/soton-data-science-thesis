## Clean environment
rm(list=ls())

### Load required packages ###
#library(fields)
#library(raster)
#library(spatial.tools)
#library(gdalUtils)
#library(rgdal)
#library(gstat)
library(automap)
source("util/my_helper.R")

### Read data ##############################################################################
epa <- readEpaDataset()
epa$sOzone <- scale(epa$Ozone) 
sites <- getSites(epa)

### Purely spatial model #####################################################################
#fm <- sOzone~Temperature+RH+Elevation
fm <- sOzone ~ Temperature+RH+logRain+sqrtWind+Dew.Point+Water.Evap+
  Geop.Height+Geop.Height.Tropo+Lat.Heat.Flux+
  Tropo.Press+Press.MSL+Vegetation+Elevation+
  Location.Setting+Doy+Dow.name+Month+Day
folds <- getFolds()
days = getDates()

## 10-fold cross validation
metrics = c()
ticToc({
  for(k in 1:1){
    print(paste("Fold", k,paste(rep("=",50),collapse = "")))
    
    epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],] 
    epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
    epa.test$sOzoneH <- NA
    
    coordinates(epa.train) <- ~UTM.X+UTM.Y
    coordinates(epa.test) <- ~UTM.X+UTM.Y
    proj4string(epa.train) <- getUTMproj()
    proj4string(epa.test) <- getUTMproj()

    out = autoKrige(fm, epa.train, epa.test)
    
    # ## Diagnostic
    # plot(out)
    # ## Test plots
    # plot(slice.test$sOzone,type="l")
    # lines(out$krige_output$var1.pred,col=2)
    # ## Cor
    # cor(slice.test$sOzone,out$krige_output$var1.pred,use="pairwise.complete.obs")
      
    epa.test$sOzoneH <- out$krige_output$var1.pred

    m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
    metrics <- rbind(metrics,m)
  }
})

## Print results
print("Metrics:")
print(metrics)
print("Cross-validation mean:")
apply(metrics,2,mean)
## Save results
saveRDS(metrics, "output/space/automap.RDS")


## Notes:
## It does not tolerate missing values