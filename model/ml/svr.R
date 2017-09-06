## Clean environment
rm(list=ls())

## Libraries
library(e1071)
library(parallel)
source("util/my_helper.R")

## Global variables
paper = setupPaper()


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

## 10-fold cross validation #################################################################################
folds <- getFolds()
cl <- makeCluster(11, outfile="")
clusterExport(cl, c("epa","sites","paper","fm"))
tryCatch({
  out <- clusterApply(cl, 1:10, function(k){
    library(e1071)
    source("util/my_helper.R")
    
    folds <- getFolds()
    print(paste("Fold",k,paste(rep("=",50),collapse = "")))
    
    ## Split data
    epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],]
    epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
    epa.test$sOzoneH <- NA  ## For the prediction
    
    ## Fit model
    ticToc({
      model <- svm(fm, epa.train)
    })
    
    ## Prediction
    a <- predict(model, epa.train)
    epa.test$sOzoneH <- predict(model, epa.test)

    print(sprintf("Fold-%d. Metrics:",k))
    print (evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH,2))
    
    return(epa.test)
  })
}, error=function(e){message(e)}, finally={stopCluster(cl)})

## Save results
epa.out <- do.call("rbind", out)
saveRDS(epa.out, "output/svr/svr.out.RDS")
#out <- readRDS("output/svr/svr.out.RDS")

## Overall metrics
print(evaluatePredictions(epa.out$sOzone,epa.out$sOzoneH,2))

## Metrics
metrics <- getMetricsByStationFromDF(epa.out)
saveRDS(metrics, "output/svr/svr.metrics.RDS")

print("Done!")
