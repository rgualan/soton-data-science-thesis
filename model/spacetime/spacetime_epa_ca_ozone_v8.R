## Ref: 
## Spatio-Temporal data in R
## https://www.r-bloggers.com/spatio-temporal-kriging-in-r/
## https://cran.r-project.org/web/packages/gstat/vignettes/st.pdf

## Clean environment
rm(list=ls())

## Libraries
library(gstat)
library(spacetime)
#library(rgdal) ## Not available in the cluster!
library(reshape2)
library(parallel)
source("util/my_helper.R")
source("model/spacetime/testStKriging.R")
source("model/spacetime/fitCovarianceModels.R")

## Global variables ##########################################################
forceRun <- T # Run heavy weight processing
debugLevel <- F # Present additional debugging information
paper <- setupPaper() # Print the plots in image format for the paper

## Read data #################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Sites
sites <- getSites(epa)
## Scale target variable #####################################################
epa$sOzone <- scale(epa$Ozone)[,1]  ## scale returns a matrix

## Assemble STFDF ############################################################
epa.st <- assembleSTFDF(epa) 
# stplot(epa.st[,"2016-01-01::2016-01-08","sOzone"])
# dim(epa.st)

## Covariance models #########################################################
fm <- sOzone~Elevation+Temperature
#tlags <- 0:7
tlags <- c(0,5,10,15)

ticToc({
  covModels <- fitCovarianceModels(epa.st, fm, tlags, paper)
})
fitSepModel <- covModels$separable
fitProdSumModel <- covModels$prodSum
fitMetricModel <- covModels$metric
fitSumMetricModel <- covModels$sumMetric
linStAni <- covModels$linStAni


## Run 10-fold CV in parallel ###############################################################
## One CPU per model
k <- 10
folds <- getFolds()
modelNames <- c("sepModel", "psModel", "metricModel", "sumMetricModel")
models <- list(fitSepModel, fitMetricModel, fitProdSumModel, fitSumMetricModel)

cl <- makeCluster(5, outfile="")
clusterExport(cl, c("epa.st","paper","k","folds","fm","modelNames","models","linStAni"))
tryCatch({
  out <- clusterApply(cl, 1:length(models), function(idx){
    library(gstat)
    library(spacetime)
    source("model/spacetime/testStKriging.R")
    print(sprintf("Process %d: %s",idx,modelNames[idx]))
    epa.st <- runFold(k, folds, fm, epa.st, modelNames[idx], models[[idx]], linStAni, nmax=50)
    return(epa.st)    
  })
}, error = function(e){
  print(e)
}, finally = {
  stopCluster(cl)
})

saveRDS(out,"output/spacetime/10cv.out.RDS")

## CV evaluation metrics ########################################################################
for(idx in 1:length(modelNames)){
  print(idx)
  print(evaluatePredictions(out[[idx]][,,modelNames[idx],drop=F]@data[[1]],
                      out[[idx]][,,"sOzone",drop=F]@data[[1]]))
}

## Assess metrics ###############################################################################
#out <- readRDS("output/spacetime/10cv.out.RDS")
for(idx in 1:length(modelNames)){
  #print(idx)
  metrics <- getMetricsByStation(out[[idx]],"sOzone",modelNames[idx])
  saveRDS(metrics, paste0("output/spacetime/",modelNames[idx],".metrics.RDS"))
}

print("Done")
