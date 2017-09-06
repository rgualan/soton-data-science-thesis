## Ref: 
## Spatio-Temporal data in R
## https://www.r-bloggers.com/spatio-temporal-kriging-in-r/
## https://cran.r-project.org/web/packages/gstat/vignettes/st.pdf

## Clean environment
rm(list=ls())

## Libraries
source("util/my_helper.R")
library(gstat)
library(spacetime)
#library(rgdal) ## Not available in the cluster!
library(reshape2)
library(parallel)
source("model/spacetime/testStKriging.R")
source("model/spacetime/fitCovarianceModels.R")

## Global variables ##########################################################
forceRun <- T # Run heavy weight processing
debugLevel <- F # Present additional debugging information
paper <- setupPaper() # Print the plots in image format for the paper

## Read data #################################################################
epa_version <- 1
epa <- readEpaDataset(epa_version)
epa <- addDateDerivedFeatures(epa)
epa <- transformFeatures(epa)
epa <- scaleTargetVariable(epa)
sites <- getSites(epa)

## Detrend the target variable ###############################################
#epa <- detrend_scaled_ozone_lm(epa)
epa <- detrend_scaled_ozone_poly(epa)

## Assemble STFDF ############################################################
epa.st <- assembleSTFDF(epa) 
# stplot(epa.st[,"2016-01-01::2016-01-08","sOzone"])
# dim(epa.st)

## Covariance models #########################################################
#fm <- sOzone~Elevation+Temperature
# fm <- sOzone ~ Temperature+RH+logRain+sqrtWind+Dew.Point+Water.Evap+
#   Geop.Height+Geop.Height.Tropo+Lat.Heat.Flux+
#   Tropo.Press+Press.MSL+Vegetation
fm <- sOzone.res~1
tlags <- 0:10
#tlags <- c(0,5,10,15)

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
runCluster <- T
k <- 10
folds <- getFolds(epa_version)
modelNames <- c("sepModel", "psModel", "metricModel", "sumMetricModel")
models <- list(fitSepModel, fitMetricModel, fitProdSumModel, fitSumMetricModel)

if(runCluster){
  cl <- makeCluster(6, outfile="")
  clusterExport(cl, c("epa.st","paper","k","folds","fm","modelNames","models","linStAni"))
  tryCatch({
    out <- clusterApply(cl, 1:length(models), function(idx){
      source("model/spacetime/testStKriging.R")
      print(sprintf("Process %d: %s",idx,modelNames[idx]))
      epa.st <- runCV(k, folds, fm, epa.st, modelNames[idx], models[[idx]], linStAni, nmax=50)
      return(epa.st)    
    })
  }, error=function(e){message(e)}, finally={stopCluster(cl)})
}else{
  out <- list()
  for(idx in 1:length(models)){
    print(modelNames[idx])
    out[[idx]] <- runCV(k, folds, fm, epa.st, modelNames[idx], models[[idx]], linStAni, nmax=50)
  }
}
saveRDS(out,"output/spacetime/st.out.RDS")
#out <- readRDS("output/spacetime/st.out.RDS")
#head(out[[1]]@data)

## CV evaluation metrics ########################################################################
for(idx in 1:length(modelNames)){
  print(modelNames[idx])
  print(evaluatePredictions(
    out[[idx]][,,modelNames[idx],drop=F]@data[[1]],
    out[[idx]][,,"sOzone.res",drop=F]@data[[1]]))
}

## ST to data.frame ########################################################################
out.df <- as.data.frame(out[[1]])
head(out.df[order(out.df$Station.Code,out.df$Date),])
for(idx in 2:length(modelNames)){
  a <- as.data.frame( out[[idx]] )
  out.df <- merge(out.df,a[,c("Station.Code","Date",modelNames[idx])])
}
head(out.df)
for(idx in 1:length(modelNames)){
  out.df[sprintf("sOzoneH.%d",idx)] <- out.df$sOzone.trend + out.df[modelNames[idx]][,1]
}
saveRDS(out.df,"output/spacetime/st.out.df.RDS")
#out.df <- readRDS("output/spacetime/st.out.df.RDS")
head(out.df)
plot(out.df$sOzone[1:1000],type="l")
lines(out.df$sOzoneH.1[1:1000],col=2)
lines(out.df$sOzoneH.2[1:1000],col=3)
lines(out.df$sOzoneH.3[1:1000],col=4)
lines(out.df$sOzoneH.4[1:1000],col=5)


## Assess metrics ###############################################################################
for(idx in 1:length(modelNames)){
  #print(idx)
  metrics <- getMetricsByStationFromDF(out.df,"sOzone",paste0("sOzoneH.",idx))
  saveRDS(metrics, paste0("output/spacetime/",modelNames[idx],".metrics.RDS"))
}

print("Done")