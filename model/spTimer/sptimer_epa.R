## Clean the environment
rm(list=ls())

## Load libraries
library(spTimer)
library(parallel)
source("util/my_helper.R")

## Global variables ###############################################################
paper <- setupPaper()

## Read data #######################################################################
epa_version <- 1
epa <- readEpaDataset(epa_version)
epa <- addDateDerivedFeatures(epa)
epa <- transformFeatures(epa)
epa <- scaleTargetVariable(epa)
sites <- getSites(epa)

## Detrend the target variable ###############################################
#epa <- detrend_scaled_ozone_lm(epa)
epa <- detrend_scaled_ozone_poly(epa)

## 10-fold cross validation #################################################################################
folds <- getFolds(epa_version)
fm <- sOzone.res ~ Temperature+Dew.Point+Water.Evap+Heat.Flux+
  Lat.Heat.Flux+Vegetation+Doy+Month

## Function that runs a fold of a CV exercise
runFold <- function(k, folds, epa){
  library(spTimer)
  source("util/my_helper.R")

  print(paste("Fold",k))
  
  ## Split data
  sites <- getSites(epa)
  epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],]
  epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
  epa.test$sOzone.res.H <- NA  ## For the prediction
  
  ## Fit model
  time.data<-spT.time(t.series=366,segments=1)
  ticToc({
    post.gp <- spT.Gibbs(formula = fm, data = epa.train, model = "GP", 
                         coords = ~UTM.X+UTM.Y, time.data=time.data,
                         newdata=epa.test, newcoords=~UTM.X+UTM.Y,
                         spatial.decay = spT.decay(distribution = Gamm(2,1), tuning = 0.1),
                         report=2)
  })
  #print(post.gp)
  summary(post.gp)
  
  ## Prediction
  epa.test$sOzone.res.H <- c(post.gp$prediction$Mean)
  ## Fold metrics
  print(sprintf("Fold-%d. Metrics:",k))
  m <- evaluatePredictions(epa.test$sOzone.res, epa.test$sOzone.res.H)
  print (round(m,2))
  
  ## Output
  ## The next line was a naive mistake 
  #epa$sOzone.res.H[epa$Station.Code %in% sites$Station.Code[folds==k]] <- epa.test$sOzone.res.H
  return(epa.test)
}

## Run Cross-validation 
out <- list()
epa$sOzone.res.H <- NA
for(k in 1:10){ 
  out[[k]] <- runFold(k, folds, epa)
}
epa.out <- do.call("rbind", out) ## Consolidate
## Save results
epa.out$sOzoneH <- epa.out$sOzone.trend+epa.out$sOzone.res.H
saveRDS(epa.out, "output/spTimer/gp.out.RDS")
#epa.out <- readRDS("output/spTimer/gp.out.RDS")

## Calculate metrics ##############################################################
## Metrics
metrics <- getMetricsByStationFromDF(epa.out,"sOzone","sOzoneH")
saveRDS(metrics, "output/spTimer/gp.metrics.RDS")
#metrics <- readRDS("output/spTimer/gp.metrics.RDS")
apply(metrics,2,function(x){round(mean(x),2)})
