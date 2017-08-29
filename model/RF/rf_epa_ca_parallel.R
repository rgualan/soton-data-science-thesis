## Clean environment
rm(list=ls())

## Load libraries
library(randomForest)
library(gstat)
library(parallel)
source("util/my_helper.R")

## Global variables
paper = setupPaper()

## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]
sites <- getSites(epa)
## Feature engineering (time dimension)
epa <- addDoyField(epa)
epa <- addDowField(epa)
epa <- addIsWeekDay(epa)
## Standardize variable 
epa$sOzone <- scale(epa$Ozone)[,1]

## Notes
## A single fold takes 42.79235 mins!

## Notes:
## http://trevorstephens.com/kaggle-titanic-tutorial/r-part-5-random-forests/
# If you were working with a larger dataset you may want to reduce the number of
# trees, at least for initial exploration, or restrict the complexity of each tree
# using nodesize as well as reduce the number of rows sampled with sampsize

fm <- sOzone ~ Temperature+RH+Rain+sqrtWind+UTM.X+UTM.Y+Elevation+Location.Setting+Doy+Dow.name+Dow.number+isWeekday
folds <- readRDS("output/folds.RDS")
metrics = c()

## 10-fold cross validation #################################################################################
cl <- makeCluster(11, outfile="")
clusterExport(cl, c("epa","sites","paper","fm","folds","metrics","ticToc"))
tryCatch({
  out <- clusterApply(cl, 1:10, function(k){
    library(randomForest)
    source("util/my_helper.R")

    print(paste("Fold",k,paste(rep("=",50),collapse = "")))

    ## Split data
    epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],]
    epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
    epa.test$sOzoneH <- NA  ## For the prediction
    
    ## Fit model
    ticToc({
      rf.fit <- randomForest(fm, epa.train[!is.na(epa.train$sOzone),], importance=TRUE)
      #saveRDS(rf.fit,paste0("output/RF/rf.fit.",k,".RDS"))
    })
    
    ## Prediction
    epa.test$sOzoneH <- predict(rf.fit, epa.test)
    # epa.test$sOzoneH <- rnorm(nrow(epa.test))
    
    m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
    print(sprintf("Fold-%d. Metrics:",k))
    print (m)

    return(epa.test)
  })
  
}, error = function(e){
  print(e)
}, finally = {
  stopCluster(cl)
})

## Save results
saveRDS(out, "output/RF/out.p.RDS")
#out <- readRDS("output/RF/out.p.RDS")

## Calculate metrics
epa.out <- do.call("rbind", out)
## Overall results
# print("Overall results:")
# print(evaluatePredictions(epa.out$sOzone,epa.out$sOzoneH))
## Metrics
metrics <- getMetricsByStationFromDF(epa.out,"sOzone","sOzoneH")
saveRDS(metrics, "output/RF/rfp.metrics.RDS")

print("Done!")
