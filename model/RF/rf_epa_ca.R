## Clean environment
rm(list=ls())

## Load libraries
library(randomForest)
library(gstat)
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
epa$sOzone <- scale(epa$Ozone)

## General ############################################################################
folds <- readRDS("output/folds.RDS")
fm <- sOzone ~ Temperature+RH+Rain+sqrtWind+UTM.X+UTM.Y+Elevation+Location.Setting+Doy+Dow.name+Dow.number+isWeekday

# ## Test (Fold-1) ####################################################################
# ## Split data for k=1
# folds <- readRDS("output/folds.RDS")
# ## Fold-1
# epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=1],] 
# epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==1],]
# epa.test$sOzoneH <- NA
# 
# ## Fit model 
# # fm <- sOzone ~ Temperature+RH+UTM.X+UTM.Y+Elevation+Location.Setting+j+wn # 0.8044
# # fm <- sOzone ~ Temperature+RH+UTM.X+UTM.Y+Elevation+Location.Setting+Doy+Dow
# fm <- sOzone ~ Temperature+RH+Rain+sqrtWind+UTM.X+UTM.Y+Elevation+Location.Setting+Doy+Dow.name+Dow.number+isWeekday
# 
# ticToc({
#   rf.fit <- randomForest(fm, epa.train[!is.na(epa.train$sOzone),], importance=TRUE)
# })
# saveRDS(rf.fit,file="output/RF/rf.fit.RDS")
# # rf.fit<-readRDS("output/RF/rf.fit.RDS")
# ## Runtime: 42.79235 mins!!!
# 
# ## Notes:
# ## http://trevorstephens.com/kaggle-titanic-tutorial/r-part-5-random-forests/
# # If you were working with a larger dataset you may want to reduce the number of
# # trees, at least for initial exploration, or restrict the complexity of each tree
# # using nodesize as well as reduce the number of rows sampled with sampsize
# 
# ## Plot covariate relevance
# printPlot(paper,"img/rf/rf_importance.jpeg",7,4,FUN=function(){
#   varImpPlot(rf.fit, main="Ozone (scaled)", bg="blue", pch=22)
# })
# 
# ## Prediction
# epa.test$sOzoneH <- predict(rf.fit, epa.test)
# ## Goodness-of-fit
# m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
# ## Print results
# print("Metrics:")
# print(m)
# ## Save results
# saveRDS(m, "output/space/metricRF.RDS")
# 
# 
# ## Check residuals
# epa.test$residuals <- epa.test$sOzone - epa.test$sOzoneH
# printPlot(T,"img/rf/rf_pred_vs_obs.jpeg",5,5,FUN=function(){
#   plot(sOzoneH~sOzone, epa.test, cex=0.5, xlab="Observations", ylab="Predictions")
#   abline(a=0,b=1,lty="dashed", col="green", lwd=2)
#   #abline(h=0,lty="dashed", col="gray"); abline(v=0,lty="dashed", col="gray")
# })
# printPlot(T,"img/rf/rf_res_vs_pred.jpeg",5,5,FUN=function(){
#   plot(residuals~sOzoneH,epa.test, cex=0.5, xlab="Residuals", ylab="Predictions")
#   abline(h=0, lty="dashed", col="gray")
# })
# 
# ## Create variogram of residuals
# epa.test.sp <- epa.test[!is.na(epa.test$residuals),]  # Ignore Nas. They create problems
# epa.test.sp <- convertDataToSp(epa.test.sp)
# v = variogram(residuals~UTM.X+UTM.Y, epa.test.sp)
# printPlot(T,"img/rf/rf_variogram_residuals.jpeg",5,5,FUN=function(){
#   print(plot(v))
# })
# # v = variogram(residuals~1, epa.test.sp)
# # plot(v)
# ## Notes:
# ## No apparently spatial correlation in residuals


# Too computationally expensive: 1 hour each fold
if(T){
  ## 10-fold cross validation
  metrics = c()
  epa$sOzoneH <- NA
  ticToc({
    for(k in 1:10){
      print(paste("Fold", k,paste(rep("=",50),collapse = "")))
  
      ## Split data
      epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],]
      epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
      #epa.test$sOzoneH <- NA  ## For the prediction
  
      ## Fit model
      rf.fit <- randomForest(fm, epa.train[!is.na(epa.train$sOzone),], importance=TRUE)
  
      ## Prediction
      epa.test$sOzoneH <- predict(rf.fit, epa.test)
      epa[epa$Station.Code %in% sites$Station.Code[folds==k],"sOzoneH"] <- epa.test$sOzoneH
  
      m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
      # plot(epa.test$sOzone,type="l")
      # lines(epa.test$OzoneTps,col=2)
      metrics <- rbind(metrics,m)
    }
  })
  ## Print results
  print("Metrics:")
  print(metrics)
  print("Cross-validation mean:")
  apply(metrics,2,mean)
  ## Save results
  saveRDS(metrics, "output/RF/rf.metrics.RDS")
  saveRDS(epa, "output/RF/rf.df.RDS")
}
