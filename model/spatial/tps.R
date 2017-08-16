## Clean environment
rm(list=ls())

### Load required packages ###
library(fields)
#library(raster)
#library(spatial.tools)
#library(gdalUtils)
#library(rgdal)
library(gstat)
source("util/my_helper.R")

### Use more computing resources ###
# sfQuickInit()

### Load data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
sites <- getSites(epa)
str(epa)
#View(epa)
## Feature engineering
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Scale variable 
epa$sOzone <- scale(epa$Ozone)
# hist(epa$Ozone)
# hist(epa$sOzone)


### Purely spatial model #####################################################################
folds <- readRDS("data/tmp/folds.RDS")
days = seq(min(epa$Date),max(epa$Date),by="days")

## 10-fold cross validation
metrics = c()
st <- Sys.time()
for(k in 1:10){
  print(paste("Fold", k,paste(rep("=",50),collapse = "")))
  
  epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k]
                   #& epa$Date<="2016-01-31"
                   ,] 
  epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k]
                  #& epa$Date<="2016-01-31"
                  ,]
  
  epa.test$sOzoneH <- NA
  for(i in 1:length(days)){
    day <- days[i]
    #print(day)
    slice.train <- epa.train[epa.train$Date==day,]
    slice.test <- epa.test[epa.test$Date==day,]
    fit <- Tps(slice.train[,c("UTM.X","UTM.Y","Elevation")],
               slice.train[,"sOzone"]) # , scale.type="range"

    # ## Output diagnostic plots ###
    # set.panel(2,2); plot(fit); set.panel(1,1)
    # surface(fit)
    # ## Train plots
    # plot(slice.train[,"sOzone"],type="l")
    # lines(fit$fitted.values, col=2)
    # ## Test plots
    # yh <- predict(fit, slice.test[,c("UTM.X","UTM.Y","Elevation")])
    # plot(slice.test$sOzone,type="l")
    # lines(yh,col=2)
    # ## Cor
    # cor(epa.test$sOzone,yh,use="pairwise.complete.obs")
    
    epa.test$sOzoneH[epa.test$Date==day] <-
      predict(fit,slice.test[,c("UTM.X","UTM.Y","Elevation")])
  }
  m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
  # plot(epa.test$sOzone,type="l")
  # lines(epa.test$OzoneTps,col=2)
  metrics <- rbind(metrics,m)
}
print("Cross-validation results:")
apply(metrics,2,mean)
print("Running time:")
Sys.time()-st

