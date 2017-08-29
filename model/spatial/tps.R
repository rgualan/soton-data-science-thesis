## Clean environment
rm(list=ls())

### Load required packages ###
library(fields)
library(gstat)
source("util/my_helper.R")

### Load data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
sites <- getSites(epa)
str(epa)
#View(epa)
## Feature engineering
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Scale variable 
epa$sOzone <- scale(epa$Ozone)


### Purely spatial model #####################################################################
folds <- readRDS("output/folds.RDS")
days = seq(min(epa$Date),max(epa$Date),by="days")

## 10-fold cross validation
metrics = c()
ticToc({
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
      predictors <- c("UTM.X","UTM.Y","Elevation")
      fit <- Tps(slice.train[,predictors],
                 slice.train[,"sOzone"]) # , scale.type="range"
      
      epa.test$sOzoneH[epa.test$Date==day] <-
        predict(fit,slice.test[,predictors])
    }
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
saveRDS(metrics, "output/space/tps.RDS")


