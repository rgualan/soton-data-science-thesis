## Clean environment
rm(list=ls())

### Load required packages ###
library(fields)
library(gstat)
source("util/my_helper.R")

### Load data ##############################################################################
epa <- readEpaDataset()
epa <- addDateDerivedFeatures(epa)
epa <- scaleTargetVariable(epa)
sites <- getSites(epa)

### Purely spatial model #####################################################################
folds <- getFolds()
days = getDates(dateOutput = F)

## 10-fold cross validation
metrics = c()
epa.out <- data.frame()
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
    epa.out <- rbind(epa.out,epa.test)
    # plot(epa.test$sOzone,type="l")
    # lines(epa.test$OzoneTps,col=2)
  }
})

metrics <- getMetricsByStationFromDF(epa.out)

## Print results
print("Metrics:")
print(metrics[1:10,])
print("Cross-validation mean:")
apply(metrics,2,mean)
## Save results
saveRDS(metrics, "output/space/tps.metrics.RDS")
saveRDS(epa.out, "output/space/tps.out.RDS")

