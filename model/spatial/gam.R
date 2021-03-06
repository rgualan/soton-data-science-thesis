## Clean environment
rm(list=ls())

### Load required packages ###
library(fields)
library(raster)
library(spatial.tools)
library(gdalUtils)
library(rgdal)
library(gstat)
library(automap)
source("util/my_helper.R")

### Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
sites <- getSites(epa)
## Scale dependent variable
epa$sOzone <- scale(epa$Ozone) 

### Purely spatial model #####################################################################
folds <- getFolds()
days = seq(min(epa$Date),max(epa$Date),by="days")

## 10-fold cross validation
metrics = c()
ticToc({
  for(k in 1:2){
    print(paste("Fold", k,paste(rep("=",50),collapse = "")))
    
    epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],] 
    epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
    epa.test$sOzoneH <- NA
    
    for(i in 1:length(days)){
      day <- days[i]
      #print(day)
      slice.train <- epa.train[epa.train$Date==day & !is.na(epa.train$sOzone),]  
      slice.test <- epa.test[epa.test$Date==day,]
      # coordinates(slice.train) <- ~UTM.X+UTM.Y
      # coordinates(slice.test) <- ~UTM.X+UTM.Y
      # proj4string(slice.train) <- getUTMproj()
      # proj4string(slice.test) <- getUTMproj()
  
      model <- gam(sOzone ~ s(Temperature) + s(RH) + s(Elevation) + s(UTM.X, UTM.Y, k = 10), 
                   data = slice.train)
      ## TODO k=???
      
      # ## Diagnostic
      # plot(out)
      # ## Test plots
      # plot(slice.test$sOzone,type="l")
      # lines(out$krige_output$var1.pred,col=2)
      # ## Cor
      # cor(slice.test$sOzone,out$krige_output$var1.pred,use="pairwise.complete.obs")
      
      epa.test$sOzoneH[epa.test$Date==day] <- out$krige_output$var1.pred
    }
    m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
    # plot(epa.test$sOzone,type="l")
    # lines(epa.test$OzoneTps,col=2)
    metrics <- rbind(metrics,m)
  }
})
print("Cross-validation results:")
apply(metrics,2,mean)

## Notes:
## It does not tolerate missing values