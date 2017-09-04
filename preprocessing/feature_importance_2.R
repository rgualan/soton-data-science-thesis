## Clean the environment
rm(list=ls())

## Load libraries
#library(randomForest)
library(Boruta)
source("util/my_helper.R")

## Global variables ###############################################################
paper <- setupPaper()

## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
sites <- getSites(epa)

## Feature engineering #############################################################
## Features derived from the Date field
epa$Month <- as.numeric(format(epa$Date, "%m"))
epa$Day <- as.numeric(format(epa$Date, "%d"))
epa <- addDoyField(epa)
epa <- addDowField(epa, sinTx = T)
epa$sqrtWind <- sqrt(epa$Wind)
epa$logRain <- log(1+epa$Rain, 10) # The additive term is to avoid -inf
epa$levelRain <- round(epa$logRain) # 6 levels derived from the log trasnformation
epa$rained <- epa$Rain>0 # Binary variable


## Rank Features By Importance ###################################################
set.seed(204)
excludedCovariates1 <- c(1,2,6,7,19,20,22,23:28,30:31)
excludedCovariates2 <- c(1,2,6,7,19,20,22,23,24,26:28,30:31)  
epaCore2 <- epa[epa$Station.Code %in% sample(unique(epa$Station.Code),50),-excludedCovariates2]
# head(epaCore2)
# dim(epaCore2)
# dim(epa)

if(F){
  ticToc({
    borImp <- Boruta(Ozone~., data = epaCore2, doTrace = 2)
  })
  saveRDS(borImp, "output/preprocessing/importanceBoruta.RDS")
}
borImp <- readRDS("output/preprocessing/importanceBoruta.RDS")
print(borImp)

printPlot(paper,"img/eda/analysis/importanceBoruta.jpeg",6,5,function(){
  par(mar=c(8,4,1,1))
  plot(borImp, xlab="", las=2, cex.axis=0.7)
})


