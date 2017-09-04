## Clean the environment
rm(list=ls())

## Load libraries
library(randomForest)
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
## Using Random Forest
set.seed(204)
excludedCovariates1 <- c(1,2,6,7,19,20,22,23:28,30:31) 
excludedCovariates2 <- c(1,2,6,7,19,20,22,23,24,26:28,30:31) 
epaCore2 <- epa[epa$Station.Code %in% sample(unique(epa$Station.Code),50),-excludedCovariates2]
# head(epaCore2)
# dim(epaCore2)
# dim(epa)

## Train RF model
if(F){
  ticToc({
    rf.model <- randomForest(Ozone~., epaCore2, importance=TRUE)
  })
  saveRDS(rf.model, "output/preprocessing/importanceRF.RDS")
}
rf.model <- readRDS("output/preprocessing/importanceRF.RDS")

rownames(rf.model$importance)[4] <- "Water.Evap" 
rownames(rf.model$importance)[6] <- "Geop.Height" 
rownames(rf.model$importance)[7] <- "Geop.Height.Tropo" 
rownames(rf.model$importance)[8] <- "Lat.Heat.Flux" 
rownames(rf.model$importance)[9] <- "Tropo.Press"
printPlot(paper,"img/eda/analysis/importanceRf.jpeg",7,4,function(){
  varImpPlot(rf.model, main=NA, bg="blue", pch=22, cex=0.8, mar=c(0,0,0,0))
})
