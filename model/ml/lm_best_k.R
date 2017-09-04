## Clean environment
rm(list=ls())

## Libraries
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Read data ########################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- addDateDerivedFeatures(epa)
epa <- transformFeatures(epa)
epa$sOzone <- scale(epa$Ozone)[,1]
sites <- getSites(epa)

## Coss-validation function #########################################################
cvByK <- function(epa, k_neighbors, fm){
  ## Add combination of closest neighbours as a predictor
  epa <- addNeighboursAverage(epa,k_neighbors)

  ## 10-fold cross validation 
  folds <- getFolds()
  metrics = c()
  ticToc({
    for(k in 1:10){
      #print(paste("Fold",k,paste(rep("=",50),collapse = "")))
      epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],] 
      epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
      
      epa.test$sOzoneH <- NA
      lmModel <- lm(fm, data=epa.train)
      epa.test$sOzoneH <- predict(lmModel,epa.test)
      m <- getMetricsByStationFromDF(epa.test, "sOzone", "sOzoneH")
      metrics <- round(rbind(metrics,m),3)
    }
  })
  
  return(metrics)
}


## Find best k number of neighbors ######################################################
fm <- sOzone ~ Temperature+Dew.Point+Water.Evaporation+
  Geopotential.Height+Geopotential.Height.Tropo+#Heat.Flux
  Tropopause.Press+Press.MSL+Longitude+Latitude+Elevation+
  Location.Setting+Doy+Neighbor

results <- data.frame()
for(k in c(4:9, seq(10,80,by=10))){
  m <- cvByK(epa,k,fm)
  results <- rbind(results, data.frame(k=k, RMSE=m[,2], R2=m[,8]))
}

## Box-plots
results$k <- as.factor(results$k)

printPlot(paper, "img/lm/bestk_boxplot_rmse.jpeg", 5, 4, function(){
  par(mar=c(4.5,4,1,1))
  plot(RMSE~k, results, xlab="Number of neighbors (k)")
})

printPlot(paper, "img/lm/bestk_boxplot_r2.jpeg", 5, 4, function(){
  par(mar=c(4.5,4,1,1))
  plot(R2~k, results, xlab="Number of neighbors (k)")
})

(aggRes <- aggregate(cbind(RMSE,R2)~k,results,mean))
aggRes[which(aggRes$RMSE == min(aggRes$RMSE)),]
aggRes[which(aggRes$R2 == max(aggRes$R2)),]

