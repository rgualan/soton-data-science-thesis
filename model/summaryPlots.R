## Clean environment
#rm(list=ls())

## Libraries
library(ggplot2)
library(xtable)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Read data ##############################################################################
## Purely spatial
mIdw <- readRDS("output/space/idw.RDS")
mK1 <- readRDS("output/space/automap.RDS")
mK2 <- readRDS("output/space/kriging.RDS")
mTps <- readRDS("output/space/tps.RDS")
## ST kriging
mKst1 <- readRDS("output/spacetime/sepModel.metrics.RDS")
mKst2 <- readRDS("output/spacetime/psModel.metrics.RDS")
mKst3 <- readRDS("output/spacetime/metricModel.metrics.RDS")
mKst4 <- readRDS("output/spacetime/sumMetricModel.metrics.RDS")
## spTimer-GP
mGp <- readRDS("output/spTimer/out.metrics.RDS")
## Machine learning
mLm <- readRDS("output/lm.metrics.RDS")
mRf <- readRDS("output/RF/rfp.metrics.RDS")

## Comparison table ########################################################################
mainTable <-rbind(apply(mK1,2,mean),
                  apply(mK2,2,mean),
                  apply(mIdw,2,mean),
                  apply(mTps,2,mean),
                  apply(mKst1,2,mean),
                  apply(mKst2,2,mean),
                  apply(mKst3,2,mean),
                  apply(mKst4,2,mean),
                  apply(mGp,2,mean),
                  apply(mRf,2,mean),
                  apply(mRf,2,mean))
rownames(mainTable) <- c("KSP-A","KSP-B","IDW","TPS","KST-A","KST-B","KST-C","KST-D","GP","LM","RF")
round(mainTable,3)
## Latex table
xtable(mainTable, digits = 3)

## Boxplot RMSE by Model ##################################################################
rmse <- rbind(data.frame(Model="IDW", RMSE=mIdw[,2]),
              data.frame(Model="KSP-A", RMSE=mK1[,2]),
              data.frame(Model="KSP-B", RMSE=mK2[,2]),
              data.frame(Model="TPS", RMSE=mTps[,2]),
              data.frame(Model="KST-SE", RMSE=mKst1[,2]),
              data.frame(Model="KST-PS", RMSE=mKst2[,2]),
              data.frame(Model="KST-ME", RMSE=mKst3[,2]),
              data.frame(Model="KST-SM",  RMSE=mKst4[,2]),
              data.frame(Model="LM",  RMSE=mLm[,2]),
              data.frame(Model="RF", RMSE=mRf[,2]))
printPlot(paper, "img/boxplot_rmse.jpeg",7,5,FUN = function(){
  plot(RMSE~Model, rmse, xlab="", las=2)
})

## Boxplot R2 by Model ######################################################################
r2 <- rbind(data.frame(Model="IDW", R2=mIdw[,8]),
            data.frame(Model="KSP-A", R2=mK1[,8]),
            data.frame(Model="KSP-B", R2=mK2[,8]),
            data.frame(Model="TPS", R2=mTps[,8]),
            data.frame(Model="KST-SE", R2=mKst1[,8]),
            data.frame(Model="KST-PS", R2=mKst2[,8]),
            data.frame(Model="KST-ME", R2=mKst3[,8]),
            data.frame(Model="KST-SM",  R2=mKst4[,8]),
            data.frame(Model="LM",  R2=mLm[,8]),
            data.frame(Model="RF", R2=mRf[,8]))
printPlot(paper, "img/boxplot_r2.jpeg",7,5,FUN = function(){
  plot(R2~Model, r2, xlab="", las=2)
})


## Time series comparison ##################################################################
## Compare the original and modeled time series of 3 arbitrarily chosen stations
outLm <- readRDS("output/lm.out.RDS")
outOri <- outLm
outRfcv <- readRDS("output/RF/out.p.RDS")
outRf <- do.call("rbind", outRfcv)

## Standardize ozone column
outOri$Ozone <- outOri$sOzone
outLm$Ozone <- outLm$sOzoneH
outRf$Ozone <- outRf$sOzoneH

set.seed(204)
#ss<-c("021-0003", "027-0101", "083-4003" ) 
ss <- sample(unique(epa$Station.Code),3)
s=ss[1]

cc <- c("Station.Code","Date","Ozone") # common columns
sdata <- rbind(
  cbind(outOri[outOri$Station.Code %in% ss, cc], Model="Original"),
  cbind(outLm[outLm$Station.Code %in% ss, cc], Model="LM"),
  cbind(outRf[outRf$Station.Code %in% ss, cc], Model="RF")
)

for(s in ss){
  printPlot(paper, paste0("img/ts_",s,".jpeg"),7,3, FUN=function(){ 
    p<-ggplot(sdata[sdata$Station.Code==s,], aes(x=Date, y=Ozone, colour=Model)) + 
      geom_line() + 
      theme(legend.justification = c("top")) + 
      labs(y="Scaled(Ozone)")
    print(p)
  })
}







