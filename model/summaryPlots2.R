## Clean environment
#rm(list=ls())

## Libraries
library(ggplot2)
source("util/my_helper.R")

## Read data ##############################################################################
## Purely spatial
mK1 <- readRDS("output/space/automap.RDS")
mK2 <- readRDS("output/space/kriging.RDS")
mIdw <- readRDS("output/space/idw.RDS")
mTps <- readRDS("output/space/tps.RDS")
## ST kriging
mKst1 <- readRDS("output/spacetime/sepModel.metrics.RDS")
mKst2 <- readRDS("output/spacetime/psModel.metrics.RDS")
mKst3 <- readRDS("output/spacetime/metricModel.metrics.RDS")
mKst4 <- readRDS("output/spacetime/sumMetricModel.metrics.RDS")
## spTimer-GP
mGp <- readRDS("output/spTimer/out.metrics.RDS")
## RF
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
                  apply(mRf,2,mean))
rownames(mainTable) <- c("KSP-A","KSP-B","IDW","TPS","KST-A","KST-B","KST-C","KST-D","GP","RF")
round(mainTable,3)

## Boxplot RMSE by Model ##################################################################
rmse <- rbind(data.frame(Model="KSP-A", RMSE=mK1[,2]),
              data.frame(Model="KSP-B", RMSE=mK2[,2]),
              data.frame(Model="IDW", RMSE=mIdw[,2]),
              data.frame(Model="TPS", RMSE=mTps[,2]),
              data.frame(Model="KST-SE", RMSE=mKst1[,2]),
              data.frame(Model="KST-PS", RMSE=mKst2[,2]),
              data.frame(Model="KST-ME", RMSE=mKst3[,2]),
              data.frame(Model="KST-SM",  RMSE=mKst4[,2]),
              data.frame(Model="RF", RMSE=mRf[,2]))
plot(RMSE~Model, rmse)

## Boxplot R2 by Model ######################################################################
r2 <- rbind(data.frame(Model="KSP-A", R2=mK1[,8]),
              data.frame(Model="KSP-B", R2=mK2[,8]),
              data.frame(Model="IDW", R2=mIdw[,8]),
              data.frame(Model="TPS", R2=mTps[,8]),
              data.frame(Model="KST-SE", R2=mKst1[,8]),
              data.frame(Model="KST-PS", R2=mKst2[,8]),
              data.frame(Model="KST-ME", R2=mKst3[,8]),
              data.frame(Model="KST-SM",  R2=mKst4[,8]),
              data.frame(Model="RF", R2=mRf[,8]))
plot(R2~Model, r2)
