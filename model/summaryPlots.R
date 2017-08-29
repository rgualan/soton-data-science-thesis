## Clean environment
rm(list=ls())

## Libraries
library(ggplot2)
source("util/my_helper.R")

### Read data ##############################################################################
mK1 <- readRDS("output/space/automap.RDS")
mK2 <- readRDS("output/space/kriging.RDS")
mIdw <- readRDS("output/space/idw.RDS")
mTps <- readRDS("output/space/tps.RDS")

## Comparison table ########################################################################
tK1 <- apply(mK1,2,mean) 
tK2 <- apply(mK2,2,mean) 
tIdw <- apply(mIdw,2,mean) 
tTps <- apply(mTps,2,mean) 
rbind(tK1,tK2,tIdw,tTps)

## RMSE boxplot ############################################################################
rmse <- rbind(data.frame(Model="SpK-A", RMSE=mK1[,2]),
              data.frame(Model="SpK-B", RMSE=mK2[,2]),
              data.frame(Model="IDW", RMSE=mIdw[,2]),
              data.frame(Model="TPS", RMSE=mTps[,2]))
plot(RMSE~Model,rmse)

## R2 boxplot ############################################################################
r2 <- rbind(data.frame(Model="SpK-A", R2=mK1[,8]),
            data.frame(Model="SpK-B", R2=mK2[,8]),
            data.frame(Model="IDW", R2=mIdw[,8]),
            data.frame(Model="TPS", R2=mTps[,8]))
plot(R2~Model,r2)
