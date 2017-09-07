## Clean environment
#rm(list=ls())

## Libraries
library(ggplot2)
library(xtable)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Summary table ##############################################################################
## Purely spatial
mIdw <- readRDS("output/space/idw.metrics.RDS")
mK1 <- readRDS("output/space/automap.metrics.RDS")
mK2 <- readRDS("output/space/kriging.metrics.RDS")
mTps <- readRDS("output/space/tps.metrics.RDS")
## ST kriging
mKst1 <- readRDS("output/spacetime/sepModel.metrics.RDS")
mKst2 <- readRDS("output/spacetime/psModel.metrics.RDS")
mKst3 <- readRDS("output/spacetime/metricModel.metrics.RDS")
mKst4 <- readRDS("output/spacetime/sumMetricModel.metrics.RDS")
## spTimer-GP
mGp <- readRDS("output/spTimer/gp.metrics.RDS")
## Machine learning
mLm <- readRDS("output/lm.metrics.RDS")
mRf <- readRDS("output/RF/rf.metrics.RDS")
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
                  apply(mLm,2,mean),
                  apply(mRf,2,mean))
rownames(mainTable) <- c("KSP-A","KSP-B","IDW","TPS","KST-SE","KST-PS","KST-ME","KST-SM","GP","LM","RF")
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
              data.frame(Model="ST-GP",  RMSE=mGp[,2]),
              data.frame(Model="LM",  RMSE=mLm[,2]),
              data.frame(Model="RF", RMSE=mRf[,2]))
printPlot(paper, "img/boxplot_rmse.jpeg",7,5,FUN = function(){
  plot(RMSE~Model, rmse, xlab="", las=2)
})

## Boxplot R2 by Model ######################################################################
r2 <- rbind(data.frame(Model="IDW", R2=mIdw[,5]),
            data.frame(Model="KSP-A", R2=mK1[,5]),
            data.frame(Model="KSP-B", R2=mK2[,5]),
            data.frame(Model="TPS", R2=mTps[,5]),
            data.frame(Model="KST-SE", R2=mKst1[,5]),
            data.frame(Model="KST-PS", R2=mKst2[,5]),
            data.frame(Model="KST-ME", R2=mKst3[,5]),
            data.frame(Model="KST-SM",  R2=mKst4[,5]),
            data.frame(Model="ST-GP",  R2=mGp[,5]),
            data.frame(Model="LM",  R2=mLm[,5]),
            data.frame(Model="RF", R2=mRf[,5]))
printPlot(paper, "img/boxplot_r2.jpeg",7,5,FUN = function(){
  plot(R2~Model, r2, xlab="", las=2)
})


## Time series comparison ##################################################################
## Compare the original and modeled time series of 3 arbitrarily chosen stations
outSTK <- readRDS("output/spacetime/st.out.df.RDS")
outSTK1 <- outSTK2 <- outSTK3 <- outSTK4 <- outSTK
#head(outSTK)
outLm <- readRDS("output/lm.out.RDS")
outOri <- outLm
outIdw <- readRDS("output/space/idw.out.RDS")
outK1 <- readRDS("output/space/automap.out.RDS")
outK2 <- readRDS("output/space/kriging.out.RDS")
outTps <- readRDS("output/space/tps.out.RDS")
outRf <- readRDS("output/RF/epa.out.RDS")
outGp <- readRDS("output/spTimer/gp.out.RDS")

## Standardize ozone column
outOri$Ozone <- outOri$sOzone
outIdw$Ozone <- outIdw$sOzoneH
outK1$Ozone <- outK1$sOzoneH
outK2$Ozone <- outK2$sOzoneH
outTps$Ozone <- outTps$sOzoneH
outSTK1$Ozone <- outSTK1$sOzoneH.1
outSTK2$Ozone <- outSTK2$sOzoneH.2
outSTK3$Ozone <- outSTK3$sOzoneH.3
outSTK4$Ozone <- outSTK4$sOzoneH.4
outGp$Ozone <- outGp$sOzoneH
outLm$Ozone <- outLm$sOzoneH
outRf$Ozone <- outRf$sOzoneH

#set.seed(204)
#ss<-c("021-0003", "027-0101", "083-4003" ) 
#(ss <- sample(unique(outGp$Station.Code),3))
ss <- c("021-0003", "071-2002", "001-0009")
s=ss[1]

cc <- c("Station.Code","Date","Ozone") # common columns
sdata <- rbind(
  cbind(outOri[outOri$Station.Code %in% ss, cc], Model="Original"),
  cbind(outIdw[outIdw$Station.Code %in% ss, cc], Model="IDW"),
  cbind(outK1[outK1$Station.Code %in% ss, cc], Model="KSP-A"),
  cbind(outK2[outK2$Station.Code %in% ss, cc], Model="KSP-B"),
  cbind(outTps[outTps$Station.Code %in% ss, cc], Model="TPS"),
  cbind(outSTK1[outSTK1$Station.Code %in% ss, cc], Model="KST-SE"),
  cbind(outSTK2[outSTK2$Station.Code %in% ss, cc], Model="KST-PS"),
  cbind(outSTK3[outSTK3$Station.Code %in% ss, cc], Model="KST-ME"),
  cbind(outSTK4[outSTK4$Station.Code %in% ss, cc], Model="KST-SM"),
  cbind(outGp[outGp$Station.Code %in% ss, cc], Model="GP"),
  cbind(outLm[outLm$Station.Code %in% ss, cc], Model="LM"),
  cbind(outRf[outRf$Station.Code %in% ss, cc], Model="RF")
)

for(s in ss){
  printPlot(paper, paste0("img/ts_",s,".jpeg"),7,3, FUN=function(){ 
    p<-ggplot(sdata[sdata$Station.Code==s,], aes(x=Date, y=Ozone, colour=Model)) + 
      geom_line(lwd=0.4, alpha=0.4) + 
      theme(legend.justification = c("top")) + 
      labs(y="Scaled(Ozone)")+
      geom_line(data=sdata[sdata$Station.Code==s & sdata$Model=="Original",], 
                aes(x=Date, y=Ozone, colour=Model)) + 
      theme(legend.justification = c("top")) + 
      labs(y="Scaled(Ozone)")
    
    print(p)
  })
}



## Hexbin plots
require(hexbin)
printPlot(paper,"img/hexbin/Idw.jpeg",5,5, FUN=function(){
  print(hexbinplot(outIdw$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
             ylab="predicted (IDW)"))
})
printPlot(paper,"img/hexbin/KspA.jpeg",5,5, FUN=function(){
  print(hexbinplot(outK1$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (KSP-A)"))
})
printPlot(paper,"img/hexbin/KspB.jpeg",5,5, FUN=function(){
  print(hexbinplot(outK2$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (KSP-B)"))
})
printPlot(paper,"img/hexbin/Tps.jpeg",5,5, FUN=function(){
  print(hexbinplot(outTps$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (TPS)"))
})
printPlot(paper,"img/hexbin/KstA.jpeg",5,5, FUN=function(){
  print(hexbinplot(outSTK1$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (KST-SE)"))
})
printPlot(paper,"img/hexbin/KstB.jpeg",5,5, FUN=function(){
  print(hexbinplot(outSTK2$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (KST-PS)"))
})
printPlot(paper,"img/hexbin/KstC.jpeg",5,5, FUN=function(){
  print(hexbinplot(outSTK3$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (KST-ME)"))
})
printPlot(paper,"img/hexbin/KstD.jpeg",5,5, FUN=function(){
  print(hexbinplot(outSTK4$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (KST-SM)"))
})
printPlot(paper,"img/hexbin/Gp.jpeg",5,5, FUN=function(){
  print(hexbinplot(outGp$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (GP)"))
})
printPlot(paper,"img/hexbin/Lm.jpeg",5,5, FUN=function(){
  print(hexbinplot(outLm$Ozone~outOri$Ozone, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (LM)"))
})
printPlot(paper,"img/hexbin/Rf.jpeg",5,5, FUN=function(){
  print(hexbinplot(sOzone~sOzoneH, data=outRf, colramp = function(n){rev(heat.colors(n))}, xlab="measured",
           ylab="predicted (RF)"))
})
