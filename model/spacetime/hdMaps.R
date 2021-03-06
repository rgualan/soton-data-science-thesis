## Ref: 
## Spatio-Temporal data in R
## https://www.r-bloggers.com/spatio-temporal-kriging-in-r/
## https://cran.r-project.org/web/packages/gstat/vignettes/st.pdf

## Clean environment
rm(list=ls())

## Libraries
library(gstat)
library(spacetime)
#library(rgdal) ## Not available in the cluster!
library(reshape2)
library(parallel)
source("util/my_helper.R")
source("model/spacetime/testStKriging.R")
source("model/spacetime/fitCovarianceModels.R")

## Global variables ##########################################################
paper <- setupPaper() # Print the plots in image format for the paper
set.seed(123)
smplDays <- sort(sample(365,8))


## Read data #################################################################
epa_version <- 1
epa <- readEpaDataset(epa_version)
epa <- addDateDerivedFeatures(epa)
epa <- transformFeatures(epa)
epa <- scaleTargetVariable(epa)
epa <- detrend_scaled_ozone_poly(epa)
sites <- getSites(epa)

## Assemble STFDF ############################################################
epa.st <- assembleSTFDF(epa) 

## Covariance models #########################################################
fm <- sOzone.res ~ 1
tlags <- 0:10

ticToc({
  covModels <- fitCovarianceModels(epa.st, fm, tlags, paper=F)
})
fitSepModel <- covModels$separable
fitProdSumModel <- covModels$prodSum
fitMetricModel <- covModels$metric
fitSumMetricModel <- covModels$sumMetric
linStAni <- covModels$linStAni

## Interpolation #################################################################################
## Build a grid over CA
gridCA <- SpatialGrid(GridTopology(epa.st@sp@bbox[,1]%/%10*10, c(10,10),
                                   cells.dim=ceiling(apply(epa.st@sp@bbox,1,diff)/10)))
proj4string(gridCA) <- getUTMproj() 
fullgrid(gridCA) <- F

## Pre-compute to avoid using the method over (unavailable library in the cluster)
ind <- over(gridCA, as(getCAmap(proj="utm"),"SpatialPolygons"))  ## TODO: There is some problem here!!!
# saveRDS(ind, file="data/others/ind.RDS")
ind<-readRDS("data/others/ind.RDS")

gridCA <- gridCA[!is.na(ind)]
# class(gridCA)
printPlot(paper,"img/spacetime/gridCA.jpeg",7,7,FUN=function(){
  print(plot(gridCA))
})

CA_pred <- STF(gridCA, epa.st@time[smplDays])
tIDS <- unique(pmax(1,pmin(as.numeric(outer(-5:5, smplDays, "+")), 365)))

ticToc({
  sepPred <- krigeST(fm, data=epa.st[,tIDS],
                     newdata=CA_pred, fitSepModel, nmax=50,
                     stAni=linStAni/24/3600,
                     progress=F) # fitMetricModel$stAni
})
ticToc({
  psPred <- krigeST(fm, data=epa.st[,tIDS],
                    newdata=CA_pred, fitProdSumModel, nmax=50,
                    stAni=linStAni/24/3600,
                    progress=F)
})
ticToc({
  metricPred <- krigeST(fm, data=epa.st[,tIDS],
                    newdata=CA_pred, fitMetricModel, nmax=50,
                    stAni=linStAni/24/3600,
                    progress=F)
})
ticToc({
  sumPred <- krigeST(fm, data=epa.st[,tIDS],
                     newdata=CA_pred, fitSumMetricModel, nmax=50,
                     stAni=fitSumMetricModel$stAni/24/3600,
                     progress=F)
})

saveRDS(sepPred, file="output/spacetime/model.sepPred.RDS")
saveRDS(psPred, file="output/spacetime/model.psPred.RDS")
saveRDS(metricPred, file="output/spacetime/model.metricPred.RDS")
saveRDS(sumPred, file="output/spacetime/model.sumPred.RDS")
## Load the last saved variables
sepPred <- readRDS("output/spacetime/sepPred.RDS")
psPred <- readRDS("output/spacetime/psPred.RDS")
metricPred <- readRDS("output/spacetime/model.metricPred.RDS")
sumPred <- readRDS("output/spacetime/sumPred.RDS")

## Spatio-temporal plots (for the paper) ############################################################
## Separated model
stpl <- stplot(sepPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
               main=NULL,
               sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
                                list("sp.points", epa.st@sp, col=gray(0.25), pch=3, cex=.5)))
printPlot(paper, "img/spacetime/pred_daily_means_ozone_sep.jpeg", 9, 6, FUN=function(){
  print(stpl)
})

## Product-sum model
stpl <- stplot(psPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
               main=NULL,
               sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
                                list("sp.points", epa.st@sp, col=gray(0.25), pch=3, cex=.5)))
printPlot(paper, "img/spacetime/pred_daily_means_ozone_ps.jpeg", 9, 6, FUN=function(){
  print(stpl)
})

## Metric model
stpl <- stplot(metricPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
               main=NULL, 
               sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
                                list("sp.points", epa.st@sp, col=gray(0.25), pch=3, cex=.5)))
printPlot(paper, "img/spacetime/pred_daily_means_ozone_metric.jpeg", 9, 6, FUN=function(){
  print(stpl)
})

# Sum meteric model
stpl <- stplot(sumPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
               main=NULL, 
               sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
                                list("sp.points", epa.st@sp, col=gray(0.25), pch=3, cex=.5)))
printPlot(paper, "img/spacetime/pred_daily_means_ozone_sum.jpeg", 9, 6, FUN=function(){
  print(stpl)
})



## Simple plots
# if(!paper){
#   stplot(sepPred, col.regions=bpy.colors, scales=list(draw=F),
#          main="spatio-temporal separable model")
#   stplot(psPred, col.regions=bpy.colors, scales=list(draw=F),
#          main="spatio-temporal product-sum model")
#   stplot(sumPred, col.regions=bpy.colors, scales=list(draw=F),
#          main="spatio-temporal sum-metric model")
# }
