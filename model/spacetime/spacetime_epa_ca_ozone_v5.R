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
source("util/my_helper.R")

## Global variables
set.seed(123)
smplDays <- sort(sample(365,8))
forceRun <- T
paper <- setupPaper()
debugLevel <- T

## Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Sites
sites <- getSites(epa)
## Scale target variable #####################################################
epa$sOzone <- scale(epa$Ozone)[,1]  ## scale returns a matrix

## Assemble STFDF ############################################################
## Space dimension
epa.sp <- getSites(epa)
rownames(epa.sp) <- epa.sp$Station.Code
head(epa.sp)
coordinates(epa.sp) <- ~UTM.X+UTM.Y
proj4string(epa.sp) <- getUTMproj()

## Time dimension
epa.tm <- sort(unique(epa$Date))
epa.tm <- as.Date(epa.tm)  ## Ignore time data?? >> Corrects x labels problem in acf
# Combine the objects spatial, data-frame and time-dim into a STIDF:
epa.st <- STFDF(epa.sp,epa.tm,epa[,c("Date","Ozone","sOzone","Temperature",
                                     "RH","Rain","logRain","Wind","sqrtWind")]) 
# summary(epa.st)
# stplot(epa.st[,"2016-01-01::2016-01-08","sOzone"])
# dim(epa.st)


## Formula
#fm <- sOzone~1
fm <- sOzone ~ Temperature
#fm <- sOzone~Temperature+RH+Elevation

## Empirical variogram
## Original empirical empiral (7 time lags)
# empVgm <- variogramST(fm, data=epa.st, tlags=0:7, na.omit=T)
# printPlot(paper,"img/spacetime/empirical_variogram_0.jpeg",6,6,FUN=function(){
#   print(plot(empVgm, map=F))
# })
## Valid variogram 
#empVgm <- variogramST(fm, data=epa.st, tlags=c(0,5,7), na.omit=T)
empVgm <- variogramST(fm, data=epa.st, na.omit=T)
printPlot(paper,"img/spacetime/empirical_variogram.jpeg",6,6,FUN=function(){
  print(plot(empVgm, map=F))
})
printPlot(paper,"img/spacetime/empirical_variogram_2.jpeg",6,6,FUN=function(){
  print(plot(empVgm, wireframe=T, zlab=NULL, xlab=list("Distance (km)", rot=30, cex=0.8),
             ylab=list("Time lag (days)", rot=-35, cex=0.8),
             scales=list(arrows=F, z = list(distance = 5), cex=0.7)))
})

## Notes:
## Ignore the oscilatory lag days
## It seems to have a weekly frequency
## If the time lags are not chosen manually, then the models do not fit the variogram ***
## They just create a plane following the trend of the oscilations
## Another problem that arises if not chosing carefully the time lags:
## When fitting: Error in vgm(1 - par[4], as.character(model$time$model[2]), par[3], par[4],  : 
## range should be positive


## (Eye) Fit spatio-temporal covariance models #######################################
(linStAni <- estiStAni(empVgm, c(10,500)))
## Lower and upper bounds
# pars.l<-c(sill.s=0,range.s=10, nugget.s=0,sill.t=0,range.t=1, nugget.t=0,sill.st=0.1,range.st=10,nugget.st=0,anis=0)
# pars.u<-c(sill.s=2,range.s=500,nugget.s=1,sill.t=2,range.t=15,nugget.t=1,sill.st=5,range.st=200,nugget.st=5,anis=150) 
# pars.c<-c(1,10,1,1,0.1,1,10)
# pars.c<-c(sill.s=0.01,range.s=10,nugget.s=0.01,sill.t=0.01,range.t=1,nugget.t=0.01,sill.st=0.01,range.st=10,nugget.st=0.01,anis=5,k.k=0.1) 

## Separable  ###################################################################################
separableModel <- vgmST("separable",
                   space = vgm(0.45,"Exp",500,0.5), ##Factor,function,range,nugget
                   time = vgm(0.45,"Exp",7,0.5), ##Factor,function,range,nugget 
                   sill=0.8) # Base
printPlot(paper,"img/spacetime/vm_separable_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, separableModel, all=T, map=F, zlab=""))
})
#separableModel <- fit.StVariogram(empVgm, separable, fit.method=0)
#attr(separable_Vgm,"MSE")
# (fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, #11? 
#                                stAni = linStAni, method = "L-BFGS-B", 
#                                #control = list(parscale=c(5,0.1,1,0.1,0.1)),
#                                lower = c(range.s=10,nugget.s=0,range.t=1,nugget.t=0,sill.st=0.1), ##"range.s"  "nugget.s" "range.t"  "nugget.t" "sill"
#                                upper = c(range.s=500,nugget.s=1,range.t=15,nugget.t=1,sill.st=2)))
#                                #lower=pars.l,
#                                #upper=pars.u))
(fitSepModel <- fit.StVariogram(empVgm, separableModel))
#extractParNames(fitSepModel)
attr(fitSepModel, "optim.output")$value
attr(fitSepModel, "MSE")
printPlot(paper,"img/spacetime/vm_separable_fit.jpeg",7,5,FUN=function(){
  print(plot(empVgm, fitSepModel, all=T, map=F))
})
#attr(fitSepModel, "MSE")
## Exp+Exp: 0.002507338
## Exp+Sph: 
## Sph+Exp: 
## Sph+Sph: 

printPlot(paper,"img/spacetime/vm_separable_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F), zlab=""))
})

## Product-sum ################################################################################
prodSumModel <- vgmST("productSum",
                      space=vgm(0.03, "Exp", 75, 0.01),
                      time= vgm(0.2, "Sph", 7,  0.003), 
                      k=70)
printPlot(paper,"img/spacetime/vm_prodsum_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, prodSumModel, all=T, map=F))
})
(fitProdSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                   stAni = linStAni, method = "L-BFGS-B", 
                                   control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                   lower = rep(0.0001,7)))
attr(fitProdSumModel, "optim.output")$value
attr(fitProdSumModel, "MSE")
printPlot(paper,"img/spacetime/vm_prodsum_fit.jpeg",7,5,FUN=function(){
  print(plot(empVgm, fitProdSumModel, all=T, map=F))
})
# Exp+Exp: 0.002605351
# Exp+Sph: 
# Sph+Exp: 
# Sph+Sph: 
printPlot(paper,"img/spacetime/vm_prodsum_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F),zlab="'"))
})

## Metric ################################################################################
metricModel <- vgmST("metric",
                     joint=vgm(0.4, "Mat", 100, 0.3, kappa = 1.5),
                     stAni=50)
printPlot(paper,"img/spacetime/vm_metric_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, metricModel, all=T, map=F))
})
(fitMetricModel <- fit.StVariogram(empVgm, metricModel, fit.method = 7,
                                   stAni = linStAni, method = "L-BFGS-B"))
attr(fitProdSumModel, "optim.output")$value
attr(fitProdSumModel, "MSE")
printPlot(paper,"img/spacetime/vm_metric_fit.jpeg",7,5,FUN=function(){
  print(plot(empVgm, fitProdSumModel, all=T, map=F))
})
# Exp+Exp: 0.002605351
# Exp+Sph: 
# Sph+Exp: 
# Sph+Sph: 
printPlot(paper,"img/spacetime/vm_metric_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F),zlab=""))
})

## Sum-metric #################################################################################
sumMetricModel <- vgmST("sumMetric",
                        space = vgm(0.15, "Sph", 75, 0.07),
                        time = vgm(0.16, "Exp", 7, 0.1),
                        joint = vgm(0.25, "Sph", 150, 0.1),
                        stAni = linStAni)
printPlot(paper,"img/spacetime/vm_summetric_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, sumMetricModel, all=T, map=F))
})
(fitSumMetricModel <- 
  fit.StVariogram(empVgm, sumMetricModel, fit.method = 7, stAni=linStAni,
                  method = "L-BFGS-B", 
                  lower = c(sill.s = 0,  range.s = 10,  nugget.s = 0,
                            sill.t = 0,  range.t = 0.1, nugget.t = 0,
                            sill.st= 0,  range.st = 10, nugget.st = 0, 
                            anis = 40),
                  upper = c(sill.s = 10, range.s = 500,  nugget.s = 5,
                            sill.t = 10, range.t = 15,   nugget.t = 5,
                            sill.st= 10, range.st = 500, nugget.st = 5,
                            anis = 500),
                  control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                 maxit=1e5)))
attr(fitSumMetricModel, "optim.output")$value
attr(fitSumMetricModel, "MSE")
printPlot(paper,"img/spacetime/vm_summetric_fit.jpeg",7,5,FUN=function(){
  print(plot(empVgm, fitSumMetricModel, all=T, map=F))
})
printPlot(paper,"img/spacetime/vm_summetric_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, fitSumMetricModel, wireframe=T, all=T, scales=list(arrows=F),zlab=""))
})
## Sph+Exp+Sph: 0.0009615683


## Interpolation #################################################################################
if(F & fm == formula("sOzone~1")){
  ## Build a grid over CA
  gridCA <- SpatialGrid(GridTopology(epa.st@sp@bbox[,1]%/%10*10, c(10,10),
                                     cells.dim=ceiling(apply(epa.st@sp@bbox,1,diff)/10)))
  proj4string(gridCA) <- getUTMproj() 
  fullgrid(gridCA) <- F
  
  ## Pre-compute to avoid using the method over (unavailable library in the cluster)
  # ind <- over(gridCA, as(getCAmap(proj="utm"),"SpatialPolygons"))
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
    sumPred <- krigeST(fm, data=epa.st[,tIDS],
                       newdata=CA_pred, fitSumMetricModel, nmax=50,
                       stAni=fitSumMetricModel$stAni/24/3600,
                       progress=F)
  })

  # saveRDS(sepPred, file="output/spacetime/sepPred.RDS")
  # saveRDS(psPred, file="output/spacetime/psPred.RDS")
  # saveRDS(sumPred, file="output/spacetime/sumPred.RDS")
  # ## Load the last saved variables
  # sepPred <- readRDS("output/spacetime/sepPred.RDS")
  # psPred <- readRDS("output/spacetime/psPred.RDS")
  # sumPred <- readRDS("output/spacetime/sumPred.RDS")

  ## Spatio-temporal model (for the paper)
  # if(forceRun){
  #   # Separated model
  #   stpl <- stplot(sepPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
  #                  main=NULL, at=seq(-1.35,0.86,(0.86+1.35)/70),
  #                  sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
  #                                   list("sp.points", epa.st@sp, col=gray(0.25), pch=3, cex=.5)))
  #   printPlot(paper, "img/spacetime/pred_daily_means_ozone_sep.jpeg", 9, 6, FUN=function(){
  #     print(stpl)  
  #   })
  # 
  #   # Product-sum model
  #   stpl <- stplot(psPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
  #                  main=NULL, at=seq(-1.35,0.86,(0.86+1.35)/70),
  #                  sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
  #                                   list("sp.points", epa.st@sp, col=gray(0.25), pch=3, cex=.5)))
  #   printPlot(paper, "img/spacetime/pred_daily_means_ozone_ps.jpeg", 9, 6, FUN=function(){
  #     print(stpl)
  #   })
  # 
  #   # Sum meteric model
  #   stpl <- stplot(sumPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
  #                  main=NULL, at=seq(-1.35,0.86,(0.86+1.35)/70),
  #                  sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
  #                                   list("sp.points", epa.st@sp, col=gray(0.25), pch=3, cex=.5)))
  #   printPlot(paper, "img/spacetime/pred_daily_means_ozone_sum.jpeg", 9, 6, FUN=function(){
  #     print(stpl)
  #   })
  #   
  # }
  
  ## Simple plots
  # if(!paper){
  #   stplot(sepPred, col.regions=bpy.colors, scales=list(draw=F),
  #          main="spatio-temporal separable model")
  #   stplot(psPred, col.regions=bpy.colors, scales=list(draw=F),
  #          main="spatio-temporal product-sum model")
  #   stplot(sumPred, col.regions=bpy.colors, scales=list(draw=F),
  #          main="spatio-temporal sum-metric model")
  # }
}


## 10-fold cross-validation ######################################################################
runFold <- function(k, folds, fm, epa.st, outputColumnName, model, linStAni, nmax=50){
  res <- matrix(NA, dim(epa.st)[1], dim(epa.st)[2])
  ticToc({
    for(i in 1:k) {
      cat("Fold", i, "\n")
      testIndices <- folds==i
      
      for(j in which(folds==i)){
        tmp <- krigeST(fm, data=epa.st[!testIndices],
                       newdata=epa.st[j,drop=F], 
                       model, nmax=nmax,
                       stAni=linStAni/24/3600,
                       progress=F)$var1.pred
        res[j, !is.na(epa.st[j,])[,"sOzone"] ] <- tmp
        # Simple partial assessment
        if (debugLevel) print(evaluatePredictions(as.vector(epa.st[j,,"sOzone"][,1]), res[j,]))
        printPlot(paper,paste0("img/spacetime/cv/site",j,".jpg"),6,5,FUN=function(){
          plot(as.numeric(epa.st[j,,"sOzone"][,1]), type="l",ylab=""); lines(res[j,],col=2)
        })
        # cat("=")
      }
      # cat("\n")
      #stop("Debug!")
    }
    epa.st@data$yh <- as.vector(res)
    names(epa.st@data)[ncol(epa.st@data)] <- outputColumnName
  })
  return(epa.st)
}

k <- 1 # TODO!
folds <- readRDS("output/folds.RDS")
if(forceRun){
  epa.st <- runFold(k, folds, fm, epa.st, "sepModel", fitSepModel, linStAni, nmax=50)
  epa.st <- runFold(k, folds, fm, epa.st, "psModel", fitProdSumModel, linStAni, nmax=50)
  epa.st <- runFold(k, folds, fm, epa.st, "sumMetricModel", fitSumMetricModel, linStAni, nmax=50) ## TODO: the var
  epa.st <- runFold(k, folds, fm, epa.st, "prodSumModel", fitProdSumModel, linStAni, nmax=50)
  #View(epa.st@data)
  saveRDS(epa.st,file="output/spacetime/epa.st.RDS")
}else{
  epa.st <- readRDS("output/spacetime/epa.st.RDS")
}

## CV evaluation metrics
rbind(
  evaluatePredictions(epa.st[,,"sepModel",drop=F]@data[[1]], 
                      epa.st[,,"sOzone",drop=F]@data[[1]]),
  evaluatePredictions(epa.st[,,"psModel",drop=F]@data[[1]], 
                      epa.st[,,"sOzone",drop=F]@data[[1]]),
  evaluatePredictions(epa.st[,,"sumMetricModel",drop=F]@data[[1]], 
                      epa.st[,,"sOzone",drop=F]@data[[1]]),
  evaluatePredictions(epa.st[,,"prodSumModel",drop=F]@data[[1]], 
                      epa.st[,,"sOzone",drop=F]@data[[1]])
)

