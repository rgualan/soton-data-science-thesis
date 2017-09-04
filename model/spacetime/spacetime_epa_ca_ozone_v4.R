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
#load("model/spacetime/ws1.RData")


## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]

## Standardize variable ############################################################
## http://machinelearningmastery.com/normalize-standardize-time-series-data-python/
# hist(epa$Ozone)
# theMean <- mean(epa$Ozone, na.rm=T)
# theSd <- sd(epa$Ozone,na.rm=T)
# epa$sOzone <- (epa$Ozone-theMean)/theSd
# hist(epa$Ozone)

## Scale and add a Bias, to avoid negative numbers
# epa$sOzone <- scale(epa$Ozone)
# epa$sOzone2 <- epa$sOzone+abs(min(epa$sOzone,na.rm=T))
epa$sOzone2 <- scalePlusBias(epa$Ozone)
# hist(epa$sOzone2)

## Normalize
# simpleNormaliza <- function(x){
#   (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
# }
# epa$nOzone <- simpleNormaliza(epa$Ozone)
# lines(epa$nOzone,col=3)
# hist(epa$nOzone)

# ## Re-scale
# epa$Ozone3 <- epa$Ozone2*theSd+theMean
# hist(epa$Ozone3)
# View(epa)


## Create a SpatialPointsDataFrame ####################################################
epa.utm <- convertDataToSp(epa,"utm")
#View(epa.utm)


## Assemble STFDF ############################################################
## Data
## HERE THE VARIABLE IS CHOSEN!!!
epa.matrix <- dcast(epa.utm@data, Date ~ Station.Code, value.var="sOzone2", fill=NA)
epa.matrix <- as.matrix(epa.matrix[,-1])
#image(epa.matrix, xlab="Stations", ylab="Date")
## Space
epa.sp <- data.frame(unique(data.frame(Station.Code=epa.utm$Station.Code,coordinates(epa.utm))))
rownames(epa.sp) <- epa.sp$Station.Code
head(epa.sp)
coordinates(epa.sp) <- ~UTM.X+UTM.Y
proj4string(epa.sp) <- getUTMproj()
## Time
epa.tm <- sort(unique(epa.utm$Date))
# Combine the objects spatial, data-frame and time-dim into a STIDF:
epa.STFDF <- STFDF(epa.sp,epa.tm,data.frame(Ozone=as.vector(epa.matrix))) 
summary(epa.STFDF)
#stplot(epa.STFDF[,"2016-01-01::2016-01-09"])
#dim(epa.STFDF)


## Empirical variogram
empVgm <- variogramST(Ozone~1, data=epa.STFDF, tlags=0:7, na.omit=T)
printPlot(paper,"img/spacetime/empirical_variogram_0.jpeg",6,6,FUN=function(){
  print(plot(empVgm, map=F))
})
empVgm <- variogramST(Ozone~1, data=epa.STFDF, tlags=c(0,5), na.omit=T)
printPlot(paper,"img/spacetime/empirical_variogram.jpeg",6,6,FUN=function(){
  print(plot(empVgm, map=F))
})
printPlot(paper,"img/spacetime/empirical_variogram_2.jpeg",6,6,FUN=function(){
  print(plot(empVgm, wireframe=T, scales=list(arrows=F)))
})

## Notes:
## Ignore the oscilatory days
## It seems to have a weekly frequency

## (Eye) Fit spatio-temporal covariance models #######################################
linStAni <- estiStAni(empVgm, c(10,200))
## Lower and upper bounds
# pars.l <- c(sill.s = 0, range.s = 10, nugget.s = 0,sill.t = 0, range.t = 1, nugget.t = 0,sill.st = 0, range.st = 10, nugget.st = 0, anis = 0)
# pars.u <- c(sill.s = 200, range.s = 1000, nugget.s = 100,sill.t = 200, range.t = 60, nugget.t = 100,sill.st = 200, range.st = 1000, nugget.st = 100,anis = 700) 
pars.l<-c(sill.s=0,range.s=10, nugget.s=0,sill.t=0,range.t=1, nugget.t=0,sill.st=0.1,range.st=10,nugget.st=0,anis=0)
pars.u<-c(sill.s=2,range.s=500,nugget.s=1,sill.t=2,range.t=15,nugget.t=1,sill.st=5,range.st=200,nugget.st=5,anis=150) 
pars.c<-c(1,10,1,1,0.1,1,10)
pars.c<-c(sill.s=0.01,range.s=10,nugget.s=0.01,sill.t=0.01,range.t=1,nugget.t=0.01,sill.st=0.01,range.st=10,nugget.st=0.01,anis=5,k.k=0.1) 

## Separable
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

# product-sum
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


## Sum-metric
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


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Interpolation #################################################################################
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Build a grid over CA
gridCA <- SpatialGrid(GridTopology(epa.STFDF@sp@bbox[,1]%/%10*10, c(10,10),
                                   cells.dim=ceiling(apply(epa.STFDF@sp@bbox,1,diff)/10)))
proj4string(gridCA) <- getUTMproj() 
fullgrid(gridCA) <- F

# ind <- over(gridCA, as(getCAmap(proj="utm"),"SpatialPolygons"))
# saveRDS(ind, file="output/spacetime/ind.RDS")
ind<-readRDS("output/spacetime/ind.RDS")

gridCA <- gridCA[!is.na(ind)]
# class(gridCA)
printPlot(paper,"img/spacetime/gridCA.jpeg",7,7,FUN=function(){
  print(plot(gridCA))
})

CA_pred <- STF(gridCA, epa.STFDF@time[smplDays])
tIDS <- unique(pmax(1,pmin(as.numeric(outer(-5:5, smplDays, "+")), 365)))

if(forceRun){
  st<-Sys.time()
  sepPred <- krigeST(Ozone~1, data=epa.STFDF[,tIDS],
                     newdata=CA_pred, fitSepModel, nmax=50,
                     stAni=linStAni/24/3600,
                     progress=F) # fitMetricModel$stAni
  Sys.time()-st
  st<-Sys.time()
  psPred <- krigeST(Ozone~1, data=epa.STFDF[,tIDS],
                    newdata=CA_pred, fitProdSumModel, nmax=50,
                    stAni=linStAni/24/3600,
                    progress=F)
  Sys.time()-st
  st<-Sys.time()
  sumPred <- krigeST(Ozone~1, data=epa.STFDF[,tIDS],
                     newdata=CA_pred, fitSumMetricModel, nmax=50,
                     stAni=fitSumMetricModel$stAni/24/3600,
                     progress=F)
  Sys.time()-st
  saveRDS(sepPred, file="output/spacetime/sepPred.RDS")
  saveRDS(psPred, file="output/spacetime/psPred.RDS")
  saveRDS(sumPred, file="output/spacetime/sumPred.RDS")
}else{
  ## Load the last saved variables
  sepPred <- readRDS("output/spacetime/sepPred.RDS")
  psPred <- readRDS("output/spacetime/psPred.RDS")
  sumPred <- readRDS("output/spacetime/sumPred.RDS")
}

# Spatio-temporal model (for the paper)
if(forceRun){
  # Separated model
  stpl <- stplot(sepPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
                 main=NULL, at=seq(-1.35,0.86,(0.86+1.35)/70),
                 sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
                                  list("sp.points", epa.STFDF@sp, col=gray(0.25), pch=3, cex=.5)))
  printPlot(paper, "img/spacetime/pred_daily_means_ozone_sep.jpeg", 9, 6, FUN=function(){
    print(stpl)  
  })

  # Product-sum model
  stpl <- stplot(psPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
                 main=NULL, at=seq(-1.35,0.86,(0.86+1.35)/70),
                 sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
                                  list("sp.points", epa.STFDF@sp, col=gray(0.25), pch=3, cex=.5)))
  printPlot(paper, "img/spacetime/pred_daily_means_ozone_ps.jpeg", 9, 6, FUN=function(){
    print(stpl)
  })

  # Sum meteric model
  stpl <- stplot(sumPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
                 main=NULL, at=seq(-1.35,0.86,(0.86+1.35)/70),
                 sp.layout = list(list("sp.polygons", getCAmap(proj="utm"), first=FALSE, col=gray(0.5)),
                                  list("sp.points", epa.STFDF@sp, col=gray(0.25), pch=3, cex=.5)))
  printPlot(paper, "img/spacetime/pred_daily_means_ozone_sum.jpeg", 9, 6, FUN=function(){
    print(stpl)
  })
  
}

## Simple plots
if(!paper){
  stplot(sepPred, col.regions=bpy.colors, scales=list(draw=F),
         main="spatio-temporal separable model")
  stplot(psPred, col.regions=bpy.colors, scales=list(draw=F),
         main="spatio-temporal product-sum model")
  stplot(sumPred, col.regions=bpy.colors, scales=list(draw=F),
         main="spatio-temporal sum-metric model")
}



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 10-fold cross-validation ######################################################################
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
k <- 10
folds <- getFolds()
#par(ask=F) # For showing TS plots en each iteration

if(forceRun){

  ## Separable model - 50 neighbours
  res <- matrix(NA, dim(epa.STFDF)[1], dim(epa.STFDF)[2])
  st<-Sys.time()
  for(i in 1:k) {
    cat("Fold", i, "\n")
    testIndices <- folds==i
  
    for(j in which(folds==i)){
      tmp <- krigeST(Ozone~1, data=epa.STFDF[!testIndices,,1],
                     newdata=epa.STFDF[j,drop=F], 
                     fitSepModel, nmax=50,
                     stAni=linStAni/24/3600,
                     progress=F)$var1.pred
      res[j, !is.na(epa.STFDF[j,])[,1] ] <- tmp
      # Simple partial assessment
      tmp2 <- cbind(epa.STFDF[j,][,1],res[j,])
      RMSE <- sqrt(mean((tmp2[,1]-tmp2[,2])^2,na.rm=T))
      COR <- cor(tmp2[,1], tmp2[,2], use = "pairwise.complete.obs")
      print(c(RMSE=RMSE, COR=COR))
      plot(tmp2[,1], type="l"); lines(tmp2[,2],col=2)
    }
  }
  epa.STFDF@data$sepModel <- as.vector(res)
  Sys.time()-st
  
  ## Re-scale TODO
  #epa.STFDF@data$sepModel10Nghbr <- epa.STFDF@data$sepModel10Nghbr*theSd+theMean
  
  ## product-sum model - 50 neighbours
  res <- matrix(NA, dim(epa.STFDF)[1], dim(epa.STFDF)[2])
  st<-Sys.time()
  for(i in 1:k) {
    cat("Fold", i, "\n")
    testIndices <- folds==i
    
    for(j in which(folds==i)){
      tmp <- krigeST(Ozone~1, data=epa.STFDF[!testIndices,,1],
                     newdata=epa.STFDF[j,drop=F], 
                     fitProdSumModel, nmax=50,
                     stAni=linStAni/24/3600)$var1.pred
      res[j, !is.na(epa.STFDF[j,])[,1] ] <- tmp
      # Simple partial assessment
      tmp2 <- cbind(epa.STFDF[j,][,1],res[j,])
      RMSE <- sqrt(mean((tmp2[,1]-tmp2[,2])^2,na.rm=T))
      COR <- cor(tmp2[,1], tmp2[,2], use = "pairwise.complete.obs")
      print(c(RMSE=RMSE, COR=COR))
      plot(tmp2[,1], type="l"); lines(tmp2[,2],col=2)
    }
  }
  epa.STFDF@data$psModel <- as.vector(res)
  Sys.time()-st
  
  ## Sum-metric model - 50 neighbours
  res <- array(NA, c(dim(epa.STFDF)[1], dim(epa.STFDF)[2], 2))
  st<-Sys.time()
  for(i in 1:k) {
    cat("Fold", i, "\n")
    testIndices <- folds==i
    
    for(j in which(folds==i)){
      tmp <- krigeST(Ozone~1, data=epa.STFDF[!testIndices,,1],
                     newdata=epa.STFDF[j,drop=F], 
                     fitSumMetricModel, nmax=50,
                     computeVar=T,
                     stAni=fitSumMetricModel$stAni/24/3600,
                     progress=F)@data[,c("var1.pred","var1.var")]
      res[j, !is.na(epa.STFDF[j,])[,1], ] <- as.matrix(tmp)
      # Simple partial assessment
      tmp2 <- cbind(epa.STFDF[j,][,1], res[j,,1])
      RMSE <- sqrt(mean((tmp2[,1]-tmp2[,2])^2,na.rm=T))
      COR <- cor(tmp2[,1], tmp2[,2], use = "pairwise.complete.obs")
      print(c(RMSE=RMSE, COR=COR))
      plot(tmp2[,1], type="l"); lines(tmp2[,2],col=2)
    }
  }
  epa.STFDF@data$sumMetricModel <- as.vector(res[,,1])
  epa.STFDF@data$sumMetricModelVar <- as.vector(res[,,2])
  epa.STFDF@data$sumMetricModel95u <- 
    apply(epa.STFDF@data, 1, function(x) {
      qnorm(0.975, x["sumMetricModel"], sqrt(x["sumMetricModelVar"]))
    })
  epa.STFDF@data$sumMetricModel95l <- 
    apply(epa.STFDF@data, 1, function(x) {
      qnorm(0.025, x["sumMetricModel"], sqrt(x["sumMetricModelVar"]))
    })
  Sys.time()-st
  
  #View(epa.STFDF@data)
  saveRDS(epa.STFDF,file="output/spacetime/epa.STFDF.RDS")
}else{
  readRDS("output/spacetime/epa.STFDF.RDS")
}

## Cross-stats
# crossStat <- function(var1, var2="Ozone", STxDF=epa.STFDF, digits=NA) {
#   evaluatePredictions(STxDF[,,var1,drop=F]@data[[1]], STxDF[,,var2,drop=F]@data[[1]])
# }

rbind(
  evaluatePredictions(epa.STFDF[,,"sepModel",drop=F]@data[[1]], 
                      epa.STFDF[,,"Ozone",drop=F]@data[[1]]),
  evaluatePredictions(epa.STFDF[,,"psModel",drop=F]@data[[1]], 
                      epa.STFDF[,,"Ozone",drop=F]@data[[1]]),
  evaluatePredictions(epa.STFDF[,,"sumMetricModel",drop=F]@data[[1]], 
                      epa.STFDF[,,"Ozone",drop=F]@data[[1]])
)

