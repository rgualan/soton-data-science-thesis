## Ref: 
## Spatio-Temporal data in R
## https://www.r-bloggers.com/spatio-temporal-kriging-in-r/
## https://cran.r-project.org/web/packages/gstat/vignettes/st.pdf

## Clean environment
rm(list=ls())
## Settings
set.seed(123)
smplDays <- sort(sample(365,8))
forceRun <- T
#load("model/spacetime/ws1.RData")
pdf("img/spacetime_ozone.pdf")

## Libraries
library(gstat)
library(spacetime)
library(rgdal)
library(reshape2)

## Pre-process data ############################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]
epa.sites <- readRDS("data/epa/epa_daily/2016/sites_ozone.RDS")
epa.sites <- epa.sites[,c("Station.Code","Longitude","Latitude")]
epa <- merge(epa, epa.sites )
names(epa) <- c("Station.Code", "Date", "Measurement", "x", "y")
#View(epa)
#View(epa.sites)

## Standardize variable ############################################################
## http://machinelearningmastery.com/normalize-standardize-time-series-data-python/
# hist(epa$Measurement)
theMean <- mean(epa$Measurement)
theSd <- sd(epa$Measurement)
epa$Measurement <- (epa$Measurement-theMean)/theSd
hist(epa$Measurement)
# ## Re-scale
# epa$Measurement3 <- epa$Measurement2*theSd+theMean
# hist(epa$Measurement3)
# View(epa)


## Create a SpatialPointsDataFrame
epa.longlat <- epa
coordinates(epa.longlat) <- ~x+y
proj4string(epa.longlat) <- "+proj=longlat +datum=WGS84"

## Transform into Mercator Projection
utm11 = CRS("+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +units=km")
epa.utm <- spTransform(epa.longlat,utm11) 
epa.utm <- epa.utm[order(epa.utm$Station.Code,epa.utm$Date),]


## California boundaries
mapUSA <- readRDS("data/maps/usa/USA_adm1.rds")
mapCA <- mapUSA[mapUSA$NAME_1=="California",]
proj4string(mapCA)
plot(mapCA)
mapCA.utm <- spTransform(mapCA,utm11) 
plot(mapCA.utm)
range(coordinates(mapCA))
range(coordinates(mapCA.utm))



## Assemble STFDF ############################################################
## Data
epa.matrix <- dcast(epa.utm@data, Date ~ Station.Code, value.var="Measurement", fill=NA)
epa.matrix <- as.matrix(epa.matrix[,-1])
image(epa.matrix, xlab="Stations", ylab="Date")
## Space
epa.sp <- data.frame(unique(data.frame(Station.Code=epa.utm$Station.Code,coordinates(epa.utm))))
rownames(epa.sp) <- epa.sp$Station.Code
head(epa.sp)
coordinates(epa.sp) <- ~x+y
proj4string(epa.sp) <- utm11
## Time
epa.tm <- sort(unique(epa.utm$Date))
# Combine the objects spatial, data-frame and time-dim into a STIDF:
epa.STFDF <- STFDF(epa.sp,epa.tm,data.frame(Measurement=as.vector(epa.matrix))) 
summary(epa.STFDF)
stplot(epa.STFDF[,"2016-01-01::2016-01-09"])
dim(epa.STFDF)


## Simple data imputation (TODO)
# sum(is.na(epa.matrix))
# epa.matrix[is.na(epa.matrix)] <- 0

## Exploratory Variogram Analysis ################################################
data(meuse)
coordinates(meuse) <- ~x+y

## Scatter plots of pairs Z(si) and Z(sj), grouped according to their 
## separation distance h
hscat(log(zinc) ~ 1, meuse, (0:9) * 100)
hscat(Measurement ~ 1, 
      epa.utm[epa.utm$Date==sample(epa.utm$Date,1),], 
      (0:9) * 100)

## variogram cloud
empVgm <- variogramST(Measurement~1, data=epa.STFDF, tlags=0:3, na.omit=T)
plot(empVgm, map=F)
plot(empVgm, wireframe=T, scales=list(arrows=F))


## Fit spatio-temporal covariance models #######################################
## Eye fit

linStAni <- estiStAni(empVgm, c(10,200))

separableModel <- vgmST("separable",
                   space = vgm(0.004,"Exp",20, 0),
                   time = vgm(0.015,"Sph",4.5, 0), sill=100) 
plot(empVgm, separableModel, all=T, map=F)
#separableModel <- fit.StVariogram(empVgm, separable, fit.method=0)
#attr(separable_Vgm,"MSE")
# fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, 
#                                stAni = linStAni, method = "L-BFGS-B")
fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, 
                               stAni = linStAni, method = "L-BFGS-B", 
                               control = list(parscale=c(100,1,10,1,100)),
                               lower = c(10,0,.1,0,0.1), 
                               upper = c(10000,1,15,1,200))
attr(fitSepModel, "optim.output")$value
attr(fitSepModel, "MSE")
plot(empVgm, fitSepModel, all=T, map=F)
#attr(fitSepModel, "MSE")
## Exp+Exp: 0.001201547
## Exp+Sph: 0.001027299
## Sph+Exp: 0.001027299
## Sph+Sph: 0.001027299

plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F))


# product-sum
prodSumModel <- vgmST("productSum",
                      space=vgm(0.04, "Exp", 500, 0.1),
                      time= vgm(0.7, "Sph",   5, 0.1), 
                      k=2.5)
plot(empVgm, prodSumModel, all=T, map=F)
fitProdSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                   stAni = linStAni, method = "L-BFGS-B", 
                                   control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                   lower = rep(0.0001,7))
attr(fitProdSumModel, "optim.output")$value
attr(fitProdSumModel, "MSE")
plot(empVgm, fitProdSumModel, all=T, map=F)
# Exp+Exp: 10.09, Exp+Sph: 6.91, Sph+Exp: 10.64, Sph+Sph: 7.59
plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F))



## Sum-metric
sumMetricModel <- vgmST("sumMetric",
                        space = vgm(0.1, "Sph", 20, 0),
                        time = vgm(0.8, "Exp", 1.5, 0),
                        joint = vgm(0.1, "Sph", 200, 0.2),
                        stAni = linStAni)
plot(empVgm, sumMetricModel, all=T, map=F)
fitSumMetricModel <- 
  fit.StVariogram(empVgm, sumMetricModel, fit.method = 7, stAni=linStAni,
                  method = "L-BFGS-B", 
                  lower = c(sill.s = 0,  range.s = 10,  nugget.s = 0,
                            sill.t = 0,  range.t = 0.1,   nugget.t = 0,
                            sill.st= 0, range.st = 10, nugget.st = 0, 
                            anis = 40),
                  upper = c(sill.s = 200,  range.s = 1E4,  nugget.s = 20,
                            sill.t = 200,  range.t = 20,   nugget.t = 20,
                            sill.st= 200, range.st = 5E3, nugget.st = 20,
                            anis = 500),
                  control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                 maxit=1e5))
attr(fitSumMetricModel, "optim.output")$value
attr(fitSumMetricModel, "MSE")
plot(empVgm, fitSumMetricModel, all=T, map=F)
plot(empVgm, fitSumMetricModel, wireframe=T, all=T, scales=list(arrows=F))



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## INTERPOLATION #################################################################################
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Build a grid over CA
gridCA <- SpatialGrid(GridTopology(epa.STFDF@sp@bbox[,1]%/%10*10, c(10,10),
                                   cells.dim=ceiling(apply(epa.STFDF@sp@bbox,1,diff)/10)))
proj4string(gridCA) <- utm11 
fullgrid(gridCA) <- F

# class(DE_NUTS1)
# class(mapCA.utm)
# plot(DE_NUTS1)
# plot(mapCA.utm)
ind <- over(gridCA, as(mapCA.utm,"SpatialPolygons"))
gridCA <- gridCA[!is.na(ind)]
# class(gridCA)
plot(gridCA)

CA_pred <- STF(gridCA, epa.STFDF@time[smplDays])
tIDS <- unique(pmax(1,pmin(as.numeric(outer(-5:5, smplDays, "+")), 365)))

if(forceRun){
  sepPred <- krigeST(Measurement~1, data=epa.STFDF[,tIDS],
                     newdata=CA_pred, fitSepModel, nmax=50,
                     stAni=linStAni/24/3600) # fitMetricModel$stAni
  psPred <- krigeST(Measurement~1, data=epa.STFDF[,tIDS],
                    newdata=CA_pred, fitProdSumModel, nmax=50,
                    stAni=linStAni/24/3600)
  sumPred <- krigeST(Measurement~1, data=epa.STFDF[,tIDS],
                     newdata=CA_pred, fitSumMetricModel, nmax=50,
                     stAni=fitSumMetricModel$stAni/24/3600)
  saveRDS(sepPred, file="output/sepPred.RDS")
  saveRDS(psPred, file="output/psPred.RDS")
  saveRDS(sumPred, file="output/sumPred.RDS")
}else{
  ## Load the last saved variables
  sepPred <- readRDS("output/sepPred.RDS")
  psPred <- readRDS("output/psPred.RDS")
  sumPred <- readRDS("output/sumPred.RDS")
}

# Spatio-temporal model (for the paper)
if(forceRun){
  stpl <- stplot(sepPred, col.regions=bpy.colors(120)[-(1:20)], scales=list(draw=F),
                 main=NULL, at=seq(-1.35,0.86,(0.86+1.35)/70),
                 sp.layout = list(list("sp.polygons", mapCA.utm, first=FALSE, col=gray(0.5)),
                                  list("sp.points", epa.STFDF@sp, col=gray(0.25), pch=3, cex=.5)))
  #png(file="img/pred_daily_means_ozone.png", width=9, height=6, "in", res=150)
  print(stpl)
  #dev.off()
}

## Simple plots
stplot(sepPred, col.regions=bpy.colors, scales=list(draw=F),
       main="spatio-temporal separable model")
stplot(psPred, col.regions=bpy.colors, scales=list(draw=F),
       main="spatio-temporal product-sum model")
stplot(sumPred, col.regions=bpy.colors, scales=list(draw=F),
       main="spatio-temporal sum-metric model")



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Cross-validation ######################################################################
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crossStat <- function(var1, var2="Measurement", STxDF=epa.STFDF, digits=NA) {
  diff <- STxDF[,,var1,drop=F]@data[[1]] - STxDF[,,var2,drop=F]@data[[1]]
  RMSE <- sqrt(mean(diff^2,na.rm=T))
  MAE <- mean(abs(diff),na.rm=T)
  ME <- mean(diff,na.rm=T)
  COR <- cor(STxDF[,,var1,drop=F]@data[[1]], STxDF[,,var2,drop=F]@data[[1]],
             use = "pairwise.complete.obs")
  res <- c(RMSE, MAE, ME, COR)
  names(res) <- c("RMSE", "MAE", "ME", "COR")
  if(is.na(digits))
    return(res)
  else
    return(round(res, digits))
}


## 10-CV
k <- 10
#folds <- cut(sample(1:dim(epa.STFDF)[1]),breaks=k,labels=F)
#saveRDS(folds, file="output/folds.RDS")
folds <- readRDS("output/folds.RDS")
par(ask=F) # For showing TS plots en each iteration

if(forceRun){

  ## Separable model - 50 neighbours
  res <- matrix(NA, dim(epa.STFDF)[1], dim(epa.STFDF)[2])
  for(i in 1:k) {
    cat("Fold", i, "\n")
    testIndices <- folds==i
  
    for(j in which(folds==i)){
      tmp <- krigeST(Measurement~1, data=epa.STFDF[!testIndices,,1],
                       newdata=epa.STFDF[j,drop=F], 
                       fitSepModel, nmax=50,
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
  epa.STFDF@data$sepModel <- as.vector(res)
  ## Re-scale TODO
  #epa.STFDF@data$sepModel10Nghbr <- epa.STFDF@data$sepModel10Nghbr*theSd+theMean
  
  ## product-sum model - 50 neighbours
  res <- matrix(NA, dim(epa.STFDF)[1], dim(epa.STFDF)[2])
  for(i in 1:k) {
    cat("Fold", i, "\n")
    testIndices <- folds==i
    
    for(j in which(folds==i)){
      tmp <- krigeST(Measurement~1, data=epa.STFDF[!testIndices,,1],
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
  
  ## Sum-metric model - 50 neighbours
  res <- array(NA, c(dim(epa.STFDF)[1], dim(epa.STFDF)[2], 2))
  for(i in 1:k) {
    cat("Fold", i, "\n")
    testIndices <- folds==i
    
    for(j in which(folds==i)){
      tmp <- krigeST(Measurement~1, data=epa.STFDF[!testIndices,,1],
                     newdata=epa.STFDF[j,drop=F], 
                     fitSumMetricModel, nmax=50,
                     computeVar=T,
                     stAni=fitSumMetricModel$stAni/24/3600)@data[,c("var1.pred","var1.var")]
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

  #View(epa.STFDF@data)
  saveRDS(epa.STFDF,file="output/epa.STFDF.RDS")
}else{
  readRDS("output/epa.STFDF.RDS")
}

## Cross-stats
rbind(
  crossStat("sepModel", digits=2),
  crossStat("psModel", digits=2),
  crossStat("sumMetricModel", digits=2)
)



