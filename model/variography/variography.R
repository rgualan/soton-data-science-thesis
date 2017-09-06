## Ref:
# Introduction to Spatio-Temporal Variography - Edzer Pebesma - 2017

## Clean environment
rm(list=ls())

## Load required libraries
library(sp)
library(gstat)
library(spacetime)
library(reshape2)
library(lattice)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Read data ##############################################################################
epa <- readEpaDataset()
epa <- addDateDerivedFeatures(epa)
epa <- transformFeatures(epa)
epa <- scaleTargetVariable(epa)
epa <- detrend_scaled_ozone_poly(epa)
sites <- getSites(epa)

# ## Simple imputation of Ozone ################################################
# source("preprocessing/gbmImpute.R")
# testStation <- "065-2002"
# ## If we are interested in filling a particular station
# ## The design matrix should be restricted to stations around that station
# kn <- getKneighbours(testStation, sites, 50, include_own=T)
# epa_sub <- epa[epa$Station.Code %in% kn$Station.Code,]
# X <- epa_sub[,c("Ozone","Temperature","Wind","Elevation","UTM.X","UTM.Y","Location.Setting","Doy","Dow.number")]
# ticToc({fitGbm <- gbmImpute(X, max.iters = 2, cv.fold = 5, verbose=T)})
# epa_sub$Ozone2 <- fitGbm$x$Ozone[,1]
# epa.test <- epa_sub[epa_sub$Station.Code==testStation,]
# plot(Ozone~Date,epa.test,type="l")
# lines(Ozone2~Date,epa.test,col=4)
# ## Notes:
# ## Biased
# ## Not very good
# 
# ## Approach trying to detrend the time series
# model <- lm(Ozone ~ poly(Doy,4), epa.test) 
# plot(Ozone~Doy,epa.test,type="l")
# trend <- predict(model,epa.test)
# lines(trend, col="deepskyblue", lwd=2)
# epa.test$residuals <- epa.test$Ozone-trend
# plot(residuals~Doy,epa.test,type="l")
# ## Weekly trend??
# wt<-aggregate(residuals~Dow.number,epa.test,mean)
# names(wt)[2] <- "WeeklyMean"
# plot(wt)
# epa.test$residualsTrend<-wt$WeeklyMean[match(epa.test$Dow.number,wt$Dow.number)]
# plot(residuals~Doy,epa.test,type="l")
# lines(residualsTrend~Doy,epa.test,col=2)
# hist(epa.test$residuals)
# 
# ## Now use DBM imputation
# epa_sub$residuals <- NA
# for(s in unique(epa_sub$Station.Code)){
#   print(s)
#   epa.test <- epa_sub[epa_sub$Station.Code==s,]
#   model <- lm(Ozone ~ poly(Doy,4), epa.test) 
#   plot(Ozone~Doy,epa.test,type="l")
#   trend <- predict(model,epa.test)
#   lines(trend, col="deepskyblue", lwd=2)
#   epa.test$residuals <- epa.test$Ozone-trend
#   #plot(residuals~Doy,epa.test,type="l")
#   epa_sub$residuals[epa_sub$Station.Code==s] <- epa.test$residuals
#   #readline("Continue?")
# }
# X <- epa_sub[,c("residuals","Temperature","Wind","Elevation","UTM.X","UTM.Y","Location.Setting","Doy","Dow.number")]
# ticToc({fitGbm <- gbmImpute(X, max.iters = 2, cv.fold = 5, verbose=T)})
# epa_sub$residuals2 <- fitGbm$x$residuals[,1]
# epa.test <- epa_sub[epa_sub$Station.Code==testStation,]
# plot(residuals~Date,epa.test,type="l")
# lines(residuals2~Date,epa.test,col=4)
# 
# ## Try GP
# library(spTimer)
# ticToc(
#   simpleGp <- spT.Gibbs(
#     formula = residuals~Temperature+Elevation+Dow.number,
#     model = "GP",
#     data = epa_sub[epa_sub$Station.Code!=testStation,], 
#     coords = ~UTM.X + UTM.Y, #scale.transform = "SQRT",
#     spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
# )
# gpRed <- predict(simpleGp, newdata=epa.test,
#                  newcoords = ~UTM.X + UTM.Y)
# epa.test$residuals3 <- gpRed$Mean
# #plot(residuals~Date,epa.test,type="l")
# lines(residuals3~Date,epa.test,col=5)
# ## Knn
# library(kknn)
# fit.kknn <- kknn(residuals ~ Temperature+Elevation+Dow.number, 
#                  epa_sub[epa_sub$Station.Code!=testStation & !is.na(epa_sub$residuals),], 
#                  epa_sub[epa_sub$Station.Code==testStation,])
# epa.test$residuals4 <- fit.kknn$fitted.values
# lines(residuals4~Date,epa.test,col=6)
# 
# ## Metrics 
# mGp <- evaluatePredictions(epa.test$residuals, epa.test$residuals3)
# mKknn <- evaluatePredictions(epa.test$residuals, epa.test$residuals4)
# metrics <- rbind(mGp,mKknn)
# rownames(metrics) <- c("gp","kknn")
# metrics
# ## GP is slightly better
# 
# ## Try with an autoregressive function


## Assemble STFDF ############################################################
epa.st <- assembleSTFDF(epa)
summary(epa.st)
stplot(epa.st[,"2016-01-01::2016-01-08","Ozone"])
dim(epa.st)

# ## Assemble STFDF only with urban stations
# ## Space dimension
# epa.sp <- getSites(epa)
# rownames(epa.sp) <- epa.sp$Station.Code
# head(epa.sp)
# coordinates(epa.sp) <- ~UTM.X+UTM.Y
# proj4string(epa.sp) <- getUTMproj()
# ## Time dimension
# epaUrban <- epa[epa$Location.Setting=="URBAN",]
# epaUrban.sp <- getSites(epaUrban)
# rownames(epaUrban.sp) <- epaUrban.sp$Station.Code
# coordinates(epaUrban.sp) <- ~UTM.X+UTM.Y
# proj4string(epaUrban.sp) <- getUTMproj()
# epaUrban.tm <- sort(unique(epaUrban$Date))
# epaUrban.tm <- as.Date(epaUrban.tm)
# # Combine the objects spatial, data-frame and time-dim into a STIDF:
# epaUrban.st <- STFDF(epaUrban.sp,epaUrban.tm,epaUrban[,c("Date","Ozone","sOzone","Temperature",
#                                      "RH","Rain","logRain","Wind","sqrtWind")]) 
# plotStations(F,epaUrban$Station.Code,6,6)


## Select 4 random stations/names (arbitrarily)
#rn <- sample(row.names(epa.st@sp),4)
rn <- row.names(epa.st@sp)[114:117] # random names
plotStations(paper, rn, "img/variogram/rnd_stations.jpeg", 6, 6, fill=2)
# Distances
ks0 <- sites[sites$Station.Code %in% rn,c("Station.Code","UTM.X","UTM.Y")]
ks <- ks0[,-1]
rownames(ks) <- ks0$Station.Code
#head(ks)
m <- as.matrix(dist(ks))
m[m==0] <- NA
min(m, na.rm=T)
max(m, na.rm=T)

## Temporal autocorrelation and cross-correlation functions ##############################################################
printPlot(paper,"img/acf/tm_acf_4s.jpeg",5,5, FUN=function(){
  par(mfrow=c(2,2))
  for(s in rn){
    acf(na.omit(epa.st[s,,"sOzone"][,1]), main = s, lag.max = 20)
  }
  par(mfrow=c(1,1))
})

## Auto- and cross correlations can be computed when a multivariate time series
## object is passed to acf
printPlot(paper,"img/acf/tm_ccf_4s.jpeg",5,5, FUN=function(){
  acf(na.omit(as(epa.st[rn,,"sOzone"], "xts")))
})

## Notes:
# The plot further more shows that for these four stations the asymmetry is not
# very strong, but that cross correlations are fairly strong and of a similar form
# of autocorrelations.
# No cross-correlation:
# and note that here we see in the last figure (DESH & DESN04) a pair of plots
# with nearly no cross correlation. This might have to do with the spatial distance
# between these two stations.

## Check distances
print(spDists(epa.st[114:117,]@sp), digits=3)


## Spatial correlation, variograms ######################################################
# In the next steps, we will sample 100 time instances randomly,
rs = sample(dim(epa.st)[2], 100) # random samples
# we select these instances as a SpatialPointsDataFrame and add a time index
# to them. After this we bind them together in a single SpatialPointsDataFrame
# which has a time index ti:
lst = lapply(rs, function(i) { x = epa.st[,i]; x$ti = i; rownames(x@coords) = NULL; x} )
pts = do.call(rbind, lst)
## Then, we can compute the pooled variogram
v = variogram(sOzone~ti, pts[!is.na(pts$Ozone),], dX=0)
#v = variogram(sOzone~1, epa.st[!is.na(epa.st$sOzone),], dX=0)
vmod = fit.variogram(v, vgm(1, "Exp", 200))
printPlot(paper,"img/variogram/vm_empirical_sOzone.jpeg",5,5,FUN=function(){
  print(plot(v, vmod))
})
vmod


## spatio-temporal variogram the usual way, by passing an object of class STFDF
#fm <- sOzone~Temperature
fm <- sOzone~Elevation+Temperature

ticToc({vv = variogram(fm, epa.st, width=20, cutoff = 200, tlags=0:7)})
ticToc({vv2 = variogram(fm, epa.st, width=20, cutoff = 200, tlags=0:3)})
plot(vv,map=F)
plot(vv2,map=F)
empVgm <- vv2 # Chose one for fitting the variogram models 
linStAni <- estiStAni(empVgm, c(10,200))

printPlot(paper,"img/variogram/emp_stvgm_1.jpeg",6,4,FUN=function(){
  print(plot(vv))
})
printPlot(paper,"img/variogram/emp_stvgm_2.jpeg",5,5,FUN=function(){
  print(plot(vv, map = FALSE))
})
printPlot(paper,"img/variogram/emp_stvgm_custom_1.jpeg",6,4,FUN=function(){
  print(plot(vv2))
})
printPlot(paper,"img/variogram/emp_stvgm_custom_2.jpeg",6,4,FUN=function(){
  print(plot(vv2, map = FALSE))
})
printPlot(paper,"img/variogram/emp_stvgm_1_wireframe.jpeg",5,5,FUN=function(){
  print(plot(vv, wireframe=T, zlab=NULL, xlab=list("Distance (km)", rot=30, cex=0.8),
             ylab=list("Time lag (days)", rot=-35, cex=0.8),
             scales=list(arrows=F, z = list(distance = 5), cex=0.7)))
})
printPlot(paper,"img/variogram/emp_stvgm_2_wireframe.jpeg",5,5,FUN=function(){
  print(plot(vv2, wireframe=T, zlab=NULL, xlab=list("Distance (km)", rot=30, cex=0.8),
             ylab=list("Time lag (days)", rot=-35, cex=0.8),
             scales=list(arrows=F, z = list(distance = 5), cex=0.7)))
})




## Fitting a spatio-temporal variogram model ########################################
## Metric model with spatio-temporal anisotropy #####################################
metricModel <- vgmST("metric",
                   joint=vgm(0.8,"Exp",50,0),
                   stAni=20)
printPlot(paper,"img/variogram/vm_metric_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, metricModel, all=T, map=F))
})
metricModel <- fit.StVariogram(empVgm, metricModel)
# As numerical criterion to judge the goodness of fit of model and sample vari-
# ogram, the root-mean-squared-difference between the surfaces can be obtained by:
attr(metricModel, "optim")$value
# The final model can be plotted with the sample variogram (Figure 5):
printPlot(paper,"img/variogram/vm_metric_fit_1.jpeg",7,5,FUN=function(){
  print(plot(empVgm, metricModel, all=T, map=F))
})
printPlot(paper,"img/variogram/vm_metric_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, metricModel, wireframe=T, all=T, scales=list(arrows=F), zlab=""))
})

## Separable model ###################################################################
sepModel <- vgmST("separable", 
                space=vgm(1,"Exp", 150, 0.1),
                time =vgm(0.9,"Exp", 5, 0.1),
                sill=100)
printPlot(paper,"img/variogram/vm_sep_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, sepModel, all=T, map=F))
})
ticToc({sepModel <- fit.StVariogram(empVgm, sepModel, method = "L-BFGS-B",
                                  lower = c(10,0,0.01,0,1),
                                  upper = c(500,1,20,1,200))})
attr(sepModel, "optim")$value
printPlot(paper,"img/variogram/vm_sep_fit_1.jpeg",7,5,FUN=function(){
  print(plot(empVgm, sepModel, all=T, map=F))
})
printPlot(paper,"img/variogram/vm_sep_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, sepModel, wireframe=T, all=T, scales=list(arrows=F), zlab=""))
})

## Product-sum ##################################################################
prodSumModel <- vgmST("productSum",
                      space=vgm(0.03, "Exp", 75, 0.01),
                      time= vgm(0.2, "Sph", 7,  0.003), 
                      k=70)
printPlot(paper,"img/variogram/vm_prodsum_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, prodSumModel, all=T, map=F))
})
(prodSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                    stAni = linStAni, method = "L-BFGS-B", 
                                    control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                    lower = rep(0.0001,7)))
attr(prodSumModel, "optim.output")$value
printPlot(paper,"img/variogram/vm_prodsum_fit_1.jpeg",7,5,FUN=function(){
  print(plot(empVgm, prodSumModel, all=T, map=F))
})
printPlot(paper,"img/variogram/vm_prodsum_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, prodSumModel, wireframe=T, all=T, scales=list(arrows=F),zlab="'"))
})

## Sum-metric ########################################################################
sumMetricModel <- vgmST("sumMetric",
                        space = vgm(0.15, "Sph", 75, 0.07),
                        time = vgm(0.16, "Exp", 7, 0.1),
                        joint = vgm(0.25, "Sph", 150, 0.1),
                        stAni = linStAni)
printPlot(paper,"img/variogram/vm_summetric_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, sumMetricModel, all=T, map=F))
})
(sumMetricModel <- 
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
attr(sumMetricModel, "optim.output")$value
printPlot(paper,"img/variogram/vm_summetric_fit.jpeg",7,5,FUN=function(){
  print(plot(empVgm, sumMetricModel, all=T, map=F))
})
printPlot(paper,"img/variogram/vm_summetric_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, sumMetricModel, wireframe=T, all=T, scales=list(arrows=F),zlab=""))
})


## All the models in one Figure ###################################################
# A wireframe (3D) plot of sample variogram and fitted variogram models can
# be obtained e.g. by
printPlot(paper,"img/variogram/all.jpeg",7,5,FUN=function(){
  print(plot(vv, list(sepModel, metricModel, sumMetricModel, prodSumModel), all=T, wireframe=T, #zlim=c(0,120),
             zlab=NULL,
             xlab=list("Distance (km)", rot=30, cex=0.6),
             ylab=list("Time lag (days)", rot=-35, cex=0.6),
             scales=list(arrows=F, z = list(distance = 5), cex=0.5)))
})


## Prediction #####################################################################
#folds <- getFolds()
# folds <- cut(sample(1:nrow(sites)),breaks=10,labels=F)
# testIndices <- folds==1
# j = which(testIndices)[1]
# a<-epa.st[!testIndices]
# ticToc({tmp <- krigeST(fm, data=a,
#                        newdata=a, 
#                        sepModel, nmax=50,
#                        stAni=linStAni/24/3600,
#                        progress=F)$var1.pred})
# dates <- getDates()
# plot(epa.st[j,,"sOzone",drop=F]@data$sOzone, type="l")
# lines(tmp,col=2)
# evaluatePredictions(epa.st[j,,"sOzone",drop=F]@data$sOzone,tmp)

