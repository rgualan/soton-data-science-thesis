## Ref: Spatio-Temporal data in R

## Clean environment
rm(list=ls())

## Libraries
library(gstat)
library(spacetime)
library(rgdal)
library(reshape2)

## Pre-process data ############################################################
epa <- readRDS("data/epa/epa_daily/2016/california.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]
epa.sites <- readRDS("data/epa/epa_daily/2016/sites.RDS")
epa.sites <- epa.sites[,c("Station.Code","Longitude","Latitude")]
epa <- merge(epa, epa.sites )
names(epa) <- c("Station.Code", "Date", "Measurement", "x", "y")
#View(epa)
#View(epa.sites)


## Create a SpatialPointsDataFrame
epa.longlat <- epa
coordinates(epa.longlat) <- ~x+y
proj4string(epa.longlat) <- "+proj=longlat +datum=WGS84"

## Transform into Mercator Projection
utm11 = CRS("+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +units=km")
epa.utm <- spTransform(epa.longlat,utm11) 
#rownames(epa.utm@coords) <- epa.utm$Station.Code
epa.utm <- epa.utm[order(epa.utm$Station.Code,epa.utm$Date),]

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
epa.STFDF <- STFDF(epa.sp,epa.tm,data.frame(PM10=as.vector(epa.matrix))) 
summary(epa.STFDF)
stplot(epa.STFDF[,"2016-01-01::2016-01-09"])
dim(epa.STFDF)

hist(log(epa$Measurement))


## Temporal autocorrelation and cross correlation ################################

## Example
data(air)
rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
#rr = rural[,"2005::2010"]
rr = rural[,4019:4138]
unsel = which(apply(as(rr, "xts"), 2, function(x) all(is.na(x))))
r5to10 = rr[-unsel,]
summary(r5to10)


## Next, we will (rather arbitrarily) select four stations, which have the following
## labels:
rn = row.names(epa.STFDF@sp)[4:7]
rn

## autocorrelation functions are computed and plotted
par(mfrow=c(2,2))
# select 4, 5, 6, 7
for(i in rn){
  acf(na.omit(epa.STFDF[i,]), main = i)
  readline("Continue?")
}
par(mfrow=c(1,1))


## This kind of plot does not work very well in layouts of e.g.  10 x 10 sub-plots;
## acf
## automatically chooses 4 x 4 as the maximum a single plot.  To try this out,
## do a 7 x 7 plot
acf(na.omit(as(r5to10[4:10,], "xts")))
acf(na.omit(as(epa.STFDF[4:10,], "xts")))

## Spatial distance between stations
print(spDists(epa.STFDF[4:10,]@sp), digits=3)


print(spDists(epa.STFDF@sp), digits=3)
acf(na.omit(as(epa.STFDF[c(9,10),], "xts")))
#946   634.55257 3629.167
#1067  642.12692 3616.405
acf(na.omit(as(r5to10[c(1,2),], "xts")))
acf(na.omit(as(epa.STFDF[9:10,], "xts")))
par(mfrow=c(2,1))
plot(as(epa.STFDF[9,], "xts"))
plot(as(epa.STFDF[10,], "xts"))
plot(as(r5to10[1,], "xts"))
plot(as(r5to10[2,], "xts"))


## Spatial correlation, variograms ########################################
rs = 1:dim(r5to10)[2]
rs = 1:ncol(epa.STFDF)

lst = lapply(rs, function(i) { x = r5to10[,i]; x$ti = i; rownames(x@coords) = NULL; x} )
lst = lapply(rs, function(i) { x = epa.STFDF[,i]; x$ti = i; rownames(x@coords) = NULL; x} )
pts = do.call(rbind, lst)

v = variogram(PM10~ti, pts[!is.na(pts$PM10),], dX=0)
vmod = fit.variogram(v, vgm(100, "Exp", 200))
plot(v, vmod)
vmod


## Source 2
v2 <- variogramST(PM10~1,data=epa.STFDF,tlags=1:5,na.omit=T)
plot(v2,map=F) 
plot(v2,map=T) 
plot(v2,wireframe=T) 




## Fit a spatio-temporal variogram the usual way, by passing an object of class
## STFDF
vv = variogram(PM10~1, epa.STFDF, width=20, cutoff = 200, tlags=0:5)
plot(vv)
plot(vv, map = FALSE)


## Fitting a spatio-temporal variogram model
metricVgm <- vgmST("metric",
                   joint=vgm(50,"Exp",100,0),
                   stAni=50)
metricVgm <- fit.StVariogram(vv, metricVgm)
attr(metricVgm, "optim")$value
plot(vv, metricVgm)

## Now, let us try to fit and plot a separable model (Figure 6):
sepVgm <- vgmST("separable",
                space=vgm(0.9,"Exp", 123, 0.1),
                time =vgm(0.9,"Exp", 2.9, 0.1),
                sill=100)
sepVgm <- fit.StVariogram(vv, sepVgm, method = "L-BFGS-B",
                          lower = c(10,0,0.01,0,1),
                          upper = c(500,1,20,1,200))
attr(sepVgm, "optim")$value
plot(vv, list(sepVgm, metricVgm))

## A wireframe (3D) plot of sample variogram and fitted variogram models can
## be obtained e.g.  by
library(lattice)
plot(vv, list(sepVgm, metricVgm), all=T, wireframe=T, zlim=c(100,400),
     zlab=NULL,
     xlab=list("distance (km)", rot=30),
     ylab=list("time lag (days)", rot=-35),
     scales=list(arrows=F, z = list(distance = 5)))


library(ggplot2)
ggplot(as.data.frame(epa.STFDF)) + 
  geom_line(aes(x=time, y=PM10, col=sp.ID), alpha=0.5, size=0.5) + 
  theme(legend.position="none")




