## Ref: Spatio-Temporal data in R

## Clean environment
rm(list=ls())

## Libraries
library(gstat)
library(spacetime)
library(rgdal)

## Pre-process data ############################################################
data("wind")
epa <- readRDS("data/epa/epa_daily/2016/california.RDS")
epa.sites <- readRDS("data/epa/epa_daily/2016/sites.RDS")
epa.sites <- epa.sites[,c("Station.Code","Longitude","Latitude")]
epa <- merge(epa, epa.sites )
names(epa) <- c("Station.Code", "Date", "Var", "x", "y")
#View(epa)
#View(epa.sites)

mapUSA <- readRDS("data/maps/usa/USA_adm1.rds")
mapCA <- mapUSA[mapUSA$NAME_1=="California",]
proj4string(mapCA)
plot(mapCA)


wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y
proj4string(wind.loc) = "+proj=longlat +datum=WGS84"
coordinates(epa.sites) = ~Longitude+Latitude
proj4string(epa.sites) = "+proj=longlat +datum=WGS84"


head(wind)
wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
wind$jday = as.numeric(format(wind$time, ' %j ' ))
epa$jday = as.numeric(format(epa$Date, ' %j ' ))
head(wind)
head(epa)

## subtract a smooth time trend of daily means ##################################
stations = 4:15
windsqrt = sqrt(0.5148 * as.matrix(wind[stations])) # knots -> m/s
epa$Var = sqrt(epa$Var) # Transform
Jday = 1:366
windsqrt = windsqrt - mean(windsqrt)
epa$Var = epa$Var - mean(epa$Var)
daymeans = sapply(split(windsqrt, wind$jday), mean)
daymeans2 = aggregate(Var~jday,epa,mean)
meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
meanVar = lowess(daymeans2$Var ~ daymeans2$jday, f = 0.1)
meanVar2 <- data.frame(jday=meanVar$x, mean=meanVar$y)
velocities = apply(windsqrt, 2, function(x) { x - meanwind })
hist(velocities)
tmp = epa[,c("jday","Var")]
head(tmp)
head(meanVar2)
tmp <- merge(tmp,meanVar2)
head(tmp)
epa$Var <- tmp$Var-tmp$mean
hist(epa$Var)
range(epa$Var)
plot(daymeans2)
lines(meanVar2$mean,col=2)
plot(Var~jday,epa)


## match the wind data to its location, by connecting station names to location
## coordinates, and create a spatial points object:
wind.loc = wind.loc[match(names(wind[4:15]), wind.loc$Code),]
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])
rownames(pts) = wind.loc$Station
pts = SpatialPoints(pts, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
pts2 = coordinates(epa.sites)
rownames(pts2) = epa.sites$Station.Code
pts2 = SpatialPoints(pts2, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
head(pts)
head(pts2)

## Project the longitude/latitude coordinates and country boundary to UTM zone 29
utm29 = CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84")
pts = spTransform(pts, utm29)
utm11 = CRS("+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84")
pts2 = spTransform(pts2, utm11)
head(pts2)

## Construct the spatio-temporal object from the space-wide table with velocities:
wind.data = stConstruct(velocities, space = list(values = 1:ncol(velocities)),
                        time = wind$time, SpatialObj = pts, interval = TRUE)
class(wind.data)
epa.data = stConstruct(epa, c("x","y"), "Date", interval = TRUE)
class(wind.data)
class(epa.data)
## TODO
proj4string(epa.data@sp) <- utm11
class(wind.data@sp)
class(epa.data@sp)




## For plotting purposes, we can obtain country boundaries from package maps:
library(maps)
library(mapdata)
library(maptools)
m = map2SpatialLines(
  #map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
  map("worldHires", xlim = c(-11.5,-6.0), ylim = c(51.3,55.0), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84"
m = spTransform(m, utm29)
plot(m)
mapCa.utm = spTransform(mapCA, utm11)
plot(mapCa.utm)

## For interpolation, we can define a grid over the area:
grd = SpatialPixels(SpatialPoints(makegrid(m, n = 300)),
                    proj4string = proj4string(m))
grd2 = SpatialPixels(SpatialPoints(makegrid(mapCa.utm, n = 300)),
                    proj4string = proj4string(mapCa.utm))
plot(m); plot(grd, add=T)
plot(mapCa.utm); plot(grd2, add=T)


## We (arbitrarily) restrict observations to those of April 1961:
wind.data = wind.data[, "1961-04"]
## and choose 10 time points from that period to form the spatio-temporal 
## prediction grid:
n = 10
library(xts)
library(RColorBrewer)
tgrd = seq(min(index(wind.data)), max(index(wind.data)), length=n)
pred.grd = STF(grd, tgrd)
tgrd2 = seq(min(index(epa.data)), max(index(epa.data)), length=n)
pred.grd2 = STF(grd2, tgrd2)

## Fit variograms



## We will interpolate with a separable exponential covariance model, 
## with ranges 750 km and 1.5 days:
v = vgmST("separable", space = vgm(1, "Exp", 750000), 
          time = vgm(1, "Exp", 1.5 * 3600), sill=0.6)
wind.ST = krigeST(values ~ 1, wind.data, pred.grd, v) # All data: ~10 min
colnames(wind.ST@data) <- "sqrt_speed"
v2 = vgmST("separable", space = vgm(1, "Exp", 750000), 
          time = vgm(1, "Exp", 1.5 * 3600), sill=0.6)
epa.ST = krigeST(Var ~ 1, epa.data, pred.grd2, v2) # All data: ~10 min
#saveRDS(epa.ST, file="data/tmp/epa.ST.RDS")
#epa.ST = readRDS("data/tmp/epa.ST.RDS")
colnames(epa.ST@data) <- "sqrt_var"

## Then creates the STFDF object with interpolated values
layout = list(list("sp.lines", m, col= ' grey ' ),
              list("sp.points", pts, first=F, cex=.5))
stplot(wind.ST, col.regions=brewer.pal(11, "RdBu")[-c(10,11)],
       at=seq(-1.375,1,by=.25),
       par.strip.text = list(cex=.7), sp.layout = layout)
stplot(epa.ST, col.regions=brewer.pal(11, "RdBu"),
       par.strip.text = list(cex=.7), sp.layout = layout)

## Space-time plots
scales=list(x=list(rot = 45))
stplot(wind.data, mode = "xt", scales = scales, xlab = NULL,
       col.regions=brewer.pal(11, "RdBu")[-c(10,11)])
stplot(epa.data, mode = "xt", scales = scales, xlab = NULL)


