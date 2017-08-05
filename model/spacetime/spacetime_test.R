## Ref: Spatio-Temporal data in R

## Clean environment
rm(list=ls())

## Libraries
library(gstat)

## Pre-process data ############################################################
data("wind")
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y
proj4string(wind.loc) = "+proj=longlat +datum=WGS84"

head(wind)
wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
wind$jday = as.numeric(format(wind$time, ' %j ' ))
head(wind)


## subtract a smooth time trend of daily means ##################################
stations = 4:15
windsqrt = sqrt(0.5148 * as.matrix(wind[stations])) # knots -> m/s
Jday = 1:366
windsqrt = windsqrt - mean(windsqrt)
daymeans = sapply(split(windsqrt, wind$jday), mean)
meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
velocities = apply(windsqrt, 2, function(x) { x - meanwind })

## match the wind data to its location, by connecting station names to location
## coordinates, and create a spatial points object:
wind.loc = wind.loc[match(names(wind[4:15]), wind.loc$Code),]
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])
rownames(pts) = wind.loc$Station
pts = SpatialPoints(pts, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

## Project the longitude/latitude coordinates and country boundary to UTM zone 29
library(rgdal)
utm29 = CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84")
pts = spTransform(pts, utm29)

## Construct the spatio-temporal object from the space-wide table with velocities:
wind.data = stConstruct(velocities, space = list(values = 1:ncol(velocities)),
                        time = wind$time, SpatialObj = pts, interval = TRUE)
class(wind.data)

## For plotting purposes, we can obtain country boundaries from package maps:
library(maps)
library(mapdata)
library(maptools)
m = map2SpatialLines(
  #map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
  map("worldHires", xlim = c(-11.5,-6.0), ylim = c(51.3,55.0), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84"
m = spTransform(m, utm29)

## For interpolation, we can define a grid over the area:
grd = SpatialPixels(SpatialPoints(makegrid(m, n = 300)),
                    proj4string = proj4string(m))

## We (arbitrarily) restrict observations to those of April 1961:
# wind.data = wind.data[, "1961-04"]
## and choose 10 time points from that period to form the spatio-temporal 
## prediction grid:
n = 10
library(xts)
library(RColorBrewer)
tgrd = seq(min(index(wind.data)), max(index(wind.data)), length=n)
pred.grd = STF(grd, tgrd)

## We will interpolate with a separable exponential covariance model, 
## with ranges 750 km and 1.5 days:
v = vgmST("separable", space = vgm(1, "Exp", 750000), 
          time = vgm(1, "Exp", 1.5 * 3600), sill=0.6)
wind.ST = krigeST(values ~ 1, wind.data, pred.grd, v) # All data: ~10 min
colnames(wind.ST@data) <- "sqrt_speed"


## Then creates the STFDF object with interpolated values
layout = list(list("sp.lines", m, col= ' grey ' ),
              list("sp.points", pts, first=F, cex=.5))
stplot(wind.ST, col.regions=brewer.pal(11, "RdBu")[-c(10,11)],
       at=seq(-1.375,1,by=.25),
       par.strip.text = list(cex=.7), sp.layout = layout)


## Space-time plots
scales=list(x=list(rot = 45))
stplot(wind.data, mode = "xt", scales = scales, xlab = NULL,
       col.regions=brewer.pal(11, "RdBu")[-c(10,11)])


