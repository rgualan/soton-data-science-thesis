# Clean environment #################################################################
rm(list=ls())

# Libraries
library(raster)
library(rasterVis)
library(maptools)
library(maps)
library(RColorBrewer)
library(lattice)
library(ncdf4)
source("util/my_helper.R")


## Interpolation test #####################################################################
if(F){

  # NEtcdf file
  ncFile <- "data/reanalysis/air.2m.2016.nc"
  nc <- nc_open(ncFile)
  print(nc)
  
  
  # read the netCDF file as a raster layer
  r <- raster(ncFile)
  #r
  #print(r)
  plot(r)
  
  ## USA map
  mapUSA <- map("state", plot=F)
  mapUSA.sp <- map2SpatialLines(mapUSA, proj4string = CRS("+proj=longlat"))
  plot(mapUSA.sp)
  mapUSA.sp <- spTransform(mapUSA.sp, proj4string(r))
  plot(mapUSA.sp)
  ## Stations 
  sites <- readRDS("data/epa/epa_daily/2016/california_ozone_sites.RDS")
  lcc <- CRS("+proj=lcc +x_0=5632642.22547 +y_0=4612545.65137 +lat_0=50 +lon_0=-107 +lat_1=50 +lat_2=50 +ellps=WGS84")
  coordinates(sites)<-~Longitude+Latitude
  proj4string(sites)<-"+proj=longlat"
  rownames(sites@coords) <- sites$Station.Code
  sites.lcc <- spTransform(sites,lcc)
  
  
  mapTheme <- rasterTheme(region = rev(brewer.pal(10, "RdBu")))
  #cutpts <- c(-2.5, -2.0, -1.5, -1, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5)
  plt <- levelplot(r, margin = F, cuts=10, pretty=TRUE, par.settings = mapTheme,
                   main="Layer from NetCDF file") #at=cutpts, 
  plt + 
    layer(sp.lines(mapUSA.sp, col = "black", lwd = 0.5)) +
    layer(sp.points(sites.lcc))
  
  # Example of data extraction
  #extract(r, matrix(c(0.125,43.875),1,2), method='bilinear', fun=mean, na.rm=TRUE)
  #extract(r, matrix(c(0.125*2,43.875+0.125),1,2), method='bilinear', fun=mean, na.rm=TRUE)
  #extract(r, matrix(c(-1.625,44.375),1,2), method='bilinear', fun=mean, na.rm=TRUE)
}



## EDA of the predictoros extracted from the NetCDF files ##############################
d <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
head(d)

# Plots
w<-7; h<-2.5
printPlot(T,"img/eda/reanalysis/narr_ozone.jpeg",w,h,FUN=function(){
  timeSeriesPlotWithInterval(d, "Ozone", "Ozone (ppb)", xLabel = F)
})
printPlot(T,"img/eda/reanalysis/narr_temperature.jpeg",w,h,FUN=function(){
  timeSeriesPlotWithInterval(d, "Temperature", "Temperature (K)", xLabel = F)
})
printPlot(T,"img/eda/reanalysis/narr_rh.jpeg",w,h,FUN=function(){
  timeSeriesPlotWithInterval(d, "RH", "Relative Humidity (%)", xLabel = F)
})
printPlot(T,"img/eda/reanalysis/narr_wind.jpeg",w,h,FUN=function(){
  timeSeriesPlotWithInterval(d, "Wind", "Wind (m/s)", xLabel = F)
})
printPlot(T,"img/eda/reanalysis/narr_rain.jpeg",w,h,FUN=function(){
  timeSeriesPlotWithInterval(d, "Rain", "Rain (Kg/m2)", xLabel = T)
})


## Histograms ###########################################################################
hist(d$Temperature)
hist(d$RH)  ## Bimodal!!!
hist(d$Rain) ## Not normal
hist(log(d$Rain))  ## Heavy tailed
hist(d$Wind)
hist(log(d$Wind)) ## Slightly more normal?


