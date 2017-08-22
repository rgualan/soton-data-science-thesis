## Clean environment 
rm(list=ls())

## Libraries
library(raster)
library(ncdf4)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Read data #########################################################################
#d <- readRDS("data/epa/epa_daily/2016/california_ozone.RDS")
d <- readRDS("data/epa/epa_daily/2016/california_ozone_2.RDS")
sites <- getSites(d)

## Relevant period
dateA <- convertStringToPOSIXct("2016-01-01")
dateB <- convertStringToPOSIXct("2016-12-31")
numDays <- as.integer((dateB-dateA)+1)

## Transform coordinates
sites.sp <- sites
lcc <- CRS("+proj=lcc +x_0=5632642.22547 +y_0=4612545.65137 +lat_0=50 +lon_0=-107 +lat_1=50 +lat_2=50 +ellps=WGS84")
coordinates(sites.sp)<-~Longitude+Latitude
proj4string(sites.sp)<-"+proj=longlat"
rownames(sites.sp@coords) <- sites.sp$Station.Code
sites.lcc <- spTransform(sites.sp,lcc)

## Extract covariates ###################################################
temp <- extractTimeSeriesFromNc("data/reanalysis/air.2m.2016.nc", 
                                dateA, dateB, sites.lcc@coords, 
                                c("Station.Code","Date","Temperature"))
rh <- extractTimeSeriesFromNc("data/reanalysis/rhum.2m.2016.nc",
                              dateA, dateB, sites.lcc@coords,
                              c("Station.Code","Date","RH"))
rain <- extractTimeSeriesFromNc("data/reanalysis/apcp.2016.nc",
                              dateA, dateB, sites.lcc@coords,
                              c("Station.Code","Date","Rain"))
uwnd <- extractTimeSeriesFromNc("data/reanalysis/uwnd.10m.2016.nc",
                                dateA, dateB, sites.lcc@coords,
                                c("Station.Code","Date","Uwnd"))
vwnd <- extractTimeSeriesFromNc("data/reanalysis/vwnd.10m.2016.nc",
                                dateA, dateB, sites.lcc@coords,
                                c("Station.Code","Date","Vwnd"))
## Calculate Wind speed (Pythagoras theorem)
wind <- merge(uwnd,vwnd)
wind$Wind <- sqrt(wind$Uwnd^2+wind$Vwnd^2)

## Transformation #######################################################
wind$sqrtWind <- sqrt(wind$Wind)
rain$logRain <- log(rain$Rain)
## rain$logRain[rain$logRain==-Inf] <- 0 ## TODO!!!

## Merge datasets #######################################################
# dim(d); dim(temp)
d2 <- merge(d,temp,all=T)
d2 <- merge(d2,rh,all=T)
d2 <- merge(d2,rain,all=T)
d2 <- merge(d2,wind[,-(3:4)],all=T)
# dim(d2)

## EDA ##################################################################

## Temperature
printPlot(paper,"img/eda/reanalysis/histogram_temperature.jpeg",5,5,FUN=function(){
  hist(d2$Temperature, main="", xlab="Air Temperature at 2m (K)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_temperature.jpeg", 6, 7, FUN=function(){
  simple_heatmap(d2, "Temperature")
})
printPlot(paper, "img/eda/reanalysis/ts_temperature.jpeg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(d2, "Temperature", "Temperature (K)")
})

## Relative Humidity
printPlot(paper,"img/eda/reanalysis/histogram_rh.jpeg",5,5,FUN=function(){
  hist(d2$RH, main="", xlab="Relative humidity (%)")
})
## Notes: bimodal distribution
## For LDA you really do not want to remove the bimodality, as it will reduce the discrimination of your analysis
## https://stats.stackexchange.com/questions/209241/what-transformation-should-i-use-for-a-bimodal-distribution
printPlot(paper, "img/eda/reanalysis/heatmap_rh.jpeg", 6, 7, FUN=function(){
  #heatmapPlusTs(d2, "RH", "Relative humidity (%)")
  simple_heatmap(d2, "RH")
})
printPlot(paper, "img/eda/reanalysis/ts_rh.jpeg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(d2, "RH", "Rel. humidity (%)")
})

## Wind speed
printPlot(paper,"img/eda/reanalysis/histogram_ws.jpeg",5,5,FUN=function(){
  hist(d2$Wind, main="", xlab="Wind speed (m/s)")
})
printPlot(paper,"img/eda/reanalysis/histogram_ws_sqrt.jpeg",5,5,FUN=function(){
  hist(sqrt(d2$Wind), main="", xlab="sqrt(Wind speed)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_ws.jpeg", 6, 7, FUN=function(){
  simple_heatmap(d2, "Wind")
})
printPlot(paper, "img/eda/reanalysis/ts_ws", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(d2, "Wind", "Wind speed (m/s)")
})

## Rain
printPlot(paper,"img/eda/reanalysis/histogram_rain.jpeg",5,5,FUN=function(){
  hist(d2$Rain, main="", xlab="Precipitation amount(kg/m^2)")
})
printPlot(paper,"img/eda/reanalysis/histogram_rain_log.jpeg",5,5,FUN=function(){
  hist(log(d2$Rain), main="", xlab="log(Precipitation)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_rain.jpeg", 6, 7, FUN=function(){
  simple_heatmap(d2, "Rain")
})
printPlot(paper, "img/eda/reanalysis/ts_rain", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(d2, "Rain", "Precipitation (Kg/m^2)")
})

## Save #################################################################
d2 <- merge(d2,sites[,c("Station.Code","Longitude","Latitude","UTM.X","UTM.Y","Elevation","Location.Setting")])
saveRDS(d2,file="data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
