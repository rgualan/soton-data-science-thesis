## Clean environment #################################################################
rm(list=ls())

## Libraries
library(raster)
library(ncdf4)
source("util/my_helper.R")

## Settings ##########################################################################

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


## Merge datasets #######################################################
# dim(d); dim(temp)
d2 <- merge(d,temp,all=T)
d2 <- merge(d2,rh,all=T)
d2 <- merge(d2,rain,all=T)
d2 <- merge(d2,wind[,-(3:4)],all=T)
# dim(d2)
d2 <- merge(d2,sites[,c("Station.Code","Longitude","Latitude","UTM.X","UTM.Y","Elevation","Location.Setting")])
saveRDS(d2,file="data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
