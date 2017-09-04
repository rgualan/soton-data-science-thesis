## Clean environment 
rm(list=ls())

## Libraries
library(raster)
library(ncdf4)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Read data #########################################################################
d <- readRDS("data.big/epa/epa_daily/2016/california_ozone_2.RDS")
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
temp <- extractTimeSeriesFromNc("data.big/reanalysis/air.2m.2016.nc", 
                                dateA, dateB, sites.lcc@coords, 
                                c("Station.Code","Date","Temperature"))
rh <- extractTimeSeriesFromNc("data.big/reanalysis/rhum.2m.2016.nc",
                              dateA, dateB, sites.lcc@coords,
                              c("Station.Code","Date","RH"))
rain <- extractTimeSeriesFromNc("data.big/reanalysis/apcp.2016.nc",
                              dateA, dateB, sites.lcc@coords,
                              c("Station.Code","Date","Rain"))
uwnd <- extractTimeSeriesFromNc("data.big/reanalysis/uwnd.10m.2016.nc",
                                dateA, dateB, sites.lcc@coords,
                                c("Station.Code","Date","Uwnd"))
vwnd <- extractTimeSeriesFromNc("data.big/reanalysis/vwnd.10m.2016.nc",
                                dateA, dateB, sites.lcc@coords,
                                c("Station.Code","Date","Vwnd"))
vwnd <- extractTimeSeriesFromNc("data.big/reanalysis/vwnd.10m.2016.nc",
                                dateA, dateB, sites.lcc@coords,
                                c("Station.Code","Date","Vwnd"))
vwnd <- extractTimeSeriesFromNc("data.big/reanalysis/vwnd.10m.2016.nc",
                                dateA, dateB, sites.lcc@coords,
                                c("Station.Code","Date","Vwnd"))
dewPoint <- extractTimeSeriesFromNc("data.big/reanalysis/dpt.2m.2016.nc",
                                dateA, dateB, sites.lcc@coords,
                                c("Station.Code","Date","Dew.Point"))
waterEvaporation <- extractTimeSeriesFromNc("data.big/reanalysis/evap.2016.nc",
                                    dateA, dateB, sites.lcc@coords,
                                    c("Station.Code","Date","Water.Evap"))
heatFlux <- extractTimeSeriesFromNc("data.big/reanalysis/gflux.2016.nc",
                                            dateA, dateB, sites.lcc@coords,
                                            c("Station.Code","Date","Heat.Flux"))
geoHeight <- extractTimeSeriesFromNc("data.big/reanalysis/hgt.sfc.nc",
                                    dateA, dateB, sites.lcc@coords,
                                    c("Station.Code","Date","Geop.Height"), singleLevel=T)
geoHeightTropo <- extractTimeSeriesFromNc("data.big/reanalysis/hgt.tropo.2016.nc",
                                     dateA, dateB, sites.lcc@coords,
                                     c("Station.Code","Date","Geop.Height.Tropo"))
latentHeatFlux <- extractTimeSeriesFromNc("data.big/reanalysis/lhtfl.2016.nc",
                                          dateA, dateB, sites.lcc@coords,
                                          c("Station.Code","Date","Lat.Heat.Flux"))
tropopausePressure <- extractTimeSeriesFromNc("data.big/reanalysis/pres.tropo.2016.nc",
                                              dateA, dateB, sites.lcc@coords,
                                              c("Station.Code","Date","Tropo.Press"))
pressureMSL <- extractTimeSeriesFromNc("data.big/reanalysis/prmsl.2016.nc",
                                       dateA, dateB, sites.lcc@coords,
                                       c("Station.Code","Date","Press.MSL"))
vegetation <- extractTimeSeriesFromNc("data.big/reanalysis/veg.2016.nc",
                                       dateA, dateB, sites.lcc@coords,
                                       c("Station.Code","Date","Vegetation"))

## Calculate Wind speed (Pythagoras theorem)
wind <- merge(uwnd,vwnd)
wind$Wind <- sqrt(wind$Uwnd^2+wind$Vwnd^2)

## Merge datasets #######################################################
# dim(d); dim(temp)
d2 <- merge(d,temp,all=T)
d2 <- merge(d2,rh,all=T)
d2 <- merge(d2,rain,all=T)
d2 <- merge(d2,wind[,-(3:4)],all=T)
d2 <- merge(d2,dewPoint,all=T)
d2 <- merge(d2,waterEvaporation,all=T)
d2 <- merge(d2,heatFlux,all=T)
d2 <- merge(d2,geoHeight,all=T)
d2 <- merge(d2,geoHeightTropo,all=T)
d2 <- merge(d2,latentHeatFlux,all=T)
d2 <- merge(d2,tropopausePressure,all=T)
d2 <- merge(d2,pressureMSL,all=T)
d2 <- merge(d2,vegetation,all=T)
# dim(d2)

## Save #################################################################
d2 <- merge(d2,sites[,c("Station.Code","Longitude","Latitude","UTM.X","UTM.Y","Elevation","Location.Setting")])
saveRDS(d2,file="data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
# d2 <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")

## TODO: Temporal solution to create california_ozone_plus_rcov_3.RDS
d3 <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")  ## Ozone missing values were replaced
d3$Ozone2 <- d3$Ozone
d4 <- merge( d2, d3[,c("Station.Code","Date","Ozone2")] )
#cor(d4$Ozone2, d4$Ozone, use="pairwise.complete.obs")
d4$Ozone <- d4$Ozone2
d4 <- d4[,c("Station.Code","Date","Ozone","Temperature",
            "RH","Rain","Wind","Dew.Point","Water.Evap",
            "Heat.Flux","Geop.Height","Geop.Height.Tropo",
            "Lat.Heat.Flux","Tropo.Press","Press.MSL",
            "Vegetation","Longitude","Latitude","UTM.X","UTM.Y",
            "Elevation","Location.Setting")]
#d4 <- d4[order(d4$Station.Code,d4$Date),]
#head(d4)

## Some stations couldnt be interpolated
#summary(d4)
MIDs <- unique(d4$Station.Code[is.na(d4$Heat.Flux)])
#plotStations(F,MIDs)
# In water?
# Bye bye to these stations
d4 <- d4[!d4$Station.Code %in% MIDs, ]
d4 <- d4[d4$Station.Code!="023-1005", ] ##TODO!!!


saveRDS(d4, "data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")  
#d4 <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")  

