# Clean environment #################################################################
rm(list=ls())

# Libraries
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgeos)
library(ggmap)
library(rgdal)

## Stations monitoring PM10  #######################################################
## Build the list of active stations starting from the data
## Relevant period:
## 2015-01-01 to 2015-04-30
## Read data
d <- read.delim("/home/ronald/projects/back_up/data/GB_AQeReporting_2013-2015/GB_5_2013-2015_aggregated_timeseries_pm10.csv", header=T)
d$DatetimeBegin <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
d$DatetimeEnd <- as.POSIXct(d$DatetimeEnd, format="%Y-%m-%d %H:%M:%S",tz="UTC")
length(unique(d$AirQualityStationEoICode)) # 74 stations


# Check what months have the most data
dcount <- d
dcount$MonthYear <- as.factor(format(d$DatetimeBegin, format="%Y-%m"))
countByMonthYear <- aggregate(AirPollutionLevel~MonthYear,data=dcount,FUN=length)
names(countByMonthYear)[2] <- "Count"
#View(countByMonthYear)
par(las=2) # make label text perpendicular to axis
barplot(countByMonthYear$Count, names.arg=countByMonthYear$MonthYear, horiz=T,
        cex.names=0.7) 


## Filter data to match relevant criteria
dateA <- as.POSIXct("2015-01-01", format="%Y-%m-%d", tz="UTC")
dateB <- as.POSIXct("2015-04-30", format="%Y-%m-%d", tz="UTC")
d2 <- d[d$DatetimeBegin>=dateA & d$DatetimeEnd<=dateB 
        & d$DataAggregationProcess == "P1D",
        c("AirQualityStationEoICode","DatetimeBegin","AirPollutionLevel")]
#View(d2)
## How many records per stations
recordsByStation <- aggregate(AirPollutionLevel~AirQualityStationEoICode,
                              data=d2,FUN=length)
names(recordsByStation)[2] <- "Count"
#View(recordsByStation)
minNumRecords <- 0.95*as.integer((dateB-dateA)+1)
## Active stations
stations <- recordsByStation[recordsByStation$Count>minNumRecords,]
#View(stations)
## Complete fields of the stations
stations_md <- read.delim("/home/ronald/projects/back_up/data/GB_AQeReporting_2013-2015/GB_2013-2015_metadata.csv", header=T)
stations_md <- unique(stations_md[stations_md$AirPollutant=="PM10",
                                  c("AirQualityStationEoICode","Projection",
                                    "Altitude","Longitude","Latitude",
                                    "MeasurementType","AirQualityStationType",
                                    "AirQualityStationArea")])
stations <- merge(stations,stations_md,all.x=T)
#View(stations)
# Remove stations outside continental territory
nrow(stations)
# Only 11!


# Plot all the active stations
spStations <- stations
coordinates(spStations) <- ~Longitude+Latitude
proj4string(spStations) <- "+proj=longlat" #TODO: +init=epsg:4979
map('worldHires', c('UK'), xlim=c(-11,2), ylim=c(49,59), mar=rep(1,4))	
points(spStations)
readline("continue?")



## Maps #########################################################
#Refs: 
#http://www.milanor.net/blog/maps-in-r-plotting-data-points-on-a-map/
#https://blog.dominodatalab.com/geographic-visualization-with-rs-ggmaps/

## Plot map of UK and the stations (coloured by type of Area)
mapGB <- get_map(location = 'Great Britain', zoom = 6, maptype = "roadmap")
#save(mapGB, file="data/maps/mapGB.Rdata")
#if(!exists("mapGB")) load("data/maps/mapGB.Rdata")

ggmap(mapGB) +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = AirQualityStationArea), 
             data = stations, alpha = .75, shape=21, size=2)
readline("Continue?")
# Notes:
# Several types of map were tested
# The most relevant types of map for the study are:
# roadmap

## Plot map of France and the stations (coloured by type of station)
ggmap(mapGB) +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = AirQualityStationType), 
             data = stations, alpha = .75, shape=21, size=2)
readline("Continue?")

