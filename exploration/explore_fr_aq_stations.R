# Clean environment #################################################################
rm(list=ls())

# Libraries
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgeos)
library(ggmap)
library(rgdal)

## Stations monitoring PM10 near Paris ###########################################
## Build the list of active stations starting from the data
## Relevant period: 4 months
## 2015-01-01 to 2015-04-30
## Read data
d <- read.delim("data/FR_AQeReporting_2013-2015/FR_5_2013-2015_aggregated_timeseries_pm10.csv", header=T)
d$DatetimeBegin <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
d$DatetimeEnd <- as.POSIXct(d$DatetimeEnd, format="%Y-%m-%d %H:%M:%S",tz="UTC")
length(unique(d$AirQualityStationEoICode))

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
stations_md <- read.delim("data/FR_AQeReporting_2013-2015/FR_2013-2015_metadata.csv", header=T)
stations_md <- unique(stations_md[stations_md$AirPollutant=="PM10",
                                  c("AirQualityStationEoICode","Projection",
                                    "Altitude","Longitude","Latitude",
                                    "MeasurementType","AirQualityStationType",
                                    "AirQualityStationArea")])
stations <- merge(stations,stations_md,all.x=T)
#View(stations)
# Remove stations outside continental territory
stations <- stations[stations$Latitude>40,] 
nrow(stations)


# Plot all the active stations
st2 <- stations
coordinates(st2) <- ~Longitude+Latitude
proj4string(st2) <- "+proj=longlat" #TODO: +init=epsg:4979
map('worldHires', c('France'), xlim=c(-4.7,8.1), ylim=c(42.4,51), mar=rep(1,4))	
points(st2)

# Highlight stations close to/inside PARIS
distanceToParis <- function(x) {
  paris <- c(2.348410,48.85800)
  point <- x[5:6]
  p <- rbind(paris, point)
  dst <- dist(p)
  return(dst)  
}

stations$distanceToParis <- apply(stations, 1, distanceToParis)
#View(stations)
#library(scales)
points(stations[stations$distanceToParis<1.5, c("Longitude","Latitude")], 
       pch=16,col=alpha("green",0.5))

nrow(stations[stations$distanceToParis<1.5,])
readline("continue?")



## Maps #########################################################
#Refs: 
#http://www.milanor.net/blog/maps-in-r-plotting-data-points-on-a-map/
#https://blog.dominodatalab.com/geographic-visualization-with-rs-ggmaps/

## Plot map of France and the stations (coloured by type of Area)
# zoom=6
#mapFrance <- get_map(location = 'France', zoom = 6, maptype = "roadmap")
#save(mapFrance, file="data/maps/mapFrance.Rdata")
if(!exists("mapFrance")) load("data/maps/mapFrance.Rdata")

ggmap(mapFrance) +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = AirQualityStationArea), 
             data = stations, alpha = .75, shape=21, size=2)
readline("Continue?")
# Notes:
# Several types of map were tested
# The most relevant types of map for the study are:
# roadmap

## Plot map of France and the stations (coloured by type of station)
ggmap(mapFrance) +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = AirQualityStationType), 
             data = stations, alpha = .75, shape=21, size=2)
readline("Continue?")

## Plot map of Paris and the stations
#mapParis <- get_map(location = 'Paris', zoom = "auto", maptype="roadmap")
#save(mapParis, file="data/maps/mapParis.Rdata")
if(!exists("mapParis")) load("data/maps/mapParis.Rdata")

ggmap(mapParis) +
  geom_point(aes(x=Longitude, y=Latitude, fill = AirQualityStationArea), 
             data = stations, alpha = .75, shape=21, size=2)
readline("Continue?")



## Plot map from a shapefile
# Ref: http://www.kevjohnson.org/making-maps-in-r/
tract <- readOGR(dsn = "data/FR_AQeReporting_2013-2015/shapefiles/regions-20170102-shp/", 
                 layer = "regions-20170102")
tract <- fortify(tract, region="insee")

ggplot() +
  geom_polygon(data = tract, aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue", size = 0.25) +
  geom_point(data = stations, 
             aes(x=Longitude, y=Latitude, fill = AirQualityStationArea), 
             alpha = .75, shape=21, size=2)
#ggsave(p, file = "map1.png", width = 6, height = 4.5, type = "cairo-png")

