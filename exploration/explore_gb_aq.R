# Clean environment #################################################################
rm(list=ls())

# Libraries
library(openair)
library(maps)       # Provides functions that let us plot the maps
library(mapdata)    # Contains the hi-resolution points that mark out the countries.

#### UK plus stations ################################################################
metadata <- read.delim("data/GB_AQeReporting_2013-2015/GB_2013-2015_metadata.csv", header = T)
metadata$ObservationDateEnd <- as.POSIXct(metadata$ObservationDateEnd,
                                          format="%Y-%m-%d %H:%M:%S",tz = "UTC")
metadata$ObservationDateBegin <- as.POSIXct(metadata$ObservationDateBegin,
                                          format="%Y-%m-%d %H:%M:%S",tz = "UTC")
head(metadata$ObservationDateEnd[!is.na(metadata$ObservationDateEnd)])
head(metadata$ObservationDateBegin[!is.na(metadata$ObservationDateBegin)])
str(metadata)

# Pollutants
p <- unique(metadata$AirPollutant)
length(p)
p
cp <- c("NO2","PM10","PM2.5","SO2","CO","NO","O3") # chosen pollutants

# All stations
all_stations <- unique(metadata[,c("AirQualityStationEoICode","Projection",
                               "Longitude","Latitude","Altitude",
                               "AirQualityStationType","AirQualityStationArea")])
coords <- as.matrix(unique(cbind(all_stations[, c("Longitude","Latitude")])))

map('worldHires',
    c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
    xlim=c(-11,2), ylim=c(50,58.6), mar=rep(1,4))	
points(coords, pch = 18, col = "lightgreen")

# Active stations
stations <- unique(metadata[
  metadata$ObservationDateBegin<="2015-01-01 00:00:00" 
  & (is.na(metadata$ObservationDateEnd) | metadata$ObservationDateEnd>="2015-12-31 12:00:00")
  & metadata$AirPollutant %in% cp, 
  c("ObservationDateBegin","ObservationDateEnd","AirPollutantCode",
    "AirQualityStationEoICode","Projection",
    "Longitude","Latitude","Altitude",
    "AirQualityStationType","AirQualityStationArea")])
stations$AirPollutantCode <- as.integer(gsub(".*/([0-9]+)$", "\\1",stations$AirPollutantCode))
#View(stations)
coords <- as.matrix(unique(cbind(stations[, c("Longitude","Latitude","AirPollutantCode")])))
map('worldHires',
    c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
    xlim=c(-11,2), ylim=c(50,58.6), mar=rep(1,4))	
#points(coords, pch = 4, col = "blue", cex=0.5)
points(coords[,1:2], pch = 4, col = coords[,3], cex=0.5)
readline("Continue?")

# Active stations 2
tmp <- as.data.frame(coords)
coords2 <- aggregate(AirPollutantCode~Longitude+Latitude,data=tmp,FUN=length)
map('worldHires',
    c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
    xlim=c(-11,2), ylim=c(50,58.6), mar=rep(1,4))	
points(coords2[,1:2], pch = as.character(coords2[,3]), col = coords2[,3])


par(.pardefault)


# Time series ###############################################################
# Read the data
#d <- read.delim("data/GB_AQeReporting_2013-2015/GB_1_2013-2015_aggregated_timeseries_so2.csv", header = T)
d <- read.delim("data/GB_AQeReporting_2013-2015/GB_5_2013-2015_aggregated_timeseries_pm10.csv", header = T)
d <- d[d$DataAggregationProcess=="P1D",]
str(d)
d$date <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
summary(d)
d$site <- as.factor(d$AirQualityStationEoICode)
summaryPlot(d[d$site %in% head(unique(d$site)), c("date","site","AirPollutionLevel")],
            pollutant = "AirPollutionLevel", type="site")
# timePlot(d[d$site %in% head(unique(d$site)), c("date","site","AirPollutionLevel")], 
#          pollutant = "AirPollutionLevel", type="site")
