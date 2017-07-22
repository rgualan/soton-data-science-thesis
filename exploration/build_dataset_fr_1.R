# Refs:
# http://geog.uoregon.edu/GeogR/topics/netcdf-to-raster.html

# Clean environment
rm(list=ls())

# Load libraries
library(sp)
library(rgdal)

## Global parameters  ##########################################################
# Time period: 
# Initial test: 4 months (3 training, 1 validation)
dateA <- as.POSIXct("2015-01-01", format="%Y-%m-%d", tz="UTC")
dateB <- as.POSIXct("2015-03-31", format="%Y-%m-%d", tz="UTC")



## Create coordinates.csv ######################################################
# Read metadata
metadata <- read.delim("data/FR_AQeReporting_2013-2015/FR_2013-2015_metadata.csv", header = T)
metadata$ObservationDateEnd <- as.POSIXct(metadata$ObservationDateEnd,
                                          format="%Y-%m-%d %H:%M:%S",tz = "UTC")
metadata$ObservationDateBegin <- as.POSIXct(metadata$ObservationDateBegin,
                                            format="%Y-%m-%d %H:%M:%S",tz = "UTC")

# Filter data to the relevant time period and pollutant
# Ignore weird stations
stations <- unique(metadata[
  metadata$ObservationDateBegin<="2015-01-01 00:00:00" 
  & (is.na(metadata$ObservationDateEnd) | metadata$ObservationDateEnd>="2015-12-31 12:00:00")
  & metadata$AirPollutant %in% 'PM10'
  & metadata$Latitude>40, 
  c("AirQualityStationEoICode",
    "Longitude","Latitude","Altitude")])
stations <- stations[order(stations$AirQualityStationEoICode),] 
names(stations)[1] <- "Station.ID"
#View(stations)

# Transform longlat to UTM
st <- stations
#plot(st$Longitude, st$Latitude)
coordinates(st) <-~Longitude+Latitude
proj4string(st) <- "+proj=longlat +datum=WGS84"
#plot(coordinates(st))

stUtm <- spTransform(st, CRS("+proj=utm +zone=31 ellps=WGS84  +units=km"))
#plot(coordinates(stUtm))
#View(stUtm)

stations$UTMX <- coordinates(stUtm)[,1]
stations$UTMY <- coordinates(stUtm)[,2]
#View(stations)


# Records by station
# Read the data
d <- read.delim("data/FR_AQeReporting_2013-2015/FR_5_2013-2015_aggregated_timeseries_pm10.csv", header = T)
d$Date <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
names(d)[5] <- "Station.ID"
# Filter
d <- d[d$DataAggregationProcess=="P1D" & d$Date>=dateA & d$Date<=dateB,]

# Insert rows for missing dates/times
fullTs <- seq(min(d$Date),max(d$Date),by="day")
fullNw <- expand.grid(unique(d$Station.ID), fullTs)
names(fullNw) <- c("Station.ID","Date")
#View(fullNw)
d <- merge(d,fullNw,all=TRUE)
#View(d2)

polByStation <- aggregate(AirPollutant~Station.ID,data=d,FUN=length)
names(polByStation)[2] <- "Count"
#View(polByStation)
#View(d[d$Station.ID=="FR01006",])

# Tolerance to missing data
stations2 <- merge(stations, polByStation)
#View(stations2)
numDays <- as.integer((dateB-dateA)+1)

# Disable stations based on the amount of missing data
stations3 <- stations2[stations2$Count>numDays*0.90,-7]
#View(stations3)

# Split stations in training and validation
stations <- stations3

set.seed(111)
stations_validation <- stations[sample(nrow(stations), 37), ]
stations_training <- stations[
  !(stations$Station.ID %in% stations_validation$Station.ID), ]

# Simple plot
plot(stations_training$Longitude, stations_training$Latitude, col="blue")
points(stations_validation$Longitude, stations_validation$Latitude, col="green")
readline("Continue?")


## Second try
## The first attemp failed in converge maybe due to the big amount of stations
## In this attemp, the same number of stations than in the paper, is used
## TODO: Randomly chose stations which are not too close
#stations_training2 <- stations_training[sample(nrow(stations_training),24),]
#stations_validation2 <- stations_validation[sample(nrow(stations_validation),10),]
# Simple plot
plot(stations_training2$Longitude, stations_training2$Latitude, col="blue")
points(stations_validation2$Longitude, stations_validation2$Latitude, col="green")
readline("Continue?")


write.csv(stations_training2, 
          file = "data/FR_AQeReporting_2013-2015/coordinates.csv", row.names = F)
write.csv(stations_validation2, 
          file = "data/FR_AQeReporting_2013-2015/coordinates_val.csv", row.names = F)



## Create data by day ####################################################################
# Read the data
#d <- read.delim("data/FR_AQeReporting_2013-2015/FR_5_2013-2015_aggregated_timeseries_pm10.csv", header = T)
#str(d)
#d$Date <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
# daily mean
#d <- d[d$DataAggregationProcess=="P1D",]

# Filter and Merge with the simple ids
d2 <- d[d$Station.ID %in% stations$Station.ID,]
d2 <- merge(d2, stations[,c("Station.ID","Altitude",
                        "Longitude", "Latitude", "UTMX", "UTMY")], 
            by="Station.ID")

# Variable
names(d2)[which(names(d2)=="AirPollutionLevel")] <- "PM10"

d2 <- d2[order(d2$Date,d2$Station.ID),]
#View(d2)


# Split data in training and testing data
# d2_training <- d2[d2$AirQualityStationEoICode %in% stations_training$AirQualityStationEoICode,]
# d2_validation <- d2[d2$AirQualityStationEoICode %in% stations_validation$AirQualityStationEoICode,]


# Write
write.csv(d2[,c("Station.ID","Date","Altitude","Longitude","Latitude",
                "UTMX","UTMY","PM10")],
          file = "data/FR_AQeReporting_2013-2015/France_data_byday_0.csv", row.names = F)




