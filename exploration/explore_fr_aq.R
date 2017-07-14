# Clean environment #################################################################
rm(list=ls())

# Libraries
library(openair)
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.

#### UK plus stations ################################################################
metadata <- read.delim("data/FR_AQeReporting_2013-2015/FR_2013-2015_metadata.csv", header = T)
metadata$ObservationDateEnd <- as.POSIXct(metadata$ObservationDateEnd,
                                          format="%Y-%m-%d %H:%M:%S",tz = "UTC")
metadata$ObservationDateBegin <- as.POSIXct(metadata$ObservationDateBegin,
                                          format="%Y-%m-%d %H:%M:%S",tz = "UTC")
#head(metadata$ObservationDateEnd[!is.na(metadata$ObservationDateEnd)])
#head(metadata$ObservationDateBegin[!is.na(metadata$ObservationDateBegin)])
#str(metadata)

# All stations
all_stations <- unique(metadata[,c("AirQualityStationEoICode","Projection",
                               "Longitude","Latitude","Altitude",
                               "AirQualityStationType","AirQualityStationArea")])
coords <- as.matrix(unique(cbind(all_stations[, c("Longitude","Latitude")])))
map('worldHires', c('France'), xlim=c(-4.7,8.1), ylim=c(42.4,51), mar=rep(1,4))	
points(coords, pch = 18, col = "lightgreen")
readline("Continue?")


# Number of stations by monitored pollutants
BP <- c("PM10","O3","NO2","SO2","CO","PM2.5") # Basic pollutants
tmp <- metadata[metadata$AirPollutant %in% BP,c("AirQualityStationEoICode","AirPollutant")]
tmp2 <- tmp[order(tmp$AirQualityStationEoICode,tmp$AirPollutant),] 
tmp3 <- aggregate(AirPollutant~AirQualityStationEoICode,data=tmp2,FUN=paste, collapse="-")
#View(tmp3)
polByStation <- aggregate(AirQualityStationEoICode~AirPollutant,data=tmp3,FUN=length)
#View(polByStation)
# The most monitored combination by station is NO2-O3-PM10 
BP2 <- c("NO2","O3","PM10")


# Active stations
stations <- unique(metadata[
  metadata$ObservationDateBegin<="2015-01-01 00:00:00" 
  & (is.na(metadata$ObservationDateEnd) | metadata$ObservationDateEnd>="2015-12-31 12:00:00")
  & metadata$AirPollutant %in% BP2, 
  c("AirPollutantCode",
    "AirQualityStationEoICode","Projection",
    "Longitude","Latitude","Altitude",
    "AirQualityStationType","AirQualityStationArea")])
stations$AirPollutantCode <- as.integer(gsub(".*/([0-9]+)$", "\\1",stations$AirPollutantCode)) # extract the code from
#View(stations)
stationsAgg <- aggregate(AirPollutantCode~AirQualityStationEoICode,data=stations,FUN=length)
stationIds <- stationsAgg$AirQualityStationEoICode[stationsAgg$AirPollutantCode==3]
stations2 <- stations[stations$AirQualityStationEoICode %in% stationIds,]
#View(stations2)

coords <- as.matrix(unique(cbind(stations2[, c("Longitude","Latitude","AirPollutantCode")])))
map('worldHires', c('France'), xlim=c(-4.7,8.1), ylim=c(42.4,51), mar=rep(1,4))	
points(coords[,1:2], pch = 4, col = "blue", cex=0.5)
readline("Continue?")


# Stations monitoring only O3
st <- tmp3$AirQualityStationEoICode[tmp3$AirPollutant=="O3"]
coords <- as.matrix(unique(cbind(stations[stations$AirQualityStationEoICode %in% st, 
                                          c("Longitude","Latitude")])))
points(coords[,1:2], pch = 4, col = "green", cex=0.5)

# Stations monitoring only NO2-PM10
st <- tmp3$AirQualityStationEoICode[tmp3$AirPollutant=="NO2-PM10"]
coords <- as.matrix(unique(cbind(stations[stations$AirQualityStationEoICode %in% st, 
                                          c("Longitude","Latitude")])))
points(coords[,1:2], pch = 4, col = "orange", cex=0.5)


par(.pardefault)
readline("Continue?")

# Time series ###############################################################
# Read the data
#d <- read.delim("data/FR_AQeReporting_2013-2015/FR_5_2013-2015_aggregated_timeseries_pm10.csv", header = T)
d <- read.delim("data/FR_AQeReporting_2013-2015/FR_7_2013-2015_aggregated_timeseries_o3.csv", header = T)
#d <- read.delim("data/FR_AQeReporting_2013-2015/FR_8_2013-2015_aggregated_timeseries_no2.csv", header = T)
#str(d)
d$date <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
d$site <- d$AirQualityStationEoICode
#simpleFilter <- c("FR01001","FR01005","FR01006","FR01009","FR01011")
simpleFilter <- head(unique(d$site))
summaryPlot(d[d$site %in% simpleFilter,
              c("date","site","AirPollutionLevel")])
#timePlot(d[d$site %in% simpleFilter,], pollutant = "AirPollutionLevel", type="site")



# # Create image files of the time series plots ##########################################
# # Create time series plots for each station monitoring PM10
# # To asess the presence of missing values
# d <- read.delim("data/FR_AQeReporting_2013-2015/FR_5_2013-2015_aggregated_timeseries_pm10.csv", header = T)
# d$date <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
# d$site <- d$AirQualityStationEoICode
# for(id in stations2$AirQualityStationEoICode){
#     print(id)
#     d2 <- d[d$site==id, c("date","site","AirPollutionLevel")]
#     if(nrow(d2)>10){      
#       jpeg(filename = paste0("data/FR_AQeReporting_2013-2015/img/PM10/",id,".jpeg"))
#       summaryPlot(d2)
#       dev.off()
#     }
#     #TODO: Check why some stations which reported to monitor PM10 in metadata
#     #return empty or less than 10 rows
# }
# 
# # Create time series plots for each station monitoring O3 
# # To asess the presence of missing values
# d <- read.delim("data/FR_AQeReporting_2013-2015/FR_7_2013-2015_aggregated_timeseries_o3.csv", header = T)
# d$date <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
# d$site <- d$AirQualityStationEoICode
# for(id in stations2$AirQualityStationEoICode){
#   print(id)
#   d2 <- d[d$site==id, c("date","site","AirPollutionLevel")]
#   if(nrow(d2)>10){      
#     jpeg(filename = paste0("data/FR_AQeReporting_2013-2015/img/O3/",id,".jpeg"))
#     summaryPlot(d2)
#     dev.off()
#   }
# }
# 
# # Create time series plots for each station monitoring NO2
# # To asess the presence of missing values
# d <- read.delim("data/FR_AQeReporting_2013-2015/FR_8_2013-2015_aggregated_timeseries_no2.csv", header = T)
# d$date <- as.POSIXct(d$DatetimeBegin, format="%Y-%m-%d %H:%M:%S",tz="UTC")
# d$site <- d$AirQualityStationEoICode
# for(id in stations2$AirQualityStationEoICode){
#   print(id)
#   d2 <- d[d$site==id, c("date","site","AirPollutionLevel")]
#   if(nrow(d2)>10){      
#     jpeg(filename = paste0("data/FR_AQeReporting_2013-2015/img/NO2/",id,".jpeg"))
#     summaryPlot(d2)
#     dev.off()
#   }
# }

# For high quality plots
# png("plot3.png", width = 4, height = 6, units = 'in', res = 300)
# summaryPlot(d2)
# dev.off()



