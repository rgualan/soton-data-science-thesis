# Clean environment #################################################################
rm(list=ls())

# Libraries
library(openair) # Provides functions for easily ploting time series
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgdal)

#### Piemonte plus stations #########################################################
border <- read.csv("data/Piemonte/Piemonte_borders.csv", header = T)
stations <- read.csv("data/Piemonte/coordinates.csv", header = T)
stations_val <- read.csv("data/Piemonte/coordinates_validation.csv", header = T)
#str(stations)
# coordinates(stations) <- c("UTMX", "UTMY")
# proj4string(stations) <- CRS("+proj=utm")  
# stations_longlat <- spTransform(stations, CRS("+proj=longlat"))
# head(cbind(coordinates(stations), coordinates(stations_longlat)))
#map('worldHires', c('Italy'), xlim=c(7.8,16.3), ylim=c(37.6,47.5), mar=rep(1,4))	
#map('worldHires', c('Italy'), xlim=c(5,20), ylim=c(30,50), mar=rep(1,4))	

plot(border, type="l")
points(stations$UTMX,stations$UTMY, pch=18, col="blue")
points(stations_val$UTMX,stations_val$UTMY, pch=18, col="green")
readline("Continue?")

# Time series ###############################################################
# Read the data  
d <- read.csv("data/Piemonte/Piemonte_data_byday.csv", header=T)
d$date <- as.POSIXct(d$Date,format="%d/%m/%y", tz="UTC")
d$site <- as.factor(d$Station.ID)
#str(d)

# Summary plots
for(p in c("WS","TEMP","HMIX","PREC","EMI","PM10")){
  print(p)  
  suppressWarnings(
    summaryPlot(d[d$site %in% 1:5,c("date","site",p)],period="months",type="site",main=p)
  )
  readline("Continue?")
}

# Time plots
for(p in c("WS","TEMP","HMIX","PREC","EMI","PM10")){
  print(p)  
  timePlot(d, pollutant=p, type="site", main=p)
  readline("Continue?")
}





