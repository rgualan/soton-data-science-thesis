# Clean environment #################################################################
rm(list=ls())

# Libraries
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgeos)
library(ggmap)
library(rgdal)
library(ggplot2)
library(openair)

## Stations ###################################################################
## Read data
aurn <- read.table("data/aurn/AURN_data_07_11.txt",header=T)
aurn$date <- as.POSIXct(sprintf("%04d-%02d-%02d", aurn$year,aurn$month,aurn$day),tz="GMT")
names(aurn)[1] <- "site"
aurn$site <- as.factor(aurn$site)
nrow(aurn) # 262944


## Plot all stations by type of station
sites <- unique(aurn[,c(1:3,6)])
map('worldHires', 'UK', xlim=c(-8.5,2), ylim=c(50,58.6), mar=rep(0,4))	
points(lat~lon,sites,pch=2,col=type)
legend("topright",legend=unique(sites$type), pch=2, col=unique(sites$type))
readline("Continue?")


## Highlight stations inside Greater London
sitesGL <- sites[sites$lon >= -0.5 & sites$lon <= 0.4
                 & sites$lat >= 51 & sites$lat <= 51.7,]
nrow(sitesGL)
map('worldHires', c('UK'), xlim=c(-8.5,2), ylim=c(50,58.6), mar=rep(0,4))	
points(lat~lon,sites,pch=".",col="blue",cex=2)
points(lat~lon,sitesGL,pch="*",col="green",cex=1)
readline("Continue?")

## Only stations inside Greater London
par(.pardefault)
plot(lat~lon,sitesGL,pch=2,col="green",cex=1)
readline("Continue?")

## Shapefile plus stations inside GL
## Ref: https://medium.com/towards-data-science/plotting-a-map-of-london-crime-data-using-r-8dcefef1c397
shpGL <- readOGR("/home/ronald/projects/aq-tngapms/data/maps/statistical-gis-boundaries-london/ESRI", 
                 layer = "London_Borough_Excluding_MHW")
proj4string(shpGL)
shpGL.wgs84 <- spTransform(shpGL, CRS("+init=epsg:4326"))

ggplot(shpGL.wgs84) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  labs(x = "Longitude", y = "Latitude", title = "Map of Greater London with the borough boundaries") +
  geom_point(data = sitesGL, aes(x = lon, y = lat, colour = type)) +
  coord_quickmap()
#+scale_colour_manual(values = rainbow(14))
readline("Continue?")


## Time series ###################################################################
aqum <- aurn[,c("site","date","aqum_no2")]
aqumSub <- aqum[aqum$site %in% sites$site[sample(nrow(sites),5)],]
summaryPlot(aqumSub, main="AQUM_NO2")
readline("Continue?")

aqum <- aurn[,c("site","date","obs_no2")]
aqumSub <- aqum[aqum$site %in% sites$site[sample(nrow(sites),5)],]
summaryPlot(aqumSub, main="OBS_NO2")
