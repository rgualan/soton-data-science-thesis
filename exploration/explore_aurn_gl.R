## AURN dataset taken from
## http://www.southampton.ac.uk/~sks/pollution_estimates/


# Clean environment #################################################################
rm(list=ls())
par(ask=T)

# Libraries
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgeos)
library(ggmap)
library(rgdal)
library(ggplot2)
library(openair)
library(lattice)
library(RColorBrewer)
library(openair)

## Stations ###################################################################
## Read data
aurn <- read.table("data/aurn/AURN_data_07_11.txt",header=T)
aurn$date <- as.POSIXct(sprintf("%04d-%02d-%02d", aurn$year,aurn$month,aurn$day),tz="GMT")
names(aurn)[1] <- "site"
aurn$site <- as.factor(aurn$site)
nrow(aurn) # 262944


## Filter to only keep data from GL and certain year ##########################
aurn0 <- aurn
#aurn <- aurn[aurn$lon>=-0.567 & aurn$lon<=0.312 &aurn$lat>=51.2 & aurn$lat<=51.8, ]
aurn <- aurn[aurn$lon>=-1.567 & aurn$lon<=1.312 &aurn$lat>=50.2 & aurn$lat<=52.8, ]


## Plot stations by type of station
sites <- unique(aurn[,c(1:3,6)])
nrow(sites) # 62
map('worldHires', 'UK', xlim=c(-8.5,2), ylim=c(50,58.6), mar=rep(0,4))	
points(lat~lon,sites,pch=2,col=type)
legend("topright",legend=unique(sites$type), pch=2, col=unique(sites$type))



## Shapefile plus stations inside GL
## Ref: https://medium.com/towards-data-science/plotting-a-map-of-london-crime-data-using-r-8dcefef1c397
shpGL <- readOGR("/home/ronald/projects/aq-tngapms/data/maps/statistical-gis-boundaries-london/ESRI", 
                 layer = "London_Borough_Excluding_MHW")
proj4string(shpGL)
shpGL.wgs84 <- spTransform(shpGL, CRS("+init=epsg:4326"))

ggplot(shpGL.wgs84) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_point(data = sites, aes(x = lon, y = lat, colour = type)) +
  labs(x = "Longitude", y = "Latitude", title = "Map of Greater London with the borough boundaries") +
  coord_quickmap()
#+scale_colour_manual(values = rainbow(14))




## Time series ###################################################################

## Which is the year having the most data?
obsByMonth <- aggregate(cbind(aurn$obs_no2,aurn$obs_pm10),by=list(ym=format(aurn$date,"%Y-%m")),function(x){sum(!is.na(x))})
barchart(ym~V1,obsByMonth)


obsByYear <- aggregate(cbind(no2=aurn$obs_no2,pm10=aurn$obs_pm10),
                       by=list(year=format(aurn$date,"%Y")),
                       function(x){sum(!is.na(x))})
barchart(year~no2+pm10,obsByYear,auto.key=T,main="Records by year")
# Answer: 2007



## Simple ts plot to assess the results from the aggregation
levelplot(obs_no2 ~ date * site, aurn[aurn$year==2009,], 
          cuts = 10, col.regions = rev(brewer.pal(11, "Spectral")),
          scales = list(y = list(cex = .7)), main="No2 - 2009")

levelplot(obs_no2 ~ date * site, aurn[aurn$year==2007,], 
          cuts = 10, col.regions = rev(brewer.pal(11, "Spectral")),
          scales = list(y = list(cex = .7)), main="No2 - 2007")

## Plot a data concentration as matrix (sites x date) ################
levelplot(obs_no2~date*site, aurn,
          cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
          scales=list(y=list(cex=.7)),main="No2 - 07-11")

## Notes:
## The difference in the amount of data between years is caused by the fact 
## that some stations only cover certain years



## How many active stations by year
tmp <- aggregate(cbind(obs_no2,obs_pm10)~site+year,aurn,function(x){sum(!is.na(x))},na.action=na.pass)
tmp$obs_no2p <- tmp$obs_no2/365
tmp$obs_pm10p <- tmp$obs_pm10/365
tmp$obs_no2.on <- tmp$obs_no2p>0.80
tmp$obs_pm10.on <- tmp$obs_pm10p>0.80
#View(tmp)
sitesByYear <- aggregate(cbind(obs_no2.on,obs_pm10.on)~year,tmp,sum)
sitesByYear$year <- as.factor(sitesByYear$year)
barchart(year~obs_no2.on+obs_pm10.on,sitesByYear,auto.key=T,
         main="Active stations by year",
         key=simpleKey(c("Active No2","Active Pm10"), points=F, rectangles=T, cex=0.7))




## Aggregated time series
dDm <- aggregate(cbind(obs_no2,obs_pm10)~date,aurn,mean) # daily mean
#View(dDm)
#summaryPlot(dDm, period="months")
timePlot(dDm, c("obs_pm10","obs_no2"))



## End
par(ask=F)