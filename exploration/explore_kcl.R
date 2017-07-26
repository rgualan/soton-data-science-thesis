# Clean environment #################################################################
rm(list=ls())

# Libraries
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgeos)
library(ggmap)
library(rgdal)

## Stations ###################################################################
## Read data
load("data/kcl/sites.RData")
nrow(sites) # 900

## Plot all stations
plot(Latitude~Longitude,sites,pch=2,col=Classification)
#plot(os_grid_y~os_grid_x,sites)
readline("Continue?")

## Plot active stations
dateA <- as.POSIXct("2015-01-01", format="%Y-%m-%d", tz="UTC")
dateB <- as.POSIXct("2015-04-30", format="%Y-%m-%d", tz="UTC")
sites2 <- sites[sites$OpeningDate<=dateA & 
                  (sites$ClosingDate>=dateB | is.na(sites$ClosingDate)),]
nrow(sites2)
plot(Latitude~Longitude,sites2,pch=2,col=Classification)
readline("Continue?")

## Plot a map with the active stations ###################################
## Refs: 
## http://www.milanor.net/blog/maps-in-r-plotting-data-points-on-a-map/
## https://blog.dominodatalab.com/geographic-visualization-with-rs-ggmaps/

## Plot map of London
#mapLondon <- get_map(location = 'London', zoom = "auto", maptype = "roadmap")
#save(mapLondon, file="data/maps/mapLondon.Rdata")
if(!exists("mapLondon")) load("data/maps/mapLondon.Rdata")

ggmap(mapLondon) +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = Classification), 
             data = sites2, alpha = .75, shape=21, size=2)
readline("Continue?")

## How many stations in London/Greater London?
sites3 <- sites2[sites2$Longitude >= -0.567 & sites2$Longitude <= 0.312
                 & sites2$Latitude >= 51.2 & sites2$Latitude <= 51.8,]
nrow(sites3)



## Explore meteorological dataset
load("data/kcl/metData.RData")
met2 <- met[met$date>=dateA & met$date<=dateB,]
nrow(met2)
#View(met2)
summaryPlot(met2, period = "months")
