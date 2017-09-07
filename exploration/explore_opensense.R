## Clean environment #################################################################
rm(list=ls())

## Libraries
library(lattice)
library(RColorBrewer)
library(ggplot2)

## Read data #########################################################################
## Variable: PM10
## Build the list of active stations starting from the data
## Relevant period: 
## 2015-01-01 to 2015-04-30
## Read data
d <- read.csv("/home/ronald/projects/back_up/data/opensense/ozon_tram1_14102011_14012012.csv", header=T, 
              stringsAsFactors = F)
#d$Datetime <- as.POSIXct(as.numeric(as.character(d$generation_time)),origin="1970-01-01",tz="GMT")
d$Datetime <- as.POSIXct(as.numeric(as.character(d$generation_time))/ 1000.0,
                         origin="1970-01-01",tz="GMT")
d$longitude <- as.double(d$longitude)
d$latitude <- as.double(d$latitude)
## Sites
sites.coords <- unique(d[,c("longitude","latitude")])
dim(d)
dim(sites.coords)

## Plot ##############################################################################
plot(latitude~longitude, sites.coords, pch=".")

