## Clean environment
rm(list=ls())

## Load libraries
library(spTimer)
library(sp)
library(raster)

##### MAAQ Background #####
## Background pollution maps at 1x1 km resolution are modelled each year under 
## Defra's Modelling of Ambient Air Quality (MAAQ) contract. 
## These maps are used to provide policy support for Defra and to fulfil the 
## UK's reporting obligations to Europe.


## Read data
d <- read.csv("data/MAAQ_background/mappm102015g.csv",header=T, 
              na.strings="MISSING", skip=5)
#View(d)

## Ugly basic visualization 
# coords <- unique(d[, 2:3]) # 281802 points!
# plot(y~x, d, pch=".", col=pm102015g)


## Convert dataframe to raster for visualization
r <- rasterFromXYZ(d[,2:4])  #Convert first two columns as lon-lat and third as value                
plot(r)

