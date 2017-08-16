## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
source("util/my_helper.R")

## Settings ##########################################################################
printPlots <- T

## Read data #########################################################################
d <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")

## Relevant period
dateA <- min(d$Date)
dateB <- max(d$Date)
numDays <- as.integer((dateB-dateA)+1)

## Read sites
sites <- readRDS("data/epa/epa_daily/2016/california_ozone_sites.RDS")

## Initial exploration ##############################################################################
hist(d$Ozone)
hist(d$Temperature)
hist(d$RH)



## Heatmap plus daily mean TS ##################################################################
levelplot(Ozone~Date*Station.Code,d,
                cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                scales=list(y=list(cex=.6)))

heatmapPlusTs(d, "Ozone", "Ozone (ppb)")
heatmapPlusTs(d, "Ozone", "Ozone (ppb)", "img/eda/heatmap_ca_ozone.jpeg")
heatmapPlusTs(d, "Temperature", "Temperature (K)")
heatmapPlusTs(d, "Temperature", "Temperature (K)", "img/eda/heatmap_ca_temperature_2.jpeg")
heatmapPlusTs(d, "RH", "RH (%)")
heatmapPlusTs(d, "RH", "RH (K)", "img/eda/heatmap_ca_rh_2.jpeg")

