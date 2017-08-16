## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
source("util/my_helper.R")

## Settings ##########################################################################
paper <- setupPaper()

## Read data #########################################################################
ozone <- getSites(readRDS("data/epa/epa_daily/2016/california_ozone_2.RDS"))
pm10 <-  getSites(readRDS("data/epa/epa_daily/2016/california_pm10.RDS"))
temperature <- getSites(readRDS("data/epa/epa_daily/2016/california_temperature.RDS"))
wind <- getSites(readRDS("data/epa/epa_daily/2016/california_wind.RDS"))
rh <- getSites(readRDS("data/epa/epa_daily/2016/california_rh.RDS"))

## Initial exploration ##############################################################################
plot(getCAmap())
points(Latitude~Longitude,ozone,cex=0.5,pch=1,col=1)
points(Latitude~Longitude,pm10,cex=0.5,pch=2,col=2)
points(Latitude~Longitude,temperature,cex=0.5,pch=3,col=3)
points(Latitude~Longitude,wind,cex=0.5,pch=4,col=4)
points(Latitude~Longitude,rh,cex=0.5,pch=5,col=5)

stations <- rbind(cbind(ozone,Parameter="Ozone"),
                  cbind(pm10,Parameter="PM10"),
                  cbind(temperature,Parameter="Temperature"),
                  cbind(wind,Parameter="WS"),
                  cbind(rh,Parameter="RH"))
printPlot(paper, "img/eda/ca_parameters.jpeg", 6, 6, FUN=function(){
  p<-ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = stations, aes(x = Longitude, y = Latitude, shape=Parameter, color=Parameter),
               alpha = .75, size=1.5) +
    scale_shape_discrete(solid=F) +
    labs(x = "Longitude", y = "Latitude") +
    coord_quickmap() +
    theme(legend.justification = c("right", "top"), legend.position = c(.95, .95),
          legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
  print(p)
})
