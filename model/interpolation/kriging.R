## Clean environment
rm(list=ls())

## Load libraries
library(sp)
library(gstat)
library(spTimer)

## Packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})

## Load datasets
data(NYdata)
glimpse(NYdata)

NYdataAgg <- aggregate(o8hrmax~Longitude+Latitude, FUN=mean, data = NYdata)
coordinates(NYdataAgg) <- ~ Longitude + Latitude

NYdataAgg %>% as.data.frame %>% 
  ggplot(aes(Longitude, Latitude)) + geom_point(aes(size=o8hrmax), color="blue", alpha=3/4) + 
  ggtitle("Temperature") + coord_equal() + theme_bw()

#### Fitting a variogram ####
NYdataAgg %>% as.data.frame %>% glimpse

lzn.vgm <- variogram(log(o8hrmax)~1, NYdataAgg) # calculates sample variogram values 
lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph")) # fit model
plot(lzn.vgm, lzn.fit) # plot the sample values, along with the fit model



#### Performing Kriging ####
# load spatial domain to interpolate over
data("meuse")
data("meuse.grid")

# to compare, recall the bubble plot above; those points were what there were values for. this is much more sparse
plot1 <- meuse %>% as.data.frame %>%
  ggplot(aes(x, y)) + geom_point(size=1) + coord_equal() + 
  ggtitle("Points with measurements")

# this is clearly gridded over the region of interest
plot2 <- meuse.grid %>% as.data.frame %>%
  ggplot(aes(x, y)) + geom_point(size=1) + coord_equal() + 
  ggtitle("Points at which to estimate")

library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)

#### Computation ####



