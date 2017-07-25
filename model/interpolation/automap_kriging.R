## Clean environment
rm(list=ls())

## Load libraries
library(automap)
library(spTimer)

## packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})

#### Very basic example ####
# loadMeuse()
# # Ordinary kriging
# kriging_result = autoKrige(zinc~1, meuse, meuse.grid)
# plot(kriging_result)
# # Universal kriging
# kriging_result = autoKrige(zinc~soil+ffreq+dist, meuse, meuse.grid)
# plot(kriging_result)


#### Automap - kriging ####
# NYdata ####
data("NYdata")
data("NYgrid")

NYdataAgg <- aggregate(cbind(o8hrmax,cMAXTMP,WDSP,RH)~s.index+Longitude+Latitude, 
                       FUN=mean, data=NYdata)
coordinates(NYdataAgg) <- ~ Longitude + Latitude
head(NYdataAgg)

NYgridAgg <- aggregate(cbind(cMAXTMP,WDSP,RH)~s.index+Longitude+Latitude,
                       FUN=mean, data=NYgrid)
coordinates(NYgridAgg) <- ~ Longitude + Latitude
head(NYgridAgg)


NYgridAgg %>% as.data.frame %>% 
  ggplot(aes(Longitude, Latitude)) + geom_point(aes(size=cMAXTMP), color="blue", alpha=3/4) + 
  ggtitle("Temperature") + coord_equal() + theme_bw()

# Kriging NYdata
# Ordinary kriging
kriging1 = autoKrige(o8hrmax~1, NYdataAgg, NYgridAgg)
plot(kriging1)
# Universal kriging
kriging2 = autoKrige(o8hrmax~cMAXTMP+WDSP+RH, NYdataAgg, NYgridAgg)
plot(kriging2)
# There really is spatial correlation?
# The range of o8hrmax is
range(NYdataAgg$o8hrmax)
max(NYdataAgg$o8hrmax) - min((NYdataAgg$o8hrmax))
# The closest samples show a very high variation (point label 9), which
# may suggest that there is not spatial correlation
# Check it!

