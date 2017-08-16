## Clean environment
rm(list=ls())

## Load libraries
library(sp)

## Load dataset
load("data/ny_ozone/NYdata.Rdata")
load("data/ny_ozone/NYcheap.Rdata")
load("data/ny_ozone/NYcheapPlusNoise.Rdata")


## Visualize error histogram
NYerror <- NYcheap[, c("s.index", "date")]
NYerror$error <- NYcheapPlusNoise$o8hrmax - NYcheap$o8hrmax
hist(NYerror$error[NYerror$s.index==100])



## Calculate BIAS
NYerrorAgg <- aggregate(cbind(NYerror$error)~NYerror$s.index
                        +NYerror@coords[,1]+NYerror@coords[,2],  
                        FUN=mean)
names(NYerrorAgg) <- c("s.index","Longitude","Latitude","mean.error")
head(NYerrorAgg)
coordinates(NYerrorAgg) <- ~Longitude+Latitude
spplot(NYerrorAgg, "mean.error")


# Calculate SD
NYerrorAgg <- aggregate(cbind(NYerror$error)~NYerror$s.index
                        +NYerror@coords[,1]+NYerror@coords[,2],  
                        FUN=sd)
names(NYerrorAgg) <- c("s.index","Longitude","Latitude","sd.error")
head(NYerrorAgg)
coordinates(NYerrorAgg) <- ~Longitude+Latitude
spplot(NYerrorAgg, "sd.error")



