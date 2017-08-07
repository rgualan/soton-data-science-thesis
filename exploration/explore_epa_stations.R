## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)

## Settings ##########################################################################
printPlots <- T

## Read data #########################################################################
sites_md <- read.csv("data/epa/sites/aqs_sites.csv", header=T)
sites_md <- sites_md[sites_md$State.Name=="California",]
sites_md$Station.Code <- as.factor(sprintf("%03d-%04d",sites_md$County.Code,sites_md$Site.Num))
#names(sites_md)
sites_md <- sites_md[,c("Station.Code","State.Code", "County.Code", "Site.Number", 
                        "Latitude", "Longitude", "Datum", "Elevation", "Location.Setting",
                        "State.Name")] # Keep relevant fields only
#View(sites_md)

## Obtain locations in kilometers using UTM 
sites_md.sp <- sites_md
coordinates(sites_md.sp) <- ~Longitude+Latitude
proj4string(sites_md.sp) <- "+proj=longlat +datum=WGS84"
utm11 = CRS("+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +units=km")
sites_md.sp <- spTransform(sites_md.sp,utm11) 

sites_md_2 <- as.data.frame(sites_md.sp)
names(sites_md_2)[9:10] <- c("x", "y")
sites_md_3 <- merge(sites_md,sites_md_2[,c("Station.Code","x","y")])
#View(sites_md_3)

## Save dataset #####################################################################
sites_md_3 <- sites_md_3[order(sites_md_3$Station.Code),]
saveRDS(sites_md_3, "data/epa/epa_daily/2016/aqs_sites.RDS")
