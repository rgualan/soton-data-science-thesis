# Clean environment #################################################################
rm(list=ls())

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


## Plot all stations by type of station
sites <- unique(aurn[,c(1:3,6)])
map('worldHires', 'UK', xlim=c(-8.5,2), ylim=c(50,58.6), mar=rep(0,4))	
points(lat~lon,sites,pch=2,col=type)
legend("topright",legend=unique(sites$type), pch=2, col=unique(sites$type))
readline("Continue?")


## Highlight stations inside Greater London
sitesGL <- sites[sites$lon >= -0.5 & sites$lon <= 0.4
                 & sites$lat >= 51 & sites$lat <= 51.7,]
nrow(sitesGL)
map('worldHires', c('UK'), xlim=c(-8.5,2), ylim=c(50,58.6), mar=rep(0,4))	
points(lat~lon,sites,pch=".",col="blue",cex=2)
points(lat~lon,sitesGL,pch="*",col="green",cex=1)
readline("Continue?")

## Only stations inside Greater London
par(.pardefault)
plot(lat~lon,sitesGL,pch=2,col="green",cex=1)
readline("Continue?")

## Shapefile plus stations inside GL
## Ref: https://medium.com/towards-data-science/plotting-a-map-of-london-crime-data-using-r-8dcefef1c397
shpGL <- readOGR("/home/ronald/projects/aq-tngapms/data/maps/statistical-gis-boundaries-london/ESRI", 
                 layer = "London_Borough_Excluding_MHW")
proj4string(shpGL)
shpGL.wgs84 <- spTransform(shpGL, CRS("+init=epsg:4326"))

ggplot(shpGL.wgs84) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  labs(x = "Longitude", y = "Latitude", title = "Map of Greater London with the borough boundaries") +
  geom_point(data = sitesGL, aes(x = lon, y = lat, colour = type)) +
  coord_quickmap()
#+scale_colour_manual(values = rainbow(14))
readline("Continue?")


## Time series ###################################################################
# aqum <- aurn[,c("site","date","aqum_no2")]
# aqumSub <- aqum[aqum$site %in% sites$site[sample(nrow(sites),5)],]
# summaryPlot(aqumSub, main="AQUM_NO2")
# readline("Continue?")
# 
# aqum <- aurn[,c("site","date","obs_no2")]
# aqumSub <- aqum[aqum$site %in% sites$site[sample(nrow(sites),5)],]
# summaryPlot(aqumSub, main="OBS_NO2")


## Plot a data concentration as matrix (sites x date) ################
st.on <- sort(unique(aurn$site))
NS<-50
for(i in 1:ceiling(length(st.on)/50)){
  print(i)
  print(levelplot(obs_no2~date*site, 
                  aurn[aurn$site
                     %in% st.on[(NS*(i-1)+1):(NS*i)],],
                  cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                  scales=list(y=list(cex=.7))))
  readline("Continue?")
}


## Aggregated time series ############################################
dailyAgg <- aggregate(obs_no2~date,aurn,mean) # daily mean
names(dailyAgg)[2] <- "mean"
ua <- aggregate(obs_no2~date,aurn,quantile,0.95) # upper
la <- aggregate(obs_no2~date,aurn,quantile,0.05) # lower
dailyAgg$upper <- ua$obs_no2
dailyAgg$lower <- la$obs_no2
#View(dailyAgg)
## Plot
ggplot() + 
  geom_smooth(aes(x=date, y=mean, ymax=upper, ymin=lower), 
              data=dailyAgg, stat='identity') +
  ggtitle("Daily average and bounds")


## Individual time series ############################################
dailyAgg <- aggregate(obs_no2~date,aurn,mean) # daily mean
names(dailyAgg)[2] <- "mean"
ua <- aggregate(obs_no2~date,aurn,quantile,0.95) # upper
la <- aggregate(obs_no2~date,aurn,quantile,0.05) # lower
dailyAgg$upper <- ua$obs_no2
dailyAgg$lower <- la$obs_no2
#View(dailyAgg)
## Plot
ggplot() + 
  geom_smooth(aes(x=date, y=mean, ymax=upper, ymin=lower), 
              data=dailyAgg, stat='identity') +
  ggtitle("Daily average and bounds")


## Spatial autocorrelation ##########################################
## Spatial autocorrelation in a variable can be exogenous (it is caused 
## by another spatially autocorrelated variable, e.g. rainfall) or endogenous 
## (it is caused by the process at play, e.g. the spread of a disease).

## Daily mean
acf(dailyAgg$mean)
## One station
acf(aurn$obs_no2[aurn$site==4], na.action = na.exclude)
acf(aurn$obs_no2[aurn$site==sample(aurn$site,1)], na.action = na.exclude)
readline("Continue?")


## Anual aggregation ##########################################
annual <- aggregate(obs_no2~site, aurn, mean)
#annual$obs_no2 <- annual$obs_no2/1000
annual <- merge(annual, sites)
#View(annual)

ggplot(annual[order(annual$obs_no2),]) +
  labs(x = "Stations", y = "Annual sum (x1000)", title = "Annual observations") +
  geom_point(aes(x = 1:nrow(annual), y = obs_no2, colour = type)) + 
  theme(legend.position="bottom")


## Annual interpolation map ###################################
## Create SpatialPointsDataFrame
dsp <- SpatialPoints(annual[,3:4], proj4string=CRS("+proj=longlat +datum=WGS84"))
dsp <- SpatialPointsDataFrame(dsp, annual)

## GB boundline shapefile 
#shpGB1 <- readOGR("data/maps/bdline_essh_gb/GB", layer = "european_region_region1")
#shpGB <- readOGR("data/maps/bdline_essh_gb/GB", layer = "european_region_region_boundary")
#shpGB <- readOGR("data/maps/strtgi_essh_gb/data", layer = "coastline")
shpGB <- readOGR("data/maps/Countries_De_2016_Clipped_Boundaries_in_GB", 
                 layer = "Countries_December_2016_Ultra_Generalised_Clipped_Boundaries_in_Great_Britain")
proj4string(shpGB)
shpGB <- spTransform(shpGB, CRS("+proj=longlat +datum=WGS84"))
ggplot(shpGB) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  labs(x = "Longitude", y = "Latitude", title = "Map of Great Britain") +
  geom_point(data = dsp@data, aes(x = lon, y = lat, colour = type)) +
  coord_quickmap()


# define groups for mapping
range(dsp$obs_no2)
cuts <- c(0,25,50,75,200)
# set up a palette of interpolated colors
#blues <- colorRampPalette(c('yellow', 'orange', 'blue', 'dark blue'))
blues <- colorRampPalette(c('lightblue', 'orange', 'red4'))
pols <- list("sp.polygons", shpGB, fill = "white")
#spplot(dsp, 'obs_no2', cuts=cuts, col.regions=blues(5), sp.layout=pols, pch=20, cex=1)
spplot(dsp, 'obs_no2', cuts=cuts, col.regions=blues(5), sp.layout=pols, 
       pch=20, cex=1)




## Proximity polygons ######################################################
## Proximity polygons can be used to interpolate categorical variables. 
## Another term for this is “nearest neighbour” interpolation.

library(dismo)
v <- voronoi(dsp)
## Loading required namespace: deldir
plot(v)

## Aggregate countries into a single boundary for clipping the Voronoi polygons
shpAggGB <- aggregate(shpGB)
plot(shpAggGB)

vclipped <- intersect(v, shpAggGB)
spplot(vclipped, 'obs_no2', col.regions=rev(get_col_regions()))

## Much better. These are polygons. We can ‘rasterize’ the results like this.
r <- raster(shpAggGB, res=0.01)  
#vr <- rasterize(vclipped, r, 'obs_no2')
vr <- rasterize(vclipped, r, 'obs_no2')
plot(vr)


## Now evaluate with 5-fold cross validation
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

set.seed(133133)
kf <- kfold(nrow(dsp))

rmse <- rep(NA, 5)
for (k in 1:5) {
  #k <- 1
  test <- dsp[kf == k, ]
  train <- dsp[kf != k, ]
  v <- voronoi(train)
  p <- extract(v, test)
  rmse[k] <- RMSE(test$obs_no2, p$obs_no2)
}
rmse
mean(rmse)
## This results are relatively good because they were obtained with aggregated data



## [1] 196.7708
1 - (mean(rmse) / null)
## [1] 0.5479875
