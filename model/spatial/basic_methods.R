## References:
## http://rspatial.org/analysis/rst/4-interpolation.html

# Clean environment #################################################################
rm(list=ls())
par(ask=T)

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

## Read data
aurn <- read.table("data/aurn/AURN_data_07_11.txt",header=T)
aurn$date <- as.POSIXct(sprintf("%04d-%02d-%02d", aurn$year,aurn$month,aurn$day),tz="GMT")
names(aurn)[1] <- "site"
aurn$site <- as.factor(aurn$site)
nrow(aurn) # 262944

sites <- unique(aurn[,c(1:3,6)])




##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
### INTERPOLATION METHODS ####
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Interpolation does not take into accoun the time dimension
## Thus, the following experiments will be executed on 
## aggregated data

## Aggregated time series (by date) ############################################
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
## One (random) station
acf(aurn$obs_no2[aurn$site==sample(aurn$site,1)], na.action = na.exclude)
acf(na.omit(aurn$obs_no2[aurn$site==sample(aurn$site,1)]))



## Anual aggregation (by station) ##########################################
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



## NULL model ######################################################
## We are going to interpolate (estimate for unsampled locations) the OBSERVED values. 
## The simplest way would be to take the mean of all observations. We can consider that a 
## “Null-model” that we can compare other approaches to. We’ll use the RMSE as evaluation 
## statistic.

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

## Get the RMSE for the Null-model
null <- RMSE(mean(dsp$obs_no2), dsp$obs_no2)

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
r <- raster(shpAggGB, res=0.05)  
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

## How does the proximity-polygon approach compare to the NULL model?
1 - (mean(rmse) / null)
## Answer: it is not a huge improvement!
## This is caused by the fact that there is not a considerable ammount of stations
## that can capture the large variability of the domain (which is considerably big) 
## Proximty polygons is not a good approach for AIR POLLUTION data, because it assumes that
## adjacent areas have the same levels 



## Nearest neighbour interpolation ###############################################################
## Here we do nearest neighbour interpolation considering multiple neighbours. 
## First make a distance matrix from all points to each other.

# control points
cp <- rasterToPoints(vr)

# distance matrix
d <- pointDistance(cp[, 1:2], dsp, lonlat=T)
nrow(dsp) ## [1] 456
nrow(cp) ## [1] 4087
# not symmetric!
dim(d) ## [1] 4087  456

## Find the nearest 5 neighbours to each point
nn <- 5
ngb <- t(apply(d, 1, function(x) order(x)[1:nn]))
## Check if this all makes sense
plot(shpGB)
points(cp[1, 1:2, drop=FALSE], col='blue', pch='x', cex=2)
points(dsp[ngb[1,], ], col='red', pch=20)
points(cp[nrow(cp), 1:2, drop=FALSE], col='blue', pch='x', cex=2)
points(dsp[ngb[nrow(cp),], ], col='red', pch=20)
## Notes:
## The stations are too far from each other

## Use this to make a map
## First make pairs
pairs <- cbind(rep(1:nrow(ngb), nn), as.vector(ngb))

## Get the values for the pairs and compute their average
values <- dsp$obs_no2[pairs[,2]]
pn <- tapply(values, pairs[,1], mean)
## And assign these to a new RasterLayer
nnr <- r
nnr[!is.na(vr)] <- pn
plot(nnr)

## Cross validate the result
rr <- r
rmse <- rep(NA, 5)
for (k in 1:5) {
  test <- dsp[kf == k, ]
  train <- dsp[kf != k, ]
  d <- pointDistance(cp[, 1:2], train, lonlat=T)
  ngb <- t(apply(d, 1, function(x) order(x)[1:nn]))
  pairs <- cbind(rep(1:nrow(ngb), nn), as.vector(ngb))
  values <- dsp$obs_no2[pairs[,2]]
  pn <- tapply(values, pairs[,1], mean)
  rr[!is.na(vr)] <- pn
  p <- extract(rr, test)
  rmse[k] <- RMSE(test$obs_no2, p)
}
rmse 
mean(rmse) 
1 - (mean(rmse) / null) ## [1] -0.0591731



## Inverse distance weighted ################################################################
## A more commonly used method is “inverse distance weighted” interpolation. 
## We can use the gstat package. 
## First fit the model. ~1 means “intercept only”. 
## In the case of spatial data, that would be only ‘x’ and ‘y’ coordinates are used.
library(gstat)
gs <- gstat(formula=obs_no2~1, locations=dsp)
idw <- interpolate(r, gs)
idwr <- mask(idw, vr)
plot(idwr)

## Cross validate. We can predict to the locations of the test points
rmse <- rep(NA, 5)
for (k in 1:5) {
  test <- dsp[kf == k, ]
  train <- dsp[kf != k, ]
  gs <- gstat(formula=obs_no2~type, locations=train)
  p <- predict(gs, test)
  rmse[k] <- RMSE(test$obs_no2, p$var1.pred)
}
rmse
mean(rmse)
1 - (mean(rmse) / null)
# RMSE
# ~1: 22.53612
# ~type: 18.48396

## Test of some extra parameters
gs2 <- gstat(formula=obs_no2~1, locations=dsp, nmax=1, set=list(idp=1))
idw2 <- interpolate(r, gs2)
idwr2 <- mask(idw2, vr)
plot(idwr2)
## Similar to the geproximity polygons???






## California air pollution data
x <- read.csv("data/airqual.csv")
x$OZDLYAV <- x$OZDLYAV * 1000
coordinates(x) <- ~LONGITUDE + LATITUDE
proj4string(x) <- CRS('+proj=longlat +datum=NAD83')
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +ellps=GRS80")
aq <- spTransform(x, TA)
## Create an template raster to interpolate
cageo <- readRDS('data/counties.rds')
ca <- spTransform(cageo, TA)
r <- raster(ca)
res(r) <- 10  # 10 km if your CRS's units are in km
g <- as(r, 'SpatialGrid')

## Fit a variogram
## Use gstat to create an emperical variogram ‘v’
gs <- gstat(formula=OZDLYAV~1, locations=aq)
v <- variogram(gs, width=20)
head(v)
plot(v)
## Fit a model variogram
fve <- fit.variogram(v, vgm(85, "Exp", 75, 20))
fve
plot(variogramLine(fve, 400), type='l', ylim=c(0,120))
points(v[,2:3], pch=20, col='red')
## Another way
plot(v, fve)


## Ordinary kriging ####################################################
k <- gstat(formula=OZDLYAV~1, locations=aq, model=fve)
## predicted values
kp <- predict(k, g)
spplot(kp)

# variance
ok <- brick(kp)
ok <- mask(ok, ca)
names(ok) <- c('prediction', 'variance')
plot(ok)



## Inverse distance weigthed parameters

f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula=OZDLYAV~1, locations=train, nmax=nmx, set=list(idp=idp))
  p <- predict(m, newdata=test, debug.level=0)$var1.pred
  RMSE(test$OZDLYAV, p)
}
set.seed(124124)
i <- sample(nrow(aq), 0.2 * nrow(aq))
tst <- aq[i,] # 20% for testing
trn <- aq[-i,] # 80% for training
opt <- optim(c(8, .5), f1, test=tst, train=trn)
opt <- optim(c(15, 1), f1, test=tst, train=trn)
opt



# The optimal IDW model
m <- gstat(formula=OZDLYAV~1, locations=aq, nmax=opt$par[1], set=list(idp=opt$par[2]))
idw <- interpolate(r, m)
idw <- mask(idw, ca)
plot(idw)



## A thin plate spline model
library(fields)
m <- Tps(coordinates(aq), aq$OZDLYAV)
tps <- interpolate(r, m)
tps <- mask(tps, idw)
plot(tps)




## Cross-validate ##########################################################
## Cross-validate the three methods (IDW, Ordinary kriging, TPS) and add 
## RMSE weighted ensemble model.

library(dismo)

nfolds <- 5
k <- kfold(aq, nfolds)

ensrmse <- tpsrmse <- krigrmse <- idwrmse <- rep(NA, 5)

for (i in 1:nfolds) {
  test <- aq[k!=i,]
  train <- aq[k==i,]
  
  ## IDW
  m <- gstat(formula=OZDLYAV~1, locations=train, nmax=opt$par[1], set=list(idp=opt$par[2]))
  p1 <- predict(m, newdata=test, debug.level=0)$var1.pred
  idwrmse[i] <-  RMSE(test$OZDLYAV, p1)
  
  ## Kriging
  m <- gstat(formula=OZDLYAV~1, locations=train, model=fve)
  p2 <- predict(m, newdata=test, debug.level=0)$var1.pred
  krigrmse[i] <- RMSE(test$OZDLYAV, p2)
  
  ## TPS 
  m <- Tps(coordinates(train), train$OZDLYAV)
  p3 <- predict(m, coordinates(test))
  tpsrmse[i] <-  RMSE(test$OZDLYAV, p3)

  ## RMSE weighted ensemble model.
  w <- c(idwrmse[i], krigrmse[i], tpsrmse[i])
  weights <- w / sum(w)
  ensemble <- p1 * weights[1] + p2 * weights[2] + p3 * weights[3]
  ensrmse[i] <-  RMSE(test$OZDLYAV, ensemble)
  
}
rmi <- mean(idwrmse)
rmk <- mean(krigrmse)
rmt <- mean(tpsrmse)
rms <- c(rmi, rmt, rmk)
rms
rme <- mean(ensrmse)
rme


## We can use the rmse scores to make a weighted ensemble. Let’s look at the maps
weights <- ( rms / sum(rms) )
s <- stack(idw, ok[[1]], tps)
ensemble <- sum(s * weights)
## Compare maps
s <- stack(idw, ok[[1]], tps, ensemble)
names(s) <- c('IDW', 'OK', 'TPS', 'Ensemble')
plot(s)


## End
par(ask=F)