## Clean environment #################################################################
rm(list=ls())


## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(mice)
library(Amelia)
library(VIM)
library(gstat)
library(raster)

## Load functions
source("util/my_helper.R")



## Read data #########################################################################
ozone <- readRDS("data/epa/epa_daily/2016/california_ozone.RDS")
temp <- readRDS("data/epa/epa_daily/2016/california_temperature.RDS")
wind <- readRDS("data/epa/epa_daily/2016/california_wind.RDS")
rh <- readRDS("data/epa/epa_daily/2016/california_rh.RDS")
ozone.sites <- readRDS("data/epa/epa_daily/2016/california_ozone_sites.RDS")
temp.sites <- readRDS("data/epa/epa_daily/2016/california_temperature_sites.RDS")
wind.sites <- readRDS("data/epa/epa_daily/2016/california_wind_sites.RDS")
rh.sites <- readRDS("data/epa/epa_daily/2016/california_rh_sites.RDS")

## Read sites
sites <- readRDS("data/epa/sites/aqs_sites.RDS")
sites <- sites[sites$State.Name=="California",]


## Intersection between stations
sum(temp.sites$Station.Code %in% ozone.sites$Station.Code)
sum(wind.sites$Station.Code %in% ozone.sites$Station.Code)
sum(rh.sites$Station.Code %in% ozone.sites$Station.Code) # Only 20!


## Merge with the other variables ########################################################
ozone <- fillMissingDates(ozone)
## Implicit filling of missing dates in the other data frames
cali <- merge(ozone, temp, all.x=T)
cali <- merge(cali, wind, all.x=T)
cali <- merge(cali, rh, all.x=T)
#View(cali)
summary(cali)
days <- seq(from=min(cali$Date), to=max(cali$Date),by='days' )



## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Convert to Spatial data ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cali.sp <- convertDataToSp(cali)

# cali.var <- variogram(cali~x+y, 
#                        na.omit(cali.sp[cali.sp$Date==sample(cali.sp$Date,1),"Ozone"])) # Random date
# plot(cali.var)



## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Analyze Missing Values ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Check missing rows
#md.pattern(cali)
aggr(cali[,-(1:2)], gap=3, cex.axis=0.8)

levelplot(Ozone~Date*Station.Code,cali,
                cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                scales=list(y=list(cex=.7)))
levelplot(Temperature~Date*Station.Code,cali,
          cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
          scales=list(y=list(cex=.7)))
levelplot(Wind.speed~Date*Station.Code,cali,
          cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
          scales=list(y=list(cex=.7)))
levelplot(RH~Date*Station.Code,cali,
          cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
          scales=list(y=list(cex=.7)))




## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Impute Temperature ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tempNasByStation <- aggregate(Temperature~Station.Code,cali,function(x){sum(is.na(x))}, na.action=na.pass)
tempNasByStation <- tempNasByStation[tempNasByStation$Temperature>0,] 
#View(tempNasByStation)
  

## Visualize Missing and Available Stations
naStations <- tempNasByStation[tempNasByStation$Temperature>100,1]
temp.sites.na <- 
  merge(data.frame(Station.Code=naStations), 
        sites[,c("Station.Code","Latitude","Longitude","Datum","Elevation",
                 "Location.Setting")], all.x=T)
temp.sites.all <- rbind(cbind(temp.sites,Type="Available"),
                        cbind(temp.sites.na,Type="Missing"))

mapUSA <- readRDS("data/maps/usa/USA_adm1.rds")
mapCA <- mapUSA[mapUSA$NAME_1=="California",]
ggplot(mapCA) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_point(data = temp.sites.all, 
             aes(x = Longitude, y = Latitude, fill=Type),
             alpha = .75, shape=21, size=2) +
  geom_text(data = ozone.sites[ozone.sites$Station.Code=="089-0009",], 
            aes(x = Longitude, y = Latitude, label=Station.Code)) +
  labs(x = "Longitude", y = "Latitude", fill = "Type") +
  coord_quickmap() +
  theme(legend.justification = c("right", "top"), legend.position = c(.95, .95), 
        legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))

print(levelplot(Temperature~Date*Station.Code,temp,
                cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                scales=list(y=list(cex=.7))))

## Observations:
## Since the amount of missing data in the "available" dataset is small, 
## a simple imputation method (Inverse distance weigthed) is used


## IDS (Inverse Distance Weigthed) Interpolation
# Compute the sample variogram; note that the f.1 trend model is one of the
# parameters passed to variogram(). This tells the function to create the 
# variogram on the de-trended data.
temp.sp <- convertDataToSp(temp)
# temp.var <- variogram(Temperature~x+y,temp.sp)
# plot(temp.var)


## Assessing IDW interpolation on different locations 
## Test 1
ss <- c(
  "083-1021", ## Optimally located between several close stations 
  "089-3003"  ## Most isolated station with 356 values for testing
)
s<- ss[2]
## Results (OS-R2, OS-BIAS)
## 1: 0.9386069, 0.789375
## 2: 0.9219895, 17.93853   ## HUGE BIAS!!!
## NOTES:
## IDW Interpolation can be used for data imputation 
## only when there is close stations
## Otherwise, the imputation is substancially biased

if(F){
  
  ggplot(mapCA) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = temp.sites.all, 
               aes(x = Longitude, y = Latitude, fill=Type),
               alpha = .75, shape=21, size=2) +
    geom_point(data = temp.sites.all[temp.sites.all$Station.Code==s,], 
               aes(x = Longitude, y = Latitude), alpha = .75, shape=21, size=3, fill="red") +
    labs(x = "Longitude", y = "Latitude", fill = "Type") +
    coord_quickmap() +
    theme(legend.justification = c("right", "top"), legend.position = c(.95, .95), 
          legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
  
  dim(temp.sp[temp.sp$Station.Code!=s,])
  dim(temp.sp[temp.sp$Station.Code==s,])
  
  ## IDW
  test <- temp[temp$Station.Code==s,]
  test <- fillMissingDates(test)
  test <- convertDataToSp(test)
  #head(coordinates(test))
  test$Temperature2 <- NA
  temp.sp$Temperature[temp.sp$Station.Code==s 
                      & temp.sp$Date>=convertStringToPOSIXct("2016-05-01")
                      & temp.sp$Date<=convertStringToPOSIXct("2016-05-31")] <- NA # Inject NAs
  #View(temp.sp[temp.sp$Station.Code==s,])
  
  for ( i in 1:length(days) ){
    print(days[i])
    a <- idw(Temperature ~ 1, 
             temp.sp[temp.sp$Date==days[i] & !is.na(temp.sp$Temperature),],  # IN-SAMPLE
             # temp.sp[temp.sp$Station.Code!=s & temp.sp$Date==days[i],],  # OUT_OF_SAMPLE
             newdata=test[test$Date==days[i],], 
             idp=2.0, debug.level=0)
    test$Temperature2[test$Date==days[i]] <- a$var1.pred
  }
  plot(Temperature~Date,test, type="l")
  lines(Temperature2~Date,test, col=2, lty="dashed")
  cor(test$Temperature,test$Temperature2,use="pairwise.complete.obs")^2
  mean(test$Temperature2-test$Temperature,na.rm=T)

  ## Kriging
  for ( i in 1:length(days) ){
    print(days[i])
    sampleVar <- variogram(Temperature~1,temp.sp[temp.sp$Date==days[i],], cutoff=500)
    #plot(sampleVar)
    modelVar <- fit.variogram(sampleVar, vgm("Exp"))
    #plot(sampleVar, modelVar)
    k <- gstat(formula=Temperature~1, loc=temp.sp[temp.sp$Date==days[i],], model=modelVar)
    kp <- predict(k, 
                  newStation[newStation$Station.Code=="083-1021" & newStation$Date==days[i],])
    #spplot(kp)
    newStation$Temperature2[newStation$Date==days[i]] <- kp$var1.pred
  }  
  lines(Temperature2~Date,newStation, lty="dashed", col=3)
}




test <- temp[temp$Station.Code==s,]
test <- fillMissingDates(test)
test <- convertDataToSp(test)
#head(coordinates(test))
test$Temperature2 <- NA
if(F){
  ## IDW (Inverse Distance Weigthed) Interpolation
  for( i in 1:length(days) ){
    #print(days[i])
    a <- idw(Temperature ~ 1, 
             temp.sp[temp.sp$Station.Code!=s & temp.sp$Date==days[i],], 
             newdata=test[test$Date==days[i],], 
             idp=2.0, debug.level=0)
    test$Temperature2[test$Station.Code==s & test$Date==days[i]] <- a$var1.pred
  }
  plot(Temperature~Date,test,type="l")      
  lines(Temperature2~Date,test,col=2,lty="dashed")      
  ## BIAS:
  mean(test$Temperature2-test$Temperature, na.rm=T)
}




### Imputation stage 1: Impute stations with less than 100 NAs by station
### using IDW interpolation
tempNasByStationS1 <- tempNasByStation[tempNasByStation$Temperature>0 & tempNasByStation$Temperature<100,] 

cali.sp$TemperatureFlag <- 0
cali.sp$TemperatureFlag[cali.sp$Station.Code %in% tempNasByStationS1$Station.Code
                         & is.na(cali.sp$Temperature)] <- 1 # Imputation method
for(station in tempNasByStationS1$Station.Code){
  print(station)
  #days <- cali.sp$Date[cali.sp$Station.Code==station & is.na(cali.sp$Temperature)]
  days <- cali.sp$Date[cali.sp$Station.Code==station 
                        & !is.na(cali.sp$TemperatureFlag) & cali.sp$TemperatureFlag==1]
  #print(days)
  
  for ( i in 1:length(days) ){
    #print(days[i])
    a <- idw(Temperature ~ 1, 
             temp.sp[temp.sp$Date==days[i],], 
             newdata=cali.sp[cali.sp$Station.Code==station & cali.sp$Date==days[i],], 
             idp=2.0, debug.level=0)
    cali.sp$Temperature[cali.sp$Station.Code==station & cali.sp$Date==days[i]] <- a$var1.pred
  }
}
#View(cali.sp[,c("Station.Code","Date","Ozone","Temperature","TemperatureIDW")])
## Test
plot(Temperature~Date,cali.sp[cali.sp$Station.Code==cali.sp$Station.Code[1],], col=0,
     ylim=range(cali.sp$Temperature,na.rm=T))
for(station in tempNasByStationS1$Station.Code[tempNasByStationS1$Temperature>0
                                             & tempNasByStationS1$Temperature<100]){
  
  #plot(Temperature~Date,cali.sp[cali.sp$Station.Code==station,], type="l")
  lines(Temperature~Date,cali.sp[cali.sp$Station.Code==station,], type="l")
  points(Temperature~Date,cali.sp[cali.sp$Station.Code==station
                                   & !is.na(cali.sp$TemperatureIM==1),], col=2, pch=18)
  #readline("Continue?")
}

ggplot(cali.sp@data) + 
  geom_line(aes(Date,Temperature, group=Station.Code, colour=(TemperatureFlag==0)), alpha=0.5) + 
  theme(legend.position="none")




## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Interpolate Unavailable Temperature stations ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Annual IDW interpolation (test)
if(F){
  annualMeanTemperature <- aggregate(Temperature~Station.Code, temp, mean)
  #View(annualMeanTemperature)
  annualMeanTemperature <- merge(annualMeanTemperature, sites[,c("Station.Code","x","y")], all.x=T)
  coordinates(annualMeanTemperature)<-~x+y
  proj4string(annualMeanTemperature) <- utm11
  
  mapCa.utm <- spTransform(mapCA, utm11)
  mapCa.utm.raster <- raster(mapCa.utm, res=10); 
  
  coordinates(temp.sites.all) <-~Longitude+Latitude
  proj4string(temp.sites.all) <- "+proj=longlat"
  temp.sites.all <- spTransform(temp.sites.all, utm11)
  
  gs <- gstat(formula=Temperature~1, locations=annualMeanTemperature)
  out <- interpolate(mapCa.utm.raster, gs)
  out <- mask(out, mapCa.utm)
  plot(out)
  points(temp.sites.all, pch=3, cex=0.5, col=temp.sites.all$Type)
}  
  

### Imputation stage 2: Interpolate unavailable stations 
### using IDW interpolation
tempNasByStationS2 <- tempNasByStation[tempNasByStation$Temperature>100,] 
#cali.sp.new <- cali.sp[cali.sp$Station.Code %in% tempNasByStationS2$Station.Code,]
days <- sort(unique(cali.sp$Date))
for ( i in 1:length(days) ){
  #print(days[i])
  a <- idw(Temperature ~ 1,
           temp.sp[temp.sp$Date==days[i],],
           newdata=cali.sp[cali.sp$Station.Code %in% tempNasByStationS2$Station.Code & cali.sp$Date==days[i],],
           idp=2.0, debug.level=0)
  cali.sp$Temperature[cali.sp$Station.Code %in% tempNasByStationS2$Station.Code & cali.sp$Date==days[i]] <- a$var1.pred
}
## Test
ggplot(cali.sp[cali.sp$Station.Code %in% tempNasByStationS2$Station.Code,]@data) + 
  geom_line(aes(Date,Temperature, colour=Station.Code), alpha=0.5) + 
  theme(legend.position="none")
levelplot(Temperature~Date*Station.Code,cali.sp@data,
          cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
          scales=list(y=list(cex=.7)))
ggplot(cali.sp@data) + 
  geom_line(aes(Date,Temperature, colour=Station.Code), alpha=0.5) + 
  theme(legend.position="none")







## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Impute Wind speed ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
varNasByStation <- aggregate(Wind.speed~Station.Code,cali,function(x){sum(is.na(x))}, na.action=na.pass)
varNasByStation <- varNasByStation[varNasByStation$Wind.speed>0,] 
#View(varNasByStation)

wind.sp <- wind
wind.sp <- merge(wind.sp, sites[,c("Station.Code","x","y")], all.x=T)
coordinates(wind.sp) <- ~x+y
proj4string(wind.sp) <- utm11


## Visualize Missing and Available Stations
naStations <- varNasByStation[varNasByStation$Wind.speed>100,1]
var.sites.na <- 
  merge(data.frame(Station.Code=naStations), 
        sites[,c("Station.Code","Latitude","Longitude","Datum","Elevation",
                 "Location.Setting")], all.x=T)
var.sites.all <- rbind(cbind(wind.sites,Type="Available"),
                       cbind(var.sites.na,Type="Missing"))

ggplot(mapCA) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_point(data = var.sites.all, 
             aes(x = Longitude, y = Latitude, fill=Type),
             alpha = .75, shape=21, size=2) +
  labs(x = "Longitude", y = "Latitude", fill = "Type") +
  coord_quickmap() +
  theme(legend.justification = c("right", "top"), legend.position = c(.95, .95), 
        legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))


# Test: 001-2005
print(levelplot(Wind.speed~Date*Station.Code,wind,
                cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                scales=list(y=list(cex=.7))))
print(levelplot(Wind.speed~Date*Station.Code,cali,
                cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                scales=list(y=list(cex=.7))))


## Test 
theStation <- "089-3003"
if(F){
  ggplot(mapCA) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = var.sites.all, 
               aes(x = Longitude, y = Latitude, fill=Type),
               alpha = .75, shape=21, size=2) +
    geom_point(data = var.sites.all[var.sites.all$Station.Code==theStation,], 
               aes(x = Longitude, y = Latitude),
               alpha = .75, shape=21, size=4, fill="red") +
    labs(x = "Longitude", y = "Latitude", fill = "Type") +
    coord_quickmap() +
    theme(legend.justification = c("right", "top"), legend.position = c(.95, .95), 
          legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
  
  ## IDW
  days <- seq(from=min(cali.sp$Date), to=max(cali.sp$Date),by='days' )
  cali.sp$Wind.speed.log2 <- NA
  cali.sp$Wind.speed.log <- log(cali.sp$Wind.speed)
  wind.sp$Wind.speed.log <- log(wind.sp$Wind.speed)
  for ( i in 1:length(days) ){
    print(days[i])
    a <- idw(Wind.speed.log ~ 1,
             wind.sp[wind.sp$Date==days[i] & wind.sp$Station.Code!=theStation 
                     & !is.na(wind.sp$Wind.speed),],
             newdata=cali.sp[cali.sp$Station.Code==theStation
                              & cali.sp$Date==days[i],],
             idp=2.0, debug.level=0)
    cali.sp$Wind.speed.log2[cali.sp$Station.Code==theStation 
                         & cali.sp$Date==days[i]] <- a$var1.pred
  }
  plot(Wind.speed.log~Date,cali.sp[cali.sp$Station.Code==theStation,], type="l")
  lines(Wind.speed.log2~Date,cali.sp[cali.sp$Station.Code==theStation,], col=2, lty="dashed")
  ## Notes: Fail!!!
  

  ## Semivariogram
  for ( i in 1:length(days) ){
    sampleVar <- variogram(Wind.speed.log~1,cali.sp[cali.sp$Date==days[i] & !is.na(cali.sp$Wind.speed),], cutoff=500)
    if(i==1) plot(gamma~dist,sampleVar,type="l") 
    lines(gamma~dist,sampleVar,col=i)
  }
  sampleVar <- variogram(Wind.speed.log~1,cali.sp[!is.na(cali.sp$Wind.speed),], cutoff=500)
  plot(sampleVar) 
  modelVar <- fit.variogram(sampleVar, vgm(150, "Exp", 400, nugget = 1))
  ## NOTE:   singular model in variogram fit

  sampleVar <- variogram(Wind.speed.log~Temperature,cali.sp[cali.sp$Date==days[1] & !is.na(cali.sp$Wind.speed.log),], cutoff=500)
  plot(sampleVar) 
  modelVar <- fit.variogram(sampleVar, vgm(0.5, "Exp", 50, nugget = 0))
  plot(sampleVar, modelVar)
  
  ## Kriging
  cali.sp$Wind.speed.log3  <- NA
  for ( i in 1:length(days) ){
    print(days[i])
    sampleVar <- variogram(Wind.speed.log~Temperature,cali.sp[cali.sp$Date==days[i] & !is.na(cali.sp$Wind.speed.log),], cutoff=500)
    #plot(sampleVar);
    modelVar <- fit.variogram(sampleVar, vgm(0.5, "Exp", 50, nugget = 0))
    #plot(sampleVar, modelVar)
    k <- gstat(formula=Wind.speed~Temperature, loc=wind.sp[cali.sp$Date==days[i] & !is.na(cali.sp$Wind.speed.log),], model=modelVar)
    kp <- predict(k, cali.sp[cali.sp$Station.Code==theStation & cali.sp$Date==days[i],])
    #spplot(kp)
    cali.sp$Wind.speed.log3[cali.sp$Date==days[i]] <- kp$var1.pred
  }  
  plot(Wind.speed.log~Date,cali.sp[cali.sp$Station.Code==theStation,], type="l")
  lines(Wind.speed.log3~Date,cali.sp[cali.sp$Station.Code==theStation,], lty="dashed", col=3)
  ## NOTE: FAIL!
  

}



### Imputation stage 1: Impute stations with less thatn 100 NAs by station
### using IDW interpolation
tempNasByStationS1 <- tempNasByStation[tempNasByStation$Temperature>0 & tempNasByStation$Temperature<100,] 

cali.sp$TemperatureFlag <- 0
cali.sp$TemperatureFlag[cali.sp$Station.Code %in% tempNasByStationS1$Station.Code
                         & is.na(cali.sp$Temperature)] <- 1 # Imputation method
for(station in tempNasByStationS1$Station.Code){
  print(station)
  #days <- cali.sp$Date[cali.sp$Station.Code==station & is.na(cali.sp$Temperature)]
  days <- cali.sp$Date[cali.sp$Station.Code==station 
                        & !is.na(cali.sp$TemperatureFlag) & cali.sp$TemperatureFlag==1]
  #print(days)
  
  for ( i in 1:length(days) ){
    #print(days[i])
    a <- idw(Temperature ~ 1, 
             temp.sp[temp.sp$Date==days[i],], 
             newdata=cali.sp[cali.sp$Station.Code==station & cali.sp$Date==days[i],], 
             idp=2.0, debug.level=0)
    cali.sp$Temperature[cali.sp$Station.Code==station & cali.sp$Date==days[i]] <- a$var1.pred
  }
}
#View(cali.sp[,c("Station.Code","Date","Ozone","Temperature","TemperatureIDW")])
## Test
plot(Temperature~Date,cali.sp[cali.sp$Station.Code==cali.sp$Station.Code[1],], col=0,
     ylim=range(cali.sp$Temperature,na.rm=T))
for(station in tempNasByStationS1$Station.Code[tempNasByStationS1$Temperature>0
                                               & tempNasByStationS1$Temperature<100]){
  
  #plot(Temperature~Date,cali.sp[cali.sp$Station.Code==station,], type="l")
  lines(Temperature~Date,cali.sp[cali.sp$Station.Code==station,], type="l")
  points(Temperature~Date,cali.sp[cali.sp$Station.Code==station
                                   & !is.na(cali.sp$TemperatureIM==1),], col=2, pch=18)
  #readline("Continue?")
}

ggplot(cali.sp@data) + 
  geom_line(aes(Date,Temperature, group=Station.Code, colour=(TemperatureFlag==0)), alpha=0.5) + 
  theme(legend.position="none")
