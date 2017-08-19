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
library(spacetime)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()


## Read data #########################################################################
ozone <- readRDS("data/epa/epa_daily/2016/california_ozone.RDS")
temp <- readRDS("data/epa/epa_daily/2016/california_temperature.RDS")
wind <- readRDS("data/epa/epa_daily/2016/california_wind.RDS")
rh <- readRDS("data/epa/epa_daily/2016/california_rh.RDS")
ozone.sites <- getSites(ozone)
temp.sites <- getSites(temp)
wind.sites <- getSites(wind)
rh.sites <- getSites(rh)
## Read sites
#sites <- readRDS("data/epa/sites/aqs_sites.RDS")
#sites <- sites[sites$State.Name=="California",]


## Merge with the other variables ########################################################
ozone <- fillMissingDates(ozone)
## Implicit filling of missing dates in the other data frames
cali <- merge(ozone, temp, all.x=T)
cali <- merge(cali, wind, all.x=T)
cali <- merge(cali, rh, all.x=T)
#View(cali)
summary(cali)
days <- seq(from=min(cali$Date), to=max(cali$Date), by='days')



## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Convert to Spatial data ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cali<-cali[order(cali$Station.Code,cali$Date),]  ## Not necessary
cali.sp <- convertDataToSp(cali)
#View(cali.sp)

## Variogram - temperature ###################################################
sp <- SpatialPoints(unique(coordinates(cali.sp)))
time <- unique(cali.sp$Date) 
cali.st <- STFDF(sp,time,cali)

v <- variogramST(Temperature~1, cali.st, cutoff = 700, tlags=0:3)
printPlot(paper,"img/preprocessing/variogram_temp.jpeg",7,5,FUN=function(){
  print(plot(v, map=F))
})




## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Analyze Missing Values ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Check missing rows
#md.pattern(cali)
printPlot(paper,"img/preprocessing/missing_epa.jpeg",7,5,FUN=function(){
  aggr(cali[,-(1:2)], gap=3, cex.axis=0.8)
})

printPlot(paper,"img/preprocessing/heatmap_ozone.jpeg",6,6,FUN=function(){
  print(levelplot(Ozone~Date*Station.Code,cali,
                  cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                  ylab="Station code",
                  scales=list(y=list(draw=FALSE))))
})
printPlot(paper,"img/preprocessing/heatmap_temperature.jpeg",6,6,FUN=function(){
  print(levelplot(Temperature~Date*Station.Code,cali,
                  cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                  ylab="Station code",
                  scales=list(y=list(draw=FALSE))))
})
printPlot(paper,"img/preprocessing/heatmap_wind.jpeg",6,6,FUN=function(){
  print(levelplot(Wind.speed~Date*Station.Code,cali,
                  cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                  ylab="Station code",
                  scales=list(y=list(draw=FALSE))))
})
printPlot(paper,"img/preprocessing/heatmap_rh.jpeg",6,6,FUN=function(){
  print(levelplot(RH~Date*Station.Code,cali,
                  cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                  ylab="Station code",
                  scales=list(y=list(draw=FALSE))))
})



## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Impute Temperature ####
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tempNasByStation <- aggregate(Temperature~Station.Code,cali,function(x){sum(is.na(x))}, na.action=na.pass)
tempNasByStation <- tempNasByStation[tempNasByStation$Temperature>0,] 
#View(tempNasByStation)
  

## Visualize Missing and Available Stations
naStations <- tempNasByStation[tempNasByStation$Temperature>100,1]
temp.sites.na <- getSites(naStations)

temp.sites.all <- rbind(cbind(temp.sites,Type="Available"),
                        cbind(temp.sites.na,Type="Missing"))

printPlot(paper,"img/preprocessing/stations_temp.jpeg",6,6,FUN=function(){
  p<-ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = temp.sites.all, 
               aes(x = Longitude, y = Latitude, fill=Type),
               alpha = .75, shape=24, size=2) +
    # geom_text(data = ozone.sites[ozone.sites$Station.Code=="089-0009",], 
    #           aes(x = Longitude, y = Latitude, label=Station.Code)) +
    labs(x = "Longitude", y = "Latitude", fill = "Type") +
    coord_quickmap() +
    theme(legend.justification = c("right", "top"), legend.position = c(.95, .95), 
          legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
  print(p)
})

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

if(T){
  printPlot(paper,paste0("img/preprocessing/station_",s,".jpeg"),6,6,FUN=function(){
    p<-ggplot(getCAmap()) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
      geom_point(data = temp.sites.all, 
                 aes(x = Longitude, y = Latitude, fill=Type),
                 alpha = .75, shape=24, size=2) +
      geom_point(data = temp.sites.all[temp.sites.all$Station.Code==s,], 
                 aes(x = Longitude, y = Latitude), alpha = .75, shape=21, size=3, fill="red") +
      labs(x = "Longitude", y = "Latitude", fill = "Type") +
      coord_quickmap() +
      theme(legend.justification = c("right", "top"), legend.position = c(.95, .95), 
            legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
    print(p)
  })
  # dim(temp.sp[temp.sp$Station.Code!=s,])
  # dim(temp.sp[temp.sp$Station.Code==s,])
  
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
    a <- idw(Temperature ~ Ozone, 
             temp.sp[temp.sp$Date==days[i] & !is.na(temp.sp$Temperature),],  # IN-SAMPLE
             # temp.sp[temp.sp$Station.Code!=s & temp.sp$Date==days[i],],  # OUT_OF_SAMPLE
             newdata=test[test$Date==days[i],], 
             idp=2.0, debug.level=0)
    test$Temperature2[test$Date==days[i]] <- a$var1.pred
  }
  printPlot(paper,paste0("img/preprocessing/station_",s,"_ts.jpeg"),5,5,FUN=function(){
    plot(Temperature~Date,test, type="l", ylim=range(c(test$Temperature,test$Temperature2),na.rm=T))
    lines(Temperature2~Date,test, col=2, lty="dashed")
  })
  evaluatePredictions(test$Temperature, test$Temperature2)
}

## Notes
## Use imputation methods for filling missing data highly depends on the position of the station



###
### FAIL
###




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
             temp.sp[temp.sp$Date==days[i] & !is.na(temp.sp$Temperature),], 
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
