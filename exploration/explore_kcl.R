# Clean environment #################################################################
rm(list=ls())

# Libraries
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgeos)
library(ggmap)
library(rgdal)


## Stations ###################################################################
## Read data
load("/home/ronald/projects/back_up/data/kcl/sites.RData")
nrow(sites) # 900

## Plot all stations
plot(Latitude~Longitude,sites,pch=2,col=Classification)
#plot(os_grid_y~os_grid_x,sites)
readline("Continue?")

## Plot active stations
dateA <- as.POSIXct("2015-01-01", format="%Y-%m-%d", tz="UTC")
dateB <- as.POSIXct("2015-04-30", format="%Y-%m-%d", tz="UTC")
sites2 <- sites[sites$OpeningDate<=dateA & 
                  (sites$ClosingDate>=dateB | is.na(sites$ClosingDate)),]
nrow(sites2)
plot(Latitude~Longitude,sites2,pch=2,col=Classification)
readline("Continue?")

## Plot a map with the active stations ###################################
## Refs: 
## http://www.milanor.net/blog/maps-in-r-plotting-data-points-on-a-map/
## https://blog.dominodatalab.com/geographic-visualization-with-rs-ggmaps/

## Plot map of London
mapLondon <- get_map(location = 'London', zoom = "auto", maptype = "roadmap")
#save(mapLondon, file="data/maps/mapLondon.Rdata")
#if(!exists("mapLondon")) load("data/maps/mapLondon.Rdata")

ggmap(mapLondon) +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = Classification), 
             data = sites2, alpha = .75, shape=21, size=2)
readline("Continue?")

## How many stations in London/Greater London?
sites3 <- sites2[sites2$Longitude >= -0.567 & sites2$Longitude <= 0.312
                 & sites2$Latitude >= 51.2 & sites2$Latitude <= 51.8,]
nrow(sites3)



## Explore meteorological dataset
load("/home/ronald/projects/back_up/data/kcl/metData.RData")
met2 <- met[met$date>=dateA & met$date<=dateB,]
nrow(met2)
#View(met2)
summaryPlot(met2, period = "months")




## Create a compact data frame ###################################################
createCompact = F

if(createCompact){
  variables <- c("site","date","nox","no2","so2","pm10_raw","pm10","pm25","o3","co","co2")
  variables2 <- c("site","date","Latitude","Longitude","Classification","nox","no2","so2","pm10_raw","pm10","pm25","o3","co","co2")
  
  files <- list.files("/home/ronald/projects/back_up/data/kcl/2015")
  for(f in files){
    print(f)
    load(paste0("data/kcl/2015/",f))
    
    #names(x)
    names(x) <- sapply(names(x),FUN=tolower)
    
    ## Remove irrelevant variables
    x <- x[ , (names(x) %in% variables)]
    #names(x)
    
    ## Add NA variables
    for(varName in variables){
      #cat(varName,"\n")
      if (!varName %in% names(x)){
        x$new <- NA
        names(x)[ncol(x)] <- varName
      }
    }
    #head(x)
    
    if(!exists("compact")){
      compact <- x
    }else{
      compact <- rbind(compact,x)  
    }
    
  }
  
  ## Merge with stations
  compact <- merge(compact,
                    sites[,c("SiteCode","Classification","Latitude","Longitude")],
                    by.x="site", by.y="SiteCode")
  
  ## Save compact dataframe
  d$site <- as.factor(d$site)
  compact <- compact[order(compact$site,compact$date),variables2]
  save(compact, file="data/kcl/compact2015.RData")
}


load("/home/ronald/projects/back_up/data/kcl/compact2015.RData")
d <- compact



## Plot data as a matrix (sites x date) ################
library(lattice)
library(RColorBrewer)

st.on <- sort(unique(d$site))
NS<-50

for(i in 1:ceiling(length(st.on)/50)){
  cat(i) #i<-1
  print(levelplot(ozone~date*site, 
                  d[d$site
                       %in% st.on[(NS*(i-1)+1):(NS*i)],],
                  cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                  scales=list(y=list(cex=.7))))
  readline("Continue?")
}




# Active stations (by pollutant)
nasBySite <- aggregate(cbind(no2,pm10)~site,d,FUN=function(x){sum(!is.na(x))}, na.action=na.pass)
View(nasBySite)
tot <- as.double(max(d$date)-min(d$date), units="hours")
nasBySite$no2p <- nasBySite$no2/tot
nasBySite$pm10p <- nasBySite$pm10/tot
summaryPlot(d[d$site=="HK6",c("site","date","no2")])
summaryPlot(d[d$site=="BT4",c("site","date","no2")])

sites.on <- nasBySite$site[nasBySite$no2p>=0.90]
sites.on <- merge(data.frame(site=sites.on), 
                  sites[,c("SiteCode","Longitude","Latitude","Classification")], 
                  by.x="site", by.y="SiteCode")
nrow(sites.on)
## 47 for pm10
## 55 for no2
plot(Latitude~Longitude,sites.on,pch="|")
points(Latitude~Longitude,sites.on,pch="2",col="red")
