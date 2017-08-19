## Clean environment
rm(list=ls())

### Load required packages ###
library(fields)
library(raster)
library(spatial.tools)
library(gdalUtils)
library(rgdal)
library(gstat)
library(automap)
source("util/my_helper.R")
source("preprocessing/gbmImpute.R")

## Global variables
paper <- setupPaper()

### Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_cov.RDS")
sites <- getSites(epa)
str(epa)

## Scale dependent variable
epa$sOzone <- scale(epa$Ozone)[,1]  ## scale returns a matrix
# hist(epa$Ozone)
# hist(epa$sOzone)
dim(epa)
head(epa)

## Count missing values per station
ozoneNasByStation <- aggregate(Ozone~Station.Code,epa,function(x){sum(is.na(x))}, na.action=na.pass)
ozoneNasByStation <- ozoneNasByStation[ozoneNasByStation$Ozone>5,] 
ozoneNasByStation <- ozoneNasByStation[order(ozoneNasByStation$Ozone,decreasing=T),]

testStation <- "065-2002"
data <- epa
data$sOzone2 <- data$sOzone
data$sOzone2[data$Station.Code==testStation & 
               data$Date>="2016-05-01" & data$Date<"2016-05-15"] <- NA ## WARNING!!!
data1 <- data[,c("sOzone2","Temperature","Wind.speed","UTM.X","UTM.Y","Elevation","Location.Setting")]
st<-Sys.time()
fit <- gbmImpute(data1, max.iters = 2, cv.fold = 5, verbose=T)
Sys.time()-st

sum(is.na(fit$x))
head(data1)
head(fit$x)

data2 <- cbind(Station.Code=data$Station.Code, Date=data$Date,
               sOzone=data$sOzone, fit$x)
tmp <- data2[data2$Station.Code==testStation,]
tmp2 <- rbind(data.frame(Date=tmp$Date,Ozone=tmp$sOzone,Type="Original"),
              data.frame(Date=tmp$Date,Ozone=tmp$sOzone2,Type="Interpolated"))
ggplot(tmp2, aes(x=Date, y=Ozone, colour=Type)) + 
    # annotate("rect",
    #          xmin=tmp2$Date[is.na(tmp2$Ozone)]-1*24*60*60,
    #          xmax=tmp2$Date[is.na(tmp2$Ozone)]+1*24*60*60,
    #          ymin=-Inf, ymax=Inf, alpha=0.75, fill="yellow") +
    geom_line() + 
    theme(legend.justification = c("top")) + 
    labs(y="Scaled(Ozone)")

##  Imputation is not Interpolation
  


## Data interpolation
metrics = c()
st <- Sys.time()
for(i in 1:nrow(ozoneNasByStation)){
  testStation <- ozoneNasByStation$Station.Code[i]
  print(paste("Fold", i, ":", testStation,paste(rep("=",50), collapse = "")))

  train <- epa.sp[epa.sp$Station.Code!=testStation,] 
  test <- epa.sp[epa.sp$Station.Code==testStation,]
  test$sOzoneH <- NA
  
  for(j in 1:length(days)){
    day <- days[j]
    #print(day)
    slice.train <- train[train$Date==day & !is.na(train$sOzone),]  
    slice.test <- test[test$Date==day,]

    fit <- idw(sOzone~Elevation, slice.train,
             newdata=slice.test,
             debug.level=0)
    
    test$sOzoneH[test$Date==day] <- fit$var1.pred
    epa.sp$sOzoneH[epa.sp$Station.Code==testStation & epa.sp$Date==day] <- fit$var1.pred
  }
  
  m <- evaluatePredictions(test$sOzone, test$sOzoneH)
  metrics <- rbind(metrics,m)
  rownames(metrics)[nrow(metrics)] <- as.character(testStation)
  
  if(m[8]<0.4 & F){
    # Diagnostic
    print(m)
    plot(test$sOzone,type="l") 
    lines(test$sOzoneH,col=2)
    readline("Continue?")
  }
}
Sys.time()-st

## Summary metrics
metricsDf <- data.frame(mean=apply(metrics,2,mean),
                        min=apply(metrics,2,min),
                        max=apply(metrics,2,max))

## Manual assessment ##################################################################
#plot(metrics[order(metrics[,8]),8])
#head(metrics[order(metrics[,8], decreasing = F),])
#View(metrics[order(metrics[,8], decreasing = F),])

ss<-c("021-0003", # >> min RMSE
      "027-0101", # >> max RMSE
      "083-4003" ) # >> worst R2

s <- ss[1] #019-2009", "083-4003"
#highlightStation(s); 
#getKneighbours(s, sites, 5)
n5 <- merge(getKneighbours(s, sites, 5, T),sites); n5[order(n5$distance),]
#getKneighboursInRadius(s, sites, 100)
nd <- merge(getKneighboursInRadius(s, sites, 75, T),sites); nd[order(nd$distance),];
range(nd$Elevation)
## Time series
tsPredictionPlot <- function(paper,epa.sp,prefix,s){
  tmp <- epa.sp[epa.sp$Station.Code==s,]
  tmp2 <- rbind(data.frame(Date=tmp$Date,Ozone=tmp$sOzone,Type="Original"),
                data.frame(Date=tmp$Date,Ozone=tmp$sOzoneH,Type="Interpolated"))
  
  printPlot(paper, paste0(prefix,s,".jpeg"),7,3, FUN=function(){ 
    p<-ggplot(tmp2, aes(x=Date, y=Ozone, colour=Type)) + 
      annotate("rect",
               xmin=tmp2$Date[is.na(tmp2$Ozone)]-1*24*60*60,
               xmax=tmp2$Date[is.na(tmp2$Ozone)]+1*24*60*60,
               ymin=-Inf, ymax=Inf, alpha=0.75, fill="yellow") +
      geom_line() + 
      theme(legend.justification = c("top")) + 
      labs(y="Scaled(Ozone)")
    print(p)
  })
}

tsPredictionPlot(paper,epa.sp,"img/preprocessing/idw/ozone_",ss[1])
tsPredictionPlot(paper,epa.sp,"img/preprocessing/idw/ozone_",ss[2])
tsPredictionPlot(paper,epa.sp,"img/preprocessing/idw/ozone_",ss[3])

## Three cases 
printPlot(paper, "img/preprocessing/idw/cases.jpeg",6,6,FUN=function(){
  p<-ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = sites[!sites$Station.Code %in% ss,], aes(x = Longitude, y = Latitude), #, fill=Location.Setting
               alpha = 0.75, shape=4, size=1) +
    geom_point(data = sites[sites$Station.Code %in% ss,], 
               aes(x = Longitude, y = Latitude, fill=Station.Code),
               alpha = 1, shape=24, size=2) +
    labs(x = "Longitude", y = "Latitude", col="Station") +
    coord_quickmap() +
    theme(legend.justification = c("right", "top"), legend.position = c(.95, .95),
          legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
  print(p)
})


# trend <- aggregate(sOzone~Date, epa.sp[epa.sp$Station.Code %in% kneighbours$Station.Code,], mean)
# lines(sOzone~Date, trend, col="blue", lwd=1)
# cor(trend$sOzone,epa.sp$sOzone[epa.sp$Station.Code==s],use="pairwise.complete.obs")^2
