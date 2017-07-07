library(automap)
library(spTimer)
# suppressPackageStartupMessages({
# })
library(ggplot2)
library(magrittr)
library(mgcv)
library(spdep)

# Load pseudo-locations of the cheap sensors
load(file="data/ny_ozone/NYdata.Rdata")

#### Fill NA data #########################################################################
# Algorithm:
# First, fill stations with small amount of NA
# Find closest station
# Build ds
# Train and apply model: o8hrmax~cMAXTMP+WDSP+RH+Month+ohrmaxSclose

# Create a flag field
NYdata$filled <- rep(0,nrow(NYdata))


# Stations with missing data
# unique(NYdata$s.index[is.na(NYdata$o8hrmax)])
# all stations
stations <- unique(NYdata[,c("s.index", "Longitude", "Latitude")])
rownames(stations) <- stations$s.index
coordinates(stations) <- ~Longitude+Latitude
# Find closest station (knn, k=1)
stations$knn1 <- knearneigh(coordinates(stations), k=1, longlat=T)$nn[,1]
# incomplete stations:
iStations <- aggregate(o8hrmax~s.index+Longitude+Latitude, data=NYdata[is.na(NYdata$o8hrmax),], 
                       FUN=function(x) sum(is.na(x)), na.action = na.pass)
rownames(iStations) <- iStations$s.index
iStations <- iStations[order(-iStations$o8hrmax), ]
iStations <- merge(iStations,stations[stations$s.index %in% iStations$s.index,]@data,by="s.index")
iStations <- iStations[order(iStations$o8hrmax),]

# Two methods are used:
# 1) Gaussian process trained using the whole dataset. It was the best "general" model
# 2) Tailored LM + the closest station. It offers very high accuracy 

# Train GP model (General)
trainDs2 <- NYdata[!(NYdata$s.index %in% unique(NYdata$s.index[is.na(NYdata$o8hrmax)])), ]
simpleGp <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = trainDs2, model = "GP", 
                      coords = ~Longitude + Latitude, scale.transform = "SQRT", 
                      spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
pred.gp <- fitted(simpleGp)^2
spT.validation(trainDs2$o8hrmax,pred.gp[,1])
#plot(trainDs2$o8hrmax, type="l")
#lines(pred.gp[,1], col=2)


# Algorithm:
# Train individual LMs
# Chose the best model GP versus LM
# And fill the missing values
# By station
# First the stations with less missing values???

#i <- 6
#i=which(iStations$s.index==20)
for(i in 1:nrow(iStations)){
  paste("Processing: ", iStations[i,]$s.index)
  
  # Linear model
  trainSx <- NYdata[NYdata$s.index==iStations[i,]$s.index, ]
  trainSx$ohrmaxSn <- NYdata$o8hrmax[NYdata$s.index==iStations[i,]$knn1]
  simpleLmSx <- lm(o8hrmax~cMAXTMP+WDSP+RH+Month+ohrmaxSn, data = trainSx)
  spT.validation(trainSx$o8hrmax, predict(simpleLmSx, newdata = trainSx))
  cor(trainSx$o8hrmax, predict(simpleLmSx, newdata = trainSx), use="pairwise.complete.obs")
  # Basic ploting
  plot(trainSx$date, trainSx$o8hrmax, type="l", xlab="", ylab="", main=iStations[i,]$s.index)
  abline(v=trainSx$date[is.na(trainSx$o8hrmax)], lty="dashed", col="gray")
  lines(trainSx$date, predict(simpleLmSx, newdata = trainSx), col="orange", lwd=1) 
  pred.gp.Sx <- predict(simpleGp, newdata=trainSx, newcoords = ~Longitude + Latitude)
  lines(trainSx$date, pred.gp.Sx$Mean, col="blue", lwd=1)
  # Model evaluation
  spT.validation(trainSx$o8hrmax, predict(simpleLmSx, newdata = trainSx))
  spT.validation(trainSx$o8hrmax, pred.gp.Sx$Mean)
  r2_lm<-cor(trainSx$o8hrmax, predict(simpleLmSx, newdata = trainSx), use = "pairwise.complete.obs")^2
  r2_gp<-cor(trainSx$o8hrmax, pred.gp.Sx$Mean, use = "pairwise.complete.obs")^2
  c(r2_lm,r2_gp)
  # Filling
  if(r2_lm>r2_gp){
    fillSx <- trainSx[is.na(trainSx$o8hrmax), ]
    filling <- predict(simpleLmSx, newdata = fillSx)
    NYdata$filled[NYdata$s.index==iStations$s.index[i] & is.na(NYdata$o8hrmax)] <- 1
    NYdata$o8hrmax[NYdata$s.index==iStations$s.index[i] & is.na(NYdata$o8hrmax)] <- filling
    lines(trainSx$date, NYdata$o8hrmax[NYdata$s.index==iStations$s.index[i]], type="l", col="pink", lwd=2)
  }else{
    pred.gp.Sx <- predict(simpleGp, newdata=trainSx, newcoords = ~Longitude + Latitude)
    filling <- pred.gp.Sx$Mean[is.na(trainSx$o8hrmax)]   
    NYdata$filled[NYdata$s.index==iStations$s.index[i] & is.na(NYdata$o8hrmax)] <- 1
    NYdata$o8hrmax[NYdata$s.index==iStations$s.index[i] & is.na(NYdata$o8hrmax)] <- filling
    lines(trainSx$date, NYdata$o8hrmax[NYdata$s.index==iStations$s.index[i]], type="l", col="purple", lwd=2)
  }
  readline("Continue?")
}

# Note:
# There was a problem in station 20
# Because the closest station was also missing one value 



#### Check the filled points #########################################################
#i <- 7
ds <- NYdata[NYdata$s.index==stations$s.index[1],]
plot(ds$date, ds$o8hrmax, type="l", xlab="", ylab="", col=1, ylim=range(NYdata$o8hrmax))
points(ds$date[ds$filled==1], ds$o8hrmax[ds$filled==1], col=1, cex=1.5, pch=19) 
points(ds$date[ds$filled==1], ds$o8hrmax[ds$filled==1], col=1, cex=1.5, pch=1) 
for(i in 2:nrow(stations)){
  ds <- NYdata[NYdata$s.index==stations$s.index[i],]
  lines(ds$date, ds$o8hrmax, col=i)
  points(ds$date[ds$filled==1], ds$o8hrmax[ds$filled==1], col=i, cex=1.5, pch=19) 
  points(ds$date[ds$filled==1], ds$o8hrmax[ds$filled==1], col=1, cex=1.5, pch=1) 
}

# Save processed data
save(NYdata, file="data/ny_ozone/NYdataFilled.Rdata")



