library(automap)
library(sp)
library(spTimer)
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


#### Create the cheap-sensor dataset ################################################
# Assumptions:
# The original stations are reference data (data collected with high-precision instruments)
# Task: 
# Generate data from cheap-sensors using interpolation from the original data
# Algorithm:
# Interpolate data into the cheap sensors
# Add noise

# Load processed data and location of the "cheap sensors"
load(file="data/ny_ozone/NYdataFilled.Rdata")
load(file="data/ny_ozone/NYcheapLoc.Rdata")

# Remove problematic station
#NYcheapLoc <- NYcheapLoc[-(length(NYcheapLoc)-1), ]
#length(NYcheapLoc)
#length(NYcheapLoc2)


# Create an interpolation for each time
theDates <- unique(NYdata$date)

# Rm previously generated dataset in case of several executions
if(exists("NYcheap")) rm("NYcheap")

# i <- 5
for(i in 1:length(theDates)){
  print(i)
  print(theDates[i])  
  # Original data (points) in a specific date
  tmp <- NYdata[NYdata$date==theDates[i],] 
  coordinates(tmp) <- ~Longitude+Latitude
  #head(tmp)
  # New data (points)
  tmp2 <- data.frame(s.index=NYcheapLoc$s.index, 
                     Longitude=NYcheapLoc$x1,
                     Latitude=NYcheapLoc$x2,
                     date=rep(theDates[i], length(NYcheapLoc)))
  coordinates(tmp2) <- ~Longitude+Latitude
  #head(tmp2)
  # Kriging
  krigingcMAXTMP = autoKrige(cMAXTMP~1, tmp, tmp2)
  krigingcWDSP = autoKrige(WDSP~1, tmp, tmp2)
  krigingcRH = autoKrige(RH~1, tmp, tmp2)
  tmp2$cMAXTMP <- krigingcMAXTMP$krige_output$var1.pred
  tmp2$WDSP <- krigingcWDSP$krige_output$var1.pred
  tmp2$RH <- krigingcRH$krige_output$var1.pred
  if(sum(is.na(tmp2$cMAXTMP))+sum(is.na(tmp2$WDSP))+sum(is.na(tmp2$RH))<20){
    krigingco8hrmax = autoKrige(o8hrmax~cMAXTMP+WDSP+RH, tmp, tmp2)
    tmp2$o8hrmax <- krigingco8hrmax$krige_output$var1.pred
  }else{
    print("Warning: missing values. Using ordinary interpolation for o8hrmax")
    krigingco8hrmax = autoKrige(o8hrmax~1, tmp, tmp2)
    tmp2$o8hrmax <- krigingco8hrmax$krige_output$var1.pred
    #plot(krigingco8hrmax)
    #readline("Continue?")
  }
  
  # Assemble the new dataset
  if(!exists("NYcheap")){
    NYcheap <- tmp2
  }else{
    NYcheap <- rbind(NYcheap, tmp2)
  }
}

# Check
naByStationDate <- aggregate(cbind(o8hrmax,cMAXTMP,WDSP,RH)~date, data=NYcheap, 
                             FUN=function(x) sum(is.na(x)), na.action = na.pass)
#naByStationDate
#readline("Continue?")

# Note:
# There is a problem with the interpolation for some dates
# Eg: i=41 -> 2006-08-10
# Warning: 50: In predict.gstat(g, newdata = newdata, block = block,  ... :
# Covariance matrix singular at location [-73.712,44.5031,0]: skipping...

## Sort by station
NYcheap <- NYcheap[order(NYcheap$s.index, NYcheap$date), ]

## Add some very simple Gaussian noise
head(NYcheap); nrow(NYcheap)
NYcheapPlusNoise <- NYcheap

noiseProportion <- 0.1
biasO8hrmax <- 5
biasCmaxtmp <- 3
biasWdsp <- 1
biasRh <- 1

NYcheapPlusNoise$o8hrmax <- NYcheap$o8hrmax + biasO8hrmax +
  rnorm(nrow(NYcheap), mean=0, sd=sd(NYcheap$o8hrmax, na.rm = T)*noiseProportion)
NYcheapPlusNoise$cMAXTMP <- NYcheap$cMAXTMP + biasCmaxtmp +
  rnorm(nrow(NYcheap), mean=0, sd=sd(NYcheap$cMAXTMP, na.rm = T)*noiseProportion)
NYcheapPlusNoise$WDSP <- NYcheap$WDSP + biasWdsp +
  rnorm(nrow(NYcheap), mean=0, sd=sd(NYcheap$WDSP, na.rm = T)*noiseProportion)
NYcheapPlusNoise$RH <- NYcheap$RH + biasRh + 
  rnorm(nrow(NYcheap), mean=0, sd=sd(NYcheap$RH, na.rm = T)*noiseProportion)

# plot(NYcheap$o8hrmax, type="l")
# lines(NYcheapPlusNoise$o8hrmax, col=2)


##### Spatial plot #################################################################

coordinates(NYdata) <- ~Longitude+Latitude
spplot(NYdata[NYdata$date==theDates[i],], "o8hrmax",
       colorkey=TRUE, main=format(theDates[i], "%Y-%m-%d"))

spplot(NYcheap[NYcheap$date==theDates[i],], "o8hrmax",
       colorkey=TRUE, main=format(theDates[i], "%Y-%m-%d"))
readline("Continue?")


##### Time series for one station ##################################################
par(mfrow=c(4,1), mar=c(0,4,0,2), oma=c(2,0,0,0))

# Ozone
plot(NYcheap$o8hrmax[NYcheap$s.index==100], type = "l", xlab = "", ylab="Ozone", axes=F)
axis(1, labels = F); axis(2)
for(index in 102:149 ){
  lines(NYcheap[NYcheap$s.index==index,]$o8hrmax, type = "l", col=index)
}
# Temperature
plot(NYcheap[NYcheap$s.index==100,]$cMAXTMP, type = "l", xlab = "", ylab="Temperature", axes=F,
     ylim=range(NYcheap$cMAXTMP, na.rm=T))
axis(1, labels = F); axis(2)
for(index in 102:149 ){
  lines(NYcheap[NYcheap$s.index==index,]$cMAXTMP, type = "l", col=index)
}
# Wind direction
plot(NYcheap[NYcheap$s.index==100,]$WDSP, type = "l", xlab = "", ylab="Wind direction", axes=F,
     ylim=range(NYcheap$WDSP, na.rm=T))
axis(1, labels = F); axis(2)
for(index in 102:149 ){
  lines(NYcheap[NYcheap$s.index==index,]$WDSP, type = "l", col=index)
}
# RH
plot(NYcheap[NYcheap$s.index==100,]$RH, type = "l", xlab = "", ylab="Relative humidity", axes=F,
     ylim=range(NYcheap$RH, na.rm=T))
axis(1); axis(2)
for(index in 102:149 ){
  lines(NYcheap[NYcheap$s.index==index,]$RH, type = "l", col=index)
}


par(.pardefault)
readline("Continue?")


# Save generated data
save(NYcheap, file="data/ny_ozone/NYcheap.Rdata")
save(NYcheapPlusNoise, file="data/ny_ozone/NYcheapPlusNoise.Rdata")


## Check a Reference station and a close Cheap station
# 26 - 130 - 101
plot(NYdata$o8hrmax[NYdata$s.index==26], type="l", xlab="", ylab="")
lines(NYcheapPlusNoise$o8hrmax[NYcheapPlusNoise$s.index==130], col=2)
#lines(NYcheapPlusNoise$o8hrmax[NYcheapPlusNoise$s.index==101], col=3)


