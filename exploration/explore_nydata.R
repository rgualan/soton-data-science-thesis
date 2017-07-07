library(maps)
library(openair)
library(spTimer)


# (Save) old graphic configuration
.pardefault <- par(no.readonly = T)


#### Read data and calculate basic fields ################################################
data("NYdata")
load(file="data/ny_ozone/NYcheapLoc.Rdata")

# New fields
NYdata$date <- as.POSIXct(strptime( sprintf("%04d-%02d-%02d",NYdata$Year,NYdata$Month,NYdata$Day),
                                    format="%Y-%m-%d", tz="GMT"))
#coordinates(NYdata) <- ~Longitude+Latitude

# Simple diagnostic plot
summaryPlot(NYdata[,c("o8hrmax","cMAXTMP","WDSP","RH","date")], period="months")
readline("Continue?")

# Save to file for future reference
save(NYdata, file="data/ny_ozone/NYdata.Rdata")


#### NY map plus stations ################################################################
coords <- as.matrix(unique(cbind(NYdata[, 2:3])))
map(database = "state", regions = "new york", mar=par("mar"))
points(coords, pch = 19, col = 3)
points(coords, pch = 1, col = 1)
legend(x = -77.5, y = 41.5, col = c(3, 4), pch = c(19, 3), cex = 0.8, legend = c("Fitted sites", 
                                                                                 "Validation sites"))
stations <- unique(NYdata[,c("s.index", "Longitude", "Latitude")])
text(stations$Longitude+0.3, stations$Latitude, stations$s.index, cex=0.7)
#with(stations[stations$s.index==7,], text(Longitude+0.3, Latitude, s.index, cex=0.7, col=2))
readline("Continue?")

par(.pardefault)

#### Only stations ######################################################
stations <- unique(NYdata[,c("s.index", "Longitude", "Latitude")])
plot(stations$Longitude, stations$Latitude, pch = 19, col = "lightgreen", cex=3, axes=F, xlab="", ylab="", asp=1)
points(stations$Longitude, stations$Latitude, pch = 1, col = 1, cex=3)
text(stations$Longitude, stations$Latitude, stations$s.index, cex=1) #pos=3
#with(stations[stations$s.index==12,], text(Longitude+0.3, Latitude, s.index, cex=0.7, col=2))
readline("Continue?")


#### References + Cheap sensors #########################################
stations <- unique(NYdata[,c("s.index", "Longitude", "Latitude")])

#jpeg('test.jpg')
plot(stations$Longitude, stations$Latitude, pch = 19, col = "lightgreen", cex=2, 
     axes=F, xlab="", ylab="", asp=1, ylim=c(min(stations$Latitude),max(stations$Latitude)+0.7))
points(stations$Longitude, stations$Latitude, pch = 1, col = 1, cex=2)
text(stations$Longitude, stations$Latitude, stations$s.index, cex=0.7) #pos=3

#Temporal IDS
#length(NYcheapLoc)
#NYcheapLoc$ID <- 1:51

points(coordinates(NYcheapLoc)[,1], coordinates(NYcheapLoc)[,2], pch = 19, col = "lightblue", cex=2)
points(coordinates(NYcheapLoc)[,1], coordinates(NYcheapLoc)[,2], pch = 1, col = 1, cex=2)
text(coordinates(NYcheapLoc)[,1], coordinates(NYcheapLoc)[,2], NYcheapLoc$s.index, cex=0.6) #pos=3

# dev.copy(png,'myplot.png')
#dev.off()

readline("Continue?")


##### Time series #####
sites <- unique(NYdata$s.index)
## Individual summary
# for(site in sites ){
#   print(site)
#   summaryPlot(NYdata[NYdata$s.index==site,c(-1:-6)], period = "months", 
#               main=paste("Station", site))
# }

## All in one
par(mfrow=c(4,1))
par(mar=c(0,4,0,2))
par(oma=c(2,0,0,0))

# Ozone
plot(NYdata[NYdata$s.index==1,]$o8hrmax, type = "l", xlab = "", ylab="Ozone", axes=F, 
     ylim=range(NYdata$o8hrmax, na.rm=T))
axis(1, labels = F)
axis(2)
for(index in 2:28 ){
  lines(NYdata[NYdata$s.index==index,]$o8hrmax, type = "l", col=index)
}
# Temperature
plot(NYdata[NYdata$s.index==1,]$cMAXTMP, type = "l", xlab = "", ylab="Temperature", axes=F,
     ylim=range(NYdata$cMAXTMP, na.rm=T))
axis(1, labels = F)
axis(2)
for(index in 2:28 ){
  lines(NYdata[NYdata$s.index==index,]$cMAXTMP, type = "l", col=index)
}
# Wind direction
plot(NYdata[NYdata$s.index==1,]$WDSP, type = "l", xlab = "", ylab="Wind direction", axes=F,
     ylim=range(NYdata$WDSP, na.rm=T))
axis(1, labels = F)
axis(2)
for(index in 2:28 ){
  lines(NYdata[NYdata$s.index==index,]$WDSP, type = "l", col=index)
}
# RH
plot(NYdata[NYdata$s.index==1,]$RH, type = "l", xlab = "", ylab="Relative humidity", axes=F,
     ylim=range(NYdata$RH, na.rm=T))
axis(1)
axis(2)
for(index in 2:28 ){
  lines(NYdata[NYdata$s.index==index,]$RH, type = "l", col=index)
}

par(.pardefault)




