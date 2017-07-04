library("spTimer")

# (Save) old graphic configuration
.pardefault <- par(no.readonly = T)

# Read data
data("NYdata")

# Simple diagnostic plot
NYdata$date <- as.POSIXct(strptime( sprintf("%04d-%02d-%02d",NYdata$Year,NYdata$Month,NYdata$Day),
                                    format="%Y-%m-%d", tz="GMT"))
summaryPlot(NYdata[,c(-1:-6)], period = "months")


# Figure 7
coords <- as.matrix(unique(cbind(NYdata[, 2:3])))
map(database = "state", regions = "new york", mar=par("mar"))
points(coords, pch = 19, col = 3)
points(coords, pch = 1, col = 1)
legend(x = -77.5, y = 41.5, col = c(3, 4), pch = c(19, 3), cex = 0.8, legend = c("Fitted sites", 
                                                                                 "Validation sites"))
par(.pardefault)


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