library("spTimer")

# These packages will be required to run the code in this file
library("akima")
library("coda")
library("spacetime")
library("fields")
library("forecast")
library("MASS")
library("mgcv")
library("spBayes")
library("colorspace")
library("maps")
library("MBA")
library("openair")

start.time <- Sys.time()
.pardefault <- par(no.readonly = T)

##### Real life Example #####
# Target variable: ground level ozone
# Predictors: cMAXTEMP (max temperature in degree Celsius), WDSP (wind speed in nautical miles) 
# and RH (percentage average relative humidity)
# Num Stations: 28 in the state of New York
# Num of stations for training: 20
# Num of stations for validation: 8

# Read data
data("NYdata")
s <- c(8, 11, 12, 14, 18, 21, 24, 28)
DataFit <- spT.subset(data = NYdata, var.name = c("s.index"), s = s, reverse = TRUE)
DataFit <- subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred <- spT.subset(data = NYdata, var.name = c("s.index"), s = s)
DataValPred <- subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 
  8)))

# Create an alternative dataset and inject noise 
NYdataNoise <- NYdata
noiseProportion <- 1/2
NYdataNoise$o8hrmax <- NYdata$o8hrmax + 
  rnorm(nrow(NYdata), mean=0, sd=sd(NYdata$o8hrmax, na.rm = T)*noiseProportion)
NYdataNoise$cMAXTMP <- NYdata$cMAXTMP + 
  rnorm(nrow(NYdata), mean=0, sd=sd(NYdata$cMAXTMP, na.rm = T)*noiseProportion)
NYdataNoise$WDSP <- NYdata$WDSP + 
  rnorm(nrow(NYdata), mean=0, sd=sd(NYdata$WDSP, na.rm = T)*noiseProportion)
NYdataNoise$RH <- NYdata$RH + 
  rnorm(nrow(NYdata), mean=0, sd=sd(NYdata$RH, na.rm = T)*noiseProportion)
plot(NYdata$o8hrmax[1:100], type="l")
lines(NYdataNoise$o8hrmax[1:100], lty="dashed", col=2)
plot(NYdata$cMAXTMP[1:100], type="l")
lines(NYdataNoise$cMAXTMP[1:100], lty="dashed", col=2)
plot(NYdata$WDSP[1:100], type="l")
lines(NYdataNoise$WDSP[1:100], lty="dashed", col=2)
plot(NYdata$RH[1:100], type="l")
lines(NYdataNoise$RH[1:100], lty="dashed", col=2)

DataFitNoise <- spT.subset(data = NYdataNoise, var.name = c("s.index"), s = s, reverse = TRUE)
DataFitNoise <- subset(DataFitNoise, with(DataFitNoise, !(Day %in% c(30, 31) & Month == 8)))
DataValPredNoise <- spT.subset(data = NYdataNoise, var.name = c("s.index"), s = s)
DataValPredNoise <- subset(DataValPredNoise, with(DataValPredNoise, !(Day %in% c(30, 31) & Month == 
                                                         8)))


# Simple diagnostic plot
NYdata$date <- as.POSIXct(strptime( sprintf("%04d-%02d-%02d",NYdata$Year,NYdata$Month,NYdata$Day),
                                    format="%Y-%m-%d", tz="GMT"))
NYdataNoise$date <- as.POSIXct(strptime( sprintf("%04d-%02d-%02d",NYdataNoise$Year,NYdataNoise$Month,NYdataNoise$Day),
                                    format="%Y-%m-%d", tz="GMT"))
summaryPlot(NYdata[,c(-1:-6)], period = "months")
summaryPlot(NYdataNoise[,c(-1:-6)], period = "months")


# Figure 7
coords <- as.matrix(unique(cbind(DataFit[, 2:3])))
pred.coords <- as.matrix(unique(cbind(DataValPred[, 2:3])))
map(database = "state", regions = "new york")
points(coords, pch = 19, col = 3)
points(coords, pch = 1, col = 1)
points(pred.coords, pch = 3, col = 4)
legend(x = -77.5, y = 41.5, col = c(3, 4), pch = c(19, 3), cex = 0.8, legend = c("Fitted sites", 
  "Validation sites"))
par(.pardefault)

## Fit GP model
set.seed(11)
post.gp <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = DataFit, model = "GP", 
                     coords = ~Longitude + Latitude, scale.transform = "SQRT", 
                     spatial.decay = spT.decay(distribution = Gamm(2,1), tuning = 0.1))
post.gp2 <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = DataFitNoise, model = "GP", 
                     coords = ~Longitude + Latitude, scale.transform = "SQRT", 
                     spatial.decay = spT.decay(distribution = Gamm(2,1), tuning = 0.1))
print(post.gp)
print(post.gp2)
summary(post.gp)
summary(post.gp2)
# The difference in the parameters is concerning!
# plot(post.gp) #It works with dev.new()
# plot(post.gp2) #It works with dev.new()

# Spatial prediction for the GP model
set.seed(11)
pred.gp <- predict(post.gp, newdata = DataValPred, newcoords = ~Longitude + Latitude)
pred.gp2 <- predict(post.gp2, newdata = DataValPredNoise, newcoords = ~Longitude + Latitude)
print(pred.gp)
names(pred.gp)
# validation criteria
spT.validation(DataValPred$o8hrmax, c(pred.gp$Median))
spT.validation(DataValPred$o8hrmax, c(pred.gp2$Median))
cor(na.omit(cbind(DataValPred$o8hrmax, c(pred.gp$Median))))^2
cor(na.omit(cbind(DataValPred$o8hrmax, c(pred.gp2$Median))))^2

# Time-series plot
s.index <- 8
#subIndex <- 1:nrow(DataValPred)
subIndex <- 1:100
plot(DataValPred$o8hrmax[subIndex], type = "o", pch = 16, ylab = "y", 
     cex = 0.9, lty = 2, lwd = 1.5)
lines(DataValPredNoise$o8hrmax[subIndex], type = "o", pch = 16,  
     cex = 0.9, lty = 2, lwd = 1.5, col=4)
lines(c(pred.gp$Median[subIndex]), type = "o", pch = 1, 
     cex = 0.9, lty = 2, lwd = 1.5, col=2)
lines(c(pred.gp2$Median[subIndex]), type = "o", pch = 1, 
      cex = 0.9, lty = 2, lwd = 1.5, col=3)
# legend("topleft", pch = c(16, 16, 1, 1), lty = rep(2, 4), col = c(1, 4, 2, 3),
#        legend = c("Observed", "Observed+noise", "Model 1", "Model 2"))

# Check input data for both models
# head(DataValPred)
# head(DataValPredNoise)








##### Comparison with GAM for real-life data #####

data("NYdata")
s <- c(8, 11, 12, 14, 18, 21, 24, 28)
DataFit <- spT.subset(data = NYdata, var.name = c("s.index"), s = s, reverse = TRUE)
DataFit <- subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred <- spT.subset(data = NYdata, var.name = c("s.index"), s = s)
DataValPred <- subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 
  8)))

# GAM model
set.seed(11)
fit.gam <- gam(o8hrmax ~ s(cMAXTMP) + s(WDSP) + s(RH) + s(Longitude, Latitude, k = 10), 
  data = DataFit)
pred.gam <- predict(fit.gam, DataValPred, interval = "prediction")
spT.validation(DataValPred$o8hrmax, pred.gam)

# spTimer model
set.seed(11)
post.gp <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = DataFit, model = "GP", 
  coords = ~Longitude + Latitude, scale.transform = "SQRT", spatial.decay = spT.decay(distribution = Gamm(2, 
    1), tuning = 0.1))
set.seed(11)
pred.gp <- predict(post.gp, newdata = DataValPred, newcoords = ~Longitude + Latitude)
spT.validation(DataValPred$o8hrmax, c(pred.gp$Median))





##### Run demo file for NY example #####

## demo("nyExample") 

print(Sys.time() - start.time)
