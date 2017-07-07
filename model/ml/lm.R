# These packages will be required to run the code in this file
library("openair")


#### Real life Example ####
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

# Simple diagnostic plots
NYdata$date <- as.POSIXct(strptime( sprintf("%04d-%02d-%02d",NYdata$Year,NYdata$Month,NYdata$Day),
                                    format="%Y-%m-%d", tz="GMT"))
summaryPlot(NYdata[,c(-1:-6)], period = "months")


#### Linear model ####
# Fit linear model
lmModel <- lm(o8hrmax~cMAXTMP+WDSP+RH, data=DataFit)
summary(lmModel)
plot(lmModel)

yh <- predict(lmModel, DataValPred)
plot(DataValPred$o8hrmax, yh)

spT.validation(DataValPred$o8hrmax, yh)
#dev.off()

# Figures 8 (a) -- (d)
dev.new()
data("NYgrid")
set.seed(11)
post.gp2 <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = NYdata, model = "GP", 
  coords = ~Longitude + Latitude, scale.transform = "SQRT", nItr = 15000, nBurn = 0, 
  spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
par(mfrow = c(2, 2))
plot(post.gp2$betap[1, ], type = "l", main = "Intercept", xlab = "Iterations", ylab = "", 
  ylim = c(-1, 6))
plot(post.gp2$betap[2, ], type = "l", main = "cMAXTMP", xlab = "Iterations", ylab = "", 
  ylim = c(0.06, 0.25))
plot(post.gp2$betap[3, ], type = "l", main = "WDSP", xlab = "Iterations", ylab = "", 
  ylim = c(-0.01, 0.25))
plot(post.gp2$betap[4, ], type = "l", main = "RH", xlab = "Iterations", ylab = "", 
  ylim = c(-0.4, 0.33))
par(mfrow = c(1, 1))
dev.off()

# Figures 9 (a) and (b) Predict on grids
set.seed(11)
post.gp2 <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = NYdata, model = "GP", 
  coords = ~Longitude + Latitude, scale.transform = "SQRT", spatial.decay = spT.decay(distribution = Gamm(2, 
    1), tuning = 0.1))
set.seed(11)
grid.pred <- predict(post.gp2, newdata = NYgrid, newcoords = ~Longitude + Latitude)

# predictive plots

# this function is used to delete values outside NY
fnc.delete.map.XYZ <- function(xyz) {
  x <- xyz$x
  y <- xyz$y
  z <- xyz$z
  xy <- expand.grid(x, y)
  eus <- (map.where(database = "state", x = xy[, 1], y = xy[, 2]))
  dummy <- rep(0, length(xy[, 1]))
  eastUS <- NULL
  eastUS <- data.frame(lon = xy[, 1], lat = xy[, 2], state = eus, dummy = dummy)
  eastUS[!is.na(eastUS[, 3]), 4] <- 1
  eastUS[eastUS[, 3] == "pennsylvania" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "new jersey" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "connecticut" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "massachusetts:main" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "new hampshire" & !is.na(eastUS[, 3]), 4] <- 0
  eastUS[eastUS[, 3] == "vermont" & !is.na(eastUS[, 3]), 4] <- 0
  a <- eastUS[, 4]
  z <- as.vector(xyz$z)
  z[!a] <- NA
  z <- matrix(z, nrow = length(xyz$x))
  xyz$z <- z
  xyz
}
## 

coords <- unique(NYdata[, c("Longitude", "Latitude")])
grid.coords <- unique(NYgrid[, c("Longitude", "Latitude")])
true.val <- matrix(NYdata$o8hrmax, 62, 28)
grid.val <- matrix(grid.pred$Median, 62, dim(grid.coords)[[1]])
grid.sd <- matrix(grid.pred$SD, 62, dim(grid.coords)[[1]])

surfplot <- function(day = 60, val, ...) {
  z <- val
  surf <- cbind(grid.coords, z[day, ])
  surf <- mba.surf(surf, 200, 200)$xyz
  surf <- fnc.delete.map.XYZ(xyz = surf)
  # map(database='state',regions='new york')
  image.plot(surf, xlab = "Longitude", ylab = "Latitude", axes = F, ...)
  contour(surf, nlevels = 10, lty = 3, add = TRUE)
  map(database = "state", regions = "new york", add = TRUE)
  axis(1)
  axis(2)
}

# Section 5: code for Figure 8(a) prediction for day 60
day <- 60
surfplot(day, val = grid.val, col = rainbow_hcl(100, start = 200, end = 0))
text(coords, labels = round(true.val[day, ], 1), cex = 0.8, col = 1)

# Section 5: code for Figure 8(b) sd for day 60
day <- 60
surfplot(day, val = grid.sd, col = diverge_hcl(100, h = c(246, 40), c = 96, l = c(65, 
  90)))
points(coords, pch = 19, cex = 1, col = 2)
points(coords, pch = 1, cex = 1, col = 1)



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
