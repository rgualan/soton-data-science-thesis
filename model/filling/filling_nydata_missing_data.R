library(openair)
library(mgcv)
library(spTimer)


# Load pseudo-locations of the cheap sensors
load(file="data/ny_ozone/NYdata.Rdata")
load(file="data/ny_ozone/lcsCoords.Rdata")


# Experiment:
# Which is the best method for filling small portions of missing data


## Deal with missing data ########################################################
sum(is.na(NYdata$o8hrmax))
trainDs <- na.omit(NYdata)
summaryPlot(trainDs[,c("o8hrmax","cMAXTMP","WDSP","RH","date")], period="months")


## Approach1: create a general model ########################################################

## Linear models
simpleLm <- lm(o8hrmax~cMAXTMP+WDSP+RH, data=trainDs)
#View(data.frame(trainDs$o8hrmax, fitted(simpleLm), predict(simpleLm)))
plot(trainDs$o8hrmax, type="l")
lines(predict(simpleLm), col=2)
spT.validation(trainDs$o8hrmax, predict(simpleLm))
cor(predict(simpleLm), trainDs$o8hrmax)^2
simpleLm2 <- lm(o8hrmax~cMAXTMP+WDSP+RH+Day, data=trainDs)
spT.validation(trainDs$o8hrmax, predict(simpleLm2))
simpleLm3 <- lm(o8hrmax~cMAXTMP+WDSP+RH+Longitude*Latitude+Month, data=trainDs)
spT.validation(trainDs$o8hrmax, predict(simpleLm3))
#summary(simpleLm3)

## GAM model
simpleGam <- gam(o8hrmax ~ s(cMAXTMP) + s(WDSP) + s(RH) + s(Longitude, Latitude, k = 10), 
               data = trainDs)
#pred.gam <- predict(fit.gam, DataValPred, interval = "prediction")
spT.validation(trainDs$o8hrmax, predict(simpleGam))
cor(trainDs$o8hrmax, predict(simpleGam))^2
plot(trainDs$o8hrmax, type="l")
lines(predict(simpleGam), col=2)

## GP model
trainDs2 <- NYdata[!(NYdata$s.index %in% unique(NYdata$s.index[is.na(NYdata$o8hrmax)])), ]
simpleGp <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = trainDs2, model = "GP", 
                     coords = ~Longitude + Latitude, scale.transform = "SQRT", 
                     spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
#pred.gp <- predict(simpleGp, newdata=NYdata, newcoords = ~Longitude + Latitude)
pred.gp <- fitted(simpleGp)^2
spT.validation(trainDs2$o8hrmax,pred.gp[,1])
plot(trainDs2$o8hrmax, type="l")
lines(pred.gp[,1], col=2)

# Comparison
a <- spT.validation(trainDs$o8hrmax, predict(simpleLm3))
b <- spT.validation(trainDs$o8hrmax, predict(simpleGam))
c <- spT.validation(trainDs2$o8hrmax,pred.gp[,1])
a; b; c;

# Grafical comparison in a station
testS7 <- NYdata[NYdata$s.index==7, ]
pred.gp.7 <- predict(simpleGp, newdata=testS7, newcoords = ~Longitude + Latitude)
spT.validation(testS7$o8hrmax, pred.gp.7$Mean)
plot(testS7$o8hrmax, type="l")
abline(v=which(is.na(testS7$o8hrmax)), col="gray", lty="dashed")

lines(pred.gp.7$Mean, col=2)
#lines(pred.gp.7$Mean+pred.gp.7$SD, col=3, lty="dashed")
#lines(pred.gp.7$Mean-pred.gp.7$SD, col=3, lty="dashed")
lines(predict(simpleLm3, newdata = testS7), col=3)
lines(predict(simpleGam, newdata = testS7), col=4)

# Observations:
# The best general model is GP, despite it does not use the data from the station with 
# missing values


# Approach 2: Create a model tailored for each station with missing values ###############

# Check R2 with neightbour stations
cor(NYdata$o8hrmax[NYdata$s.index==7],NYdata$o8hrmax[NYdata$s.index==8],use="pairwise.complete.obs")^2
cor(NYdata$o8hrmax[NYdata$s.index==7],NYdata$o8hrmax[NYdata$s.index==11],use="pairwise.complete.obs")^2
cor(NYdata$o8hrmax[NYdata$s.index==7],NYdata$o8hrmax[NYdata$s.index==10],use="pairwise.complete.obs")^2
cor(NYdata$o8hrmax[NYdata$s.index==7],NYdata$o8hrmax[NYdata$s.index==9],use="pairwise.complete.obs")^2
# plot(NYdata$o8hrmax[NYdata$s.index==7], type="l")
# lines(NYdata$o8hrmax[NYdata$s.index==8], col=2)
# lines(NYdata$o8hrmax[NYdata$s.index==9], col=3)
# lines(NYdata$o8hrmax[NYdata$s.index==10], col=4)
# lines(NYdata$o8hrmax[NYdata$s.index==11], col=4)
# There is a clear correlation between close stations

# Check next case
# cor(NYdata$o8hrmax[NYdata$s.index==12],NYdata$o8hrmax[NYdata$s.index==14],use="pairwise.complete.obs")
# cor(NYdata$o8hrmax[NYdata$s.index==12],NYdata$o8hrmax[NYdata$s.index==26],use="pairwise.complete.obs")
# cor(NYdata$o8hrmax[NYdata$s.index==12],NYdata$o8hrmax[NYdata$s.index==24],use="pairwise.complete.obs")
# cor(NYdata$o8hrmax[NYdata$s.index==12],NYdata$o8hrmax[NYdata$s.index==4],use="pairwise.complete.obs")
# cor(NYdata$o8hrmax[NYdata$s.index==12],NYdata$o8hrmax[NYdata$s.index==21],use="pairwise.complete.obs")
# plot(NYdata$o8hrmax[NYdata$s.index==12], type="l")
# lines(NYdata$o8hrmax[NYdata$s.index==14], col=2)
# lines(NYdata$o8hrmax[NYdata$s.index==26], col=3)
# lines(NYdata$o8hrmax[NYdata$s.index==24], col=4)
# lines(NYdata$o8hrmax[NYdata$s.index==4], col=5)


# Linear model
trainDsS7 <- NYdata[NYdata$s.index==7, ]
trainDsS7$ohrmaxS8 <- NYdata$o8hrmax[NYdata$s.index==8]
simpleLmS7 <- lm(o8hrmax~cMAXTMP+WDSP+RH+Month+ohrmaxS8, data = trainDsS7)
spT.validation(trainDsS7$o8hrmax, predict(simpleLmS7, newdata = trainDsS7))
cor(trainDsS7$o8hrmax, predict(simpleLmS7, newdata = trainDsS7), use="pairwise.complete.obs")
#plot(trainDsS7$o8hrmax, type="l")
lines(predict(simpleLmS7, newdata = trainDsS7), col="orange", lwd=2)

# # Linear model2
# trainDsS7$ohrmaxS9 <- NYdata$o8hrmax[NYdata$s.index==9]
# trainDsS7$ohrmaxS10 <- NYdata$o8hrmax[NYdata$s.index==10]
# trainDsS7$ohrmaxS11 <- NYdata$o8hrmax[NYdata$s.index==11]
# simpleLmS7b <- lm(o8hrmax~cMAXTMP+WDSP+RH+Month+Day+ohrmaxS8+ohrmaxS9+ohrmaxS10+ohrmaxS11, 
#                   data = trainDsS7)
# spT.validation(trainDsS7$o8hrmax, predict(simpleLmS7b, newdata = trainDsS7))
# cor(trainDsS7$o8hrmax, predict(simpleLmS7b, newdata = trainDsS7), use="pairwise.complete.obs")
# #plot(trainDsS7$o8hrmax, type="l")
# lines(predict(simpleLmS7b, newdata = trainDsS7), col="pink", lwd=2)
# 
# # Linear model 3
# #trainDsS7$ohrmaxSo <- rowMeans(cbind(trainDsS7$ohrmaxS8,trainDsS7$ohrmaxS9,trainDsS7$ohrmaxS10,trainDsS7$ohrmaxS11))
# trainDsS7$ohrmaxSo <- 
#   trainDsS7$ohrmaxS8*.7+trainDsS7$ohrmaxS9*0.1+trainDsS7$ohrmaxS10*0.1+trainDsS7$ohrmaxS11*0.1
# simpleLmS7c <- lm(o8hrmax~cMAXTMP+WDSP+RH+Month+Day+ohrmaxSo, 
#                   data = trainDsS7)
# spT.validation(trainDsS7$o8hrmax, predict(simpleLmS7c, newdata = trainDsS7))
# cor(trainDsS7$o8hrmax, predict(simpleLmS7c, newdata = trainDsS7), use="pairwise.complete.obs")
# #plot(trainDsS7$o8hrmax, type="l")
# lines(predict(simpleLmS7c, newdata = trainDsS7), col="purple", lwd=2)

# GP
trainS7Gp <- NYdata[NYdata$s.index %in% c(8,9,10,11), ]
modelGpS7 <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = trainS7Gp, model = "GP", 
                      coords = ~Longitude + Latitude, scale.transform = "SQRT", 
                      spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
pred.gp.S7 <- fitted(modelGpS7)^2
spT.validation(trainS7Gp$o8hrmax,pred.gp.S7[,1])
#plot(trainDs2$o8hrmax, type="l")
lines(pred.gp[,1], col=5)


# Observations:
# For some cases, the second approach leads to highest precision

