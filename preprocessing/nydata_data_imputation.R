## Clean environment
rm(list=ls())

## Load libraries
library(openair)
library(mgcv)
library(spTimer)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Load pseudo-locations of the cheap sensors
#load(file="data/ny_ozone/NYdata.Rdata")
#load(file="data/ny_ozone/NYcheapLoc.Rdata")
data("NYdata")
NYdata$date <- as.POSIXct(strptime( sprintf("%04d-%02d-%02d",NYdata$Year,NYdata$Month,NYdata$Day),
                                    format="%Y-%m-%d", tz="GMT"))


## Deal with missing data ########################################################
## Which is the best method for filling small portions of missing data 

## Missing data by station
sum(is.na(NYdata$o8hrmax))
nasByStation <- aggregate(o8hrmax~s.index,data=NYdata,
                          FUN=function(x){sum(is.na(x))}, na.action = na.pass)
#View(nasByStation)
nasByStation <- nasByStation[order(nasByStation$o8hrmax,decreasing=T),]
nasByStation[nasByStation$o8hrmax>0,]



## Define train and validation datasets
trainDs <- na.omit(NYdata)
summaryPlot(trainDs[,c("o8hrmax","cMAXTMP","WDSP","RH","date")], period="months")



## Approach1: create a general model ########################################################

## Linear models
simpleLm <- lm(o8hrmax~cMAXTMP+WDSP+RH, data=trainDs)
simpleLm2 <- lm(o8hrmax~cMAXTMP+WDSP+RH+Day, data=trainDs)
simpleLm3 <- lm(o8hrmax~cMAXTMP+WDSP+RH+Longitude*Latitude+Month, data=trainDs)

printPlot(paper, "img/nydata/linear_model.jpeg", 5, 5, FUN=function(){
  plot(trainDs$o8hrmax, type="l", ylab="Ozone")
  lines(predict(simpleLm), col=2, cex=0.5)
  lines(predict(simpleLm2), col=3)
  lines(predict(simpleLm3), col=4)
})
# summary(simpleLm3)
# Some metrics
vlm1 <- spT.validation(trainDs$o8hrmax, predict(simpleLm))
vlm2 <- spT.validation(trainDs$o8hrmax, predict(simpleLm2))
vlm3 <- spT.validation(trainDs$o8hrmax, predict(simpleLm3))
r2lm1 <- cor(predict(simpleLm), trainDs$o8hrmax)^2
r2lm2 <- cor(predict(simpleLm2), trainDs$o8hrmax)^2
r2lm3 <- cor(predict(simpleLm3), trainDs$o8hrmax)^2
vlm <- rbind(vlm1,vlm2,vlm3)
vlm <- cbind( vlm, R2=c(r2lm1,r2lm2,r2lm3) )
vlm


## Simple plot to compare the linear models
pds <- NYdata[NYdata$s.index==7,] # prediction dataset
plot(pds$o8hrmax, type="l", ylab="Ozone", main="Ozone - Station 7", ylim=c(20,75))
lines(predict(simpleLm,newdata=pds), col=2, cex=0.5)
lines(predict(simpleLm2,newdata=pds), col=3, cex=0.5)
lines(predict(simpleLm3,newdata=pds), col=4, cex=0.5)



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

## Comparison (all data)
va <- spT.validation(trainDs$o8hrmax, predict(simpleLm3))
vb <- spT.validation(trainDs$o8hrmax, predict(simpleGam))
vc <- spT.validation(trainDs2$o8hrmax,pred.gp[,1])
vom <- rbind(va,vb,vc)
vom


## Comparison (Station 7)
testS7 <- NYdata[NYdata$s.index==7, ]
pred.gp.7 <- predict(simpleGp, newdata=testS7, newcoords=~Longitude+Latitude)
target <- testS7$o8hrmax[testS7$s.index==7]
va <- spT.validation(target, predict(simpleLm3, newdata=testS7))
vb <- spT.validation(target, predict(simpleGam, newdata=testS7))
vc <- spT.validation(target, pred.gp.7$Mean)
vom7 <- rbind(va,vb,vc)
vom7


## Grafical comparison in a station
spT.validation(testS7$o8hrmax, pred.gp.7$Mean)
plot(testS7$o8hrmax, type="l", main="Data imputation models for station 7", ylab="Ozone")
abline(v=which(is.na(testS7$o8hrmax)), col="gray", lty="dashed")
lines(predict(simpleLm3, newdata = testS7), col=2)
lines(predict(simpleGam, newdata = testS7), col=3)
lines(pred.gp.7$Mean, col=4)
#lines(pred.gp.7$Mean+pred.gp.7$SD, col=3, lty="dashed")
#lines(pred.gp.7$Mean-pred.gp.7$SD, col=3, lty="dashed")
legend("topright", legend=c("Obs.","GAM","GP"), 
       lty=rep(1,3), lwd=rep(1,3), col=c(1:3)) # gives the legend lines the correct color and width


## Observations:
## The best general model is GP, despite it does not use the data from the station with 
## missing values



## Approach 2: Create a model tailored for each station with missing values ###############

## Check R2 with neightbour stations
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

## Notes:
## Posible problem of adding several neighbors: multicollinearity

## GAM
# trainDsS7_gam <- NYdata[NYdata$s.index %in% c(7,8,9,10), ]
# simpleGamS7 <- gam(o8hrmax ~ s(cMAXTMP) + s(WDSP) + s(RH) + s(Longitude, Latitude, k = 10), 
#                  data = trainDsS7_gam)
## Error with GAM
## A term has fewer unique covariate combinations than specified maximum degrees of freedom


## GP
trainS7Gp <- NYdata[NYdata$s.index %in% c(8,9,10,11), ]
modelGpS7 <- spT.Gibbs(formula = o8hrmax ~ cMAXTMP + WDSP + RH, data = trainS7Gp, model = "GP", 
                      coords = ~Longitude + Latitude, scale.transform = "SQRT", 
                      spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
pred.gp.S7 <- fitted(modelGpS7)^2

## Metrics
t3ValLm <- spT.validation(trainDsS7$o8hrmax, predict(simpleLmS7, newdata = trainDsS7))
#cor(trainDsS7$o8hrmax, predict(simpleLmS7, newdata = trainDsS7), use="pairwise.complete.obs")
t3ValGp <- spT.validation(trainS7Gp$o8hrmax,pred.gp.S7[,1])
t3Val <- rbind(t3ValLm, t3ValGp)
t3Val


## Plot 
plot(trainDsS7$o8hrmax, type="l", main="Station 7", ylab="Ozone", xlim=c(0,70), ylim=c(25,80))
abline(v=which(is.na(testS7$o8hrmax)), col="gray", lty="dashed")
lines(predict(simpleLmS7, newdata = trainDsS7), col="orange", lwd=1)
lines(pred.gp[,1], col="cyan")
legend("topright", legend=c("Obs.","LM","GP"), 
       lty=rep(1,3), lwd=rep(1,3), col=c("black", "orange","cyan")) # gives the legend lines the correct color and width
#


# Observations:
# For some cases, the second approach leads to highest precision

