## Clean environment
rm(list=ls())

## Load libraries
library(spBayes)
source("util/my_helper.R")

## Read data #######################################################################
#epa <- readRDS("data/epa/epa_daily/2016/california_ozone.RDS")
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]
sites <- getSites(epa)
## Test
epa <- addJfield(epa)
epa$wn <- weekdays(epa$Date)
## Standardize variable 
epa$sOzone <- scale(epa$Ozone)

## Split data
folds <- readRDS("data/tmp/folds.RDS")
## Fold(1)
epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=1],] 
epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==1],]
epa.test$sOzoneH <- NA

## Fit GP model ######################################################################
model <- spLM(sOzone~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)
m.2 <- spLM(y~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.2, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)
burn.in <- 0.5*n.samples
##recover beta and spatial random effects
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)
m.2 <- spRecover(m.2, start=burn.in, verbose=FALSE)
round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.2$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.2$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
m.1.w.summary <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]
m.2.w.summary <- summary(mcmc(t(m.2$p.w.recover.samples)))$quantiles[,c(3,1,5)]
plot(w, m.1.w.summary[,1], xlab="Observed w", ylab="Fitted w",
     xlim=range(w), ylim=range(m.1.w.summary), main="Spatial random effects")
arrows(w, m.1.w.summary[,1], w, m.1.w.summary[,2], length=0.02, angle=90)
arrows(w, m.1.w.summary[,1], w, m.1.w.summary[,3], length=0.02, angle=90)
lines(range(w), range(w))
points(w, m.2.w.summary[,1], col="blue", pch=19, cex=0.5)
arrows(w, m.2.w.summary[,1], w, col="blue", m.2.w.summary[,2], length=0.02, angle=90)
arrows(w, m.2.w.summary[,1], w, col="blue", m.2.w.summary[,3], length=0.02, angle=90)
###########################
##Predictive process model
###########################
m.1 <- spLM(y~X-1, coords=coords, knots=c(6,6,0.1), starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)
m.2 <- spLM(y~X-1, coords=coords, knots=c(6,6,0.1), starting=starting,
            tuning=tuning, priors=priors.2, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)
burn.in <- 0.5*n.samples
round(summary(window(m.1$p.beta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.2$p.beta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.1$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.2$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)

