rm(list=ls())
library(spBayes)
library(MBA)
library(fields)
library(xtable)
library(colorspace)
library(maps)

###################################################
### code chunk number 17: NYOzoneData
###################################################
data("NYOzone.dat")

N.t <- 62  ## Number of days (cols)
n <- 28  ## Number of stations (rows)

##hold 31 days out of a stations 1, 5, and 10
hold.out <- c(1,5,10)

par(mar=c(0,0,0,0))
map(database="state",regions="new york")
points(NYOzone.dat[,c("Longitude","Latitude")], cex=2)
points(NYOzone.dat[hold.out,c("Longitude","Latitude")], pch=19, cex=2)

## Inject missing data in the hold.out stations at 31 random days
missing.obs <- sample(1:N.t, 31)
NYOzone.dat[,paste("O3.8HRMAX.",1:N.t,sep="")][hold.out,missing.obs] <- NA


###################################################
### code chunk number 18: NYOzoneModel
###################################################
## A formula for each day???
## Eg: O3.8HRMAX.62 ~ cMAXTMP.62 + WDSP.62 + RH.62
mods <- lapply(paste("O3.8HRMAX.",1:N.t, "~cMAXTMP.",1:N.t, 
                     "+WDSP.",1:N.t, "+RH.",1:N.t, sep=""), as.formula)

p <- 4 ##number of predictors

coords <- NYOzone.dat[,c("X.UTM","Y.UTM")]/1000
max.d <- max(iDist(coords))

starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
                 "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                 "sigma.eta"=diag(rep(0.01, p)))

tuning <- list("phi"=rep(2, N.t)) 

priors <- list("beta.0.Norm"=list(rep(0,p), diag(100000,p)),
               "phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
               "sigma.sq.IG"=list(rep(2,N.t), rep(25,N.t)),
               "tau.sq.IG"=list(rep(2,N.t), rep(25,N.t)),
               "sigma.eta.IW"=list(2, diag(0.001,p)))

n.samples <- 5000

m.i <- spDynLM(mods, data=NYOzone.dat, coords=as.matrix(coords), 
               starting=starting, tuning=tuning, priors=priors, get.fitted=TRUE,
               cov.model="exponential", n.samples=n.samples, n.report=2500)


###################################################
### code chunk number 19: plotFunction
###################################################
ts.params.plot <- function(x, xlab="", ylab="", cex.lab=2.5, cex.axis=2.5, prob=c(0.5, 0.025, 0.975),
                           mar=c(5,5.5,0.1,1), zero.line=TRUE){
  N.t <- ncol(x)
  q <- apply(x, 2, function(x){quantile(x, prob=prob)})
  
  par(mar=mar)##bottom, left, top and right
  plot(1:N.t, q[1,], pch=19, cex=0.5, xlab=parse(text=xlab), ylab=parse(text=ylab), ylim=range(q), cex.lab=cex.lab, cex.axis=cex.axis)
  arrows(1:N.t, q[1,], 1:N.t, q[3,], length=0.02, angle=90)
  arrows(1:N.t, q[1,], 1:N.t, q[2,], length=0.02, angle=90)
  if(zero.line){abline(h=0, col="blue")}
}


###################################################
### code chunk number 20: NYOzoneBeta
###################################################
burn.in <- floor(0.75*n.samples)
beta <- m.i$p.beta.samples[burn.in:n.samples,]

par(mfrow=c(4,1))
beta.0 <- beta[,grep("Intercept", colnames(beta))]           
ts.params.plot(beta.0, ylab="beta[0]")

beta.1 <- beta[,grep("cMAXTMP", colnames(beta))]
ts.params.plot(beta.1, ylab="beta[cMAXTMP]")

beta.2 <- beta[,grep("WDSP", colnames(beta))]
ts.params.plot(beta.2, ylab="beta[WDSP]")

beta.3 <- beta[,grep("RH", colnames(beta))]
ts.params.plot(beta.3, ylab="beta[RH]", xlab="Days")


###################################################
### code chunk number 21: NYOzoneTheta
###################################################
theta <- m.i$p.theta.samples[burn.in:n.samples,]

par(mfrow=c(3,1))
sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
ts.params.plot(sigma.sq, ylab="sigma^2", zero.line=FALSE)

tau.sq <- theta[,grep("tau.sq", colnames(theta))]
ts.params.plot(tau.sq, ylab="tau^2", zero.line=FALSE)

phi <- theta[,grep("phi", colnames(theta))]
ts.params.plot(3/phi, ylab="3/phi (km)", xlab="Days", zero.line=FALSE)


###################################################
### code chunk number 22: NYOzonePred
###################################################
y.hat <- apply(m.i$p.y.samples[,burn.in:n.samples], 1, quants)

data(NYOzone.dat)

obs.O3 <- NYOzone.dat[,paste("O3.8HRMAX.",1:N.t,sep="")]

ylim <- c(15,75)
ylab <- "O3.8HRMAX"

par(mfrow=c(3,1),mar=c(5,5.5,0.1,1))
station.1 <- y.hat[,seq(hold.out[1],ncol(y.hat),n)]

plot(1:N.t, obs.O3[hold.out[1],], ylim=ylim, ylab=ylab, cex.lab=2.5, cex.axis=2.5, xlab="")
lines(1:N.t, station.1[1,1:N.t])
lines(1:N.t, station.1[2,1:N.t], col="blue", lty=3)
lines(1:N.t, station.1[3,1:N.t], col="blue", lty=3)
points(missing.obs, obs.O3[hold.out[1],missing.obs], pch=19)

station.2 <- y.hat[,seq(hold.out[2],ncol(y.hat),n)]

plot(1:N.t, obs.O3[hold.out[2],], ylim=ylim, ylab=ylab, cex.lab=2.5, cex.axis=2.5, xlab="")
lines(1:N.t, station.2[1,1:N.t])
lines(1:N.t, station.2[2,1:N.t], col="blue", lty=3)
lines(1:N.t, station.2[3,1:N.t], col="blue", lty=3)
points(missing.obs, obs.O3[hold.out[2],missing.obs], pch=19)

station.3 <- y.hat[,seq(hold.out[3],ncol(y.hat),n)]

plot(1:N.t, obs.O3[hold.out[3],], ylim=ylim, ylab=ylab, cex.lab=2.5, cex.axis=2.5, xlab="Days")
lines(1:N.t, station.3[1,1:N.t])
lines(1:N.t, station.3[2,1:N.t], col="blue", lty=3)
lines(1:N.t, station.3[3,1:N.t], col="blue", lty=3)
points(missing.obs, obs.O3[hold.out[3],missing.obs], pch=19)


