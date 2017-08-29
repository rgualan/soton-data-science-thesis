### R code from vignette source 'FBG13.Rnw'

###################################################
### code chunk number 1: mkData1
###################################################
rm(list=ls())

options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)

library(spBayes)
library(MBA)
library(fields)
library(xtable)
library(colorspace)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- cbind(1, rnorm(n))

B <- as.matrix(c(1,5))
p <- length(B)

sigma.sq <- 2
tau.sq <- 1
phi <- 6

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X%*%B + w, sqrt(tau.sq))

res <- 100
col <- rev(heat_hcl(33, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3)))

w.surf <- mba.surf(cbind(coords, w), no.X=res, no.Y=res, extend=TRUE)$xyz.est

par(mar=c(5,5,0.2,0.2), cex.lab=2, cex.axis=2)
image.plot(w.surf, xlab="Easting", ylab="Northing", xaxs="r", yaxs="r", col = col)
points(coords)


###################################################
### code chunk number 2: spLM1
###################################################
n.samples <- 5000

starting <- list("tau.sq"=1, "sigma.sq"=1, "phi"=6)

tuning <- list("tau.sq"=0.01, "sigma.sq"=0.01, "phi"=0.1)

priors <- list("beta.Flat", "tau.sq.IG"=c(2, 1),
               "sigma.sq.IG"=c(2, 1), "phi.Unif"=c(3, 30))

m.i <- spLM(y~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors, cov.model="exponential",
            n.samples=n.samples, n.report=2500)

burn.in <- floor(0.75*n.samples)

round(summary(window(m.i$p.theta.samples, 
                     start=burn.in))$quantiles[,c(3,1,5)],2)


###################################################
### code chunk number 3: spRecover1
###################################################
m.i <- spRecover(m.i, start=burn.in, thin=5, n.report=100)

round(summary(m.i$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)


###################################################
### code chunk number 4: wHatSurf1
###################################################
w.hat <- apply(m.i$p.w.recover.samples, 1, median)

w.hat.surf <- mba.surf(cbind(coords, w.hat), 
                       no.X=res, no.Y=res, extend=TRUE)$xyz.est

par(mar=c(5,5,0.2,0.2), cex.lab=2, cex.axis=2)
image.plot(w.hat.surf, xlab="Easting", ylab="Northing", xaxs="r", yaxs="r", col=col)


###################################################
### code chunk number 5: runTime1
###################################################
run.time <- round(m.i$run.time[3]/60,3)


###################################################
### code chunk number 6: ppData1
###################################################
n <- 3000
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X%*%B + w, sqrt(tau.sq))

##subset for prediction
hold.out <- 1:1000
coords.ho <- coords[hold.out,]
X.ho <- X[hold.out,]
w.ho <- w[hold.out]
y.ho <- y[hold.out]

coords <- coords[-hold.out,]
X <- X[-hold.out,]
w <- w[-hold.out]
y <- y[-hold.out]


###################################################
### code chunk number 7: ppModI
###################################################
m.i <- spLM(y~X-1, coords=coords, knots=c(5, 5, 0), starting=starting,
            tuning=tuning, priors=priors, cov.model="exponential",
            modified.pp=FALSE, n.samples=n.samples, n.report=2500)


###################################################
### code chunk number 8: ppModII-IV
###################################################
m.ii <- spLM(y~X-1, coords=coords, knots=c(5, 5, 0), starting=starting,
            tuning=tuning, priors=priors, cov.model="exponential",
            modified.pp=TRUE, n.samples=n.samples, verbose=FALSE)

m.iii <- spLM(y~X-1, coords=coords, knots=c(10, 10, 0), starting=starting,
            tuning=tuning, priors=priors, cov.model="exponential",
            modified.pp=FALSE, n.samples=n.samples, verbose=FALSE)

m.iv <- spLM(y~X-1, coords=coords, knots=c(10, 10, 0), starting=starting,
            tuning=tuning, priors=priors, cov.model="exponential",
            modified.pp=TRUE, n.samples=n.samples, verbose=FALSE)

quants <- function(x){
  quantile(x, prob=c(0.5, 0.025, 0.975))
}

ci.print <- function(x, digits=2){
  apply(round(apply(x, 2, quants), digits), 2, function(x){paste(x[1]," (",x[2],", ",x[3],")",  sep="")})
}

m.i <- spRecover(m.i, start=burn.in, thin=5, verbose=FALSE)
m.ii <- spRecover(m.ii, start=burn.in, thin=5, verbose=FALSE)
m.iii <- spRecover(m.iii, start=burn.in, thin=5, verbose=FALSE)
m.iv <- spRecover(m.iv, start=burn.in, thin=5, verbose=FALSE)

sub <- burn.in:n.samples

ests <- cbind(c(B[1], B[2], sigma.sq, tau.sq, phi), 
              c(ci.print(m.i$p.beta.samples[sub,]), ci.print(m.i$p.theta.samples[sub,])),
              c(ci.print(m.ii$p.beta.samples[sub,]), ci.print(m.ii$p.theta.samples[sub,])),
              c(ci.print(m.iii$p.beta.samples[sub,]), ci.print(m.iii$p.theta.samples[sub,])),
              c(ci.print(m.iv$p.beta.samples[sub,]), ci.print(m.iv$p.theta.samples[sub,])))

run.times <- c("", round(c(m.i$run.time[3]/60, m.ii$run.time[3]/60, m.iii$run.time[3]/60,  m.iv$run.time[3]/60), 2))

rel.run.times <- c("", round(c(m.i$run.time[3]/m.iv$run.time[3], m.ii$run.time[3]/m.iv$run.time[3], m.iii$run.time[3]/m.iv$run.time[3],  m.iv$run.time[3]/m.iv$run.time[3]), 2))

tab <- rbind(ests, run.times, rel.run.times)

rownames(tab) <- c("$\\beta_0$","$\\beta_1$","$\\sigma^2$", "$\\tau^2$", "$\\phi$", "Time", "Rel. time")

colnames(tab) <- c("true", "i", "ii", "iii", "iv")


###################################################
### code chunk number 9: tab1
###################################################
print(xtable(tab, 
             caption="Candidate predictive process models' parameter estimates, run-time (wall time) in minutes, and run-time relative to model $iv$. Parameter posterior summary 50 (2.5, 97.5) percentiles.", label="tab1", align="cccccc"),
      table.placement="!ht",  caption.placement="bottom", sanitize.text.function=function(x){x})


###################################################
### code chunk number 10: wSurf2
###################################################
w.hat.i <- apply(m.i$p.w.recover.samples, 1, median)
w.hat.ii <- apply(m.ii$p.w.recover.samples, 1, median)
w.hat.iii <- apply(m.iii$p.w.recover.samples, 1, median)
w.hat.iv <- apply(m.iv$p.w.recover.samples, 1, median)

w.surf <- mba.surf(cbind(coords, w), no.X=res, no.Y=res, extend=TRUE)$xyz.est

par(mar=c(5,5,0.2,0.2), cex.lab=2, cex.axis=2)
image.plot(w.surf, xlab="Easting", ylab="Northing", xaxs="r", yaxs="r", col=col)


###################################################
### code chunk number 11: wHatSurfI
###################################################
w.hat.i.surf <- mba.surf(cbind(coords, w.hat.i), no.X=res, no.Y=res, extend=TRUE)$xyz.est

par(mar=c(5,5,0.2,0.2), cex.lab=2, cex.axis=2)
image.plot(w.hat.i.surf, xlab="Easting", ylab="Northing", xaxs="r", yaxs="r", col=col)
points(m.ii$knot.coords, pch=19, cex=1)


###################################################
### code chunk number 12: wHatSurfII
###################################################
w.hat.ii.surf <- mba.surf(cbind(coords, w.hat.ii), no.X=res, no.Y=res, extend=TRUE)$xyz.est

par(mar=c(5,5,0.2,0.2), cex.lab=2, cex.axis=2)
image.plot(w.hat.ii.surf, xlab="Easting", ylab="Northing", xaxs="r", yaxs="r", col=col)
points(m.ii$knot.coords, pch=19, cex=1)


###################################################
### code chunk number 13: wHatSurfIII
###################################################
w.hat.iii.surf <- mba.surf(cbind(coords, w.hat.iii), no.X=res, no.Y=res, extend=TRUE)$xyz.est

par(mar=c(5,5,0.2,0.2), cex.lab=2, cex.axis=2)
image.plot(w.hat.iii.surf, xlab="Easting", ylab="Northing", xaxs="r", yaxs="r", col=col)
points(m.iii$knot.coords, pch=19, cex=1)


###################################################
### code chunk number 14: wHatSurfIV
###################################################
w.hat.iv.surf <- mba.surf(cbind(coords, w.hat.iv), no.X=res, no.Y=res, extend=TRUE)$xyz.est

par(mar=c(5,5,0.2,0.2), cex.lab=2, cex.axis=2)
image.plot(w.hat.iv.surf, xlab="Easting", ylab="Northing", xaxs="r", yaxs="r", col=col)
points(m.iv$knot.coords, pch=19, cex=1)


###################################################
### code chunk number 15: prediction
###################################################
m.iv.pred <- spPredict(m.iv, start=burn.in, thin=2, pred.covars=X.ho, 
                       pred.coords=coords.ho, verbose=FALSE)


###################################################
### code chunk number 16: yHatScatter
###################################################
y.hat <- apply(m.iv.pred$p.y.predictive.samples, 1, quants)
par(mar=c(5,5,5,5))
plot(y.ho, y.hat[1,], pch=19, cex=0.5, xlab="Observed y", ylab="Predicted y", 
     ylim=range(y.hat), xlim=range(y.hat), cex.lab=2, cex.axis=2)
arrows(y.ho, y.hat[1,], y.ho, y.hat[2,], angle=90, length=0.05)
arrows(y.ho, y.hat[1,], y.ho, y.hat[3,], angle=90, length=0.05)
lines(-20:20,-20:20, col="blue")


###################################################
### code chunk number 17: NYOzoneData
###################################################
library(maps)
data("NYOzone.dat")

N.t <- 62
n <- 28

##hold 31 days out of a stations 1, 5, and 10
hold.out <- c(1,5,10)

par(mar=c(0,0,0,0))
map(database="state",regions="new york")
points(NYOzone.dat[,c("Longitude","Latitude")], cex=2)
points(NYOzone.dat[hold.out,c("Longitude","Latitude")], pch=19, cex=2)

missing.obs <- sample(1:N.t, 31)
NYOzone.dat[,paste("O3.8HRMAX.",1:N.t,sep="")][hold.out,missing.obs] <- NA


###################################################
### code chunk number 18: NYOzoneModel
###################################################
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


