##################################
## Using the SPDE approach and the INLA algorithm for spatio-temporal
## modelling and mapping
##################################

## Clean environment
rm(list=ls())

## User-overrideable configuration:
## Set any of these variables to alternate values before sourcing this file.

## To force a rerun of inla, set force.rerun=TRUE
## If "result" does not exist, inla will be rerun regardless of this setting.
if (!exists("force.rerun")) { force.rerun=FALSE }
## For which day do you want the prediction?
if (!exists("i_day")) { i_day = 89 } # The last day of the training period 
## Do you want to print (to png/pdf files)?
if (!exists("do.print")) { do.print = FALSE }
## Path for the Covariates directory
if (!exists("fr.path")) { fr.path = "data/FR_AQeReporting_2013-2015/" }

library(INLA)
library(fields)
library(abind)
library(lattice)

graphics.off()
#if (force.rerun || !exists("result")) {

##--- Source the function for selecting the covariates for a given day
source(paste0(fr.path,"Covariates/Covariates_selector.R"))

## ################################
## Load the data
## ################################
##--- for the 24 stations and 182 days
Piemonte_data <-read.csv(paste0(fr.path,"France_data_byday.csv"),header=TRUE)
coords <- read.csv(paste0(fr.path,"coordinates.csv"),header=TRUE)
coords <- coords[coords$Station.ID %in% unique(Piemonte_data$Station.ID),]

##--- for the 10 validation stations and 182 days
Piemonte_data_validation <-read.csv(
  paste0(fr.path,"France_data_byday_val.csv"),header=TRUE)
coords_validation <-read.csv(paste0(fr.path,"coordinates_val.csv"),header=TRUE)
coords_validation <- coords_validation[coords_validation$Station.ID %in% unique(Piemonte_data_validation$Station.ID),]

rownames(coords) = coords[,"Station.ID"]
rownames(coords_validation) = coords_validation[,"Station.ID"]

which_date = unique(Piemonte_data$date)[i_day]
print(paste("**---- You will get a prediction for", which_date, "---**"))

## ################################
## Work out how many stations and days there are
## ################################
n_stations <- length(coords$Station.ID)
n_stations_val <- length(coords_validation$Station.ID)
n_data <- length(Piemonte_data$Station.ID)
n_days <- as.integer(n_data/n_stations)

##--- Check that the data is OK
if (n_data %% n_stations != 0) {
  print("The number of data points needs to be an integer multiple of the number of stations!")
  #return
  exit
}

##--- Check the the validation data is OK
if (length(Piemonte_data_validation$Station.ID) != n_stations_val*n_days) {
  print("Something is wrong with the dimension of the validation data!")
  #return
  exit
}

## ################################
##Standardize covariates and calculate LOG of PM10
## ################################
##--- The covariates are standardised using the mean and std.dev. for
##--- the 24 data sites.
mean_covariates = apply(Piemonte_data[,3:10],2,mean)
sd_covariates = apply(Piemonte_data[,3:10],2,sd)

Piemonte_data[,3:10] =
  scale(Piemonte_data[,3:10],
        mean_covariates, sd_covariates)
Piemonte_data_validation[,3:10] =
  scale(Piemonte_data_validation[,3:10],
        mean_covariates, sd_covariates)

Piemonte_data$logPM10 = log(Piemonte_data$PM10)
Piemonte_data_validation$logPM10 = log(Piemonte_data_validation$PM10)

Piemonte_data$time = rep(1:n_days,each=n_stations)
Piemonte_data_validation$time = rep(1:n_days,each = n_stations_val)

## ################################
## Estimation
## ################################
map('worldHires', c('France'), xlim=c(-4.7,8.1), ylim=c(42.4,51), mar=rep(1,4))	
points(coords$Longitude, coords$Latitude, pch=16, col="green")
points(coords_validation$Longitude, coords_validation$Latitude)
points(coords_validation$Longitude, coords_validation$Latitude, pch=16, col="cyan")
legend("topleft", c("Training", "Validation"), pch=c(16,16), col=c("green","cyan"))

## ################################
## Triangulation using borders
## ################################
mesh =
  inla.mesh.2d(loc=cbind(coords$UTMX,
                         coords$UTMY),
               #loc.domain=borders, # TODO: Create borders???
               #offset=c(10, 140),
               #max.edge=c(50, 1000),
               #min.angle=c(26, 21),
               cutoff=0#,
               #plot.delay=10 #NULL
               ,max.n=-1
  )

##--- Plot the triangulation
if (do.print)
  pdf("triangulation_piemonte.pdf")
plot(mesh)
#lines(borders, lwd=3)
#points(coords2$UTMX, coords2$UTMY, pch=20, cex=2, col=2)
points(coords$UTMX, coords$UTMY, pch=20, cex=2, col=2)
points(coords_validation$UTMX, coords_validation$UTMY, pch=20, cex=2, col=4)
if (do.print) { dev.off() } else { inla.dev.new() }





## ################################
## Make the SPDE object and the formula
## ################################

##--- Construct the SPDE object
spde = inla.spde2.matern(mesh=mesh, alpha=2)

##--- Observation structure for estimation data
A.est =
  inla.spde.make.A(mesh,
                   loc=
                     as.matrix(coords[Piemonte_data$Station.ID,
                                      c("UTMX","UTMY")]),
                   group=Piemonte_data$time,
                   n.group=n_days
  )
##--- Observation structure for validation data
A.val =
  inla.spde.make.A(mesh,
                   loc=
                     as.matrix(coords_validation[Piemonte_data_validation$Station.ID,
                                                 c("UTMX","UTMY")]),
                   group=Piemonte_data_validation$time,
                   n.group=n_days
  )
##--- Observation structure for field prediction
A.pred =
  inla.spde.make.A(mesh, group=i_day, n.group=n_days)

field.indices =
  inla.spde.make.index("field",
                       n.spde=spde$n.spde,
                       n.group=n_days)
stack.est =
  inla.stack(data=list(logPM10=Piemonte_data$logPM10),
             A=list(A.est, 1),
             effects=
               list(c(field.indices,
                      list(Intercept=1)),
                    list(Piemonte_data[,c(3,6:10)])),
             tag="est")
stack.val =
  inla.stack(data=list(logPM10=NA),
             A=list(A.val, 1),
             effects=
               list(c(field.indices,
                      list(Intercept=1)),
                    list(Piemonte_data_validation[,c(3,6:10)])),
             tag="val")
scaled.mesh.loc =
  list(UTMX=(rep(scale(mesh$loc[,1],
                       mean_covariates["UTMX"],
                       sd_covariates["UTMX"]),
                 n_days)),
       UTMY=(rep(scale(mesh$loc[,2],
                       mean_covariates["UTMY"],
                       sd_covariates["UTMY"]),
                 n_days)))
stack.pred =
  inla.stack(data=list(logPM10=NA),
             A=list(A.pred),
             effects=
               list(c(field.indices,
                      scaled.mesh.loc,
                      list(Intercept=1)
               )),
             tag="pred")
stack = inla.stack(stack.est, stack.val, stack.pred)

formula <- (logPM10 ~ -1 + Intercept + Altitude + UTMX + UTMY + TMP + RAIN + PRESS + f(field, model=spde, group=field.group, control.group=list(model="ar1")))


## ################################
## Call INLA and get results
## ################################
## Attention!!! Set do.remote=FALSE to disable inla.remote
## reordering="metis" is needed to prevent crashes on some systems.

st <- Sys.time()
result =
  inla(formula,
       data=inla.stack.data(stack, spde=spde),
       family="gaussian",
       control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
       control.compute=list(cpo=FALSE),
       control.inla = list(reordering = "metis"),
       keep=FALSE, verbose=TRUE)
print(summary(result))
#First try. Started: 17:56. Ended: ~19:15. Result: Failed 
print(Sys.time()-st)

print(summary(result))
#force.rerun = FALSE


## ############################
## Extract results
## ############################

##--- Results for fixed effects - covariate coeffs -
beta = result$summary.fixed[,"mean"]
beta_sd = result$summary.fixed[,"sd"]

## Calculate validation data information
validation0 =
  list(p = rep(NA, length(Piemonte_data$logPM10)))
index = inla.stack.index(stack,"est")$data
tmp.mean = result$summary.linear.predictor[index,"mean"]
tmp.sd = result$summary.linear.predictor[index,"sd"]
validation0$res = Piemonte_data$logPM10 - tmp.mean
validation0$res.std = validation0$res /
  sqrt(tmp.sd^2 + 1/result$summary.hyperpar[1,"mean"])
validation0$p = pnorm(validation0$res.std)

validation = list()
index = inla.stack.index(stack,"val")$data
tmp.mean = result$summary.linear.predictor[index,"mean"]
tmp.sd = result$summary.linear.predictor[index,"sd"]
validation$res = Piemonte_data_validation$logPM10 - tmp.mean
validation$res.std = validation$res /
  sqrt(tmp.sd^2 + 1/result$summary.hyperpar[1,"mean"])
validation$p = pnorm(validation$res.std)
validation$rmse = sqrt(mean(validation$res^2, na.rm=TRUE))
validation$cor = cor(Piemonte_data_validation$logPM10, tmp.mean,
                     use="pairwise.complete.obs", method="pearson")
validation$cover = mean((validation$p>0.025)&(validation$p<0.975), na.rm=TRUE)

## Prediction interval coverage:

#Import predictions from Cameletti et al 2011 (using MCMC)
# predMCMC = read.table("pred_MCMC.csv",sep=";", header=T)
# predMCMC = predMCMC$x
# 
# mcmc = list()
# mcmc$rmse = sqrt(mean((predMCMC-Piemonte_data_validation$logPM10)^2,na.rm=TRUE))
# mcmc$cor = cor(Piemonte_data_validation$logPM10 , predMCMC, use="pairwise.complete.obs",method="pearson")

#plots
par(mfrow=c(3,1))
plot(Piemonte_data_validation$logPM10,t="l",ylim=c(1,5.5))
#plot(predMCMC, col=2,t="l",ylim=c(1,5.5))
plot(tmp.mean, col=3,t="l",ylim=c(1,5.5))
#inla.dev.new()


#RMSE
sqrt(mean((validation$res)^2,na.rm=T)) #spde
#sqrt(mean((predMCMC-Piemonte_data_validation$logPM10)^2,na.rm=T)) #mcmc
#CORR
cor(Piemonte_data_validation$logPM10 , tmp.mean, use="pairwise.complete.obs",method="pearson")
#cor(Piemonte_data_validation$logPM10 , predMCMC, use="pairwise.complete.obs",method="pearson")


##--- Extract SPDE results
result.field = inla.spde.result(result, "field", spde, do.transform=TRUE)


field_mean = matrix(result.field$summary.values$mean, mesh$n, n_days)
field_sd   = matrix(result.field$summary.values$sd, mesh$n, n_days)
field_pred_mean =
  result$summary.linear.predictor[inla.stack.index(stack,"pred")$data, "mean"]
field_pred_sd =
  result$summary.linear.predictor[inla.stack.index(stack,"pred")$data, "sd"]

##--- ar1 parameter
result$summary.hyperpar["GroupRho for field",]

##--- sigma2eps (1/precision)
sigma2eps_marg =
  inla.tmarginal(function(x) 1/x,
                 result$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2eps_m1 = inla.emarginal(function(x) x, sigma2eps_marg)
sigma2eps_m2 = inla.emarginal(function(x) x^2, sigma2eps_marg)
sigma2eps_stdev = sqrt(sigma2eps_m2 - sigma2eps_m1^2)
sigma2eps_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2eps_marg)

var.nom.marg = result.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)

range.nom.marg = result.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)

approx.hyperpar =
  rbind(obs.var=c(sigma2eps_m1, sigma2eps_stdev, sigma2eps_quantiles),
        spde.var.nom=c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
        spde.range.nom=c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
        AR.rho=result$summary.hyperpar["GroupRho for field",1:5])
print(approx.hyperpar)



# Notes:
# The results were terribly inprecise for the initial experiment



