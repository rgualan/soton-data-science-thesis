## Clean the environment
rm(list=ls())

## Load libraries
source("util/my_helper.R")

## Global variables ###############################################################
paper <- setupPaper()

## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- addDateDerivedFeatures(epa)
epa <- transformFeatures(epa)
epa$sOzone <- scale(epa$Ozone)[,1]
sites <- getSites(epa)

## Outlier detection ####################################################
boxplot(epa$Rain, main="Precipitation", boxwex=0.1)
outliers <- boxplot.stats(epa$Rain)$out  # outlier values
## To many outliers. It is going to be too time consuming.

## Transformation #######################################################
## Many procedures that assume normality of residuals are robust to modest 
## violations of normality of residuals
## https://stats.stackexchange.com/questions/1601/what-other-normalizing-transformations-are-commonly-used-beyond-the-common-ones

## Wind
hist(epa$Wind)
hist(sqrt(epa$Wind))
epa$sqrtWind <- sqrt(epa$Wind)
## Rain
hist(epa$Rain)
# library(LambertW)
# epa$Rain2 <- Gaussianize(epa$Rain)
# hist(epa$Rain2)
plot(epa$Rain,type="l")
epa$logRain <- log(1+epa$Rain, 10) # The additive term is to avoid -inf
epa$levelRain <- round(epa$logRain) # 6 levels derived from the log trasnformation
epa$rained <- epa$Rain>0 # Binary variable
par(mfrow=c(3,1))
hist(epa$Rain)
hist(epa$logRain)
hist(epa$levelRain, breaks=0:5)
par(mfrow=c(1,1))
## Rainy days by month?
## There should be a trend
tmp <- epa
tmp$Month <- as.numeric(format(tmp$Date, "%m"))
rainyDaysByMonth <- aggregate(Rain~Month, tmp, FUN=sum)
plot(Rain~Month,rainyDaysByMonth)
## Less rain in summer as expected

## Dew Point
hist(epa$Dew.Point)
plot(Dew.Point~Temperature, epa)
cor(epa$Temperature,epa$Dew.Point)
## Water evaporation
hist(epa$Water.Evap)
## Heat flux
hist(epa$Heat.Flux) 
plot(Heat.Flux~Temperature, epa)
cor(epa$Heat.Flux, epa$Temperature) # Negative correlated to Temperature
## Geop.Height
hist(epa$Geop.Height)
hist(log(epa$Geop.Height)) ## Log!
## Geop.Height.Tropo
hist(epa$Geop.Height.Tropo)
## Lat.Heat.Flux
hist(epa$Lat.Heat.Flux)
## Tropo.Press
hist(epa$Tropo.Press)
## Press.MSL
hist(epa$Press.MSL)
## Vegetation
hist(epa$Vegetation)


## EDA ##################################################################
## Ozone
printPlot(paper,"img/eda/histogram_ozone.jpeg",5,5,FUN=function(){
  hist(epa$Temperature, main="", xlab="Ozone (ppm)")
})
printPlot(paper, "img/eda/heatmap_ozone.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Ozone")
})
printPlot(paper, "img/eda/ts_ozone.jpeg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Ozone", "Ozone (ppm)")
})
## Temperature
printPlot(paper,"img/eda/reanalysis/histogram_temperature.jpeg",5,5,FUN=function(){
  hist(epa$Temperature, main="", xlab="Air Temperature at 2m (K)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_temperature.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Temperature")
})
printPlot(paper, "img/eda/reanalysis/ts_temperature.jpeg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Temperature", "Temperature (K)")
})

## Relative Humidity
printPlot(paper,"img/eda/reanalysis/histogram_rh.jpeg",5,5,FUN=function(){
  hist(epa$RH, main="", xlab="Relative humidity (%)")
})
## Notes: bimodal distribution
## For LDA you really do not want to remove the bimodality, as it will reduce the discrimination of your analysis
## https://stats.stackexchange.com/questions/209241/what-transformation-should-i-use-for-a-bimodal-distribution
printPlot(paper, "img/eda/reanalysis/heatmap_rh.jpeg", 6, 7, FUN=function(){
  #heatmapPlusTs(epa, "RH", "Relative humidity (%)")
  simple_heatmap(epa, "RH")
})
printPlot(paper, "img/eda/reanalysis/ts_rh.jpeg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "RH", "Rel. humidity (%)")
})

## Wind speed
printPlot(paper,"img/eda/reanalysis/histogram_ws.jpeg",5,5,FUN=function(){
  hist(epa$Wind, main="", xlab="Wind speed (m/s)")
})
printPlot(paper,"img/eda/reanalysis/histogram_ws_sqrt.jpeg",5,5,FUN=function(){
  hist(sqrt(epa$Wind), main="", xlab="sqrt(Wind speed)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_ws.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Wind")
})
printPlot(paper, "img/eda/reanalysis/ts_ws.jpeg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Wind", "Wind speed (m/s)")
})

## Rain
printPlot(paper,"img/eda/reanalysis/histogram_rain.jpeg",5,5,FUN=function(){
  hist(epa$Rain, main="", xlab="Precipitation (kg/m^2)")
})
printPlot(paper,"img/eda/reanalysis/histogram_rain_log.jpeg",5,5,FUN=function(){
  hist(log(epa$Rain), main="", xlab="log(Precipitation)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_rain.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Rain")
})
printPlot(paper, "img/eda/reanalysis/ts_rain.jpeg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Rain", "Precipitation (Kg/m^2)")
})
## Dew.Point
printPlot(paper,"img/eda/reanalysis/histogram_dp.jpeg",5,5,FUN=function(){
  hist(epa$Dew.Point, main="", xlab="Dew point (K)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_dp.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Dew.Point")
})
printPlot(paper, "img/eda/reanalysis/ts_dp", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Dew.Point", "Dew point (K)")
})
## Water.Evaporation
printPlot(paper,"img/eda/reanalysis/histogram_we.jpeg",5,5,FUN=function(){
  hist(epa$Water.Evap, main="", xlab="Water evaporation (kg/m^2)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_we.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Water.Evap")
})
printPlot(paper, "img/eda/reanalysis/ts_we", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Water.Evap", "Water Evaporation (kg/m^2)")
})
## Heat.Flux
printPlot(paper,"img/eda/reanalysis/histogram_hf.jpeg",5,5,FUN=function(){
  hist(epa$Heat.Flux, main="", xlab="Ground Heat flux (W/m^2)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_hf.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Heat.Flux")
})
printPlot(paper, "img/eda/reanalysis/ts_hf", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Heat.Flux", "Ground Heat Flux (W/m^2)")
})
## Geop.Height
printPlot(paper,"img/eda/reanalysis/histogram_gh.jpeg",5,5,FUN=function(){
  hist(epa$Geop.Height, main="", xlab="Surface Geopotential Height (m)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_gh.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Geop.Height")
})
printPlot(paper, "img/eda/reanalysis/ts_gh", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Geop.Height", "Surface Geopotential Height (m)")
})
## Geop.Height.Tropo
printPlot(paper,"img/eda/reanalysis/histogram_ght.jpeg",5,5,FUN=function(){
  hist(epa$Geop.Height.Tropo, main="", xlab="Tropopause Geopotential Height (m)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_ght.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Geop.Height.Tropo")
})
printPlot(paper, "img/eda/reanalysis/ts_ght", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Geop.Height.Tropo", "Tropopause Geopotential Height (m)")
})
## Lat.Heat.Flux
printPlot(paper,"img/eda/reanalysis/histogram_lhf.jpeg",5,5,FUN=function(){
  hist(epa$Lat.Heat.Flux, main="", xlab="Latent Heat Flux (W/m^2)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_lhf.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Lat.Heat.Flux")
})
printPlot(paper, "img/eda/reanalysis/ts_lhf", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Lat.Heat.Flux", "Latent Heat Flux (W/m^2)")
})
## Tropo.Press
printPlot(paper,"img/eda/reanalysis/histogram_tp.tp",5,5,FUN=function(){
  hist(epa$Tropo.Press, main="", xlab="Tropopause Pressure (Pa)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_tp.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Tropo.Press")
})
printPlot(paper, "img/eda/reanalysis/ts_tp", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Tropo.Press", "Tropopause Pressure (Pa)")
})
## Press.MSL
printPlot(paper,"img/eda/reanalysis/histogram_pmsl.tp",5,5,FUN=function(){
  hist(epa$Press.MSL, main="", xlab="Pressure at MSL (Pa)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_pmsl.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Press.MSL")
})
printPlot(paper, "img/eda/reanalysis/ts_pmsl", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Press.MSL", "Pressure at MSL (Pa)")
})
## Vegetation"
printPlot(paper,"img/eda/reanalysis/histogram_veg.tp",5,5,FUN=function(){
  hist(epa$Vegetation, main="", xlab="Vegetation (%)")
})
printPlot(paper, "img/eda/reanalysis/heatmap_veg.jpeg", 6, 7, FUN=function(){
  simple_heatmap(epa, "Vegetation")
})
printPlot(paper, "img/eda/reanalysis/ts_veg", 5, 3, FUN=function(){
  timeSeriesPlotWithInterval(epa, "Vegetation", "Vegetation (%)")
})


## Correlations ##################################################################
library(corrplot)
names(epa)
excludedCovariates1 <- c(1,2,6,7,19,20,22,23:28,31:32) 
excludedCovariates2 <- c(1,2,6,7,19,20,22,23,24,26:28,31:32) 
sf <- # selected features
  c("Ozone","Temperature","RH","Dew.Point","Water.Evap","Heat.Flux",
    "Geop.Height","Geop.Height.Tropo","Lat.Heat.Flux","Tropo.Press","Press.MSL",
    "Vegetation","Longitude","Latitude","Elevation","Doy","sqrtWind","logRain")
epaCore <- epa[,sf]
head(epaCore)
M <- cor(epaCore)
printPlot(paper, "img/eda/analysis/correlogram.jpeg", 6, 6, FUN=function(){
  corrplot(M, type="lower", tl.col="black", tl.srt=45,tl.cex=0.8,order="hclust", mar=c(0,0,0,0))
})

## Ozone versus date fields
printPlot(paper, "img/eda/analysis/OzoneVsMonths.jpeg", 6, 5, FUN=function(){
  boxplot(Ozone~Month, epa, xlab="Month", ylab="Ozone")  # clear pattern is noticeable
})
printPlot(paper, "img/eda/analysis/OzoneVsDow.jpeg", 6, 5, FUN=function(){
  boxplot(Ozone~Dow.number, epa, xlab="Day of week", ylab="Ozone")  # this may not be significant, as day of week variable is a subset of the month var.
})

## Ozone versus metereological variables
printPlot(paper, "img/eda/analysis/OzoneVsTemp.jpeg", 6, 5, FUN=function(){
  boxplot(Ozone~round(Temperature), epa, xlab="discretized(Temperature)", ylab="Ozone")  
  ## Non-linear relationship
})
printPlot(paper, "img/eda/analysis/OzoneVsRH.jpeg", 6, 5, FUN=function(){
  boxplot(Ozone~round(RH/2), epa, xlab="discretized(RH)", ylab="Ozone")  
})
boxplot(Ozone~round(Rain), epa, main="Ozone")  #Nothing
boxplot(Ozone~round(Wind), epa, main="Ozone")  
boxplot(Ozone~round(Dew.Point), epa, main="Ozone")  
boxplot(Ozone~round(Water.Evap*10), epa, main="Ozone")  
boxplot(Ozone~round(Heat.Flux), epa, main="Ozone")  
boxplot(Ozone~round(Geop.Height/10), epa, main="Ozone")  
boxplot(Ozone~round(Geop.Height.Tropo/100), epa, main="Ozone")  
boxplot(Ozone~round(Lat.Heat.Flux/10), epa, main="Ozone")  
boxplot(Ozone~round(Tropo.Press/1000), epa, main="Ozone")  
printPlot(paper, "img/eda/analysis/OzoneVsPressMSL.jpeg", 6, 5, FUN=function(){
  boxplot(Ozone~round(Press.MSL/100), epa, xlab="discretized(Pressure at MSL)", ylab="Ozone")  
})
boxplot(Ozone~round(Vegetation), epa, main="Ozone")  ## Nothing



