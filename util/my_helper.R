## Libraries
library(ggplot2)
library(sp)
library(lattice)
library(RColorBrewer)
library(grid)
library(gridExtra)


## Setup printing figures configuration
setupPaper <- function(){
  if(exists(".paper")){paper<-.paper}else{paper<-T}  
  return(paper)
}


## Convert regular date (yyyy-mm-dd) to  POSIXct 
## Fixing certain tz (TODO: consider to avoid this)
convertStringToPOSIXct <- function(a){
  return( as.POSIXct(a,format="%Y-%m-%d", tz="GMT") )
}

## Add Day-Of-Year field
addDoyField <- function(d){
  # x <- 1:366
  # y <- cos((x-1)*(2*pi)/(366))
  # plot(y, type="l")

  d$Doy <- as.integer(format(d$Date,"%j"))
  d$Doy <- cos((d$Doy-1)*(2*pi)/(366))
  return(d)
}

## Get folds for a 10-fold CV
getFolds <- function(){
  folds <- readRDS("output/folds.RDS")
  return(folds)
}

## Add Day-Of-Week field
addDowField <- function(d, sinTx=F){
  d$Dow.name <- weekdays(d$Date)
  ## Build mapping data.frame
  Dow.name <- c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")
  Dow.number <- 1:7
  if(sinTx)
    Dow.number <- sin(2*pi*(Dow.number-1)/6)
  dfDow <- data.frame(Dow.name=Dow.name,Dow.number=Dow.number)
  ## Combine
  d$Dow.number <- match(d$Dow.name, dfDow$Dow.name)
  head(d,n=14)
  ## Change str to factor
  d$Dow.name <- as.factor(d$Dow.name)
  return(d)
}

## Add Day-Of-Week field
addIsWeekDay <- function(d){
  d$isWeekday <- T
  d$isWeekday[d$Dow.name %in% c("Saturday","Sunday")] <- F
  d$isWeekday <- as.factor(d$isWeekday)
  return(d)
}

## Transformation
addIsWeekDay <- function(d){
  d$isWeekday <- T
  d$isWeekday[d$Dow.name %in% c("Saturday","Sunday")] <- F
  d$isWeekday <- as.factor(d$isWeekday)
  return(d)
}

## Attributes derived from Date
addDateDerivedFeatures <- function(epa){
  epa <- addDoyField(epa)
  epa <- addDowField(epa)
  epa <- addIsWeekDay(epa)
  epa$Month <- as.numeric(format(epa$Date, "%m"))
  epa$Day <- as.numeric(format(epa$Date, "%d"))
  return(epa)
}

## Transform features
## sqrtWind
## logRain
transformFeatures <- function(epa){
  epa$sqrtWind <- sqrt(epa$Wind)
  epa$logRain <- log(1+epa$Rain, 10) # The additive term is to avoid -inf
  return(epa)
}

## Scale variable and add a Bias that avoids negative values
scalePlusBias <- function(column){
  column2 <- scale(column)
  column2 <- column2+abs(min(column2,na.rm=T))
  return(column2)    
}

## Get days of the year in Date format
getDates <- function(dateA="2016-01-01", dateB="2016-12-31"){
  da <- convertStringToPOSIXct(dateA)
  db <- convertStringToPOSIXct(dateB)
  days <- seq(from=da, to=db, by='days')
  days2 <- as.Date(days)
  return(days2)
}

## Add combination of closest neighbours as a predictor
## New column: Neighbor
addNeighboursAverage <- function(epa, k=5, plotTs=F){
  sites <- getSites(epa)
  epa$Neighbor <- NA
  for(s in sites$Station.Code){
    #print(s)
    kn <- getKneighbours(s,sites,k,include_own = F)
    ## Most of the cases there will be a strong correlation, between the current station
    ## and the kneighbor stations
    original <- epa$Ozone[epa$Station.Code==s]
    
    ## Inverse distance weighted average
    w <- kn$distance/sum(kn$distance)
    A <- matrix(NA,nrow=length(original), ncol=k); colnames(A) <- kn$Station.Code
    for(j in 1:nrow(kn)){
      A[,j] <- epa$Ozone[epa$Station.Code==kn$Station.Code[j]]
    }
    b <- A%*%w
    epa$Neighbor[epa$Station.Code==s] <- b
    
    ## Debug:
    #C <- cbind(original,A,combined=b); cor(C)
    ## Plot:
    if(cor(original,b)<0.5 && plotTs){
      plotStations(F,kn$Station.Code,redIds = s); readline("Continue?")
      plot(A[,1],main=cor(original,b),type="l",lty="dashed",lwd=0.5)
      for(i in 2:k){lines(A[,i],lty="dashed",lwd=0.5,col=i)}
      lines(b,col=2,lwd=2); lines(original,col=3,lwd=2)
      readline("Continue?") 
    }
  }
  ## Problem:
  ## Station 023-1005
  ## It is isolated in the northwest R:-0.18!!!
  ## Luckily this station is the only problem
  # epa$Neighbor <- scale(epa$Neighbor)[,1]
  # cor(epa$sOzone,epa$Neighbor)
  ## Notes:
  ## k was chosen based on the best R2
  return(epa)
}

scaleTargetVariable <- function(epa){
  epa$sOzone <- scale(epa$Ozone)[,1]
  return(epa)
}

## Fill missing dates in the original ds 
fillMissingDates <- function(d, begin="2016-01-01", end="2016-12-31"){
  da <- convertStringToPOSIXct(begin)
  db <- convertStringToPOSIXct(end)
  days <- seq(from=da, to=db, by='days' )
  
  ## initial column names
  colNames <- names(d)
  
  # Fill mising dates for each station
  complete <- data.frame()
  if("Station.Code" %in% colNames){
    for(s in unique(d$Station.Code)){
      d2 <- d[d$Station.Code==s,]
      d2 <- merge(d2,data.frame(Date=days), all.y=T)
      d2$Station.Code[is.na(d2$Station.Code)] <- s
      complete <- rbind(complete,d2)
    }
    complete <- complete[,colNames]
  }else{
    complete <- merge(d,data.frame(Date=days),all.y=T)
  }
  #View(complete)

  return(complete)
}


## Convert a "Regular" time series data frame 
## Example:
#      Station.Code       Date Temperature
# 1       005-0002 2016-01-01    33.57083
# 141     005-0002 2016-01-02    39.91667
# 282     005-0002 2016-01-03    45.08333
## from one of the EPA's stations
## and convert it into an SpatialPointDataFrame(sp)
convertDataToSp <- function(d, projPar="utm") {
  if(projPar=="utm"){
    p = CRS("+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +units=km")
    #newColumns = c("Station.Code","UTM.X","UTM.Y","Elevation")
    newColumns = c("Station.Code","UTM.X","UTM.Y")
    locFormula = ~UTM.X+UTM.Y
    #if(z) locFormula = ~UTM.X+UTM.Y+Elevation
  }else if(proj=="longlat"){
    p = CRS("+proj=longlat +datum=WGS84")
    #newColumns = c("Station.Code","Longitude","Latitude","Elevation")
    newColumns = c("Station.Code","Longitude","Latitude")
    locFormula = ~Longitude+Latitude
    #if(z) locFormula = ~Longitude+Latitude+Elevation
  }else{
    stop("Wrong projection parameter")
  }
  
  # if (!exists("sites")){
  #   sites <- readRDS("data/epa/sites/aqs_sites.RDS")
  # }
  tmpSites <- getSites(d)
  tmpSites <- tmpSites[, newColumns]
  names(tmpSites)[2:3] <- c("x","y") 

  #d2 <- merge(d, tmpSites[,newColumns]) # , all.x=T # Implicitly ignore stations
  d2 <- merge(d, tmpSites) # , all.x=T # Implicitly ignore stations
  coordinates(d2) <- ~x+y
  proj4string(d2) <- p
  
  return(d2)
}

getUTMproj <- function() {
  return(CRS("+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +units=km"))
}



## Creates a heatmat of a variable (Station.Code versus Date) and
## a mean daily time series with a confidence interval.
## Data:
# Station.Code       Date Measurement
# 1601     001-2005 2016-01-01    0.026882
# 1543     001-2005 2016-01-02    0.020412
# 1864     001-2005 2016-01-03    0.024588
## Daily aggregate:
# Date       mean      upper      lower
# 1 2016-01-01 0.02790952 0.04233570 0.01515330
# 2 2016-01-02 0.02745810 0.04363825 0.01100570
# 3 2016-01-03 0.02888175 0.04411945 0.01431195

heatmapPlusTs <- function(d, hm.z.name="Ozone", ts.y.lab="Ozone (ppb)", file.name=NA){
  ## Daily mean
  f=formula(paste0(hm.z.name,"~Date"))
  dailyAgg <- aggregate(f,d,mean) # daily mean
  names(dailyAgg)[2] <- "mean"
  ua <- aggregate(f,d,quantile,0.95) # upper
  la <- aggregate(f,d,quantile,0.05) # lower
  dailyAgg$upper <- ua[,hm.z.name]
  dailyAgg$lower <- la[,hm.z.name]
  #View(dailyAgg)
  
  # heatmap.z.name
  p1<-ggplot(d[d$Station.Code %in% sample(unique(d$Station.Code),50),], 
             aes_string(x="Date", y="Station.Code", fill=hm.z.name)) +
    geom_raster()  +
    scale_fill_distiller(palette = "Spectral") +
    scale_x_datetime(expand=c(0,0)) +
    labs(y = "Station Code") + #, fill="Ozone (ppm)"
    theme_bw() + 
    theme(legend.position = "none", panel.grid=element_blank(), 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size=7),  
          plot.margin=unit(c(0.5,0.5,0.1,0.5),"cm")) 
  
  p2<-ggplot(dailyAgg, aes(Date, mean)) + 
    geom_smooth(aes(ymax=upper, ymin=lower), stat='identity') +
    scale_x_datetime(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x = "Date", y = ts.y.lab) +
    theme_bw() + 
    theme(panel.grid=element_blank(),
          plot.margin=unit(c(0.1,0.5,0.5,0.5),"cm")) 
  
  grid.newpage()
  grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first"))
  
  gA=ggplot_gtable(ggplot_build(p1))
  gB=ggplot_gtable(ggplot_build(p2))
  maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
  gA$widths[2:3] <- as.list(maxWidth)
  gB$widths[2:3] <- as.list(maxWidth)
  grid.newpage()

  if(!is.na(file.name)) jpeg(file.name, 6, 7, "in", bg="white", res=150)
  grid.arrange(
    arrangeGrob(gA,gB,nrow=2,heights=c(.8,.3))
  )
  if(!is.na(file.name)) dev.off()
  
}

simple_heatmap <- function(d, hm.z.name="Ozone"){
  # heatmap.z.name
  p1<-ggplot(d[d$Station.Code %in% sample(unique(d$Station.Code),50),], 
             aes_string(x="Date", y="Station.Code", fill=hm.z.name)) +
    geom_raster()  +
    scale_fill_distiller(palette = "Spectral") +
    scale_x_datetime(expand=c(0,0)) +
    labs(x="Date", y = "Station Code") +
    theme_bw() + 
    theme(legend.position = "none", panel.grid=element_blank()) 
  
  print(p1)
}

# dtmp<-addCountyColumn(d)
# #dtmp<-dtmp[dtmp$County.Code %in% head(unique(dtmp$County.Code)),]
# ggplot(dtmp, aes_string(x="Date", y="Station.Code", fill="Measurement")) +
#   geom_raster()  +
#   scale_fill_distiller(palette = "Spectral") +
#   scale_x_datetime(expand=c(0,0)) +
#   labs(y = "Station Code") + #, fill="Ozone (ppm)"
#   theme_bw() + 
#   facet_grid(County.Code~., scales="free", space = "free") + #, switch = "y"
#   theme(legend.position = "none", 
#         axis.title.x= element_blank(),
#         axis.text.x= element_blank(),
#         axis.text.y= element_text(size=7),  
#         plot.margin= unit(c(0.5,0.5,0.1,0.5),"cm"),
#         #plot.margin= unit(c(0.0,0.0,0.0,0.0),"cm"),
#         panel.spacing= unit(0, "lines"), 
#         panel.border = element_blank(),
#         panel.grid = element_blank(),
#         strip.background= element_blank() )


## Creates a time series plot
## with a confidence interval of 0.05 to 0.95 percetiles
timeSeriesPlotWithInterval <- function(d, var.name="Ozone", label="Ozone (ppb)", file.name=NA, xLabel=T){
  ## Daily mean
  f=formula(paste0(var.name,"~Date"))
  dailyAgg <- aggregate(f,d,mean) # daily mean
  names(dailyAgg)[2] <- "mean"
  ua <- aggregate(f,d,quantile,0.95) # upper
  la <- aggregate(f,d,quantile,0.05) # lower
  dailyAgg$upper <- ua[,var.name]
  dailyAgg$lower <- la[,var.name]
  #View(dailyAgg)
  
  p2<-ggplot(dailyAgg, aes(Date, mean)) + 
    geom_smooth(aes(ymax=upper, ymin=lower), stat='identity') +
    labs(y = label)
  if(xLabel){
    p2 <- p2 + labs(x = "Date")
  }else{
    p2 <- p2 + theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank() #,axis.ticks.x=element_blank()
                    )
  }
  
  #if(!is.na(file.name)) jpeg(file.name, 6, 4, "in", bg="white", res=150)
  print(p2)
  #if(!is.na(file.name)) dev.off()
}

## Extracts time series from a nc file
## nfile: an ncFile
## dateA: initial date
## dateB: end date
## stations: stations, with coordinates to extract the values. Eg:
extractTimeSeriesFromNc <- function(ncFile, dateA, dateB, stations, colNames, singleLevel=F){
  # read the netCDF file as a raster layer
  r <- raster(ncFile)
  #r; print(r)
  
  covariate <- data.frame( date=c(),Station.ID=c(),tmp=c())
  #i=1
  if(!singleLevel){
    for(i in 1:nbands(r)){
      #for(i in 1:5){
      print(i)  
      r <- raster(ncFile, band=i)
      z <- convertStringToPOSIXct(getZ(r))
      if(z>=dateA & z<=dateB){
        print("Processing")
        vals <- extract(r, as.matrix(stations), #stations[,2:3]
                        method='bilinear', fun=mean, na.rm=TRUE)
        tmp <- data.frame( Station.ID=rownames(stations), #stations[,1]
                           Date=rep(z,length(vals)),
                           var=vals)
        covariate <- rbind(covariate, tmp)
      }
    }
  }else{
    r <- raster(ncFile)
    vals <- extract(r, as.matrix(stations), 
                    method='bilinear', fun=mean, na.rm=TRUE)
    theDates <- getDates()
    covariate <- data.frame( Station.ID=rep(rownames(stations), length(theDates)), 
                       Date=rep(theDates, each=nrow(stations)),
                       var=rep(vals,length(theDates)))
  }
  #names(covariate)[which(names(covariate)=="var")] <- colName 
  names(covariate) <- colNames
  return(covariate)  
}


## Get all sites
getAllSites <- function(){
  sitesMd <- readRDS("data/epa/sites/aqs_sites.RDS")
  return(sitesMd)
}

## Get sites from the data 
getSites <- function(d){
  sitesMd <- readRDS("data/epa/sites/aqs_sites.RDS")
  
  if(class(d)=="data.frame"){
    ss <- sort(unique(d$Station.Code))  
  }else{
    ss <- sort(d)
  }
  
  sites <- merge(data.frame(Station.Code=ss), 
                 sitesMd[,c("Station.Code","Latitude","Longitude","Elevation",
                            "Location.Setting","UTM.X","UTM.Y")])
  return(sites)
}

## Add County column
addCountyColumn <- function(d){
  sitesMd <- getAllSites()
  d2 <- merge(d, sitesMd[,c("Station.Code","County.Code")])
  #d2$County.Code <- d2$County.Code
  d2$County.Code <- as.factor(d2$County.Code)
  return(d2)
}


## Get USA map
getUSAmap <- function(){
  mapUSA <- readRDS("data/maps/usa/USA_adm1.rds")
  #proj4string(mapUSA)
  mapUSA <- mapUSA[!mapUSA$NAME_1 %in% c("Alaska","Hawaii"),]
  return(mapUSA)
}

## Get California map
getCAmap <- function(proj="longlat"){
  mapUSA <- getUSAmap()
  mapCA <- mapUSA[mapUSA$NAME_1=="California",]
  
  if(proj=="utm"){
    #mapCA <- spTransform(mapCA,getUTMproj())
    #saveRDS(mapCA, file="data/maps/mapCA.RDS")
    mapCA <- readRDS("data/maps/mapCA.RDS")
  }

  return(mapCA)
}

## Get days of the year 2016 in Posixct format
getStudyDays <- function(){
  # numDays <- as.integer((convertStringToPOSIXct("2016-12-31")
  #                        -convertStringToPOSIXct("2016-01-01"))+1)
  days <- seq(from=convertStringToPOSIXct("2016-01-01"), 
              to=convertStringToPOSIXct("2016-12-31"), by='days' )
  return(days)
}


## Calculates goodness of fit
evaluatePredictions <- function (z, zhat, n=3) 
{
  VMSE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sum(u^2)/length(u), n)
  }
  RMSE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sqrt(sum(u^2)/length(u)), n)
  }
  MAE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- abs(c(zhat - z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), n)
  }
  MAPE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- abs(c(zhat - z))/z
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u) * 100, n)
  }
  BIAS <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), n)
  }
  rBIAS <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u) * mean(z, na.rm = TRUE)), n)
  }
  rMSEP <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat, na.rm = TRUE) - z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), n)
  }
  r2 <- function(z, zhat) {
    a <- cor(z,zhat,use="pairwise.complete.obs")^2
    round(a, n)
  }
  #cat("##\n Mean Squared Error (MSE) \n Root Mean Squared Error (RMSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Bias (BIAS) \n Relative Bias (rBIAS) \n Relative Mean Separation (rMSEP)\n##\n")
  out <- NULL
  out$MSE <- VMSE(c(z), c(zhat))
  out$RMSE <- RMSE(c(z), c(zhat))
  out$MAE <- MAE(c(z), c(zhat))
  out$MAPE <- MAPE(c(z), c(zhat))
  out$BIAS <- BIAS(c(z), c(zhat))
  out$rBIAS <- rBIAS(c(z), c(zhat))
  out$rMSEP <- rMSEP(c(z), c(zhat))
  out$r2 <- r2(c(z), c(zhat))
  out <- unlist(out)
  return(round(out, digits=n))
}


## Wraper to the sequence of instructions for printing a plot
## in a jpeg file
printPlot <- function(paper=T,file,width=6,height=6,FUN){
  if(paper) jpeg(file, width, height, "in", bg="white", res=150)

  # Evalute the desired series of expressions inside of tryCatch
  tryCatch({
    FUN()
  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("THERE WAS AN ERROR: ",err))
  }, finally = {
    if(paper) dev.off()    
  }) # END tryCatch

}

## Simple timer
ticToc <- function(expr){
  st<-Sys.time()
  expr
  print(Sys.time()-st)
}


## Highlight individual station
highlightStation <- function(s){
  
  ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = getAllSites(), aes(x = Longitude, y = Latitude, col=Location.Setting), #, fill=Location.Setting
               alpha = 0.75, shape=17, size=1) +
    geom_point(data = sites[sites$Station.Code==s,], 
               aes(x = Longitude, y = Latitude),
               alpha = 0.75, shape=24, size=2, fill="red") +
    labs(x = "Longitude", y = "Latitude", fill="Type") +
    coord_quickmap() +
    theme(legend.justification = c("right", "top"), legend.position = c(.95, .95),
          legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
  
}


getKneighbours<-function(s, sites, k, include_own=F){
  ## Calculate distance matrix
  ks <- sites[,c("Station.Code","UTM.X","UTM.Y")] # known sites
  m <- as.matrix(dist(ks[,-1]))
  #dim(m)
  #m[1:5,1:5]
  
  ## Get nearest neighbours
  actualRow <- m[which(s == sites$Station.Code),]
  kNb <-sort(actualRow)[1:(k+1)]

  includeIndex <- 1:length(kNb)
  if(!include_own)
    includeIndex <- 2:length(kNb)

  output <- data.frame(
    Station.Code=sites$Station.Code[as.numeric(names(kNb[includeIndex]))], 
    distance=kNb[includeIndex])
  
  return(output)
}

getKneighboursInRadius<-function(s, sites, radius, include_own=F){
  ## Calculate distance matrix
  ks0 <- sites[,c("Station.Code","UTM.X","UTM.Y")] # known sites
  ks <- ks0[,-1]
  rownames(ks) <- ks0$Station.Code
  #head(ks)
  m <- as.matrix(dist(ks))
  #dim(m)
  #m[1:5,1:5]
  
  ## Get nearest neighbourS
  actualRow <- m[s,]
  #kIds <- names(sort(actualRow))[1:(k+1)]
  #kNb <-sort(actualRow)[1:(k+1)]
  kNb <-sort(actualRow[actualRow<radius])
  
  #return(kIds[-1])
  if(!include_own)
    output <- data.frame(Station.Code=names(kNb[-1]), distance=kNb[-1])
  else
    output <- data.frame(Station.Code=names(kNb), distance=kNb)
  return(output)
}

## Fil:
## 1: Location.type
## 2: Station.Code
plotStations <- function(paper=T, IDS, fileName, width, height, redIds=NA, fill=1){
  sites <- getAllSites()
  if(!anyNA(redIds)){
    regularPlusRed <- unique(rbind(data.frame(sc=IDS),
                                   data.frame(sc=redIds)))
  }else{
    regularPlusRed <- unique(data.frame(sc=IDS)) 
  } 

  sites <- sites[sites$Station.Code %in% regularPlusRed$sc,]
  
  if(!anyNA(redIds)){
    IDS <- IDS[!IDS %in% redIds]
  }
  
  printPlot(paper,fileName,width,height,FUN= function(){
    p <- ggplot(getCAmap()) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
      geom_point(data = sites[sites$Station.Code %in% IDS,], aes(x = Longitude, y = Latitude, fill= if(fill==1) Location.Setting else Station.Code),
                 alpha = 0.75, shape=24, size=2) +
      labs(x = "Longitude", y = "Latitude", fill="Type") +
      coord_quickmap() +
      theme(legend.justification = c("right", "top"), 
            legend.position = c(.95, .95),
            legend.box.background = element_rect()#, 
            #legend.box.margin = margin(6, 6, 6, 6)
      )
    if(!anyNA(redIds)){
      p <- p + geom_point(data = sites[sites$Station.Code %in% redIds,], 
                          aes(x = Longitude, y = Latitude),
                          fill="red", alpha = 1, shape=24, size=2.5)
    }
    
    plot(p)
  })
}


## Assemble STFDF 
assembleSTFDF <- function(epa){
  ## Space dimension
  epa.sp <- getSites(epa)
  epa.sp <- epa.sp[,c("Station.Code","Latitude","Longitude","Location.Setting","UTM.X","UTM.Y")]
  rownames(epa.sp) <- epa.sp$Station.Code
  head(epa.sp)
  coordinates(epa.sp) <- ~UTM.X+UTM.Y
  proj4string(epa.sp) <- getUTMproj()
  
  ## Time dimension
  epa.tm <- sort(unique(epa$Date))
  epa.tm <- as.Date(epa.tm)  ## Ignore time data?? >> Corrects x labels problem in acf
  
  ## Data
  epa_2 <- epa[,c("Date","Ozone","sOzone","Temperature",
                  "RH","Rain","logRain","Wind","sqrtWind","Elevation")]
  ## TODO: Sort?
  
  # Combine the objects spatial, data-frame and time-dim into a STIDF:
  epa.st <- STFDF(epa.sp,epa.tm,epa_2) 
  # summary(epa.st)
  # stplot(epa.st[,"2016-01-01::2016-01-08","sOzone"])
  # dim(epa.st)

  return(epa.st)  
}

## Function to get the goodness-of-fit metrics by station 
getMetricsByStation <- function(epa.st,colName1,colName2){
  metrics = c()
  for(s in 1:dim(epa.st)[1]){
    m <- evaluatePredictions(as.vector(epa.st[s,,colName1][,1]), 
                             as.vector(epa.st[s,,colName2][,1]))
    metrics <- rbind(metrics,m)
    rownames(metrics)[nrow(metrics)] <- as.character(epa.st@sp[s,]$Station.Code)
  }
  return(metrics)
}

## Function to get the goodness-of-fit metrics by station (from a data.frame) 
getMetricsByStationFromDF <- function(epa,colName1,colName2){
  metrics = c()
  for(s in unique(epa$Station.Code)){
    m <- evaluatePredictions(epa[epa$Station.Code==s,colName1], 
                             epa[epa$Station.Code==s,colName2])
    metrics <- rbind(metrics,m)
    rownames(metrics)[nrow(metrics)] <- as.character(as.character(s))
  }
  return(metrics)
}
