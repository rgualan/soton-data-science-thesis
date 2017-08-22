## Libraries
library(ggplot2)
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
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
  d$Doy <- as.integer(format(d$Date,"%j"))
  return(d)
}

## Add Day-Of-Week field
addDowField <- function(d){
  d$Dow.name <- weekdays(d$Date)
  ## Build mapping data.frame
  Dow.name <- c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")
  Dow.number <- 1:7
  dfDow <- data.frame(Dow.name=Dow.name,Dow.number=Dow.number)
  ## Combine
  d$Dow.number <- match(d$Dow.name, dfDow$Dow.name)
  head(d,n=14)
  ## Change str to factor
  d$Dow.name <- as.factor(d$Dow.name)
  return(d)
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
timeSeriesPlotWithInterval <- function(d, var.name="Ozone", label="Ozone (ppb)", file.name=NA){
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
    #scale_x_datetime(expand=c(0,0)) +
    #scale_y_continuous(expand=c(0,0)) +
    labs(x = "Date", y = label) #+
    # theme_bw() + 
    # theme(panel.grid=element_blank(),
    #       plot.margin=unit(c(0.1,0.5,0.5,0.5),"cm")) 

  if(!is.na(file.name)) jpeg(file.name, 6, 4, "in", bg="white", res=150)
  print(p2)
  if(!is.na(file.name)) dev.off()
  
}




## Extracts time series from a nc file
## nfile: an ncFile
## dateA: initial date
## dateB: end date
## stations: stations, with coordinates to extract the values. Eg:
extractTimeSeriesFromNc <- function(ncFile, dateA, dateB, stations, colNames){
  # read the netCDF file as a raster layer
  r <- raster(ncFile)
  #r; print(r)
  
  covariate <- data.frame( date=c(),Station.ID=c(),tmp=c())
  #i=1
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
evaluatePredictions <- function (z, zhat) 
{
  VMSE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sum(u^2)/length(u), 4)
  }
  RMSE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sqrt(sum(u^2)/length(u)), 4)
  }
  MAE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- abs(c(zhat - z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  MAPE <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- abs(c(zhat - z))/z
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u) * 100, 4)
  }
  BIAS <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  rBIAS <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u) * mean(z, na.rm = TRUE)), 4)
  }
  rMSEP <- function(z, zhat) {
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat, na.rm = TRUE) - z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), 4)
  }
  r2 <- function(z, zhat) {
    a <- cor(z,zhat,use="pairwise.complete.obs")^2
    round(a, 4)
  }
  r2(1:10,2:11)
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
  unlist(out)
}


## Wraper to the sequence of instructions for printing a plot
## in a jpeg file
printPlot <- function(paper=T,file,width=6,height=6,res=150,FUN){
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
  Sys.time()-st
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
  kNb <-sort(actualRow)[1:(k+1)]

  #return(kIds[-1])
  if(!include_own)
    output <- data.frame(Station.Code=names(kNb[-1]), distance=kNb[-1])
  else
    output <- data.frame(Station.Code=names(kNb), distance=kNb)
  
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


plotStations <- function(paper=T, IDS, fileName,width,height, redIds=NA){
  sites <- getAllSites()
  sites <- sites[sites$Station.Code %in% IDS,]

  if(!is.na(redIds)){
    IDS <- IDS[!IDS %in% redIds]
  }
    
  printPlot(paper,fileName,width,height,FUN= function(){
    p <- ggplot(getCAmap()) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
      geom_point(data = sites[sites$Station.Code %in% IDS,], aes(x = Longitude, y = Latitude, fill=Location.Setting),
                 alpha = 0.75, shape=24, size=2) +
      labs(x = "Longitude", y = "Latitude", fill="Type") +
      coord_quickmap() +
      theme(legend.justification = c("right", "top"), legend.position = c(.95, .95),
             legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
    if(!is.na(redIds)){
      p <- p + geom_point(data = sites[sites$Station.Code %in% c(redIds),], 
                          aes(x = Longitude, y = Latitude),
                          fill="red", alpha = 0.75, shape=24, size=4)
    }
      
    plot(p)
  })
}


