## Clean environment
rm(list=ls())

## Load required libraries
library(sp)
library(gstat)
library(spacetime)
library(reshape2)
library(lattice)
library(spTimer)
source("util/my_helper.R")
source("preprocessing/gbmImpute.R")

## Global variables
paper <- setupPaper()

## Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
epa <- addIsWeekDay(epa)
## Sites
sites <- getSites(epa)
#str(epa)


## Simple imputation of Ozone ################################################
covByStation <- aggregate(cbind(Ozone,Temperature)~Station.Code,epa,length)
covByStation <- covByStation[order(covByStation$Ozone),]
testStation <- "065-2002" # The one with more missing data (35 NAs)


## Idea:
## Since there is a clear annual trend, and according to previous experiments the predictions
## are biased, the idea is to decompose the time series in:
## 1) a trend (4th order)
## 2) a random variable
## A weekly trend was explored, but it was not very clear the presence of a weekly trend
## Additionally, since we are interested in filling a particular station
## The design matrix should be restricted to N stations around that station
imputeNasByStation <- function(testStation, epa, fn_predictors="Temperature+Elevation+Dow.number+Doy", 
                               N=50){
  #fn <- paste("Ozone.res~Temperature+Rain+Elevation+Dow.number+Doy")
  fn <- formula(paste("Ozone.res~", fn_predictors))
  
  sites <- getSites(epa)
  kn <- getKneighbours(testStation, sites, N, include_own=T)
  epa_sub <- epa[epa$Station.Code %in% kn$Station.Code,] ## Zone of interest (near the target variable)
  ## Detrend the time series
  epa_sub$Ozone.trend <- NA
  epa_sub$Ozone.res <- NA
  for(s in unique(epa_sub$Station.Code)){
    #print(s)
    epa.test <- epa_sub[epa_sub$Station.Code==s,]
    model <- lm(Ozone ~ poly(Doy,4), epa.test) 
    trend <- predict(model,epa.test)
    epa.test$Ozone.res <- epa.test$Ozone-trend
    epa_sub$Ozone.trend[epa_sub$Station.Code==s] <- trend
    epa_sub$Ozone.res[epa_sub$Station.Code==s] <- epa.test$Ozone.res
    #plot(Ozone~Doy,epa.test,type="l")
    #lines(trend, col="deepskyblue", lwd=2)
    #plot(residuals~Doy,epa.test,type="l")
    #readline("Continue?")
  }
  # plot(Ozone~1,epa_sub, type="l")
  # lines(Ozone.trend~1,epa_sub, col=2)
  # plot(Ozone.res~1,epa_sub, type="l")
  
  ## Now use GP to model the residuals
  epa_sub$Ozone.res100 <- epa_sub$Ozone.res*100  ## Re-scale target variable!!!
  # plot(epa_sub$Ozone.res2, type="l")
  # hist(epa_sub$Ozone.res2)
  ticToc(
    simpleGp <- spT.Gibbs(
      formula = Ozone.res100~Temperature+Elevation+Dow.number+Doy,
      #formula = fn,
      model = "GP",
      data = epa_sub[epa_sub$Station.Code!=testStation,], 
      coords = ~UTM.X + UTM.Y,
      spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
  )
  epa.test <- epa_sub[epa_sub$Station.Code==testStation,]
  gpPred <- predict(simpleGp, newdata=epa.test, newcoords=~UTM.X+UTM.Y)
  epa.test$Ozone.res2 <- gpPred$Mean
  # plot(epa.test$Ozone.res100, type="l")
  # lines(epa.test$Ozone.res2, col=2)
  # evaluatePredictions(epa.test$Ozone.res2,epa.test$Ozone.res100)

  epa.test$Ozone.res2 <- epa.test$Ozone.res2/100## Re-scale target variable
    
  epa.test$Ozone2 <- epa.test$Ozone.trend + epa.test$Ozone.res2
  ## Plot
  # plot(Ozone~Date,epa.test,type="l")
  # lines(Ozone.trend~Date,epa.test,col=2)
  # lines(Ozone2~Date,epa.test,col=3)

  ## Metrics 
  mGp <- evaluatePredictions(epa.test$Ozone, epa.test$Ozone2)
  
  ## output
  output <- list(Station.Code=testStation,
                 prediction=epa.test$Ozone2,
                 trend=epa.test$Ozone.trend,
                 metrics=mGp)
  return(output)
}

## Test in one station
if(F){
  ## Inject missing data
  set.seed(1345)
  epa$Ozone.backup<-epa$Ozone
  randomDay <- sample(getStudyDays(), 1)
  epa$Ozone[epa$Station.Code == testStation 
            & epa$Date>=randomDay & epa$Date<=(randomDay+15*24*60*60)] <- NA
  #plot(epa$Ozone2[epa$Station.Code==testStation], type="l")
  rds <- sample(getStudyDays(), 10)
  for(i in 1:length(rds)){
    epa$Ozone[epa$Station.Code == testStation &
                epa$Date>=rds[i] & epa$Date<=(rds[i]+2*24*60*60)] <- NA
  }
  #sum(is.na(epa$Ozone2[epa$Station.Code==testStation]))

  ## Try different combinations of predictors
  a <- imputeNasByStation(testStation, epa, "Temperature+Elevation+Dow.number+Doy")
  a$metrics
  # b <- imputeNasByStation(testStation, epa, "Temperature+Rain+Wind+Elevation+Dow.number+Doy")
  # b$metrics
  # c <- imputeNasByStation(testStation, epa, "Temperature+Rain+sqrtWind+Elevation+Dow.number+Doy")
  # c$metrics
  # d <- imputeNasByStation(testStation, epa, "Temperature+Elevation+Dow.number+Doy")
  # d$metrics
  ## Dow.name >> Failed to find determinant: Matrix not positive definite !
  # a$metrics
  # b$metrics
  # c$metrics
}


## Results
# MSE    RMSE     MAE    MAPE    BIAS   rBIAS   rMSEP      r2 
# a) 0.0001  0.0079  0.0062 19.3290  0.0005  0.0120  0.2350  0.7764 ***
# b) 0.0001  0.0079  0.0063 19.6739  0.0005  0.0129  0.2371  0.7632 
# c) 0.0001  0.0079  0.0063 19.8662  0.0006  0.0160  0.2354  0.7654 



## Count the number of missing data blocks #######################################
for(s in covByStation$Station.Code[covByStation$Ozone<366]){
  tot <- covByStation$Ozone[covByStation$Station.Code==s]
  cat(sprintf("Count gaps of missing data: %s. Total: %d. ",s, tot))
  ## Only two kinds of gap:
  ## Gap of one / Gap of more than one
  counter <- c(0,0)
  tmp <- 0
  epa.test <- epa[epa$Station.Code==s,]
  for(v in epa.test$Ozone){
    if(is.na(v)){
      tmp<-tmp+1
    }else if(tmp>0){
      if(tmp==1){
        counter[1] <- counter[1]+1
      } else {
        counter[2] <- counter[2]+1 
      }
      tmp<-0
    }
  }
  print(counter)
}
## Notes:
## The amount of individual is small and do not motivate 
## creating a first module integrated 
## Only the stations with one missing value will be processed using this method



## Fill stations with one missing value ##########################################
epa$Flag<-0
for(s in covByStation$Station.Code[covByStation$Ozone<366]){
  #s="073-1016"
  tot <- covByStation$Ozone[covByStation$Station.Code==s]
  print(sprintf("Count gaps of missing data: %s. Total: %d",s, tot))
  epa.test <- epa[epa$Station.Code==s,]
  values <- epa.test$Ozone[order(epa.test$Doy)]  
  flags <- rep(0, length(values))
  ii <- which(is.na(values))
  counter <- 0
  for(i in ii){
    if(i>1 & i<366){
      if(!is.na(values[i-1]) & !is.na(values[i+1]) ){
        m <- mean(values[(i-1):(i+1)],na.rm=T)
        values[i] <- m
        flags[i] <- 1
        counter <- counter + 1
      }
      
    }      
  } 
  print(sprintf("Individual replacements: %d", counter))
  # plot(values, type="l")
  # points(which(!is.na(flags)), values[!is.na(flags)], col=2)
  # Replace values
  epa.test$Ozone[order(epa.test$Doy)] <- values 
  epa.test$Flag[order(epa.test$Doy)] <- flags 
  epa$Ozone[epa$Station.Code==s] <- epa.test$Ozone
  epa$Flag[epa$Station.Code==s] <- epa.test$Flag
  # plot(Ozone~Date,epa[epa$Station.Code==s,], type="l")
  # points(Ozone~Date,epa[epa$Station.Code==s & epa$Flag==1,], col=2)
  # readline("Continue?")
}


## Function to plot imputations
plotImputations <- function(paper, s, epa){
  test <- epa[epa$Station.Code==s,]
  # plot(Ozone~Date,epa.test,type="l")
  # lines(Ozone.trend~Date,epa.test,col=2)
  # lines(Ozone2~Date,epa.test,col=3)
  # readline("Continue?")
  
  ## Assemble temporal data.frame
  tmp <- rbind(data.frame(Date=test$Date,Ozone=test$Ozone,Type="Original",Flag=test$Flag==1),
               data.frame(Date=test$Date,Ozone=test$Ozone2,Type="Imputed",Flag=is.na(test$Ozone)),
               data.frame(Date=test$Date,Ozone=test$Ozone.trend,Type="Trend", Flag=F))
  ## Plot
  printPlot(paper, sprintf("img/preprocessing/ozonei/ts_ozone_%s.jpeg",s),7,4,FUN=function(){ 
    p<-ggplot(tmp) +
      annotate("rect",
               xmin=tmp$Date[tmp$Flag & tmp$Type=="Imputed"]-1*24*60*60,
               xmax=tmp$Date[tmp$Flag & tmp$Type=="Imputed"]+1*24*60*60,
               ymin=-Inf, ymax=Inf, alpha=0.75, fill="lightyellow") +
      geom_line(aes(x=Date, y=Ozone, col=Type)) +
      geom_point(data=tmp[tmp$Flag & tmp$Type=="Original",], 
                 aes(x=Date, y=Ozone), col="red", size=2) +#shape=21, stroke=1 +
      theme(legend.position = "top")
    print(p)
  })
}

## Apply the imputation method to all the stations
if(F){
  covByStation <- aggregate(cbind(Ozone)~Station.Code,epa,length)
  covByStation <- covByStation[covByStation$Ozone<366,]
  covByStation <- covByStation[order(covByStation$Ozone),]
  #View(covByStation)
  
  epa$Ozone2 <- NA
  for(s in covByStation$Station.Code){
    print(sprintf("Impute %s",s))
    out <- imputeNasByStation(s, epa, "Temperature+Elevation+Dow.number+Doy")
    #out$metrics
    epa$Ozone2[epa$Station.Code==s] <- out$prediction
    epa$Ozone.trend[epa$Station.Code==s] <- out$trend
    ## Plot  
    plotImputations(paper,s,epa)
  }
  ## Save imputed dataset
  saveRDS(epa, "data/epa/epa_daily/2016/california_ozone_plus_rcov_2.RDS")
}

## Print results (all at onece)
if(F){
  for(s in covByStation$Station.Code){
    plotImputations(paper, s,epa)
  }
}


## Metrics
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_2.RDS")
epa
sum(is.na(epa$Ozone2))
head(epa)

nasByStation <- aggregate(Ozone~Station.Code,epa,length)
nasByStation <- nasByStation[nasByStation$Ozone<366,]
nrow(nasByStation)
epa_sub <- epa[epa$Station.Code %in% nasByStation$Station.Code,]
#View(epa_sub)
apply(epa_sub[,c("Ozone","Ozone2")],2,function(x){sum(is.na(x))})
nrow(epa_sub)

(metrics <- evaluatePredictions(epa_sub$Ozone,epa_sub$Ozone2,3))


## Replace missing values in Ozone ######################################################
#epa$Ozone3 <- epa$Ozone
ii <- is.na(epa$Ozone)
epa$Ozone[ii] <- epa$Ozone2[ii]
#View(epa[epa$Station.Code %in% nasByStation$Station.Code,c("Ozone","Ozone2","Ozone3")])
#apply(epa,2,function(x){sum(is.na(x))})
epa <- epa[,!names(epa) %in% c("Flag","Ozone2","Ozone.trend")]
saveRDS(epa, "data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")

