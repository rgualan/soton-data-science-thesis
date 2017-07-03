library(openair)
library(reshape2)

# Read the data
co <- read.csv('data/epa_hourly/2017/hourly_42101_2017_co.csv',header=T)
pm10 <- read.csv('data/epa_hourly/2017/hourly_81102_2017_pm10.csv',header=T)
so2 <- read.csv('data/epa_hourly/2017/hourly_42401_2017_so2.csv',header=T)
no2 <- read.csv('data/epa_hourly/2017/hourly_42602_2017_no2.csv',header=T)
ozone <- read.csv('data/epa_hourly/2017/hourly_44201_2017_ozone.csv',header=T)
temp <- read.csv('data/epa_hourly/2017/hourly_TEMP_2017.csv',header=T)
wind <- read.csv('data/epa_hourly/2017/hourly_WIND_2017.csv',header=T)


# Parse datetime
parseDate <- function(df){ 
  dateColumn <- as.POSIXct(strptime( paste(df$Date.Local,df$Time.Local),
                       format="%Y-%m-%d %H:%M", tz="GMT"))
  return(dateColumn)
}

co$date <- parseDate(co)
pm10$date <- parseDate(pm10)
so2$date <- parseDate(so2)
no2$date <- parseDate(no2)
ozone$date <- parseDate(ozone)
temp$date <- parseDate(temp)
wind$date <- parseDate(wind)

# Filter relevant columns only
FILTER = c(-5, -9, -10:-13, -15:-21, -24)
co_2 = co[,FILTER]
pm10_2 = pm10[,FILTER]
so2_2 = so2[,FILTER]
no2_2 <- no2[,FILTER]
ozone_2 <- ozone[,FILTER]
temp_2 <- temp[,FILTER]
wind_2 <- wind[,FILTER]


# Concatenate by rows
tmp = rbind(co_2, pm10_2, so2_2, no2_2, ozone_2, temp_2, wind_2)
# Create unified site id
tmp$site <- sprintf("%02d-%03d-%04d", tmp$State.Code, tmp$County.Code, tmp$Site.Num)


## Test non-unique dates by state
# for( state in unique(tmp$State.Code)){
#   print(state)
#   a <- dcast(tmp[tmp$State.Code==state,],  ...~Parameter.Code, value.var="Sample.Measurement")
# } 
# 
# a <- tmp[tmp$State.Code==15,]
# x <- aggregate(Sample.Measurement~site+date+Parameter.Code, data = a, FUN=length)
# x[x$Sample.Measurement>1,]
# 
# tmp2 = rbind(co, pm10, so2)
# tmp2$site <- sprintf("%02d-%03d-%04d", tmp2$State.Code, tmp2$County.Code, tmp2$Site.Num)
# x <- tmp2[tmp2$site=="15-003-0010" & tmp2$Parameter.Code=="42101",]
# head(x)
# unique(x$POC)
# a <- x[x$POC==2,]
# b <- x[x$POC==3,]
# plot(a$date, a$Sample.Measurement, type="l")
# lines(b$date, b$Sample.Measurement, col=2)
# lines(a$date, (a$Sample.Measurement + b$Sample.Measurement)/2, col=3)


# Reshape data from single column to a column by variable
# TODO: Check the effect of the aggregation for variables with two or more POC
df <- dcast(tmp,  ...~Parameter.Code, value.var="Sample.Measurement", 
            fun.aggregate = mean, na.rm=TRUE) 
names(df)[which(names(df)=="42101")] <- "co"  
names(df)[which(names(df)=="81102")] <- "pm10" 
names(df)[which(names(df)=="42401")] <- "so2" 
names(df)[which(names(df)=="44201")] <- "ozone"
names(df)[which(names(df)=="42602")] <- "no2"
names(df)[which(names(df)=="62101")] <- "temp"
names(df)[which(names(df)=="61103")] <- "ws"
names(df)[which(names(df)=="61104")] <- "wd"

df <- df[df$date<"2017-05-01 00:00:00 GMT",]
epa_complete <- df
save(epa_complete, file="data/epa_hourly/2017/epa_complete.Rda")


# Plots
summaryPlot(df[,c(-1:-8)], period = "months")

# Count available data by site 
# dataBySite  <- aggregate(cbind(co,pm10,so2,ozone,no2,temp,ws,wd)~site, data=df, FUN=
#                            function(x) sum(!is.na(x)), na.action = na.pass)
dataBySite  <- aggregate(df[c('co','pm10','so2','ozone','no2','temp','ws','wd')], 
                         by=df['site'], 
                         FUN=function(x) sum(!is.na(x)))
View(dataBySite)

colSums(dataBySite[-1]>0)

# Check one good station
complete <- df[df$site=="22-033-0009",] #  
summaryPlot(complete[,c(9:18)], period = 'months')
# Check an incomplete station
incomplete <- df[df$site=="04-001-1003",]
summaryPlot(incomplete[,c(9:18)], period = 'months')
# Check one good station
other <- df[df$site=="04-013-9997",] # 4133002 
summaryPlot(other[,c(9:18)], period = 'months')


# Filter dataset to containing a minimum of data per variable 
MIN_ROWS = 1500#1300
chosenSites <- dataBySite[dataBySite$pm10>1500
                          & dataBySite$so2>MIN_ROWS
                          & dataBySite$ozone>MIN_ROWS
                          & dataBySite$no2>MIN_ROWS
                          & dataBySite$temp>MIN_ROWS
                          & dataBySite$ws>MIN_ROWS
                          & dataBySite$wd>MIN_ROWS, 1]
filteredDf <- df[df$site %in% chosenSites, ]
epa <- filteredDf

save(epa, file="data/epa_hourly/2017/epa_filtered.Rda")

#data("NYdata")
load("data/epa_hourly/2017/epa_filtered.Rda")
load("data/epa_hourly/2017/epa_complete.Rda")

s <- unique(epa$site)
DataFit <- spT.subset(data = epa, var.name = c("site"), s = s, reverse = FALSE)
#DataFit <- subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
#DataValPred <- spT.subset(data = NYdata, var.name = c("s.index"), s = s)
#DataValPred <- subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 
#                                                         8)))


# Create map with all the stations
allCoords <- as.matrix(unique(cbind(epa_complete$Longitude, epa_complete$Latitude)))
map('usa')
points(allCoords, pch = 1, col = 1)


# Create map with stations
coords <- as.matrix(unique(cbind(DataFit$Longitude, DataFit$Latitude)))
#pred.coords <- as.matrix(unique(cbind(DataValPred[, 2:3])))
# dev.new()
#map(database = "state", regions = c("arizona"))
map('usa')
#map("state", ".*")
points(coords, pch = 19, col = 3)
points(coords, pch = 1, col = 1)
#points(pred.coords, pch = 3, col = 4)
# legend(x = -77.5, y = 41.5, col = c(3, 4), pch = c(19, 3), cex = 0.8, legend = c("Fitted sites", 
#                                                                                  "Validation sites"))
# dev.off()


# Fit GP model
set.seed(11)
post.gp <- spT.Gibbs(formula = pm10 ~ co + so2, data = DataFit, model = "GP", 
                     coords = ~Longitude + Latitude, scale.transform = "SQRT", 
                     spatial.decay = spT.decay(distribution = Gamm(2,1), tuning = 0.1))
print(post.gp)
summary(post.gp)
# plot(post.gp)




