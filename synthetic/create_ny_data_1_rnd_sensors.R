library(sp)
library(maptools)
library(spTimer)
library(maps)


############################################################
# Location of the reference stations
data("NYdata")
coords <- as.matrix(unique(cbind(NYdata[, 2:3])))
map(database = "state", regions = "new york")
points(coords, pch = 19, col = 3)
points(coords, pch = 1, col = 1)
par(.pardefault)
readline("Continue")

############################################################
# (Random) location of the cheap-stations
#dev.new()
nyMap <- map(database = "state", regions = "new york", fill = TRUE)
#dev.off()
IDs <- sapply(strsplit(nyMap$names, ":"), function(x) x[1])
nyPoly <- map2SpatialPolygons(nyMap, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
#summary(nyPoly)
#plot(nyPoly)

nPoints <- 50
# plot(nyPoly)
# points(spsample(nyPoly, n = nPoints, "regular"), pch = 3)
# plot(nyPoly)
# points(spsample(nyPoly, n = nPoints, "random"), pch = 3)
# plot(nyPoly)
# points(spsample(nyPoly, n = nPoints, "stratified"), pch = 3)
# plot(nyPoly)
# points(spsample(nyPoly, n = nPoints, "nonaligned"), pch = 3)

set.seed(999)
plot(nyPoly)
lcsCoords <- spsample(nyPoly, n = nPoints, "nonaligned")
points(lcsCoords, pch = 3)
save(lcsCoords, file="data/ny_ozone/lcsCoords.Rdata")

# References (again)
points(coords, pch = 19, col = 3)
points(coords, pch = 1, col = 1)

