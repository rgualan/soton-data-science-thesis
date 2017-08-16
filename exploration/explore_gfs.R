# Clean environment #################################################################
rm(list=ls())

# Libraries
library(raster)
library(rasterVis)
library(maptools)
library(maps)
library(RColorBrewer)
library(lattice)
library(ncdf4)


# NEtcdf file
nc <- nc_open("data/e_obs_grid/tg_0.25deg_reg_v15.0-2015-box.nc")
print(nc)


# read the netCDF file as a raster layer
r <- raster("data/e_obs_grid/tg_0.25deg_reg_v15.0-2015-box.nc")
#r
#print(r)

# map the data
world.outlines <- map("world", regions="France", plot=F)
world.outlines.sp <- map2SpatialLines(world.outlines, proj4string = CRS("+proj=longlat"))

mapTheme <- rasterTheme(region = rev(brewer.pal(10, "RdBu")))
#cutpts <- c(-2.5, -2.0, -1.5, -1, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5)
plt <- levelplot(r, margin = F, cuts=10, pretty=TRUE, par.settings = mapTheme,
                 main="Layer from NetCDF file") #at=cutpts, 
plt + layer(sp.lines(world.outlines.sp, col = "black", lwd = 0.5))


# Example of data extraction
#extract(r, matrix(c(0.125,43.875),1,2), method='bilinear', fun=mean, na.rm=TRUE)
#extract(r, matrix(c(0.125*2,43.875+0.125),1,2), method='bilinear', fun=mean, na.rm=TRUE)
#extract(r, matrix(c(-1.625,44.375),1,2), method='bilinear', fun=mean, na.rm=TRUE)
