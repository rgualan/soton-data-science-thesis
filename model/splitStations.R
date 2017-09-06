## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(maps)
source("util/my_helper.R")

## Settings ##########################################################################
paper <- F

## Read data and split #########################################################################
epa <- readEpaDataset(1)
sites <- getSites(epa)
set.seed(204)
folds <- cut(sample(1:nrow(sites)),breaks=10,labels=F)
saveRDS(folds, file="data/folds.RDS")

## Read data and split #########################################################################
epa <- readEpaDataset(2)
sites <- getSites(epa)
set.seed(204)
folds <- cut(sample(1:nrow(sites)),breaks=10,labels=F)
saveRDS(folds, file="data/folds2.RDS")

## Split ##################################################################
#set.seed(204)
#folds <- cut(sample(1:nrow(sites)),breaks=10,labels=F)
#saveRDS(folds, file="output/folds.RDS")
#saveRDS(folds, file="output/folds2.RDS")
folds <- getFolds()
length(folds)

## Simple splitting ##################################################################
sites$Test <- "Training"
sites$Test[folds==1] <- "Validation"

if(paper) jpeg("img/benchmark_splitting.jpeg", 6, 6, "in", bg="white", res=150)
ggplot(sites) +
  geom_polygon(data=getCAmap(), aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_point(aes(x = Longitude, y = Latitude, colour=Test), 
             alpha = .75, size=2, shape=17) + 
  labs(x = "Longitude", y = "Latitude", colour="Type") +
  coord_quickmap() + 
  theme(legend.justification = c("right", "top"), legend.position = c(.95, .95),
        legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
if(paper) dev.off()


