library(openair)

# Output folder
FOLDER <- "data/aurn/"

# Only the IDs of the sites cotaining data between
# 2006 - 2015
site_ids <- read.csv("data/aurn/site_ids.csv") 

# Function that downloads/imports the data from AURN
get_data_for_site <- function(id){
  df <- importAURN(site = id, year = 2006:2015)
  write.csv(df, file=paste(FOLDER, id, ".csv", sep = ""), row.names = FALSE)
  #rm(df)
  #return(df)
}

#Test
#a <- get_data_for_site("ABD7")
#head(a)

apply(site_ids, 1, get_data_for_site)

