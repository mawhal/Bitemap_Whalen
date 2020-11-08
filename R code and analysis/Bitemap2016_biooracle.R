###################################################################################
# EXTRACT BIO-ORACLE DATA                                                         #
###################################################################################

# This script prepares Bio-ORACLE and WorldClim data for input sites,
# then merges these data with data from the input
# http://www.oracle.ugent.be/download.html (Shift-click to follow link)
# downloaded .rar data were extracted in Ubuntu using unrar function
# Note that BioOracle offers many useful oceanographic predictors
library(raster)
# function to extract raster data within a radius defined by buffer
buff <- function(rast,buffer=10000) {
  unlist(lapply(extract( rast, input, buffer=buffer ),mean,na.rm=T))
}

########### Bio-ORACLE ##################
# get a list of files to work with
files <- list.files( "../Global Databases/BioOracle Data/",pattern = ".asc")
# read in the raster data files (will assume lat,long for projection)
r <- lapply( paste0("../Global Databases/BioOracle Data/",files), raster )
# crop all rasters to same extent, keeping southern hemisphere
e <- extent(-180,180,-70,70)
r2 <- lapply( r, function(rast) crop(rast,e) )

########### WorldClim Precipitation ###############
p <- raster( "C:/Users/mawha/Dropbox/Global Databases/WorldClim/WorldClim_precip_2-5.tif" )

################## input data ######################
# extract lat and long from input
siteGPS <- read.csv( "Data/raw/Bitemap_sites.csv" )
input <- siteGPS[,c("meanLong","meanLat")]


############ extract environmental data #############
# average values in all raster cells within a given radius from the GPS pointS
# Use a "buffer" or radius over which to look for raster cells surrounding each GPS point
# that has data (note there is no data on land or freshwater for Bio-ORACLE, and no ocean data for WorldClim
# This will take a while for many variables
buffer <- 10000  # this is in meters if the map is projected correctly
# if raster coordinate reference system (CRS) is undefined, it will assume lat/long, which is correct in this case

# use the buff function (defined above) to all rasters
oracle <- lapply( r2, buff ) # this takes a very long time, average Lat and Long by site first
precip <- buff( p )

# combine all of these into a data.frame and give them names
Environmentals <- data.frame( do.call( cbind, oracle ) )
names(Environmentals) <- unlist( strsplit( files, ".asc") )
# Environmentals$precip <- precip

sites <- cbind(siteGPS,Environmentals)
# factor for hemisphere
sites$hemi <- factor(ifelse( sites$meanLat>0,"North","South"))

# write Environmentals to disk
write.csv( sites, "Data/processed/Bitemap_BioORACLE_20201107.csv", row.names = FALSE )
