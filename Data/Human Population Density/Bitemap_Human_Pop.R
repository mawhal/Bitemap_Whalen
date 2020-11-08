#############################################################################
## Bitemap
## Population density estimates (proxy for human impact) from NASA
## http://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev10
## 22 Feb 2019
## by Matt Whalen
#############################################################################

# load libraries
library( raster )
library( rgdal )

# read geotiff info
GDALinfo( "gpw_v4_population_density_rev10_2010_15_min.tif" )

# read the geotiff
hp <- raster( "gpw_v4_population_density_rev10_2010_15_min.tif" )



# read site data to get Lats and Longs
sites <- read.csv( "../Data/Bitemap_sites.csv" )
xy <- sites[, c("meanLong","meanLat")]
xy <- SpatialPointsDataFrame(xy, xy)
crs(xy) <- crs(hp)


# get the cell values for each Lat Long -- might need a buffer around this
sites$pop   <- extract( hp, xy ) 
pop2        <- extract( hp, xy, buffer=30000 ) # buffer of 30 km selects multiple areas
sites$pop30 <- unlist(lapply( pop2, mean, na.rm=T ))
hist(log10(sites$pop))
hist(log10(sites$pop30))

# write to disk
write.csv( sites, "NASA_human_pop.csv", row.names=FALSE)
