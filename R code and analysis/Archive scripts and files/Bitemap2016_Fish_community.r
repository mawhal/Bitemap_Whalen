#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Fish seines - animal community analysis
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2017.10.05
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 
# Squid pops provide estimates of predation intensity inside and outside seagrass beds
# These data are, for most sites, paired with seines
# Usually, replicate squid pops and replicate (n=3) seines were conducted

#####################################################################

## UPDATES
# 2017.10.05 - modify prelimary analysis script to create the current script
# 2017.10.05 - fish + inverts community anlaysis using NMDS



###############################
## Analysis goals/strategies ##
###############################

# In this script, focus on seining data, analyze communites as multivariate data
# Need to ensure consistency in naming of species across sites??



#############################################################
#################################################
#######################################
###############################
###############################
#######################################
#################################################
#############################################################
########################################################################
#############################################################
#################################################
#######################################
###############################



# load libraries
# data.frame manipulation and joining
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
# geospatial data
library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(rfishbase)
library(taxize) # to look up taxa and generate classification tables
library(vegan)



###########################################
## DATA ON PREDATION AND FISH COMMUNITIES #
###########################################

# read in data

# FISH SEINING DATA
seines <- read.csv( '../Data/Bitemap_Seine_ALL-DATA_20170925.csv', strip.white = TRUE)
names(seines)[4] <- "habitat"
# convert all vegetated and unvegetated sites to common categories
seines$habitat[ seines$habitat %in% c("Seagrass " )] <- "Seagrass"
seines$habitat[ seines$habitat %in% c("unveg","Unveg","Unvegetated" )] <- "Unveg"
seines <- droplevels(seines)
# for all organisms from France, UNC2, Wales, QLD2, QLD3, multiple lengths by 10 to convert from cm to mm
seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] <- 
  seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] * 10
# combine genus and species names
seines$SPECIES_NAME <- with(seines, paste(Genus, Species) )



# several seines yielded zero animals
# for now change these numbers to zero, which will be reflected in the mean numbers for each site/country
sort(unique(seines$SPECIES_NAME))
seines[ which( is.na(seines$Genus) ), ]





###############################
#######################################
#################################################
#############################################################
########################################################################
#############################################################
#################################################
#######################################
###############################


###################################################################################
# EXTRACT BIO-ORACLE DATA                                                         #
###################################################################################

# This part of the script prepares Bio-ORACLE and WorldClim data for input sites,
# then merges these data with data from the input
# http://www.oracle.ugent.be/download.html (Shift-click to follow link)
# downloaded .rar data were extracted in Ubuntu using unrar function
# Note that BioOracle offers many useful oceanographic predictors

# # function to extract raster data within a radius defined by buffer
# buff <- function(rast,buffer=10000) {
#   unlist(lapply(extract( rast, input, buffer=buffer ),mean,na.rm=T))
# }


########### Bio-ORACLE ##################
# get a list of files to work with
files <- list.files( "C:/Users/mawha/Dropbox/Global Databases/BioOracle Data/",pattern = ".asc")
# # read in the raster data files (will assume lat,long for projection)
# r <- lapply( paste0("C:/Users/mawha/Dropbox/Global Databases/BioOracle Data/",files), raster )
# # crop all rasters to same extent, keeping southern hemisphere
# e <- extent(-180,180,-70,70)
# r2 <- lapply( r, function(rast) crop(rast,e) )

########### WorldClim Precipitation ###############
# p <- raster( "C:/Users/mawha/Dropbox/Global Databases/WorldClim/WorldClim_precip_2-5.tif" )

################## input data ######################
# extract lat and long from input
# input <- siteGPS[,c("meanLong","meanLat")]


############ extract environmental data #############
# # average values in all raster cells within a given radius from the GPS pointS
# # Use a "buffer" or radius over which to look for raster cells surrounding each GPS point
# # that has data (note there is no data on land or freshwater for Bio-ORACLE, and no ocean data for WorldClim
# # This will take a while for many variables
# buffer <- 10000  # this is in meters if the map is projected correctly
# # if raster coordinate reference system (CRS) is undefined, it will assume lat/long, which is correct in this case


# use the buff function (defined above) to all rasters
# oracle <- lapply( r2, buff ) # this takes a very long time, average Lat and Long by site first
#   write.csv( oracle, "Output Data/Bitemap_BioOracle.csv")
oracle <- read.csv( "Output Data/Bitemap_BioOracle.csv" )[,-1] 
# precip <- buff( p )
#   write.csv( precip, "Output Data/Bitemap_WorlClimPrecip.csv")
precip <- read.csv( "Output Data/Bitemap_WorlClimPrecip.csv" )[,-1]

# MISSING VALUES FROM Bio-ORACLE
# 7, 110, 130
# GPSjoin[ c(7,110,130), ]
# no data for 
# New South Wales, NSW1.1
# University of Delaware
# Croatia

# combine all of these into a data.frame and give them names
Environmentals <- data.frame( do.call( cbind, oracle ) )
names(Environmentals) <- unlist( strsplit( files, ".asc") )
Environmentals$precip <- precip

sites <- cbind(siteGPS,Environmentals)








#############################################################
### Fish and invertebrate communities
#############################################################
names(seines)

# convert each element of the list to a community matrix (rows as site, species as columns)
d2comm <- function(df){
  melted <- melt( df, id.vars = c(1:4,17,21), measure.vars = 16 )
  comm <- dcast( melted, Site.Name+Country+Lat+Long+habitat~SPECIES_NAME, 
                 fun.aggregate=sum, na.rm=TRUE  )

  return(comm)
}


scomm <- d2comm(seines)
# fill NA with zero
scomm[is.na(scomm)] <- 0

# mean of all seines within a country
mscomm <- aggregate( scomm[,c(3,4,6:ncol(scomm))], 
                     by=list(Country=scomm$Country,habitat=scomm$habitat), mean )


# remove columns with zero abundance
comm <- mscomm[ , which(colSums(mscomm[,-c(1:4)])>0)+4 ]
# # remove empty rows
# mscomm <- mscomm[ which(rowSums(mscomm[,-c(1:4)])>0)+4, ]

env  <- mscomm[, 1:4]

# NMDS
md <- metaMDS( comm, distance = 'bray' )
plot(md, type='t', display='sites')
plot(md, type='t', display='sites',
     xlim=c(610,611),ylim=c(-5,5)) # 38,10
plot(md, type='t', display='sites',
     xlim=c(-4,0),ylim=c(-0.5,0.5)) # all others

# Check out sites for outliers
env[c(10,38),]


