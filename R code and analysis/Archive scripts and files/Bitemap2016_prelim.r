#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Preliminary data analysis and exploration
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2017.09.21
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################

## UPDATES
# 2017.09.21: incorporate fish biomass coefficients to calculate mass from length
# 2017.09.26: finished first draft of code to compile length-weight relationships and calculate biomass
# 2017.09.27: construct fish size distribution figures


##################################
## Quetions, hypotheses, theory ##
##################################

# Big questions:
# How does predation in marine systems vary on a global scale?
# What controls (or at least is correlated with) predation rates across the world ocean?
# Does the presence (and quality??) of structured habitat influence the predator community and predation rates?

# Hypotheses for patterns:
# stronger predation at lower latitudes because of metabolism, patterns and history of diversification, ...
    # perhaps expectation is a normal-ish distribution of predation intensity across latitude, regardless of hemisphere?

# Latitudinal gradient in fish diversity known, but we expect to see this, too. 
# Does latitudinal diversity gradient hold at size classes and relevant to mesograzer control in seagrass communities?
# Similarly, is there a latitudinal gradient when traits are considered (e.g. only consider mesopredators)?

# Hypotheses for controls:
# environment: biogeography (includes evolutionary history, climate history), modern climate, biodiversity

# Interactions between predation rate and fish diversity. Possible that these are correlated, but can we know
# whether this potential relationship is causal or that both are responses to a different, shared driver?

########################################################################


###############################
## Analysis goals/strategies ##
###############################

# consider two datasets, one for predation assay (squidpops) and other for fish community observations (seines)
# some analysis on each independently becuase conducted at different scales, times, reps
# obtain summaries for each habitat type and place and merge the two datasets together 
  # look for statistical relationships between predation and fish community (and environment)

## Environmental data: for each site (mean of all latitudes/longitutde from each site)
# Source: Bio-ORACLE, 
# temperature: mean, max, min, range
# salinity, precipitation
# nutrients, productivity, sunlight

## Seine Data
# calculate different summaries of fish community that can be related to predation intensity
# summaries: total abundance, functional group abundance, richness, diversity metrics
#            fish biomass (need to calculate based on size and taxonomy)




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
library(ggplot2)
# geospatial data
library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(rfishbase)



###########################################
## DATA ON PREDATION AND FISH COMMUNITIES #
###########################################

# read in data
# SQUID POP PREDATION ASSAYS
pops   <- read.csv( '../Data/OceanBitemap_Squidpop_Data_20170904_MW.csv' )
# rename pops columns
names(pops) <- c("timeStamp","name","email","Lat","Long","dateRetrieved","timeRetrieved",
                 "dateDeployed","timeDeployed","habitat","habitatDescription","numberDeployed","missing1",
                 "numberRetrieved","missing24","Country","Institution","prop1","prop24","notes")
# convert vegetated and unvegetated sites to common categories
pops$habitat[ pops$habitat %in% c("Seagrass","Seagrass Meadow")] <- "Seagrass"
pops$habitat[ pops$habitat %in% c("Muddy Bottom","Sandy Bottom", "unvegetated ")] <- "unvegetated"
pops$habitat[ pops$habitat %in% c("Artificial Habitat (dock, breakwater, weir, etc.)","Rocky Reef")] <- NA
pops <- pops[ !is.na(pops$habitat), ]
pops <- droplevels(pops)
levels(pops$habitat) <- c("Seagrass","Unveg")
# coerce proportion after 24 hours to be a numbers
pops$prop24[ pops$prop24=="#DIV/0!" ] <- NA
pops$prop24 <- as.numeric(as.character(pops$prop24))

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


# merging data
# this is difficult because the replicate seines do not pair up with squidpop reps
# so, in order to merge the two we need to use summaries at the site level
# try to match things based on Lat/Long info
popGPS   <- unique( pops[,c('Lat','Long','Institution','Country')])
seineGPS <- unique( seines[,c('Lat','Long','Country','Site.Name','habitat')])
# compare sites used in seining and predation assays
GPSjoin <- full_join( popGPS,seineGPS )
GPSjoin <- GPSjoin[ with(GPSjoin, order(Lat,Long)), ]
# write.csv( GPSjoin, "Output Data/pops_sienes_GPSmatch.csv", row.names=F )
# to complete merge, we need to average predation and fish metrics at site level

# combine all squidpop and seine replicates from a given location (e.g. NSW2)
siteGPS <- ddply( GPSjoin, .(Country), 
                  summarise, meanLat=mean(Lat), meanLong=mean(Long) )








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


## COULDN'T GET THIS TO WORK
# # Instead of extract with buffer, use velox extract, which requires spatial polygons as the input
# # make spatial polygons as circles of a given radius around a point
# library(dismo) # circles function
# ?circles
# inputc <- circles( input, 10000, lonlat=TRUE, dissolve=TRUE )
# inputc@polygons
# # list of circles
# inputcl <- list()
# for(i in 1:nrow(input) ){
#   inputcl[[i]] <- gBuffer( input[i,], lonlat=TRUE ) 
# }
#   apply( input, 1, function(z) circles( (z),10000,lonlat=TRUE) )
# class(inputc)
# geometry(inputc)
# plot(geometry(inputc))
# 
# v <- velox(p)
# it messes up at this stage, Error stating "Extent is non-positive in at least one dimension
# v$extract(sp=inputcl[[1]]@polygons, fun=mean) 








####################################
# Predation patterns and summaries #
####################################



# combine pops and sites
pop.site <- left_join( pops, sites )

# patterns with latitude
ggplot( pop.site, aes(x=meanLat,y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,formula = y~poly(x,2),method.args=list(family=quasibinomial))
ggplot( pop.site, aes(x=abs(meanLat),y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='lm',se=F)
ggplot( pop.site, aes(x=abs(meanLat),y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,method.args=list(family=quasibinomial))
ggplot( pop.site, aes(x=abs(meanLat),y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,method.args=list(family=quasibinomial)) + facet_wrap(~habitat)

# sites at higher Longitude tend to be at lower latitude
ggplot( pop.site, aes(x=(meanLong),y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,method.args=list(family=quasibinomial)) + facet_wrap(~habitat)


# temperature and latitude
ggplot( pop.site, aes(x=abs(meanLat),y=sstmean,col=habitat) ) + geom_point() + 
  geom_smooth(method='lm',se=T) + facet_wrap(~habitat)
ggplot( pop.site, aes(x=abs(meanLat),y=sstrange,col=habitat) ) + geom_point() + 
  geom_smooth(method='lm',se=T) + facet_wrap(~habitat)

# predation and temperature
ggplot( pop.site, aes(x=sstmean,y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,method.args=list(family=quasibinomial)) + facet_wrap(~habitat)

ggplot( pop.site, aes(x=sstrange,y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,method.args=list(family=quasibinomial)) + facet_wrap(~habitat)

ggplot( pop.site, aes(x=salinity,y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,method.args=list(family=quasibinomial)) + facet_wrap(~habitat)



# split by ocean basin


# average predation rate by location
pop.mean <- ddply( pop.site, .(Country,meanLat,meanLong,habitat,sstmean), 
                   summarise, prop24=mean(prop24) )
ggplot( pop.mean, aes(x=abs(meanLat),y=prop24,col=habitat) ) + geom_point() + 
  geom_smooth(method='glm',se=T,method.args=list(family=quasibinomial)) + facet_wrap(~habitat)


# 










###########################################################################################################
######                                       FISH SUMMARIES                                           #####
###########################################################################################################



##########################
# Fish size distribution #
##########################


# prob first need to omit rows lacking a size, then use rep() to make a vector of lengths for for each count


# remove inverts
fish <- seines[ seines$Phylum=="Chordata" & !is.na(seines$Phylum), ]

# only include rows with a length estimate (fish biomass assumed average lengths for fishes not measured...see below)
fish.size <- fish[ !is.na(fish$Length), ]

# use ggplot to make violin plots for each site x habitat combination
ggplot( fish.size, aes(x=habitat, y=log10(Length) )) + geom_boxplot(width=0.5) + facet_wrap(~Country)
ggplot( fish.size, aes(x=habitat, y=log10(Length) )) + geom_violin() + geom_boxplot(width=0.5)+ facet_wrap(~Country)


# fish length by latitude
ggplot( fish.size, aes(x=Lat, y=log10(Length), group=Country )) + geom_boxplot(width=0.5)
ggplot( fish.size, aes(x=abs(Lat), y=log10(Length), col=habitat, group=Country )) + geom_boxplot(width=0.5)


##################
# Fish diversity #
##################


# combine seines and sites
seine.site <- left_join( seines, siteGPS )

# summarize seine data based on counts of each taxon
seine.abund <- ddply( seine.site, .(Country,meanLat,meanLong,habitat,Phylum,Genus,Species), 
                      summarise, Abundance=sum(Abundance) )

# for species diversity metrics, we need counts of each species in each habitat
# omit invertebrates
fish.abund <- seine.abund[ seine.abund$Phylum=="Chordata", ]

# consider using different diversity metrics (e.g. ENSpie as in Chase & Knight 2013)
ENSPIE <- function(prop)   1 / sum(prop^2)
prop <- c(0.5,0.4,0.1)
ENSPIE(prop)
# ENSPIE is less scale dependent than many other diversity metrics (arguable important here), 
# but tends to be sensitive to aggregation (which can be affected by spatial scale, 
# and likely an issue for schooling fish and seining)

# convert abundance to proportion
# calculate total abundance for each site (lump all replicate seines together)
fish.totals <- ddply( fish.abund, .(Country,meanLat,meanLong,habitat),
       summarize, Total = sum(Abundance))
# get rid of NA rows
fish.totals <- fish.totals[!is.na(fish.totals$Country),]

# include totals as a column in fish.abund
fish.comm <- left_join(fish.abund,fish.totals)

# calculate the proportion of each species in each site
fish.comm$prop <- fish.comm$Abundance/fish.comm$Total

# calculate ENSPIE for each site
ENSPIE.site <- ddply( fish.comm, .(Country,meanLat,meanLong,habitat), 
       summarize, ENSPIE = ENSPIE(prop))
# get rid of NA rows
ENSPIE.site <- ENSPIE.site[!is.na(ENSPIE.site$Country),]

# calculate effect size as difference between ENSPIE estimates
ENSPIE.site <- ENSPIE.site[ with(ENSPIE.site,order(Country,habitat)), ]

# make all Unveg sites negative
ENSPIE.site$ENSPIEdiff <- ENSPIE.site$ENSPIE
ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] <- -(ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] )
ENSPIE.diff <- ddply( ENSPIE.site, .(Country,meanLat,meanLong), summarise, ENSPIE=sum(ENSPIEdiff) )



ggplot( ENSPIE.site, aes(x=habitat,y=ENSPIE)) + geom_point() + 
  geom_line(aes(group=Country)) + facet_wrap(~Country)
ggplot( ENSPIE.site, aes(x=habitat,y=ENSPIE)) + geom_boxplot(notch=F)
ggplot( ENSPIE.site, aes(x=abs(meanLat),y=ENSPIE)) + geom_point() + facet_wrap(~habitat, ncol=1)
ggplot( ENSPIE.diff, aes(x=abs(meanLat),y=ENSPIE)) + geom_point()



# merge average squidpops and seines by location and habitat
popseine <- left_join(pop.mean,ENSPIE.site)
popseine <- left_join(popseine,fish.totals) 
ggplot( popseine, aes(x=ENSPIE, y=prop24)) + geom_point()
ggplot( popseine, aes(x=Total, y=prop24)) + geom_point()





###################################################################################
# FISH BIOMASS                                                                    #
###################################################################################


# remove inverts
fish <- seines[ seines$Phylum=="Chordata" & !is.na(seines$Phylum), ]

# What to do with rows where we have fish abundance but no length info?
# options: ignore these rows (underestimate biomass, potentially substantially)
        #  use an average size based on other length estimates for biomass (unclear if under- or overestimate)
        #  use largest or smallest lengths
# for now, use average (median) size for all other length estimates

# find rows with NA for size, but omit rows where nothing was found
fishNA <- fish[ is.na(fish$Length) & fish$Abundance>0, ]

# calculate median sizes for all other fish species from each seine
fishL  <- fish[ !is.na(fish$Length) & fish$Abundance>0, ]
aveLengths <- ddply( fishL, .(Site.Name,habitat,Country,Date,Genus,Species), 
                         summarize, medLength=median(Length), meanLength=mean(Length) )
# remove NA values for unknown Genera
aveLengths[is.na(aveLengths$Genus),]
aveLengths <- aveLengths[!is.na(aveLengths$Genus),]

# replace NA length values with calculated median
# rename fishNA as fishEST for adding length ESTimates
fishEST <- fishNA
# loop over all rows and get median lengths
for( i in 1:nrow(fishEST) ){
  fishEST$Length[i] <- aveLengths$meanLength[ aveLengths$Country==fishEST$Country[i] &
                                              aveLengths$habitat==fishEST$habitat[i] &  
                                              aveLengths$Site.Name==fishEST$Site.Name[i] & 
                                              aveLengths$Date==fishEST$Date[i] &
                                              aveLengths$Genus==fishEST$Genus[i] &
                                              aveLengths$Species==fishEST$Species[i] ]
}
### errors from this step that were 'corrected'
## Syngnathus from NC, no other pipefish were caught but this one not measured, 
  # so use median of all syngnathus captured from that site
## Ambassis jacksoniensis from NSW2 20170223 + 201702024 not measured to maximize survival. 
  # Notes state that all were between 15 and 60 mm, so using mean of these values (37.5mm)

# merge size estimates back with original fish data
fish.clean <- full_join( fishL, fishEST )

# combine genus and species with name to match length-weight coefficients (see below)
fish.clean$SPECIES_NAME <- with(fish.clean, paste(Genus, Species) )



# read biomass coefficients
coef <- read.csv("../Data/Fish Biomass/20160710_RLS_biomass_coefs_20170921.csv")
# omit rows from coef for which we have no estimates
coef <- coef[ !is.na(coef$A), ]
# remove duplicates from coef
coef <- coef[ !duplicated(coef), ]



# identify species not in character list
lookup <- sort(unique(fish.clean$SPECIES_NAME[!(fish.clean$SPECIES_NAME %in% coef$SPECIES_NAME)]))
sort(unique(fish.clean$SPECIES_NAME[(fish.clean$SPECIES_NAME %in% coef$SPECIES_NAME)]))

# # use rfishbase to look up taxa not represented in the reef life survey biomass conversion data.frame (coef)
# fishbaseLW <- length_weight(lookup)
# sort(unique(fishbaseLW$sciname)) # added 67 taxa
# unrepresented <- sort(unique(fish.clean$SPECIES_NAME[!(fish.clean$SPECIES_NAME %in% c(fishbaseLW$sciname, as.character(coef$SPECIES_NAME)) )]))
# write.csv( fishbaseLW, "../Data/Fish Biomass/fishbaseLW.csv", row.names = FALSE)
# write.csv( data.frame(unrepresented), "../Data/Fish Biomass/fishbaseMISSING.csv", row.names = FALSE)

# read the length-weight relationships looked up in fishbase
fishbaseLW <- read.csv( "../Data/Fish Biomass/fishbaseLW_20170921.csv" )
# only keep rows that were chosen from the database (another option would be to average all estimates by taxon)
fishbaseLW <- fishbaseLW[ fishbaseLW$Keep == "Y", ]
# average by taxon
fishbaseLW <- ddply( fishbaseLW, .(sciname), 
                     summarize, a=mean(a),aTL=mean(aTL),b=mean(b) )
# accept the estimates for a (intercept) that account for standard vs total length
for(i in 1:nrow(fishbaseLW) ){
  if( !is.na(fishbaseLW$aTL[i]) ) fishbaseLW$a[i] <- fishbaseLW$aTL[i]
}
# rename sciname to Species
names(fishbaseLW) <- c( "SPECIES_NAME", "A", "aTL", "B" )

# read in table for taxa that were looked up manually in fishbase
fishbaseManual <- read.csv( "../Data/Fish Biomass/Bitemap_FishbaseLW_manual.csv" )[,1:3]
names(fishbaseManual) <- c( "SPECIES_NAME", "A", "B" )

# combine all data from fishbase, including those provided by Reef Life Survey
fishbaseFull <- full_join( fishbaseLW, fishbaseManual )
biomassCoef <- full_join(  coef, fishbaseFull )
biomassCoef <- biomassCoef %>%
  dplyr::select( SPECIES_NAME, A, B )
length(unique(biomassCoef$SPECIES_NAME))
# omit duplicated rows (separate estimates for a taxon), default to use Reef Life Survey estimates for consistency
biomassCoef <- biomassCoef[ !duplicated(biomassCoef$SPECIES_NAME), ]

# match seine data with biomass conversion coefficients
fish_biom <- left_join(fish.clean, biomassCoef, by="SPECIES_NAME" )

# remove extraneous columns
fish_biom <- fish_biom %>%
  dplyr::select(Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, Q.PA, Phylum, Genus, Species, 
         Length, Abundance, Country, SPECIES_NAME, A, B )


# which taxa don't have length-weight regression estimates
miss <- fish_biom[ is.na(fish_biom$A), ]
sort(unique(miss$SPECIES_NAME))

# create biomass column -- note that I've divided all lengths by 10 to convert mm to cm, which is the unit 
fish_biom$biomass <- fish_biom$Abundance*(fish_biom$A*(fish_biom$Length/10)^fish_biom$B)


# make biomass into kg
fish_biom$biomass <- fish_biom$biomass*0.001

# Calculate total fish biomass in each seine
biomass_total <- ddply( fish_biom, .(Site.Name,Country,Date,Lat,Long,habitat,Depth), 
                        summarize, biomass=sum(biomass, na.rm=TRUE) )



# graph biomass by habitat type
ggplot(biomass_total, aes(x=habitat, y=log10(biomass), fill = habitat)) + geom_boxplot(notch=TRUE) +
  geom_point(alpha=0.1) +
  # facet_grid(.~Year) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Biomass of All Fish By Habitat Type",
       x = "Type of Habitat", y="Biomass") 

# total biomass by Latitude
ggplot(biomass_total, aes(x=Lat, y=log10(biomass))) + geom_point() + facet_wrap(~habitat,ncol=1)
ggplot(biomass_total, aes(x=abs(Lat), y=log10(biomass))) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='lm')


# calculate average biomass by seine and habitat
biomass_mean <- ddply( fish_biom, .(Site.Name,Country,Date,Lat,Long,habitat,Depth), 
                       summarize, biomass=mean(biomass, na.rm=TRUE) )
ggplot(biomass_mean, aes(x=abs(Lat), y=log10(biomass))) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='lm')


# calculate average biomass by site and habitat
biomass_site <- ddply( fish_biom, .(Country,habitat), 
                       summarize, biomass=mean(biomass, na.rm=TRUE) )









#########################
# Fish Trophic Levels   #
#########################

fishSpecies <- sort( unique( fish.clean$SPECIES_NAME ) )

fishEcology <- ecology( fishSpecies, fields = c("SpecCode", "FeedingType", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph") )

# sometimes DietTroph (based on whole range of food items) unavailable but FoodTroph (based on individual food items) more often is available
fishEcology$Trophic <- fishEcology$DietTroph
for( i in 1:nrow(fishEcology) ) {
  if( is.na(fishEcology$DietTroph[i]) )   fishEcology$Trophic[i] <- fishEcology$FoodTroph[i]
}

# which taxa are missing trophic levels?
trophMiss <- fishEcology[ which(is.na(fishEcology$Trophic)), ]
# write.csv( trophMiss, "../Data/Fish Biomass/Bitemap_2017_fishBase_missingTrophic.csv", row.names = FALSE )
# read in data.frame with info added from fishbase
trophFill <- read.csv( "../Data/Fish Biomass/Bitemap_2017_fishBase_missingTrophic.csv" )
fishEcology <- full_join( fishEcology, trophFill[,c('sciname','fishbaseTrophEst')] )
for( i in 1:nrow(fishEcology) ) {
  if( is.na(fishEcology$Trophic[i]) )   fishEcology$Trophic[i] <- fishEcology$fishbaseTrophEst[i]
}


fishSpecies[ !(fishSpecies %in% fishEcology$sciname)]
fishbaseMiss <- write.csv( data.frame(fishSpecies), "../Data/Fish Biomass/Bitemap_2017_fishBase_noEcology.csv", row.names=FALSE )

#













############################################################
###       COMBINE PREDATION DATA WITH FISH DATA          ###
############################################################


## combine squid pops, fish diversity, fish biomass
popseine <- left_join(popseine,biomass_site) 

# Summaries of Countries, Latitudes, etc.
# How many total countries?
sort( unique(popseine$Country) )
# range of Latitudes
diff(range( popseine$meanLat ))

# Which sites conducted both squidpops and seining?
popseine_both <- popseine[ which( apply(!is.na(popseine[, c("prop24","Total") ]),1,all) ), ]
sort( unique(popseine_both$Country) )

# FIGURES

ggplot(popseine, aes(x=abs(meanLat), y=sstmean)) + geom_point() + geom_smooth(method='lm')

ggplot(popseine, aes(x=abs(meanLat), y=log10(biomass))) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='lm')

ggplot(popseine, aes(x=abs(meanLat), y=log10(biomass))) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='lm')

ggplot(popseine, aes(x=abs(meanLat), y=prop24 )) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='glm', formula=y~poly(x,2), method.args=list(family=binomial) ) + ylab('Predation intensity (prop squid missing)') +
  xlab('Latitude (absolute value)')

ggplot(popseine, aes(x=sstmean, y=prop24 )) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='glm', method.args=list(family=binomial) ) + ylab('Predation intensity (prop squid missing)') +
  xlab('Annual mean SST')

ggplot(popseine, aes(x=sstmean, y=prop24 )) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='glm', formula=y~poly(x,2),method.args=list(family=binomial) ) + ylab('Predation intensity (prop squid missing)') +
  xlab('Annual mean SST')

# with Country labels
# make country labels smaller
country_split <- strsplit( popseine$Country, split = '[()]' )
popseine$site_label <- unlist(lapply( country_split, function(z){
  if(length(z)==2) z[2] else z[1]
}
))

ggplot(popseine, aes(x=sstmean, y=prop24, label=site_label )) + geom_text(size=2.5) + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='glm', formula=y~poly(x,2),method.args=list(family=binomial) ) + ylab('Predation intensity (prop squid missing)') +
  xlab('Annual mean SST')

ggplot(popseine, aes(x=log10(biomass), y=prop24 )) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='glm', method.args=list(family=binomial) ) + ylab('Predation intensity (prop squid missing)') +
  xlab('log10(fish biomass)')

ggplot(popseine, aes(x=log10(Total), y=prop24 )) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='glm', method.args=list(family=binomial) ) + ylab('Predation intensity (prop squid missing)') +
  xlab('log10(total fish abundance)')

ggplot(popseine, aes(x=ENSPIE, y=prop24 )) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='glm', method.args=list(family=binomial) ) + ylab('Predation intensity (prop squid missing)') +
  xlab('Fish diversity (Effective no. species PIE)')






###########################################################################
#####    STATISTICAL MODELS
#########################################################################
library(bbmle)

glm1 <- glm( prop24~habitat, popseine, family='quasibinomial' )
summary(glm1)

glm2 <- glm( prop24~abs(meanLat), popseine, family='quasibinomial' )
summary(glm2)

glm2 <- glm( prop24~sstmean, popseine, family='quasibinomial' )
summary(glm2)

glm2.5 <- glm( prop24~sstmean+I(sstmean^2), popseine, family='quasibinomial' )
summary(glm2.5)

  # glm2 <- glm( prop24~sstmean, popseine, family='binomial' )
  # glm2.5 <- glm( prop24~sstmean+I(sstmean^2), popseine, family='binomial' )
  # AICctab( glm2, glm2.5, nobs=nrow(popseine)-4 )

glm3 <- glm( prop24~ENSPIE*habitat, popseine, family='quasibinomial' )
summary(glm3)

glm4 <- glm( prop24~log10(Total), popseine, family='quasibinomial' )
summary(glm4)

glm5 <- glm( prop24~sstmean+log10(Total), popseine, family='quasibinomial' )
summary(glm5)

glm6 <- glm( prop24~sstmean+ENSPIE*habitat, popseine, family='quasibinomial' )
summary(glm6)

#############################################################################################################
# NOTES + QUESTIONS
# which seagrass species are we dealing with? Potentially a useful predictor?
# many sites did not seine so some comparisons will use a smaller subset of data
# to strengthen analysis of latitudinal trends on predation (in seagrass only) can we include ZEN data?
#############################################################################################################
