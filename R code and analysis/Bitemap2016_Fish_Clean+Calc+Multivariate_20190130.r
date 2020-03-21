#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Cleaning fish community data, calculating summaries
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2019.02.19
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################

## UPDATES FOR THIS SCRIPT
# 2018.10.15: Matt adds Bitemap project to GitHub. Futher changes tracked there, and old changes tracked in archived script files


## Seine Data
# calculate different summaries of fish community that can be related to predation intensity
# summaries: total abundance, functional group abundance, richness, diversity metrics
#            fish biomass (need to calculate based on size and taxonomy)

########################################################################################
###  Seines were pulled for different distances. We take this into account,
###  but we still lack information on how wide seines were!
########################################################################################


## Video Data
# only use video data to identify which fish families were present


########################################################################################
########################################################################################





# load libraries
# data.frame manipulation and joining
library(tidyverse)
library(plyr)
library(ggrepel) # for plotting with text
library(reshape2)
library(cowplot) # for arranging multiple plots
# geospatial data
# library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(rfishbase)
library(taxize)



##################################
## DATA ON CONSUMER COMMUNITIES #
################################

# FISH SEINING DATA
seines.raw <- read_csv( '../Data/Bitemap_Seine_ALL-DATA_20200104.csv' )
# Which participants did we have at each site?
seines.raw %>% 
  select( Country, Participant.List ) %>%
  mutate( Participant.List = substr(Participant.List,1,20) ) %>%
  distinct() 
# Contact for seine widths
# Holger+Macreadie, Alistair, Kelaher, Thiel+Pino, Guca, Paul York, Diskin, Hereu, Yeager+Hovel,
# Fodrie, Silliman, Kun-seop Lee, O'Leary, BrentHughes, Lane Johnston, Hanley+RandallHughes,
# Galloway+Bree, Claudia, Hultgren, Harvell+Olivia, Mathieu, Zach, MaxRobinson+Rich, 
# Rod Connolly, Francesca, PG Ross

# 3 Jan 2018 - widths added to seine data
# data entered as Site, width of transect, width of seine
# default width of transect (how far apart seine stretched before pulling) is 10m as per protocol
seine.width <- data.frame( do.call( rbind,
  list( c( "USA (CA)", 9.144, 9.144 ),
        c( "USA (CA3)", "default", 10 ),
        c( "Chile", "default", 10 ),
        c( "Korea", "default", 15 ),
        c( "Wales", "default", 12.2 ),
        c( "France", "default", 60),
        c( "Mexico (BN)", 9, 9 ),
        c( "Brazil", "default", 10 ),
        c( "USA (OR)", "default", 10 ),
        c( "Australia (QLD)", "default", 16 ),
        c( "Canada (QC)", "default", 19 ),
        c( "Croatia", "default", 13 ),
        c( "USA (MA)", 9, 9 ),
        c( "USA (CA2)", "default", 10.2 ),
        c( "Australia (VIC)", "default", 10 ),
        c( "USA (NC)", 7, 7 ),
        c( "Canada (BC)", "default", 11),
        c( "USA (WA2)", "default", 10),
        c( "USA (WA)", "default", 10))
), stringsAsFactors = F )
names(seine.width) <- c("Country","width.transect","width.seine")
seine.width$width.seine <- as.numeric(seine.width$width.seine)
seine.width$width.transect[ seine.width$width.transect == "default" ] <- 10 
# seine.width$width.transect[ seine.width$width.transect == "circle" ]  <- seine.width$width.seine[ seine.width$width.transect == "circle" ]  / pi
seine.width$width.transect <- as.numeric(seine.width$width.transect)
# merge seine width metadata
seines.raw <- left_join( seines.raw, seine.width )


# for other sites, assume default conditions (10m wide transect, pull to shore in 100cm)
seines.raw %>% 
  select( Country ) %>%
  distinct() 
seines.raw$width.transect[ is.na(seines.raw$width.transect)] <- 10

# SITE LOCATIONS
siteGPS <- read.csv( '../Data/Bitemap_sites.csv')
# calculate distances between all points
siteGPS[,3:4]
# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)
rads    <- deg2rad(siteGPS[,3:4])


# Predation Rates from other script
# rate.env <- read.csv( "Output Data/Bitemap_DATA_analysis_20180825.csv", stringsAsFactors = FALSE )
rate.env <- read_csv( "Output Data/Bitemap_rate.env.20200222.csv" )
# below, we deal with site and habitat level data, not at the same scale
# summarize to average predation rates at site and habitat level -- need to be careful with habitat because this has problems with traits
rate.mean <- rate.env %>%
  dplyr::group_by(Country,  habitat) %>%
  dplyr::summarise( rate=mean( rate,na.rm=TRUE ) )
rate.mean2 <- rate.env %>%
  dplyr::group_by(Country) %>%
  dplyr::summarise( rate=mean( rate,na.rm=TRUE ) )

###########################################################################################################
######                                    Cleaning the data                                           #####
###########################################################################################################

# clean up  data
seines <- seines.raw
names(seines)[4] <- "habitat"
# convert all vegetated and unvegetated sites to common categories
seines$habitat[ seines$habitat %in% c("Seagrass ","seagrass" )] <- "Seagrass"
seines$habitat[ seines$habitat %in% c("unveg","Unveg","Unvegetated" )] <- "Unvegetated"
seines <- droplevels(seines)
# for all organisms from France, UNC2, Wales, QLD2, QLD3, multiple lengths by 10 to convert from cm to mm
seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] <- 
  seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] * 10


# select relevant columns
seines <- seines %>%
  dplyr::select( Site.Name='Site Name', Lat, Long, habitat, Date, Time, Depth, Distance, Q.PA='Q/PA', Phylum,
                 Genus, Species, Length, Abundance, Country, width.transect, width.seine )

# fix mispelling of several site names in Virginia
seines$Site.Name[ seines$Site.Name == "PungSand1" ] <- "PungSAV1"
seines$Site.Name[ seines$Site.Name == "PungSand2" ] <- "PungSAV2"
seines$Site.Name[ seines$Site.Name == "PungSand3" ] <- "PungSAV3"

# 
# # DON'T REMOVE INVERTS AT THIS STAGE
# # which sites sampled crabs?
# pods <- seines[ seines$Phylum=="Arthropoda" & !is.na(seines$Phylum), ]
# # remove inverts
# fish <- seines[ seines$Phylum=="Chordata" & !is.na(seines$Phylum), ]


### Fish length info
# What to do with rows where we have fish abundance but no length info?
# options: ignore these rows (underestimate biomass, potentially substantially)
#  use an average size based on other length estimates for biomass (unclear if under- or overestimate)
#  use largest or smallest lengths
# for now, use average (median) size for all other length estimates

# find rows with NA for size, but omit rows where nothing was found
seinesNA <- seines[ is.na(seines$Length) & seines$Abundance>0, ]

# calculate median and mean sizes for all other seines species from each seine
# first separate out the data for which we have a length esimtate
seinesL  <- seines[ !is.na(seines$Length) & seines$Abundance>0, ]
# add a flag so we know which were estimated later
seinesL$length.measured <- 'yes'
# calculate median and mean
aveLengths <- ddply( seinesL, .(Site.Name,habitat,Country,Date,Genus,Species), 
                     summarize, medLength=median(Length), meanLength=mean(Length) )
# remove NA values for unknown Genera
aveLengths[is.na(aveLengths$Genus),]
aveLengths <- aveLengths[!is.na(aveLengths$Genus),]

# replace NA length values with calculated median or mean
# rename seinesNA as seinesEST for adding length ESTimates
# Note this will only work for fishes
seinesEST <- seinesNA
# loop over all rows and get mean lengths
for( i in 1:nrow(seinesEST) ){
  if(seinesEST$Phylum[i] == "Chordata" ) seinesEST$Length[i] <- aveLengths$meanLength[ aveLengths$Country==seinesEST$Country[i] &
                                                aveLengths$habitat==seinesEST$habitat[i] &  
                                                aveLengths$Site.Name==seinesEST$Site.Name[i] & 
                                                aveLengths$Date==seinesEST$Date[i] &
                                                aveLengths$Genus==seinesEST$Genus[i] &
                                                aveLengths$Species==seinesEST$Species[i] ]
}
### errors from this step that were 'corrected'
## Syngnathus from NC, no other pipeseines were caught but this one not measured, 
# so use median of all syngnathus captured from that site
## Ambassis jacksoniensis from NSW2 20170223 + 201702024 not measured to maximize survival. 
# Notes state that all were between 15 and 60 mm, so using mean of these values (37.5mm)

# merge size estimates back with original seines data
seines.clean <- full_join( seinesL, seinesEST )

# combine genus and species with name to match length-weight coefficients (see below)
seines.clean$SPECIES_NAME <- with(seines.clean, paste(Genus, Species) )






## include data from partners (using trait information) to assess likelihood of interaction with squidpops
# Note this data is qualitative in nature. It is based on expert opinion, known feeding traits, but not yet to direct observation
traits <- read.csv( '../Data/Bitemap_REQUEST_trait_siene+video_EAT_SQUID.csv', 
                    stringsAsFactors = FALSE, strip.white = TRUE )
# 2018.12.27 - edited traits file to change family of Caesio sp. from Lutjanidae to Caesionidae
# remove rows with NA for family
traits <- traits[ !is.na(traits$family),]
# fix spelling error of site names in Viriginia
traits$siteName[ traits$siteName == "PungSand1" ]

str(traits)
# select columns
traits <- traits %>%
  dplyr::select( Country, Site.Name=siteName, sciName, totalCount, family, feeding=feedingTypeFishbase, trophic=trophicGroupRLS,
          waterColumn=waterColumnRLS, diel=dielActivityRLS, habitatRLS, eat.squid )
# which families are missing squid eating info?
# add phylum
# tax <- classification( sort(unique(traits$family)), db="worms" )
# taxbind <- rbind(tax)
# phyla <- taxbind %>% filter( rank=="Phylum" ) %>% 
#   dplyr::select( phylum=name, family=query )
# write.csv( phyla, "Output Data/predator_families+phyla.csv", row.names = FALSE )
phyla <- read.csv( "Output Data/predator_families+phyla.csv", stringsAsFactors = FALSE )
traits <- left_join( traits, phyla )

# Label sites with whether data are based on seines or video or transect
# start with a character vector labelled seine
traits$method <- "seine"
# overwrite with new categories based on representation on method of data collection
# which sites in traits are not in seines?
noseine <- unique(traits$Country)[ !(unique(traits$Country) %in% unique(seines.clean$Country)) ]
# label all of these sites unique to traits as "video"
traits$method[ traits$Country %in% noseine ] <- "video"
# for India, data are actually based on a diver transect
traits$method[ traits$Country == "India" ] <- "transect"

# Note some were left ambiguous because a clear decision could not be made
sort( unique( traits$family[is.na(traits$eat.squid)] ) )

# decision for different taxa
# not likely to eat squid
no.squid <- c( "Pleuronectidae", "Paralichthyidae", "Rhombosoleidae", "Hippolytidae", 
               "Gasterosteidae", "Agonidae","Pholidae",
               "Gobiidae", "Paguridae",  # blennies kept as yes.squid
               "Stichaeidae","Triglidae","Batrachoididae","Platycephalidae","Synodontidae" )#,  # these five added 28 December 2018
               # "Cottidae" ) # added 07 March 2018
               # likely to eat squid
yes.squid <- c( "Labridae", "Penaeidae", "Monodactylidae", 
                "Hemiramphidae", "Acanthuridae", "Blenniidae","Scyliorhinidae", # last frou added 28 Janurary 2019 
                "Embiotocidae" ) # re-added on 19 March 2019
# note some labrids already marked as not eating squid because they are herbivorous, but they might also eat squid
# Monodactylidae added to yes.squid on 6 Jan 2018

# impose decisions on data
traits$eat.squid[ traits$family %in% no.squid ] <- 0
traits$eat.squid[ traits$family %in% yes.squid ] <- 1    # is.na(traits$eat.squid) & # previously restrictred to non-NA entries

traits[ traits$family %in% c("Pholidae","Agonidae"), ]
traits[ traits$family %in% c("Sparidae"), ]
traits[ traits$family %in% c("Labridae") & traits$Country == "Australia (QLD3)", ]
traits[ traits$family %in% c("Monodactylidae"), ]
traits[ traits$family %in% c("Cottidae"), ]

# write to disk
# write_csv( traits, "Output Data/consumer_eat_squid.csv" )

## SENSITIVITY TEST
# include all species
# traits$eat.squid <- 1    


# Pick a Site
traits %>%
  filter( siteName=="Noosa river") %>%
  select( Country, sciName, family, feeding, eat.squid, totalCount )

# write a data.frame of unique families and their squid "trait"
famlist <- traits %>%
  select( family, eat.squid ) %>%
  distinct() 

write_csv( famlist, "Output Data/Bitemap_predator_families.csv" )


# # traits$omnivory[ is.na(traits$omnivory) ] <- 0
# # traits$herbivory[ is.na(traits$herbivory) ] <- 0
# # Acropomatidae  -- no
# # Balistidae     -- triggerfish, no?
# # Carcharhinidae -- no
# # Chaetodontidae -- chaetodon auriga - yes to omnivory
# # Echeneidae     -- no
# # Elopidae       -- no
# # Hippolytidae   -- yes to omnivory and herbivory
# # Caesionidae    -- no
# # Lysmatidae     -- yes
# # Muraenidae     -- no
# # Nemipteridae   -- no
# # Panopeidae     -- no
# # Pomacentridae  -- yes to omnivory and herbivory
# # Rhinobatidae   -- no
# # Salangidae     -- no
# # Syliorhinidae  -- no
# # Trachinidae    -- no
# # Urotrygonidae  -- no
# # Scaridae       -- yes to omni and herbivory
# traits$omnivory[ traits$family %in% c("Chaetodontidae","Hippolytidae","Lysmatidae","Pomacentridea",
#                                       "Scaridae") ]  <- 1
# traits$herbivory[ traits$family %in% c("Hippolytidae","Pomacentridea","Scaridae") ]  <- 1


# add site information
sites <- read.csv( "Output Data/Bitemap_BioORACLE_20190107.csv", stringsAsFactors = FALSE )
sites <- sites[ ,c(1:31)]
# sites <- left_join( sites,siteGPS[,1:2], by=c("Country") )

# # add ocean basin information here
# sites$Country
# sites$basin <- c( 1,1,1,1,1,1,
#                   1,1,1,2,2,1,
#                   2,1,3,3,4,2,
#                   3,1,1,2,2,2,
#                   2,1,1,1,1,2,
#                   2,2,2,2,1,2,
#                   2,2,1,1,2)
# sites$basin <- factor( sites$basin, levels=1:4, labels=c("Pacific","Atlantic", "Mediterranean", "Indian") )
# sites$coast <- c( 1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,2,2,2,2,1,2,1,1,2,1,2,2,2,2,1,1,1,1,1,2,1,1,1,2,2,2)
# sites$coast <- factor( sites$coast, levels=1:2, labels=c("West","East") )

# merge sites with traits
traits <- left_join( traits, sites, by=c("Country") )

# select columns
traits.clean <- traits %>%
  dplyr::select( SPECIES_NAME=sciName, phylum, family, eat.squid, sstmean ) %>%  #, omnivory, herbivory
  distinct()


# merge trait information into seine data
seine <- left_join( seines.clean, traits.clean, by=c("SPECIES_NAME") )


#### filter predators by likelihod of interacting with squidpops
predators <- seine[ seine$eat.squid==1, ]
# remove NA values
predators <-  predators[ !is.na(predators$habitat), ]

# How certain are we that we have a good estimate of species richness at each site?
# If a species is labeled with a genus or family designation, do we count one species? YES
# If a species is labelled as one potential member of a genus, but other known taxa
#      are present, should we count all species??
# aggregate species list for each site, then interrogate each list
by( predators, predators$Site.Name, function(z) sort(unique( z$SPECIES_NAME )) )
by( predators, predators$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# these all look good
# Sparidae
by( predators[predators$family=="Sparidae",], 
    predators[predators$family=="Sparidae",]$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# Cottidae
by( predators[predators$family=="Cottidae",], 
    predators[predators$family=="Cottidae",]$Country, function(z) sort(unique( z$SPECIES_NAME )) )

# standardize predator abundance by seined area and volume
predators <- predators %>%
  # calculate area and approximate volume seined
  # assuming the seined volume is a wedge-shaped object (triangular in cross-section) because seines pulled from depth to shore, hence the division by 2
  dplyr::mutate( area = Distance*width.transect, volume=(Distance*Depth)/2 * width.transect  ) %>%
  # standardize abundance by area and volume (Catch per unit effort) units are meters squared and cubed
  dplyr::mutate( cpua=Abundance/area, cpuv=Abundance/volume )

# get total abundance by Site and seine 
abundance.all <- predators %>%
  dplyr::group_by( Site.Name, Country, Date, Time, Lat, Long, Distance, area, volume, habitat) %>%
  dplyr::summarise( total.cpua=sum(cpua, na.rm=T), total.cpuv=sum(cpuv, na.rm=T))

# NEED TO DO THIS FOR FISH SEPARATELY?
abundance.fish <- predators %>%
  filter( phylum=="Chordata" ) %>%
  dplyr::group_by( Site.Name, Country, Date, Time, Lat, Long, Distance, area, volume, habitat) %>%
  dplyr::summarise( total.cpua.fish=sum(cpua, na.rm=T), total.cpuv.fish=sum(cpuv, na.rm=T))

# join these
abundance <- left_join( abundance.all, abundance.fish )
abundance$total.cpua.fish[ is.na(abundance$total.cpua.fish) ] <- 0
abundance$total.cpuv.fish[ is.na(abundance$total.cpuv.fish) ] <- 0

# Build family level abundance (+presence-absence) matrices for families thought to eat squid
# For each seine
# Note: totalCount from data.frame traits is not valid and must be standardized
# use predators
families <- predators %>%
  dplyr::group_by( Country, Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, width.transect,
            Phylum, family ) %>%
  dplyr::summarise( cpua = sum(cpua), cpuv = sum(cpuv) ) 
  
# add some trait information back in
traits.clean2 <- traits %>%
  select( Country, method, meanLat, sstmean, hemi, Site, basin, coast, family, eat.squid ) %>%
  distinct() %>%
  arrange( Country )
  
fam.trait <- left_join( families, traits.clean2 )

# a few families to consider
fam.trait[fam.trait$family=="Fundulidae",]
fam.trait[fam.trait$family=="Palaemonidae",]
fam.trait[fam.trait$family=="Portunidae",]
fam.trait[fam.trait$family=="Cancridae",]
fam.trait[fam.trait$family=="Cottidae",]
fam.trait[fam.trait$family=="Haemulidae",]
fam.trait[fam.trait$family=="Lutjanidae",]

# use gather and spread instead of melt and cast
# first get average abundance across all seines within a particular habitat at a site
fam.cpua.spread <- fam.trait %>%
  dplyr::group_by( Country, Site.Name, Lat, Long, habitat, family ) %>%
  dplyr::summarise( mean.cpua = mean(cpua,na.rm=T) ) %>%
  ungroup()
  # # then spread out the families
  # spread( family, mean.cpua )
# convert NA to 0
fam.cpua.spread[ is.na(fam.cpua.spread) ] <- 0
  




# expand these data to show where families were NOT found
fam.cpua.expand <- fam.cpua.spread %>%
  select( Country, Site.Name, habitat, family, mean.cpua ) %>%
  complete(  nesting(Country, habitat), family, fill=list(mean.cpua = 0)  )

# join with rate.mean (note again that individual seines DO NOT match up with predation assays)
fam.cpua.rate <- left_join( fam.cpua.spread, rate.mean )
#
# look family by family at the relationship between abundance and predation rate
# windows()
ggplot( fam.cpua.rate, aes(x=log(mean.cpua),y=rate, col=habitat) ) + facet_wrap(~family, scales="free") +
  geom_point() + geom_smooth(method='lm')






# calculate family richness (how many families at a site)
famrich   <- traits %>%
  dplyr::group_by(Site,method) %>%
  filter( eat.squid==1 ) %>%
  dplyr::summarise( famrich=length(unique(family)) )
  
 
# famrich.squid <- ddply( family.squid, .(Site, method), summarize, 
#                         famrich.squid = length(unique(family)) )

# famrich <- full_join( famrich.all, famrich.squid )
famrich.geo <- left_join( famrich, sites )

# number of predator families (likely to eat squidpops) as a function of distance from equator
ggplot( data=famrich.geo, aes(x=abs(meanLat),y=famrich)) + geom_point() # India is low, but based on videos
# windows(4.5,3.5)
# by method
ggplot( data=famrich.geo, aes(x=abs(meanLat),y=famrich, col=method)) + 
  geom_smooth( ) + #geom_smooth( aes(group=1) ) + 
  geom_point(size=3) + guides( fill = guide_legend(override.aes = list(linetype = 0,fill=NA)),
                               color = guide_legend(override.aes = list(linetype = 0,fill=NA)) ) +
  ylab("Number of predator families") + xlab("Degrees from equator")
# by hemisphere
ggplot( data=famrich.geo, aes(x=abs(meanLat),y=famrich, col=hemi, lty=hemi)) + 
  geom_smooth( data=famrich.geo[famrich.geo$method == "seine",],se=T ) + #geom_smooth( aes(group=1) ) + 
  geom_point( size=3 ) + guides( fill = guide_legend(override.aes = list(fill=NA)),
                                 color = guide_legend(override.aes = list(fill=NA)) ) +
  ylab("Number of predator families") + xlab("Degrees from equator") + ylim(c(0,20))

# presence-absence matrix for families (all and thought to eat squid)
# include all data from seining, video, and transects
trait.rate <- left_join( traits, rate.mean2 )
fam.pres <- trait.rate %>%
  filter( eat.squid==1 ) %>%
  select( Site,method,hemi,basin,coast,meanLat,sstmean,rate,family,totalCount ) %>%  #,omnivory,herbivory
  mutate( totalCount = ifelse(totalCount>0,1,0) ) %>%
  distinct()

# spread out the families
fam.spread <- fam.pres %>%
  # group_by(Site) %>%
  spread( family, totalCount )
# fam.pa <- fam.spread
fam.spread[ is.na(fam.spread) ] <- 0



# separate community data from site data
group <- fam.spread
fam.meta <- group[,1:8]
fam.data <- group[,-c(1:8)]


# site and family summaries
colSums(fam.data)
rowSums(fam.data)


# add omnivory and herbivory?

# add summarized trophic trait information from FishBase script
trophic <- read.csv( "Output Data/omnivory_site.csv", stringsAsFactors = FALSE )
trophic$meanprop <- apply( trophic[,c("omniprop","herbprop")], 1, mean )
trophsite <- left_join(sites,trophic)
fam.meta2 <- left_join( fam.meta, trophsite )


# which sites only had 2 families
group[ rowSums(fam.data) <= 2, ] 












##### MULTIVARIATE SECTION

# NMDS on families
library(vegan)
library(viridis)

# # test a number of dimensions
# mds2 <- metaMDS( fam.data, "raup", k = 2, binary=TRUE, trymax = 100 )
# mds3 <- metaMDS( fam.data, "raup", k = 3, binary=TRUE, trymax = 100 )
# mds4 <- metaMDS( fam.data, "raup", k = 4, binary=TRUE, trymax = 100 )
# mds5 <- metaMDS( fam.data, "raup", k = 5, binary=TRUE, trymax = 100 )
# mds6 <- metaMDS( fam.data, "raup", k = 6, binary=TRUE, trymax = 100 )
# mds7 <- metaMDS( fam.data, "raup", k = 7, binary=TRUE, trymax = 100 )
# mds8 <- metaMDS( fam.data, "raup", k = 8, binary=TRUE, trymax = 100 )
# mds9 <- metaMDS( fam.data, "raup", k = 9, binary=TRUE, trymax = 100 )
# mds10 <- metaMDS( fam.data, "raup", k = 10, binary=TRUE, trymax = 100 )
# mds11 <- metaMDS( fam.data, "raup", k = 11, binary=TRUE, trymax = 100 )
# mds12 <- metaMDS( fam.data, "raup", k = 12, binary=TRUE, trymax = 100 )
# mds13 <- metaMDS( fam.data, "raup", k = 13, binary=TRUE, trymax = 100 )
# mds.comp <- list( mds2, mds3, mds4, mds5, mds6, mds7, mds8, mds9, mds10, mds11, mds12, mds13 )
# stress <- unlist( lapply( mds.comp, function(z) z$stress ) )
# plot( stress ) # no break point

# windows(7,4)
mds <- metaMDS( fam.data, "raup", k = 3, binary=TRUE, trymax = 100 )
mds

stressplot( mds )


# capscale
row.names(fam.data) <- fam.meta$Site
# reduce number of families to most abundant
y <- fam.data[-10, which( colSums(fam.data) > 3 )] # remove Delaware. No SST data
x <- group[-10,]

fam.cap <- capscale( y ~ meanLat + poly(sstmean,2,raw = F) , x, dist="raup")
plot(fam.cap)
points( fam.cap, pch=21, cex=x$rate*10)

fam.rda <- dbrda( y ~ poly(sstmean,2), x, dist="raup")
fam.bioenv <- bioenv( y ~ sstmean, x,index="raup")
# pull out constrained axes and look at relationship with SST
cap1 <- fam.cap$CCA$u
cap2 <- fam.cap$Ybar
plot(x$rate ~ cap1[,1] )
plot(x$rate ~ cap2[,1] )
rda1

## rotate the NMDS to mean annual SST
mds.rot <- MDSrotate( mds, group$sstmean, na.rm=TRUE )

## rotate the original NMDS to mean predation rate
mds.rot.rate <- MDSrotate( mds, group$rate, na.rm=TRUE )

mds.use <- mds

plot( mds.rot.rate, "sites", cex=(fam.meta$rate*3)+1, xlim=c(-1,1) )
vec.fam <- envfit( mds.rot.rate$points, fam.data, perm=1000 )
# add to plot
plot( vec.fam, p.max=0.05, col="blue" )

vec.fam.df <- as.data.frame(vec.fam$vectors$arrows*sqrt(vec.fam$vectors$r))
vec.fam.df$family <- rownames(vec.fam.df)
vec.fam.df$p      <- vec.fam$vectors$pvals

# summarize species richness of families associated with colder or warmer temps
vec.pos <- vec.fam.df %>% filter( MDS1>0 )
vec.neg <- vec.fam.df %>% filter( MDS1<0 )

# significant 
vec.pos[ vec.pos$p < 0.05, ]
# strongly in direction of increasing predation rate
vec.pos.strong <- vec.pos[ vec.pos$MDS1 > 0.25, ]
# write to disk
# write.csv( vec.pos.strong, "Output Data/NMDS_rateRotate_positive.csv", row.names = FALSE )

# plot them all together
mdsdf <- cbind( group, mds.rot.rate$points )
windows(3,3)
ggplot( mdsdf, aes(x=MDS1,y=MDS2) ) + 
  geom_point(  aes(size=rate),shape=1, col='slateblue') +
  geom_segment(data=vec.pos.strong,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.2, "cm")),colour="black", alpha=0.5) + 
  geom_text_repel(data=vec.pos.strong,aes(x=MDS1,y=MDS2,label=family),
                  point.padding=0.9, box.padding = 0.1,
                  size=4, col="black", segment.color="slategray", segment.alpha = 0.5) +
  # xlim(xlimits) + ylim(ylimits) +
  theme_classic() + theme( axis.line = element_blank() )



plot( mds.use, "sites", cex=(fam.meta$rate*3)+1, xlim=c(-1,1) )
vec.fam <- envfit( mds.rot.rate$points, fam.data, perm=1000 )
# add to plot
plot( vec.fam, p.max=0.05, col="blue" )

vec.fam.df <- as.data.frame(vec.fam$vectors$arrows*sqrt(vec.fam$vectors$r))
vec.fam.df$family <- rownames(vec.fam.df)
vec.fam.df$p      <- vec.fam$vectors$pvals
#write to disk
# write.csv( vec.fam.df, "Output Data/family_vectors.csv", row.names = FALSE )

fam.sel <- rownames(vec.fam$vectors[[1]])[ vec.fam$vectors$pvals < 0.05 ]

# plot species richness within each families
by( predators, predators$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# only select relevant columns of predator data.frame
pred.sel <- predators %>%
  select( Country, Site.Name, Lat, Long, habitat, SPECIES_NAME, family ) %>%
  distinct( )
# aggregate to count the number of rows in each family at each site
pred.richness <- ddply( pred.sel, .(Country, Lat, Long, habitat, family), 
                        summarize, S = length(SPECIES_NAME) )
summary(pred.richness)
# tend to be very few species within families
# look at species list for sites with high species richness within particular families
pred.richness[ pred.richness$S >2, ]
predators[ predators$Country=="Canada (BC)" & predators$habitat=="Seagrass" & predators$family=="Cottidae",]
predators[ predators$Country=="USA (CA)" & predators$habitat=="Seagrass" & predators$family=="Serranidae",]
predators[ predators$Country=="USA (TX)" & predators$habitat=="Seagrass" & predators$family=="Sciaenidae",]
predators[ predators$family=="Diodontidae",]
unique( predators$Country[ predators$family=="Sillaginidae" ] )
# merge with the environmental data
pred.rich.env <- full_join( pred.richness, sites )
# select influential families
pre.sel <- pred.rich.env[ pred.rich.env$family %in% fam.sel, ]


# pad with zeros
S <- 0
zeros <- expand.grid( Site = fam.meta$Site, family=unique(phyla$family), S=S )
# pre.pad <- full_join( pre.sel,zeros, by=c("Site","family","S") )

# gather the families
group.sel <- group  # cbind( fam.meta, group[,fam.sel] ) to pick families with "signifcant" loadings on NMDS
group.sel.melt <- melt( group.sel, id.vars=1:8, value.name = "S", variable.name="family"  ) 
# add phylum back in
group.sel.melt <- left_join( group.sel.melt, phyla, by="family" )

# 
# # windows( 6,3 )
# ggplot( group.sel.melt, aes(x=sstmean, y=S ) ) + facet_wrap(phylum~family, ncol=5) +
#   geom_smooth(se=F) + 
#   geom_point( pch=1 ) +
#   ylab( "Abundance of each family" ) + xlab("Mean Annual SST (?C)") 
# 
# ggplot( pre.pad, aes(x=sstmean, y=S ) ) + facet_wrap(~family, ncol=10) +
#   geom_smooth(se=F) + 
#   geom_point( pch=1 ) +
#   ylab( "Species richness" ) + xlab("Mean Annual SST (?C)")

# dev.off()


# summarize species richness of families associated with colder or warmer temps
vec.pos <- vec.fam.df %>% filter( MDS1>0 )
vec.neg <- vec.fam.df %>% filter( MDS1<0 )

# aggregate and combine
pos.richness <- plyr::ddply( group.sel.melt[ group.sel.melt$family %in% vec.pos$family, ], .(sstmean), 
                             summarize, S=sum(S) )
pos.richness$thermal <- "S.warm"
neg.richness <- plyr::ddply( group.sel.melt[ group.sel.melt$family %in% vec.neg$family, ], .(sstmean), 
                             summarize, S=sum(S) )
neg.richness$thermal <- "S.cold"
total.richness <- plyr::ddply( group.sel.melt, .(sstmean), summarize, S=sum(S) )
total.richness$thermal <- "S.all"
richness.groups <- rbind(pos.richness,neg.richness,total.richness)

# merge richness and sites information
# first melt and cast richness.groups
richness.melt <- melt( richness.groups, id=c(1,3) )
richness.cast <- dcast( richness.melt, sstmean~thermal )
sites.rich <- full_join( sites, richness.cast )



# plot them
windows( 3.5,3 )
ggplot( richness.groups, aes(x=sstmean,y=S) ) + facet_wrap( ~thermal ) +
  # geom_smooth( se=F, col='black' ) + 
  geom_point( pch=1, col='slateblue', alpha=0.7 ) +
#  stat_smooth(method="glm", method.args=list(family="quasipoisson"), 
 #               formula = y ~ splines::ns(x, 2), col="black") +
  ylab( "Species richness" ) + xlab("Mean Annual SST (?C)")



#
adonis2( fam.data ~ coast:(hemi+basin), data=group )
adonis2( fam.data ~ rate, data=group )
adonis2( fam.data ~ coast:(hemi+basin)+rate, data=group )

# create interaction dummies
group$hemicoast  <- with( group, paste(hemi,coast) )
group$coastbasin <- with( group, paste(coast,basin) )

# combine data and NMDS results
mdsfd <- cbind( group, mds.use$points )
# get rid of some columns
mdsfd <- mdsfd %>%
  select( -sstmean )

# now merge in the rotated nmds results
sites.nmds <- full_join( sites.rich, mdsfd, by=c("Site","meanLat","hemi","basin","coast" ) )
sites.nmds <- sites.nmds[ !is.na(sites.nmds$Country), ]
# strong correlation with mean annual SST
plot( MDS1 ~ sstmean, sites.nmds )



### Ellipses for groups
plot( mds, "sites" )
ord<-ordiellipse( mds, group$coastbasin, display = "sites", 
                  kind = "se", conf = 0.95, label = T )
## new data.frame to show ellipses
# https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
# define function for ellipse...DELETED

## unwrap ordisurf object
# https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/
#ordisurf:
ordi<-ordisurf( mds.use, group$sstmean, display = "sites", 
               kind = "se", conf = 0.95, label = T )
ordi.grid <- ordi$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.pred <- expand.grid(x = c(ordi.grid$x), y = c(ordi.grid$y)) #get x and ys
ordi.pred$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.na <- data.frame(na.omit(ordi.pred)) #gets rid of the nas
# ordi.na #looks ready for plotting!
# # repeat with omnivory
# ordi2<-ordisurf( mds.use, fam.meta2$herbprop, display = "sites", 
#                 kind = "se", conf = 0.95, label = T )
# ordi.grid2 <- ordi2$grid #extracts the ordisurf object
# str(ordi.grid2) #it's a list though - cannot be plotted as is
# ordi.pred2 <- expand.grid(x = c(ordi.grid2$x), y = c(ordi.grid2$y)) #get x and ys
# ordi.pred2$z <- as.vector(ordi.grid2$z) #unravel the matrix for the z scores
# ordi.na2 <- data.frame(na.omit(ordi.pred2)) #gets rid of the nas



# the plots
# windows(6,4)
xlimits <- c(-1.5,1.5)
ylimits <- c(-1.0,1.0)
labelsize <- 12
( sstsurf <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
  stat_contour( data = ordi.na, aes(x = x, y = y, z = z, colour = ..level..),
               binwidth = 1.85, size=0.8 ) + #can change the binwidth depending on how many contours you want
  geom_text_repel( aes(label=Site), point.padding = 0.1 ) +
  geom_point( aes(size=group$rate) ) +
  scale_color_viridis( name = "mean\nannual SST" ) +
  scale_size_continuous( name = "predation\nrate" ) +
  xlim(xlimits) + ylim(ylimits) +
  theme_classic() + theme( axis.line = element_blank(), 
                           axis.text=element_text(size=labelsize),
                           axis.title=element_text(size=labelsize+2) ) ) 

# # omnivory
# windows(6,4)
# xlimits <- c(-1.5,1.5)
# ylimits <- c(-1.0,1.0)
# labelsize <- 12
# ( sstsurf <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
#     stat_contour( data = ordi.na2, aes(x = x, y = y, z = z*100, colour = ..level..),
#                   binwidth = 1, size=0.8 ) + #can change the binwidth depending on how many contours you want
#     geom_text_repel( aes(label=Site), point.padding = 0.1 ) +
#     geom_point( aes(size=rate, col=fam.meta2$herbprop*100) ) +
#     scale_color_viridis( name = "Number omnivorous taxa" ) +
#     scale_size_continuous( name = "predation\nrate" ) +
#     xlim(xlimits) + ylim(ylimits) +
#     theme_classic() + theme( axis.line = element_blank(), 
#                              axis.text=element_text(size=labelsize),
#                              axis.title=element_text(size=labelsize+2) ) ) 


with( mdsfd, cor.test( MDS1, sstmean ) )
cor.test( mds$points[,1], mdsfd$sstmean )
cor.test( mds.rot$points[,1], mdsfd$sstmean )

plot( x=fam.meta$propomn, y=group$rate )

# 
vec.fam.sel <- vec.fam.df[vec.fam.df$p < 0.05,]

(vectors <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
  geom_point(  col='slateblue', size=2 ) +
  geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.2, "cm")),colour="slateblue", alpha=0.5) + 
  geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family),
                  point.padding=0.01, box.padding = 1,
                  size=4, col="black", segment.color="slategray", segment.alpha = 0.5) +
  xlim(xlimits) + ylim(ylimits) +
  theme_classic() + theme( axis.line = element_blank() ))

(centroids <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
    geom_point(  col='slateblue', size=2, alpha=0.3 ) +
    geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2),colour="black", alpha=0.3) +
    geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family),
                    point.padding=0.05, box.padding = 0,
                    size=3.5, col="black") +
    xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank() ))

(centroids2 <- ggplot( mdsdf, aes(x=MDS1,y=MDS2) ) + 
    # geom_point(  col='slateblue', size=2, alpha=0.3 ) +
    geom_point( aes(size=rate), col='slateblue', alpha=0.3 ) +
    geom_segment(data=vec.pos.strong,aes(x=0,xend=MDS1,y=0,yend=MDS2),
                 colour="black", alpha=0.8) +
    geom_text_repel(data=vec.pos.strong,aes(x=MDS1,y=MDS2,label=family),
                    point.padding=0.5, box.padding = 0.5,
                    size=4.5, col="slateblue4") +
    scale_size_continuous( name = "predation\nrate" ) +
    xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank(), 
                             axis.text=element_text(size=labelsize),
                             axis.title=element_text(size=labelsize+2),
                             legend.position = "none") ) 


vec.fam.sel <- vec.fam.df[ vec.fam.df$p < 0.05, ]
percentage <- 0.15
percentage <- 0.2  # for UNC IMS interview
mds.move <- diff(range(mds.rot.rate$points[,1])) * percentage
vec.fam.sel <- vec.fam.df[ vec.fam.df$MDS1 < -mds.move | vec.fam.df$MDS1 > mds.move, ]
diff(range(mds.rot.rate$points[,1]))
# add omnivory from Fishbase script
troph2 <- read.csv( "Output Data/trophic_family_strong_omni.csv" )
vec.fam.sel <- left_join(vec.fam.sel,troph2)

# make a color for centrouds
vec.fam.sel$color <- factor( ifelse( vec.fam.sel$MDS1 >0, "black", "firebrick" ) )
windows(4,4)
(centroids3 <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
    # geom_point(  col='slateblue', size=2, alpha=0.3 ) +
    geom_point( aes(size=rate), col='slateblue', alpha=0.3 ) +
    # geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2,col=color),
    #             alpha=0.75, arrow = arrow(length = unit(0.1,"cm")), linejoin='mitre' ) +
    # geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family, col=color),
    #                 point.padding=0.1, box.padding = 0.3,
    #                 size=4, #col="sla5teblue4",
    #                 segment.color="slateblue1", segment.size=0.3, alpha=0.8) +
    scale_size_continuous( name = "predation\nrate" ) +
    scale_color_manual( values=c("black","firebrick")) +
    xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank(), 
                             axis.text=element_text(size=labelsize),
                             axis.title=element_text(size=labelsize+2),
                             legend.position = "none") ) 

   windows(10,4)
# plot_grid(centroids,sstsurf, labels=c("A","B"), ncol = 2, nrow =1,label_size=18, rel_widths = c(1,1.25) )
plot_grid(sstsurf,centroids3, labels=c("A","B"), ncol = 2, nrow =1,label_size=18, rel_widths = c(1.25,1) )


# supplementary figure to show all predator families
vec.fam.sel <- vec.fam.df[vec.fam.df$p < 1,]  # to show all families
windows( 6,6 )
(vectors <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
    geom_point(  col='slateblue', size=1.5 ) +
    geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2),
                 arrow = arrow(length = unit(0.2, "cm")),colour="slateblue", alpha=0.5) + 
    geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family),
                    point.padding = 1,
                    size=3, col="black", segment.color="black", segment.alpha = 0.5,
                    force=6) +
     xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank() ))
# dev.off()


# what's driving communities?
# SIMPER? Bubble Plots

# try PCA, best fits??
# RDA? Correspondence analysis


## Redundant Descriminant analysis
rda1 <- rda(  fam.data ~ sstmean + meanLat, group, na.action=na.exclude )
plot(rda1)

## Canonical Correspondence analysis
cca1 <- cca(  fam.data ~ sstmean, group, na.action=na.exclude )
cca1 <- cca(  fam.data ~ sstmean, group, na.action=na.exclude )
plot(cca1)

# get rid of NA
fam.na <- na.omit( fam.pa )
# separate community data from site data
fam.meta <- fam.na[,1:6]
fam.data <- fam.na[,-c(1:6)]
## environmental data matrix
vare.cca <- cca( fam.data ~ sstmean + basin, fam.meta )
vare.cca
plot(vare.cca)

## RDA
rda1 <- rda( fam.data ~ sstmean + hemi*basin*coast, fam.meta )
anova(rda1)
summary(rda1)
plot(rda1) 

rda1$Ybar


################################
## END OF MULTIVARIATE SECTION
##
################################





### add a switch to only allow taxa associated with warm waters in the calculation
pred.all  <- predators
pred.warm <- predators[ predators$family %in% vec.pos$family, ]

# decide which predators set to use
pred.use <- pred.all

# only use the subset of sites that have both seagrass and unvegetated habitats
# which sites have both
habtab <- with(pred.use, table(Country,habitat))
habboth <- which(apply(habtab,1,function(z) all(z>0) ))
pred.both <- pred.use[ pred.use$Country %in% names(habboth),]

# pick the dataset
pred.use <- pred.both

# merge with site data
pred.sites      <- left_join( pred.use, sites )
pred.sites.both <- left_join( pred.both, sites )

##############################
# Predator size distribution #
##############################

# How to deal with all of the fish for which we have abundance but not size?
# For biomass, I use the average length
# but here, I only include data with measured length
#
# only include rows with a length estimate (pred biomass assumed average lengths for preds not measured...see below)
pred.size <- pred.sites.both[ !is.na(pred.sites$length.measured), ]


# use ggplot to make violin plots for each site x habitat combination
windows()
ggplot( pred.size, aes(x=habitat, y=log10(Length) )) + facet_wrap(~Country) + 
  geom_boxplot(width=0.5)   # Note that USA (VA) and USA (NC2) lack data from unvegetated sediments
ggplot( pred.size, aes(x=habitat, y=log10(Length) )) + geom_violin() + facet_wrap(~Country)


# pred length by latitude
ggplot( pred.size, aes(x=abs(Lat), y=log10(Length) )) + facet_grid(hemi~habitat) +
  geom_point() +    geom_smooth(method='lm')
ggplot( pred.size, aes(x=abs(Lat), y=log10(Length) )) + geom_point() + facet_wrap(~habitat) + 
  geom_smooth(  )

# quantile regression
ggplot( pred.size, aes(x=abs(Lat), y=log10(Length) )) + geom_point() + facet_wrap(~habitat) + 
  geom_quantile()
ggplot( pred.size, aes(x=abs(Lat), y=log10(Length) )) + geom_point() + facet_wrap(~habitat) + 
  geom_quantile(quantiles=seq(0.05, 0.95, by = 0.05))
ggplot( pred.size, aes(x=abs(Lat), y=(Length) )) + geom_point() + facet_wrap(~habitat) + 
  geom_quantile(quantiles=c(0.9,0.8,0.5,0.2,0.1))
ggplot( pred.size, aes(x=abs(Lat), y=(Length) )) + geom_point() + facet_grid(hemi~habitat) + 
  geom_quantile(quantiles=c(0.9,0.8,0.5,0.2,0.1))

# calculate the median size
pred.size.median <- ddply( pred.size, .(Country,Lat,habitat,hemi), summarize, Length=median(Length,na.rm=T) )
ggplot( pred.size.median, aes(x=abs(Lat), y=Length, col=habitat)) + 
  geom_smooth(method='lm',se=T) + geom_point()  # little difference between habitats in fish size
# pred.size.median <- ddply( pred.size, .(Country,Lat), summarize, Length=median(Length,na.rm=T) )
# ggplot( pred.size.median, aes(x=abs(Lat), y=Length)) + 
#   geom_smooth() + geom_point() # + geom_smooth(method='lm',se=T) 

# calculate median size for each seine
pred.med.len <- ddply( pred.size, .(Country,Site.Name,Lat,sstmean,Date,Time),
                       summarize, medLength=median(Length,na.rm=T) )
windows(4,4)
ggplot( pred.med.len, aes(x=abs(Lat), y=medLength )) + geom_point() + geom_smooth() +
  xlab('degrees latitude from equator') + ylab('median length of fish per seine')
ggplot( pred.med.len, aes(x=sstmean, y=medLength )) + geom_point() + geom_smooth() +
  xlab('degrees latitude from equator') + ylab('median length of fish per seine')



# consider that the relationship could change with hemisphere
ggplot( pred.size.median, aes(x=abs(Lat), y=Length, col=hemi)) + 
  geom_smooth() + geom_point() 

# temperature might be a better indicator anyway
ggplot( pred.sites, aes(x=sstmean, y=log10(Length) )) + facet_grid(~habitat) +
  geom_point() +    geom_smooth(method='lm') # slight increase in size with temperature?? 
ggplot( pred.sites, aes(x=sstmean, y=(Length) )) + facet_grid(~habitat) +
  geom_point() +    geom_quantile()

# look at size of top XXth percentile of size
# calculate size of 80th percentile 
pred.top20 <- pred.size %>%
  dplyr::group_by( Country,Site.Name,Distance,Date,Time,Lat,Long,habitat ) %>%
  dplyr::summarize( top20=quantile(Length,probs=0.8), meanL=mean(Length,na.rm=T) )



# #########################
# # Fish Trophic Levels   #
# #########################
# library(rfishbase)
# fishSpecies <- sort( unique( fish.clean$SPECIES_NAME ) )
# 
# fishEcology <- ecology( fishSpecies, fields = c("SpecCode", "FeedingType", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph") )
# 
# # sometimes DietTroph (based on whole range of food items) unavailable but FoodTroph (based on individual food items) more often is available
# fishEcology$Trophic.Level <- fishEcology$DietTroph
# for( i in 1:nrow(fishEcology) ) {
#   if( is.na(fishEcology$DietTroph[i]) )   fishEcology$Trophic.Level[i] <- fishEcology$FoodTroph[i]
# }
# 
# # which taxa are missing trophic levels?
# trophMiss <- fishEcology[ which(is.na(fishEcology$Trophic.Level)), ]
# write.csv( trophMiss, "../Data/Fish Biomass + Traits/Bitemap_2017_fishBase_missingTrophic.csv", row.names = FALSE )
# # read in data.frame with info added from fishbase
# trophFill <- read.csv( "../Data/Fish Biomass + Traits/Bitemap_2017_fishBase_missingTrophic.csv" )
# trophFill$Trophic.Level <- trophFill$Trophic
# trophFill$fishbaseTrophEst <- c(2.9,3.2,3.3,3.3,3.4,3.4,3.3,3.2,3.2,3.4,3.7)
# fishEcology <- full_join( fishEcology, trophFill[,c('sciname','fishbaseTrophEst')] )
# for( i in 1:nrow(fishEcology) ) {
#   if( is.na(fishEcology$Trophic[i]) )   fishEcology$Trophic.Level[i] <- fishEcology$fishbaseTrophEst[i]
# }
# 
# 
# # 
# # some trophic level info added manually
# trophManual <- read.csv( "../Data/Fish Biomass + Traits/Bitemap_2017_fishBase_noEcology_ecologyAdded.csv" )
# names(trophManual) <- c("sciname","Trophic.Level")
# #
# fishINFO <- full_join( fishEcology, trophManual )
# fishSpecies[ !(fishSpecies %in% fishINFO$sciname)] # several species not represented
# 
# ## join Bitemap fish data with trait data from Reef Life Survey
# # read in fish trait data from (source = RLS?)
# fishtraits <- read.csv("../Data/Fish Biomass + Traits/Traits_all-species_edit.csv")
# # create column for joining with fishEcology
# fishtraits$sciname <- fishtraits$CURRENT_TAXONOMIC_NAME
# 
# 
# # first, make sure all bitemap species are represented
# fishINFO <- full_join(fishINFO,data.frame(sciname=fishSpecies))
# # now, left_join so we don't include taxa that aren't represented in Bitemap
# fishINFO <- left_join( fishINFO, fishtraits ) # lots of gaps
# # select relevant columns
# fishINFO <- fishINFO %>%
#   dplyr::select( sciname, FeedingType, Trophic.group, Water.column, Diel.Activity, Habitat, Trophic.Level )
# 
# 
# # (fishbaseMiss <- fishSpecies[ !(fishSpecies %in% fishINFO$sciname)])
# #  write.csv( data.frame(fishbaseMiss), "../Data/Fish Biomass/Bitemap_2017_fishBase_noEcology.csv", row.names=FALSE )
# 
# 
# 
# 
# 
# 
# #################
# # Fish Taxonomy #
# #################
# # add taxonomy for fishes
# library(taxize)
# # strip sp from unknown species before running classification
# scisplit <- strsplit( fishINFO$sciname, " ", fixed=TRUE )
# which(lapply(scisplit,length)==1)
# # scisplit[[1460]] <- c("Clupeidae","sp.")
# genusSpecies <- lapply(scisplit,function(z) data.frame(genus=z[1],species=z[2]))
# genusSpecies <- do.call(rbind,genusSpecies)
# fishINFO <- cbind(fishINFO,genusSpecies)
# 
# # omit sp designations
# for(i in 1:nrow(fishINFO)){
#   if(fishINFO$species[i] %in% c("sp.","spp.","sp. [Belize]","sp. [bottomei]")) fishINFO$species[i] <- ""
# }
# 
# # recombine genus and species without sp. designations (spaces should not matter here)
# fishINFO$sciname2 <- with(fishINFO, paste(genus,species) )
# fishINFO <- fishINFO[ order(fishINFO$sciname2), ]
# # get rid of duplicate entries
# fishINFO <- fishINFO[ !duplicated(fishINFO$sciname2), ]
# 
# # run classification on all fish taxa
# fishtax <- classification(fishINFO$sciname2,db="ncbi")
# 
# 
# # only keep non-NA ones
# which( unlist(lapply( fishtax, function(z) all(!is.na(z)) )) )
# cf <- fishtax[which( unlist(lapply( fishtax, function(z) all(!is.na(z)) )) )]
# 
# # combine the results
# # cannot cbind with NA values, so rbind
# cr <- do.call( rbind,cf )
# cr$unique <- rownames(cr)
# cr$sciname2 <- unlist(lapply( strsplit( cr$unique, split = ".", fixed=TRUE ), function(z) z[1] ))
# 
# # melt and recast
# cm <- melt(cr)
# # omit no rank
# cm <- cm[cm$rank!="no rank",]
# ccast <- dcast( cm, sciname2~rank, value.var = "name" )
# 
# # merge taxonomy with traits
# fishes <- full_join( fishINFO, ccast, "sciname2" )
# head(fishes)

# write this to disk
# write.csv( fishes, "Output Data/fish_traits+taxonomy_20180122.csv", row.names=FALSE )




##################
# Fish diversity #
##################


# combine seines and sites
pred.site <- left_join( pred.all, siteGPS )


# summarize seine data based on counts of each taxon
pred.abund <- ddply( pred.site, .(Country,meanLat,meanLong,habitat,phylum,Genus,Species), 
                      summarise, Abundance=sum(Abundance) )
# summarize seine data based on each indivdual seine (use site.name, date, time)
pred.abund.each <- ddply( pred.site, .(Country,Site.Name,Distance, Depth, width.transect,Date,Time,Lat,Long,habitat,phylum,Genus,Species), 
                      summarise, Abundance=sum(Abundance) )



# for species diversity metrics, we need counts of each species in each habitat
# # omit invertebrates
# pred.abund <- seine.abund[ seine.abund$Phylum=="Chordata", ]
# pred.abund.each <- seine.abund.each[ seine.abund.each$Phylum=="Chordata", ]

# consider using different diversity metrics (e.g. ENSpie as in Chase & Knight 2013)
ENSPIE <- function(prop){
  ifelse( sum(prop,na.rm=T)>0, 1 / sum(prop^2, na.rm=T), NA ) 
}   
prop <- c(0.5,0.4,0.1,NA)
prop <- c(0.5,0.4,0.1,0)
prop <- c(0.5,0.4,0.1)
ENSPIE(prop)
ENSPIE(0)
# ENSPIE is less scale dependent than many other diversity metrics (arguable important here), 
# but tends to be sensitive to aggregation (which can be affected by spatial scale, 
# and likely an issue for schooling predators and seining)

# convert abundance to proportion
# calculate ENSPIE for each seine
# calculate total abundance for each seine
pred.totals.each <- ddply( pred.abund.each, .(Country,Site.Name,Distance,Depth,width.transect,Date,Time,Lat,Long,habitat),
                      summarize, Total = sum(Abundance), TotalFish = sum(Abundance[phylum=="Chordata"]) )
head(pred.comm.each[ with(pred.comm.each,order(Country,Site.Name,Distance,Date,Time,habitat)), ])

pred.comm.each <- left_join(pred.abund.each,pred.totals.each)
# remove fish total information for taxa that are not fish
pred.comm.each$TotalFish[ pred.comm.each$phylum != "Chordata" ] <- NA
# calculate the proportion of each species in each site
pred.comm.each$prop <- pred.comm.each$Abundance/pred.comm.each$Total
pred.comm.each$propFish <- pred.comm.each$Abundance/pred.comm.each$TotalFish
# calculate ENSPIE and species richness for each seine
# for richness, inlcude all taxa, regardless of whether species known (this will inflate richess at some sites)
# exceptions -- Lethrinus sp. (QLD2), Cottoidea sp.
pred.comm.each$species <- with(pred.comm.each, paste(Genus, Species) )
ENSPIE.seine <- pred.comm.each %>%
  dplyr::group_by( Country,Site.Name,Distance,Date,Time,Lat,Long,habitat ) %>%
  dplyr::summarise( ENSPIE = ENSPIE(prop), ENSPIEfish = ENSPIE(propFish),
                    richness = length(Total), richness.fish = length(TotalFish))


  
  
# calculate total abundance for each site (lump all replicate seines together)
pred.totals <- ddply( pred.abund, .(Country,meanLat,meanLong,habitat),
       summarize, Total = sum(Abundance), TotalFish = sum(Abundance[phylum=="Chordata"]),
       richness=length(abundance), richness.fish=length(Abundance[phylum=="Chordata"]) )
# get rid of NA rows
pred.totals <- pred.totals[!is.na(pred.totals$Country),]
# include totals as a column in pred.abund
pred.comm <- left_join(pred.abund.each,pred.totals)
# calculate the proportion of each species in each site
pred.comm$prop     <- pred.comm$Abundance/pred.comm$Total
pred.comm$propFish <- pred.comm$Abundance/pred.comm$TotalFish
# calculate ENSPIE for each site
ENSPIE.site <- ddply( pred.comm, .(Country,meanLat,meanLong,habitat), 
       summarize, ENSPIE = ENSPIE(prop), ENSPIEfish = ENSPIE(propFish), 
       richness=mean(richness), richness.fish=mean(richness.fish) )
# get rid of NA rows
ENSPIE.site <- ENSPIE.site[!is.na(ENSPIE.site$Country),]

# calculate effect size as difference between ENSPIE estimates
# make all Unveg sites negative
ENSPIE.site$ENSPIEdiff <- ENSPIE.site$ENSPIE
ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] <- -(ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] )
ENSPIE.diff <- ddply( ENSPIE.site, .(Country,meanLat,meanLong), summarise, ENSPIEdiff=sum(ENSPIEdiff) )


ggplot( ENSPIE.site, aes(x=habitat,y=ENSPIE)) + geom_point() + 
  geom_line(aes(group=Country)) + facet_wrap(~Country)
ggplot( ENSPIE.site, aes(x=abs(meanLat),y=ENSPIE)) + geom_point() + facet_wrap(~habitat, ncol=2) + geom_smooth()
ggplot( ENSPIE.site, aes(x=abs(meanLat),y=richness)) + geom_point() + facet_wrap(~habitat, ncol=2) + geom_smooth()
ggplot( ENSPIE.diff, aes(x=abs(meanLat),y=ENSPIEdiff)) + geom_point() + geom_smooth()

ggplot(ENSPIE.site, aes(x=habitat, y=ENSPIE, fill = habitat)) + geom_boxplot(width=0.5, notch=FALSE) +
  geom_point(alpha=0.1) +
  # facet_grid(.~Year) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  labs(title="Diveristy of Select Predators By Habitat Type",
       x = "Type of Habitat", y="Diversity\n(effective number of species)") 









###################################################################################
# PREDATOR BIOMASS                                                                #
###################################################################################

# read biomass coefficients
coef <- read.csv("../Data/Fish Biomass + Traits/20160710_RLS_biomass_coefs_20170921.csv")
# omit rows from coef for which we have no estimates
coef <- coef[ !is.na(coef$A), ]
# remove duplicates from coef
coef <- coef[ !duplicated(coef), ]



# identify species not in character list
library(rfishbase)
dlookup <- fish.clean
lookup <- sort(unique(dlookup$SPECIES_NAME[!(dlookup$SPECIES_NAME %in% coef$SPECIES_NAME)]))
sort(unique(fish.clean$SPECIES_NAME[(fish.clean$SPECIES_NAME %in% coef$SPECIES_NAME)]))

# # use rfishbase to look up taxa not represented in the reef life survey biomass conversion data.frame (coef)
# fishbaseLW <- length_weight(lookup)
# sort(unique(fishbaseLW$sciname)) # added 67 taxa
# unrepresented <- sort(unique(fish.clean$SPECIES_NAME[!(fish.clean$SPECIES_NAME %in% c(fishbaseLW$sciname, as.character(coef$SPECIES_NAME)) )]))
# write.csv( fishbaseLW, "../Data/Fish Biomass/fishbaseLW.csv", row.names = FALSE)
# write.csv( data.frame(unrepresented), "../Data/Fish Biomass/fishbaseMISSING.csv", row.names = FALSE)

# read the length-weight relationships looked up in fishbase
fishbaseLW <- read.csv( "../Data/Fish Biomass + Traits/fishbaseLW_20170921.csv" )
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
fishbaseManual <- read.csv( "../Data/Fish Biomass + Traits/Bitemap_FishbaseLW_manual.csv" )[,1:3]
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
pred_biom <- left_join(pred.site, biomassCoef, by="SPECIES_NAME" )

# remove extraneous columns
pred_biom <- pred_biom %>%
  dplyr::select(Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, width.transect, area, volume, Q.PA, Phylum, family, Genus, Species, 
         Length, Abundance, Country, SPECIES_NAME, A, B )


# which taxa don't have length-weight regression estimates
miss <- pred_biom[ is.na(pred_biom$A), ]
sort(unique(miss$SPECIES_NAME))
# should be able to get Aphanius fasciatus??

# create biomass column -- note that I've divided all lengths by 10 to convert mm to cm, which is the unit 
pred_biom$biomass <- with(pred_biom, Abundance*(A*(Length/10)^B))

# which ones have no estimate
pred_biom[ is.na(pred_biom$biomass), ]

# standardize biomass by seine distance
pred_biom <- pred_biom %>%
  mutate( biomass.area = biomass/area, biomass.vol = biomass/volume )


pred_biom %>%
  filter( Genus=="Lagodon") %>%
  dplyr::group_by( Country, Genus, Species) %>% 
  dplyr::summarise( sum=sum(Abundance))
pred_biom %>%
  filter( Country=="USA (NC2)") %>%
  dplyr::group_by( Genus, Species) %>% 
  dplyr::summarise( sum=sum(Abundance))
pred_biom %>%
  filter( Country=="USA (NC)") %>%
  dplyr::group_by( Genus, Species) %>% 
  dplyr::summarise( sum=sum(Abundance))

# Calculate total fish biomass in each seine
biomass_total <- ddply( pred_biom, .(Site.Name,Country,Date,Time,Lat,Long,Distance,area, volume,habitat), 
                        summarize, biomass.area=sum(biomass.area, na.rm=TRUE) )
sort(unique(biomass_total$Country))
sort(unique(predators$Country))


# graph biomass by habitat type
ggplot(biomass_total, aes(x=habitat, y=log10(biomass.area+0.01), fill = habitat)) + geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  # facet_grid(.~Year) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  labs(title="Fish biomass by habitat type",
       x = "Type of Habitat", y="Biomass") 
  

# total biomass by Latitude
ggplot(biomass_total, aes(x=Lat, y=log10(biomass.area+0.01))) + geom_point() + facet_wrap(~habitat,ncol=1)
ggplot(biomass_total, aes(x=abs(Lat), y=log10(biomass.area+0.01))) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='lm')

ggplot( biomass_total, aes(x=abs(Lat),y=log10(biomass))) + geom_point() + facet_wrap(~habitat, ncol=2) + geom_smooth()

# calculate average biomass by seine and habitat
biomass_mean <- ddply( biomass_total, .(Site.Name,Country,Date,Lat,Long,Distance,habitat), 
                       summarize, biomass.area=mean(biomass.area, na.rm=TRUE), 
                       biomass.vol=mean(biomass.vol,na.rm=TRUE) )
ggplot(biomass_mean, aes(x=abs(Lat), y=log10(biomass.area+0.01))) + geom_point() + facet_wrap(~habitat,ncol=2) +
  geom_smooth(method='lm')


# calculate average biomass by site and habitat
biomass_site <- ddply( biomass_total, .(Country,habitat), 
                       summarize, meanLat=mean(Lat), 
                       biomass.area=mean(biomass.area,na.rm=TRUE) )

ggplot(biomass_site, aes(x=abs(meanLat), y=log10(biomass.area+0.01))) + geom_point() + facet_wrap(~habitat,ncol=2) +
  geom_smooth(se=F) + geom_smooth(method='lm')

# graph biomass by habitat type
ggplot(biomass_site, aes(x=habitat, y=log10(biomass.area+0.01), fill = habitat)) + 
  geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  # facet_grid(.~Year) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  labs(title="Fish biomass by habitat type",
       x = "Type of Habitat", y="Biomass") 




# get family level biomass for each site
# first get total biomass within each family for each seine
fam.biomass <- pred_biom %>%
  dplyr::group_by( Country, Site.Name, Lat, Long, habitat, Time, Depth, Distance, width.transect, family ) %>%
  dplyr::summarise( biomass.area = sum(biomass.area,na.rm=T) )

# now get average biomass at the level of seining sites
fam.mass.mean <- fam.biomass %>%
  dplyr::group_by( Country, Site.Name, Lat, Long, habitat, family ) %>%
  dplyr::summarise( mean.biomass = mean(biomass.area,na.rm=T) ) %>%
  ungroup()

# join with rate.mean (note again that individual seines DO NOT match up with predation assays)
fam.mass.rate <- left_join( fam.mass.mean, rate.mean )
#
# look family by family at the relationship between abundance and predation rate
windows(15,12)
ggplot( fam.mass.rate, aes(x=log(mean.biomass),y=rate, col=habitat) ) + facet_wrap(~family, scales="free") +
  geom_point() + geom_smooth(method='lm', aes(group=1))

# pick a few taxa with more than two data points across the study
# Blenniidae
# Clinidae
# Cottidae
# Fundulidae
# Gadidae
# Gerreidae
# Haemulidae
# Hemiramphidae
# Hexagrammidae
# Labridae
# Lutjanidae
# Mullidae
# Sciaenidae
# Sebastidae
# Serranidae
# Sillaginidae
# Sparidae
# Sphyraenidae
# Terapontidae
# Tetraodontidae
# Tetrarogidae

selected <- c( 'Blenniidae', 'Clinidae', 'Cottidae', 'Embiotocidae', 'Fundulidae', 'Gadidae', 'Gerreidae', 'Haemulidae', 'Hemiramphidae', 'Hexagrammidae',
   'Labridae', 'Lutjanidae', 'Mullidae', 'Sciaenidae', 'Sebastidae', 'Serranidae', 'Sillaginidae', 'Sparidae', 
   'Sphyraenidae', 'Terapontidae', 'Tetraodontidae', 'Tetrarogidae',
   'Mugilidae','Atherinidae','Gasterosteidae'     #  added 07 March 2019
   )

fam.mass.rate.sel <- fam.mass.rate %>%
  filter( family %in% selected)

# look family by family at the relationship between abundance and predation rate
windows(10,7)
seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
unveg.color    <- "#d8b365"
ggplot( fam.mass.rate.sel, aes(x=log(mean.biomass),y=rate, col=habitat) ) + facet_wrap(~family, scales="free_x") +
  # geom_smooth(method='glm', aes(group=1), method.args=list(family=quasibinomial), col='black', se=T) +
  stat_ellipse( aes(group=1), col='black' ) +
  # geom_smooth( aes(group=1), se=F, col='black') +
  # geom_smooth( method='lm', aes(group=1), se=F, col='black' ) +
  geom_point(size=3) + ylim(c(0,1)) + 
  xlab( expression(paste(log[10],"( mean biomass (g/",m^2,") )")) ) +
  ylab( "Predation rate" ) +
  scale_color_manual( values=c(seagrass.color,unveg.color)) +
  theme( text = element_text(size=12),
         axis.text=element_text(size=10) )


# OTHERS with not enough data but right direction
# Albulidae
# Carangidae
# Cyclopteridae
# Monodactylidae
# Ostraciidae (boxfish)
# Phycidae (hakes)

# individual glms for each family (note, this is multiple testing...but haven't accounted for this)

fam.glm <- fam.mass.rate.sel %>%
  group_by( family ) %>%
  do( fitFam = lme4::glmer(rate ~ log10(mean.biomass+0.001) + (1|Country), 
                           family="binomial", data=.) )
# get coefficients and model summaries
library(broom)
fam.coef <- tidy(fam.glm, fitFam)
# write to disk
write.csv( fam.coef, "Output Data/family_biomass_rate_models.csv", row.names = FALSE )
augment(fam.glm, fitFam)
glance(fam.glm, fitFam)

lapply( fam.glm$fitFam, summary )


fam.corr <- fam.mass.rate.sel %>%
  group_by( family ) %>%
  do( corrFam = cor.test( ~ rate + log10(mean.biomass+0.001), data=., method = "pearson" ) )
fam.pearson <- tidy(fam.corr, corrFam)
write.csv( fam.pearson, "Output Data/family_biomass_rate_corr.csv", row.names = FALSE )

# # sample size
# View(fam.mass.rate.sel %>%
#   dplyr::group_by(family) %>%
#   dplyr::summarize( n=length(rate) ))

################################################
## Combine data summaries and save to disk    ##
################################################

# summarize these at the highest level of detail (smalles scale) possible
# by country, habitat, replicate seine (lat, long, date): not resolvable for squidpops at every site
# only means by habitat and country will be directly comparable to squidpops at every site



# check abundance for missing values
abundance

# length frequency
pred.top20

# diversity (ENSPIE)
diversity <- ENSPIE.seine

# Biomass
biomass <- biomass_total 


## merge
dim(abundance); dim(pred.top20); dim(diversity); dim(biomass)
predator2 <- left_join(abundance,diversity)
predator1  <- left_join(predator2,biomass)
predator  <- left_join(predator1,pred.top20)
# when different ways of counting species are used, some sites are dropped
# merge with a bigger data set, the select relevant column
pred.all <- left_join( predators,predator )

left_join( pred.all, sites )


head(pred.all)
unique( pred.all$Country[ is.na(pred.all$biomass) ] )
# # convert NA to 0 because there is no abundance, diversity, or biomass if enough species filtered out
# pred.all[ is.na(pred.all$biomass), 
#           c("Total","TotalFish","ENSPIE","ENSPIEfish","biomass") ] <- 0

biomass[biomass$Country=="USA (NC2)",]


# relationships between predator summaries
ggplot( pred.all, aes(x=log10(cpua),y=log10(biomass.area+0.01))) + geom_point() +
  geom_smooth()
ggplot( pred.all[pred.all$biomass>0,], aes(x=ENSPIEfish,y=log10(biomass.area+0.01))) + geom_point() +
  geom_smooth( ) 
ggplot( pred.all, aes(x=log10(cpua),y=ENSPIEfish)) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=log10(cpua),y=ENSPIE)) + geom_point() +geom_smooth(method='lm')

# latitudinal gradients
ggplot( pred.all, aes(x=abs(Lat),y=ENSPIEfish)) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=(cpua))) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=log10(cpua))) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=log10(biomass.area))) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=log10(biomass.area+0.01))) + geom_point() +geom_smooth(method='lm')

# boxplots of abundance, diversity, and biomass by habitat
# point colors
seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
unveg.color    <- "#d8b365"
# seagrass.color <- "#4fc179"    # see http://www.color-hex.com/color-palette/31668
# unveg.color    <- "blanchedalmond"

# get site means
pred.graph <- left_join( predator, sites )
pred.means <- pred.graph %>%
  dplyr::group_by( Country, habitat ) %>%
  dplyr::summarize( sstmean=mean(sstmean,na.rm=T), total.cpua.fish=mean(total.cpua.fish),
             ENSPIE=mean(ENSPIE,na.rm=T), ENSPIEfish=mean(ENSPIEfish,na.rm=T),
             richness=mean(richness,na.rm=T), richness.fish=mean(richness.fish, na.rm=T),
             biomass.area=mean(biomass.area,na.rm=T), top20=mean(top20,na.rm=T),meanL=mean(meanL,na.rm=T) )

# axis label size needs to be smaller
labsize = 12

# graph fish abundance by habitat type
a <- ggplot(pred.means, aes(x=sstmean, y=(total.cpua.fish+0.01), col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  scale_y_log10(breaks=c(0.001,0.01,0.1,1)) +
  # scale_x_discrete( labels=c("Seagrass", "Unvegetated")) +
  # scale_fill_manual(values=c(seagrass.color,unveg.color)) + 
  labs( x = "Mean annual SST (?C)", y=expression(paste("Catch per ", m^2)) ) 

# graph predator diversity by habitat type
b <- ggplot(pred.means, aes(x=sstmean, y=ENSPIEfish, col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Mean annual SST (?C)", y=expression(paste("Effective number species")) ) 

# graph fish size (80th percentile) by habitat type
c <- ggplot(pred.means, aes(x=sstmean, y=top20, col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Mean annual SST (?C)", y=expression(paste("Fish total length (cm)")) ) 

# graph fish biomass by habitat type
d <- ggplot(pred.means, aes(x=sstmean, y=(biomass.area+0.01), col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  scale_y_log10(breaks=c(0.1,1,10)) +
  labs( x = "Mean annual SST (?C)", y=expression(paste("Fish biomass (g per ",m^2,")")) ) 

# graph fish richness by habitat type
e <- ggplot(pred.means, aes(x=sstmean, y=richness.fish, col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Mean annual SST (?C)", y=expression(paste("Species richness")) ) 

# plot all of it together
windows(15,3)
cowplot::plot_grid( a,e,b,c,d, labels="AUTO", ncol=5,
                    align="v",
                    label_x = c(0.05,0.15,0.15,0.15,0.15) )


### repeat the above graphing, but instead of showing trends with SST just show habitat differences
# graph fish abundance by habitat type
a <- ggplot(pred.means, aes(x=habitat, y=(total.cpua.fish+0.01), fill = habitat)) + 
  geom_boxplot( width=0.5 ) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  scale_y_log10(breaks=c(0.001,0.01,0.1,1)) +
  labs( x = "Habitat", y=expression(paste("Fish catch per ", m^2)) ) 
# graph fish species richenss by habitat type
b <- ggplot(pred.means, aes(x=habitat, y=richness.fish, fill = habitat)) + 
  geom_boxplot( width=0.5 ) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Habitat", y="Fish species richness" )
# graph fish diversity by habitat type
c <- ggplot(pred.means, aes(x=habitat, y=ENSPIEfish, fill = habitat)) + 
  geom_boxplot( width=0.5 ) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Habitat", y="Effective number of species" )
# graph fish size by habitat type
d <- ggplot(pred.means, aes(x=habitat, y=top20, fill = habitat)) + 
  geom_boxplot( width=0.5 ) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Habitat", y="Total length (cm)" )
# graph fish biomass by habitat type
e <- ggplot(pred.means, aes(x=habitat, y=(biomass.area+0.01), fill = habitat)) + 
  geom_boxplot( width=0.5 ) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  scale_y_log10(breaks=c(0.1,1,10)) +
  labs( x = "Habitat", y=expression(paste("Fish biomass (g per ",m^2,")")) )

# plot all of it together
windows(9,3)
cowplot::plot_grid( a,b,e, labels=NA, ncol=3,
                    align="v",
                    label_x = c(0.05,0.15,0.15,0.15,0.15) )


plot( richness ~ log10(total.cpua.fish), pred.means)
plot( richness ~ ENSPIE, pred.means)

# show bivariate correlations
library(psych)
corr.graph <- pred.means %>%
  mutate( biomass = log10(biomass.area+0.01), cpua = log10(total.cpua.fish+0.01)) %>%#,
          # richness.fish=log10(richness.fish),ENSPIEfish=log10(ENSPIEfish)) %>%
  select( biomass, cpua, top20, 
          richness.fish, ENSPIEfish) 
windows()
pairs.panels( corr.graph[,-1] )

## run models for habitat differences for all of these
# consider also running GAMS for this
library( lme4 )
library( lmerTest )
# abundance
glm.abun <- lmer( log10(total.cpua.fish+0.001) ~  habitat + (1|Country), data=pred.means )
summary(glm.abun)

# diversity
glm.div <- lmer( ENSPIEfish ~  habitat + (1|Country), data=pred.means )
summary(glm.div)

# richness
glm.div <- lmer( richness.fish ~  habitat + (1|Country), data=pred.means )
summary(glm.div)

# length
glm.len <- lmer( top20 ~ habitat + (1|Country), data=pred.means )
summary(glm.len)

# biomass
glm.mass <- lmer( log10(biomass.area+0.01) ~  habitat + (1|Country), data=pred.means )
summary(glm.mass)


############################################################################################
# relationship between biomass and diversity and temperature
# merge predator summaries with environmental data
pred.env <- left_join( pred.all, sites, by=c("Country") )
# select relevant columns and remove duplicate entries
pred.env <- pred.env %>%
  select( Site, Site.Name, Country, Lat, Long, habitat, 
          cpua, ENSPIEfish, biomass.area, sstmean, hemi ) %>%
  distinct()
# get rid of rows without sstmean estimates
pred.env <- pred.env[ !is.na(pred.env$sstmean), ]
  

# plot biomass ~ diversity
ggplot( data=pred.env, aes( x=ENSPIEfish, y=log10(biomass.area+0.01) )) +
  geom_point() + geom_smooth()

# two bins for temperature (choose a cutoff based on differences in fish composition or predation rate...
#    current cutoff is 15.6 C)
ggplot( data=pred.env, aes( x=ENSPIEfish, y=log10(biomass.area+0.01), 
                            col=cut(sstmean,breaks=2), group=cut(sstmean,breaks=2) )) +
  geom_point( size=3, alpha=0.5 ) + geom_smooth( method="lm", se=F) +
  ylab( expression(paste(log[10],"(biomass)") )) + xlab( "Effective Number of Species (PIE)" ) +
  scale_color_manual( name="SST bins (?C)", values=c("blue","red") ) +
  theme(legend.justification=c(1,0), legend.position=c(1,0.77))

# maybe not surprising, since lots of cold water species are removed...what if we look at all species?
library(lmerTest)
pred.env$tempbin <- cut(pred.env$sstmean,breaks=2)
pred.env$tempbin <- factor( pred.env$tempbin, levels=c( "(15.6,26]","(5.3,15.6]" ) )
summary( lmer( log10(biomass.std+1) ~ ENSPIEfish * tempbin + (1|Country), data=pred.env ))




# sites with high abundance
pred.all[pred.all$cpua>3,]  # Texas, Urunga Lagoon

## write predator summaries to disk
write.csv(  pred.all, "Output Data/Bitemap_SEINE_summaries_20190322.csv", row.names=FALSE )
write.csv(  sites.nmds, "Output Data/Bitemap_BioORACLE_20190322.csv", row.names = FALSE )
