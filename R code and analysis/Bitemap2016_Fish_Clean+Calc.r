#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Cleaning fish community data, calculating summaries
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2019.10.10 -- separate multivariate analysis
#                       use this script only for calculations
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
# calculate different summaries of fish community that can be related to consumption intensity
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
seines.raw <- read.csv( '../Data/Bitemap_Seine_ALL-DATA_20191010.csv', stringsAsFactors = FALSE, strip.white = TRUE)
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


# consumption Rates from other script
rate.env <- read.csv( "Output Data/Bitemap_DATA_analysis_20180825.csv", stringsAsFactors = FALSE )
# below, we deal with site and habitat level data, not at the same scale
# summarize to average consumption rates at site and habitat level -- need to be careful with habitat because this has problems with traits
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
seines$Length[ seines$Country %in% c('France', 'USA (CA2)', 'USA (NC2)', 'Wales', 
                                     'Australia (QLD2)', 'Australia (QLD3)')] <- 
    seines$Length[ seines$Country %in% c('France', 'USA (CA2)', 'USA (NC2)', 'Wales', 
                                         'Australia (QLD2)', 'Australia (QLD3)')] * 10


# select relevant columns
seines <- seines %>%
  dplyr::select( Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, Q.PA, Phylum,
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
write_csv( traits, "Output Data/consumer_eat_squid.csv" )


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
  dplyr::select( SPECIES_NAME=sciName, phylum, family, eat.squid ) %>%  #, omnivory, herbivory
  distinct()


# merge trait information into seine data
seine <- left_join( seines.clean, traits.clean, by=c("SPECIES_NAME") )
# get rid of scallop
seine <- seine[seine$SPECIES_NAME != "Argopecten irradians",]
seine <- seine[seine$SPECIES_NAME != "Aplysia californica",]
# misidentified crab? change to Cancer irroratus
seine$SPECIES_NAME[ seine$SPECIES_NAME == "Crab n/a" ] <- "Cancer irroratus"
seine$phylum[ seine$SPECIES_NAME == "Cancer irroratus" ] <- "Arthropoda"
# some missing family information
seine$family[seine$Genus == 'Romaleon'] <- 'Cancridae'
seine$family[seine$Genus == 'Hyporhamphus'] <- 'Hemiramphidae'
seine$eat.squid[seine$Genus == 'Hyporhamphus'] <- 1
seine$family[seine$Genus == 'Arripis'] <- 'Arripidae'
seine$eat.squid[seine$Genus == 'Arripis'] <- 1
seine$family[seine$Genus == 'Gadus'] <- 'Gadidae'
seine$eat.squid[seine$Genus == 'Gadus'] <- 1
seine$family[seine$Genus == 'Tozeuma'] <- 'Hippolytidae'
seine$eat.squid[seine$Genus == 'Tozeuma'] <- 0
seine$family[seine$SPECIES_NAME == 'Cancer irroratus'] <- 'Cancridae'
seine$SPECIES_NAME[is.na(seine$family)] 

# write to disk
write_csv( seine, "Output Data/Bitemap_SEINE_clean.csv" )






# new data.frame, SWITCH for all consumers or just ones with eat.squid == 1
consumers <- seine
# consumers <- seine[ seine$eat.squid==1, ]

# How certain are we that we have a good estimate of species richness at each site?
# If a species is labeled with a genus or family designation, do we count one species? YES
# If a species is labelled as one potential member of a genus, but other known taxa
#      are present, should we count all species??
# aggregate species list for each site, then interrogate each list
by( consumers, consumers$Site.Name, function(z) sort(unique( z$SPECIES_NAME )) )
by( consumers, consumers$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# these all look good
# Sparidae
by( consumers[consumers$family=="Sparidae",], 
    consumers[consumers$family=="Sparidae",]$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# Cottidae
by( consumers[consumers$family=="Cottidae",], 
    consumers[consumers$family=="Cottidae",]$Country, function(z) sort(unique( z$SPECIES_NAME )) )

# standardize predator abundance by seined area and volume
consumers <- consumers %>%
  # calculate area and approximate volume seined
  # assuming the seined volume is a wedge-shaped object (triangular in cross-section) because seines pulled from depth to shore, hence the division by 2
  dplyr::mutate( area = Distance*width.transect, volume=(Distance*Depth)/2 * width.transect  ) %>%
  # standardize abundance by area and volume (Catch per unit effort) units are meters squared and cubed
  dplyr::mutate( cpua=Abundance/area, cpuv=Abundance/volume )

# get total abundance by Site and seine 
abundance.all <- consumers %>%
  dplyr::group_by( Site.Name, Country, Date, Time, Lat, Long, Distance, area, volume, habitat) %>%
  dplyr::summarise( total.cpua=sum(cpua, na.rm=T), total.cpuv=sum(cpuv, na.rm=T))

# NEED TO DO THIS FOR FISH SEPARATELY?
abundance.fish <- consumers %>%
  filter( phylum=="Chordata" ) %>%
  dplyr::group_by( Site.Name, Country, Date, Time, Lat, Long, Distance, area, volume, habitat) %>%
  dplyr::summarise( total.cpua.fish=sum(cpua, na.rm=T), total.cpuv.fish=sum(cpuv, na.rm=T))

# join these
abundance <- left_join( abundance.all, abundance.fish )
abundance$total.cpua.fish[ is.na(abundance$total.cpua.fish) ] <- 0
abundance$total.cpuv.fish[ is.na(abundance$total.cpuv.fish) ] <- 0

# Build family level abundance (+presence-absence) matrices for families 
# Note: totalCount from data.frame traits is not valid and must be standardized
# use consumers
families <- consumers %>%
  dplyr::group_by( Country, Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, width.transect,
            phylum, family ) %>%
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

# join with rate.mean (note again that individual seines DO NOT match up with consumption assays)
fam.cpua.rate <- left_join( fam.cpua.spread, rate.mean )
#




# calculate family richness (how many families at a site)
famrich   <- traits %>%
  dplyr::group_by(Site,method) %>%
  # filter( eat.squid==1 ) %>%
  dplyr::summarise( famrich=length(unique(family)) )
  
 
# famrich.squid <- ddply( family.squid, .(Site, method), summarize, 
#                         famrich.squid = length(unique(family)) )

# famrich <- full_join( famrich.all, famrich.squid )
famrich.geo <- left_join( famrich, sites )

# number of predator families (likely to eat squidpops) as a function of distance from equator
ggplot( data=famrich.geo, aes(x=abs(meanLat),y=famrich)) + geom_point() # India is low, but based on videos
# windows(4.5,3.5)
# by method
# ggplot( data=famrich.geo, aes(x=abs(meanLat),y=famrich, col=method)) + 
#   geom_smooth( ) + #geom_smooth( aes(group=1) ) + 
#   geom_point(size=3) + guides( fill = guide_legend(override.aes = list(linetype = 0,fill=NA)),
#                                color = guide_legend(override.aes = list(linetype = 0,fill=NA)) ) +
#   ylab("Number of families") + xlab("Degrees from equator")
# # by hemisphere
# ggplot( data=famrich.geo, aes(x=abs(meanLat),y=famrich, col=hemi, lty=hemi)) + 
#   geom_smooth( data=famrich.geo[famrich.geo$method == "seine",],se=T ) + #geom_smooth( aes(group=1) ) + 
#   geom_point( size=3 ) + guides( fill = guide_legend(override.aes = list(fill=NA)),
#                                  color = guide_legend(override.aes = list(fill=NA)) ) +
#   ylab("Number of families") + xlab("Degrees from equator") + ylim(c(0,20))

# presence-absence matrix for families (all and thought to eat squid)
# include all data from seining, video, and transects
trait.rate <- left_join( traits, rate.mean2 )
fam.pres <- trait.rate %>%
  # filter( eat.squid==1 ) %>%
  select( Site,method,hemi,basin,coast,meanLat,sstmean,rate,family,totalCount ) %>%  #,omnivory,herbivory
  mutate( totalCount = ifelse(totalCount>0,1,0) ) %>%
  distinct()

# spread out the families
fam.spread <- fam.pres %>%
  # group_by(Site) %>%
  spread( family, totalCount, fill=0 )


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
group[ rowSums(fam.data) <= 2, ]   # NC2 - two fish species recovered from seines




# merge with site data
consumer.sites      <- left_join( consumers, sites )






###############################
# Consumer size distribution #
#############################

# How to deal with all of the fish for which we have abundance but not size?
# For biomass, I use the average length
# but here, I only include data with measured length
#
# only include rows with a length estimate (pred biomass assumed average lengths for preds not measured...see below)
consumer.size <- consumer.sites[ !is.na(consumer.sites$length.measured), ]


# # use ggplot to make violin plots for each site x habitat combination
# ggplot( consumer.size, aes(x=habitat, y=log10(Length) )) + facet_wrap(~Country) + 
#   geom_boxplot(width=0.5)   # Note that USA (VA) and USA (NC2) lack data from unvegetated sediments
# ggplot( consumer.size, aes(x=habitat, y=log10(Length) )) + geom_violin() + facet_wrap(~Country)
# ggsave( "Figs/Length_habitat_site_violin.png", width = 10, height = 7, dpi = 600 )

# # pred length by latitude
# ggplot( consumer.size, aes(x=abs(Lat), y=log10(Length) )) + facet_grid(hemi~habitat) +
#   geom_point() +    geom_smooth(method='lm')
# ggplot( consumer.size, aes(x=abs(Lat), y=log10(Length) )) + geom_point() + facet_wrap(~habitat) + 
#   geom_smooth(  )
# 
# # quantile regression
# ggplot( consumer.size, aes(x=abs(Lat), y=log10(Length) )) + geom_point() + facet_wrap(~habitat) + 
#   geom_quantile()
# ggplot( consumer.size, aes(x=abs(Lat), y=log10(Length) )) + geom_point() + facet_wrap(~habitat) + 
#   geom_quantile(quantiles=seq(0.05, 0.95, by = 0.05))
# ggplot( consumer.size, aes(x=abs(Lat), y=(Length) )) + geom_point() + facet_wrap(~habitat) + 
#   geom_quantile(quantiles=c(0.9,0.8,0.5,0.2,0.1))
# ggplot( consumer.size, aes(x=abs(Lat), y=(Length) )) + geom_point() + facet_grid(hemi~habitat) + 
#   geom_quantile(quantiles=c(0.9,0.8,0.5,0.2,0.1))

# calculate the median size
consumer.size.median <- consumer.size %>% 
  dplyr::group_by(Country,Lat,habitat,hemi) %>% 
  dplyr::summarize( Length=median(Length,na.rm=T) )
ggplot( consumer.size.median, aes(x=abs(Lat), y=Length, col=habitat)) + 
  # geom_smooth(method='lm',se=T) + 
  geom_smooth() +
  geom_point()

# little difference between habitats in fish size
# consumer.size.median <- ddply( consumer.size, .(Country,Lat), summarize, Length=median(Length,na.rm=T) )
# ggplot( consumer.size.median, aes(x=abs(Lat), y=Length)) + 
#   geom_smooth() + geom_point() # + geom_smooth(method='lm',se=T) 

# calculate median size for each seine
consumer.med.len <- ddply( consumer.size, .(Country,Site.Name,Lat,sstmean,Date,Time),
                       summarize, medLength=median(Length,na.rm=T) )
# windows(4,4)
# ggplot( consumer.med.len, aes(x=abs(Lat), y=medLength )) + geom_point() + geom_smooth() +
#   xlab('degrees latitude from equator') + ylab('median length of fish per seine')
# ggplot( consumer.med.len, aes(x=sstmean, y=medLength )) + geom_point() + geom_smooth() +
#   xlab('SST') + ylab('median length of fish per seine')


# 
# # consider that the relationship could change with hemisphere
# ggplot( consumer.size.median, aes(x=abs(Lat), y=Length, col=hemi)) + 
#   geom_smooth() + geom_point() 
# 
# # temperature might be a better indicator anyway
# ggplot( consumer.sites, aes(x=sstmean, y=log10(Length) )) + facet_grid(~habitat) +
#   geom_point() +    geom_smooth(method='lm') # slight increase in size with temperature?? 
# ggplot( consumer.sites, aes(x=sstmean, y=(Length) )) + facet_grid(~habitat) +
#   geom_point() +    geom_quantile()
# 
# # look at size of top XXth percentile of size
# # calculate size of 80th percentile 
# consumer.top20 <- consumer.size %>%
#   dplyr::group_by( Country,Site.Name,Distance,Date,Time,Lat,Long,habitat ) %>%
#   dplyr::summarize( top20=quantile(Length,probs=0.8), meanL=mean(Length,na.rm=T) )



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




#######################
# Consumer diversity #
#####################


# combine seines and sites
consumer.site <- left_join( consumers, siteGPS )


# summarize seine data based on counts of each taxon
consumer.abund <- ddply( consumer.site, .(Country,meanLat,meanLong,habitat,phylum,Genus,Species), 
                      summarise, Abundance=sum(Abundance) )
# summarize seine data based on each indivdual seine (use site.name, date, time)
consumer.abund.each <- ddply( consumer.site, .(Country,Site.Name,Distance, Depth, width.transect,Date,Time,Lat,Long,habitat,phylum,Genus,Species), 
                      summarise, Abundance=sum(Abundance) )



# for species diversity metrics, we need counts of each species in each habitat
# # omit invertebrates
# consumer.abund <- seine.abund[ seine.abund$Phylum=="Chordata", ]
# consumer.abund.each <- seine.abund.each[ seine.abund.each$Phylum=="Chordata", ]

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
consumer.totals.each <- ddply( consumer.abund.each, .(Country,Site.Name,Distance,Depth,width.transect,Date,Time,Lat,Long,habitat),
                      summarize, Total = sum(Abundance), TotalFish = sum(Abundance[phylum=="Chordata"]) )

consumer.comm.each <- left_join(consumer.abund.each,consumer.totals.each)
head(consumer.comm.each[ with(consumer.comm.each,order(Country,Site.Name,Distance,Date,Time,habitat)), ])

# remove fish total information for taxa that are not fish
consumer.comm.each$TotalFish[ consumer.comm.each$phylum != "Chordata" ] <- NA
# calculate the proportion of each species in each site
consumer.comm.each$prop <- consumer.comm.each$Abundance/consumer.comm.each$Total
consumer.comm.each$propFish <- consumer.comm.each$Abundance/consumer.comm.each$TotalFish
# calculate ENSPIE and species richness for each seine
# for richness, inlcude all taxa, regardless of whether species known (this will inflate richess at some sites)
# exceptions -- Lethrinus sp. (QLD2), Cottoidea sp.
consumer.comm.each$species <- with(consumer.comm.each, paste(Genus, Species) )
ENSPIE.seine <- consumer.comm.each %>%
  dplyr::group_by( Country,Site.Name,Distance,Date,Time,Lat,Long,habitat ) %>%
  dplyr::summarise( ENSPIE = ENSPIE(prop), ENSPIEfish = ENSPIE(propFish),
                    richness = length(Total), richness.fish = length(TotalFish))


  
  
# calculate total abundance for each site (lump all replicate seines together)
consumer.totals <- ddply( consumer.abund, .(Country,meanLat,meanLong,habitat),
       summarize, Total = sum(Abundance), TotalFish = sum(Abundance[phylum=="Chordata"]),
       richness=length(abundance), richness.fish=length(Abundance[phylum=="Chordata"]) )
# get rid of NA rows
consumer.totals <- consumer.totals[!is.na(consumer.totals$Country),]
# include totals as a column in consumer.abund
consumer.comm <- left_join(consumer.abund.each,consumer.totals)
# calculate the proportion of each species in each site
consumer.comm$prop     <- consumer.comm$Abundance/consumer.comm$Total
consumer.comm$propFish <- consumer.comm$Abundance/consumer.comm$TotalFish
# calculate ENSPIE for each site
ENSPIE.site <- ddply( consumer.comm, .(Country,meanLat,meanLong,habitat), 
       summarize, ENSPIE = ENSPIE(prop), ENSPIEfish = ENSPIE(propFish), 
       richness=mean(richness), richness.fish=mean(richness.fish) )
# get rid of NA rows
ENSPIE.site <- ENSPIE.site[!is.na(ENSPIE.site$Country),]

# calculate effect size as difference between ENSPIE estimates
# make all Unveg sites negative
ENSPIE.site$ENSPIEdiff <- ENSPIE.site$ENSPIE
ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] <- -(ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] )
ENSPIE.diff <- ddply( ENSPIE.site, .(Country,meanLat,meanLong), summarise, ENSPIEdiff=sum(ENSPIEdiff) )


# ggplot( ENSPIE.site, aes(x=habitat,y=ENSPIE)) + geom_point() + 
#   geom_line(aes(group=Country)) + facet_wrap(~Country)
# ggplot( ENSPIE.site, aes(x=abs(meanLat),y=ENSPIE)) + geom_point() + facet_wrap(~habitat, ncol=2) + geom_smooth()
# ggplot( ENSPIE.diff, aes(x=abs(meanLat),y=ENSPIEdiff)) + geom_point() + geom_smooth()



# ggplot(ENSPIE.site, aes(x=habitat, y=ENSPIE, fill = habitat)) + geom_boxplot(width=0.5, notch=FALSE) +
#   geom_point(alpha=0.1) +
#   # facet_grid(.~Year) +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14),
#         legend.position = "none") +
#   labs(title="Diveristy of Select Predators By Habitat Type",
#        x = "Type of Habitat", y="Diversity\n(effective number of species)") 









####################################################################################
# CONSUMER BIOMASS                                                                #
##################################################################################

# read biomass coefficients
coef <- read_csv("../Data/Fish Biomass + Traits/20160710_RLS_biomass_coefs_20170921.csv")
# omit rows from coef for which we have no estimates
coef <- coef[ !is.na(coef$A), ]
# remove duplicates from coef
coef <- coef[ !duplicated(coef), ]



# # identify species not in character list
# library(rfishbase)
# dlookup <- sort(unique(consumers$SPECIES_NAME))
# lookup <- dlookup[ !(dlookup %in% coef$SPECIES_NAME) ]
# 
# # use rfishbase to look up taxa not represented in the reef life survey biomass conversion data.frame (coef)
# fishbaseLW <- length_weight(lookup)
# sort(unique(fishbaseLW$Species)) 
# unrepresented <- sort(unique( dlookup[!(dlookup %in% c(fishbaseLW$Species, as.character(coef$SPECIES_NAME)) )]))
# sealifebaseLW <- length_weight(unrepresented,server = "sealifebase")
# baseLW <- bind_rows( fishbaseLW, sealifebaseLW )
# 
# 
# write_csv( baseLW, "../Data/Fish Biomass + Traits/baseLW.csv" )
# # write_csv( data.frame(unrepresented), "../Data/Fish Biomass/fishbaseMISSING.csv" )

# read the length-weight relationships looked up in fishbase
fishbaseLWread <- read_csv( "../Data/Fish Biomass + Traits/fishbaseLW_20170921.csv" )


baseLW <- read_csv( "../Data/Fish Biomass + Traits/baseLW.csv" )


# average by taxon
meanbaseLW <- baseLW %>% 
  dplyr::group_by(Species) %>% 
  dplyr::summarize(  a=mean(a,na.rm=T),aTL=mean(aTL,na.rm=T),b=mean(b,na.rm=T) )
  
# accept the estimates for a (intercept) that account for standard vs total length
for(i in 1:nrow(meanbaseLW) ){
  if( !is.na(meanbaseLW$aTL[i]) ) meanbaseLW$a[i] <- meanbaseLW$aTL[i]
}
# rename sciname to Species
names(meanbaseLW) <- c( "SPECIES_NAME", "A", "aTL", "B" )

# read in table for taxa that were looked up manually in fishbase
fishbaseManual <- read_csv( "../Data/Fish Biomass + Traits/Bitemap_FishbaseLW_manual.csv" )[,1:3]
names(fishbaseManual) <- c( "SPECIES_NAME", "A", "B" )

# combine all data from fishbase, including those provided by Reef Life Survey
fishbaseFull <- full_join( meanbaseLW, fishbaseManual )
biomassCoef <- full_join(  coef, fishbaseFull )
biomassCoef <- biomassCoef %>%
  dplyr::select( SPECIES_NAME, A, B ) %>% 
  dplyr::arrange( SPECIES_NAME )
length(unique(biomassCoef$SPECIES_NAME))
# omit duplicated rows (separate estimates for a taxon), default to use Reef Life Survey estimates for consistency
biomassCoef <- biomassCoef %>% 
  distinct()

biomassCoef <- biomassCoef %>% 
  dplyr::group_by( SPECIES_NAME ) %>% 
  dplyr::summarize( A=mean(A,na.rm=T), B=mean(B,na.rm=T) )




# match seine data with biomass conversion coefficients
consumer_biomass <- left_join(consumer.site, biomassCoef, by="SPECIES_NAME" )

# which taxa don't have length-weight regression estimates
miss <- consumer_biomass[ is.na(consumer_biomass$A), ]
sort(unique(miss$SPECIES_NAME))
sort(unique(miss$Genus))

# use existing L-W relationships in the dataset to estimate missing relationships using averages at the genus level
missGenus <- consumer_biomass %>% 
  dplyr::filter(Genus %in% sort(unique(miss$Genus))) %>% select(Genus,A,B) %>% distinct() %>%
  dplyr::group_by( Genus ) %>% 
  dplyr::summarize( A=mean(A,na.rm=T), B=mean(B,na.rm=T) )
missGenus <- na.omit(missGenus)
# replace missing values with new genus-level averages where possible
for(i in missGenus$Genus ){
  consumer_biomass$A[ consumer_biomass$Genus == i & is.na(consumer_biomass$A) ] <- missGenus$A[missGenus$Genus==i]
  consumer_biomass$B[ consumer_biomass$Genus == i & is.na(consumer_biomass$B) ] <- missGenus$B[missGenus$Genus==i]
}

# repeat the genus-level exercise at the family level
miss <- consumer_biomass[ is.na(consumer_biomass$A), ]
missFam <- consumer_biomass %>% 
  dplyr::filter( family %in% sort(unique(miss$family)) ) %>% 
  select(family,A,B) %>% 
  distinct() %>%
  dplyr::group_by( family ) %>% 
  dplyr::summarize( A=mean(A,na.rm=T), B=mean(B,na.rm=T) )
missGenus <- na.omit(missGenus)
# replace missing values with new family-level averages where possible
for(i in missFam$family ){
  consumer_biomass$A[ consumer_biomass$family == i & is.na(consumer_biomass$A) ] <- missFam$A[missFam$family==i]
  consumer_biomass$B[ consumer_biomass$family == i & is.na(consumer_biomass$B) ] <- missFam$B[missFam$family==i]
}

consumer_biomass %>% 
  filter( SPECIES_NAME %in% sort(unique(miss$SPECIES_NAME))) %>% 
  dplyr::group_by( SPECIES_NAME, Site ) %>% 
  dplyr::summarize( Total = sum(Abundance) )




# remove extraneous columns
consumer_biomass <- consumer_biomass %>%
  dplyr::select(Site.Name, Lat, Long, habitat, Date, Time, Depth, 
                Distance, width.transect, area, volume, Q.PA, phylum, 
                family, Genus, Species, Length, Abundance, Country, 
                SPECIES_NAME, cpua, A, B )


# rules for taxa without L-W regressions (ignore for now)
# rules for taxa with only counts and not lengths (e.g., some shrimp taxa)





## Calculate Biomass from Length-Weight regressions
# create biomass column -- note that I've divided all lengths by 10 to convert mm to cm, which is the unit 
consumer_biomass$biomass <- with(consumer_biomass, Abundance*(A*(Length/10)^B))

# which ones have no estimate
consumer_biomass[ is.na(consumer_biomass$biomass), ]
consumer_biomass %>% 
  dplyr::filter( is.na(biomass) ) %>% 
  dplyr::group_by( phylum, family, SPECIES_NAME ) %>% 
  dplyr::summarize( N=length(Abundance), size=mean(Length,na.rm=T) )

# standardize biomass by seine distance
consumer_biomass <- consumer_biomass %>%
  mutate( biomass.area = biomass/area, biomass.vol = biomass/volume )

## write to disk
write_csv( consumer_biomass, "Output Data/Bitemap_seine_abundance_biomass.csv" )



consumer_biomass %>%
  filter( Genus=="Lagodon") %>%
  dplyr::group_by( Country, Genus, Species) %>% 
  dplyr::summarise( sum=sum(Abundance))
consumer_biomass %>%
  filter( Country=="USA (NC2)") %>%
  dplyr::group_by( Genus, Species) %>% 
  dplyr::summarise( sum=sum(Abundance))
consumer_biomass %>%
  filter( Country=="USA (NC)") %>%
  dplyr::group_by( Genus, Species) %>% 
  dplyr::summarise( sum=sum(Abundance))

# Calculate total fish biomass in each seine
biomass_total <- ddply( consumer_biomass, .(Site.Name,Country,Date,Time,Lat,Long,Distance,area, volume,habitat), 
                        summarize, biomass.area=sum(biomass.area, na.rm=TRUE) )
sort(unique(biomass_total$Country))
sort(unique(predators$Country))


# # graph biomass by habitat type
# ggplot(biomass_total, aes(x=habitat, y=log10(biomass.area+0.01), fill = habitat)) + geom_boxplot(width=0.5,notch=TRUE) +
#   geom_point(alpha=0.1) +
#   # facet_grid(.~Year) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14),
#         legend.position = "none") +
#   labs(title="Fish biomass by habitat type",
#        x = "Type of Habitat", y="Biomass") 
#   
# 
# # total biomass by Latitude
# ggplot(biomass_total, aes(x=Lat, y=log10(biomass.area+0.01))) + geom_point() + facet_wrap(~habitat,ncol=1)
# ggplot(biomass_total, aes(x=abs(Lat), y=log10(biomass.area+0.01))) + geom_point() + facet_wrap(~habitat,ncol=1) +
#   geom_smooth(method='lm')
# 
# ggplot( biomass_total, aes(x=abs(Lat),y=log10(biomass))) + geom_point() + facet_wrap(~habitat, ncol=2) + geom_smooth()
# 
# # calculate average biomass by seine and habitat
# biomass_mean <- ddply( biomass_total, .(Site.Name,Country,Date,Lat,Long,Distance,habitat), 
#                        summarize, biomass.area=mean(biomass.area, na.rm=TRUE), 
#                        biomass.vol=mean(biomass.vol,na.rm=TRUE) )
# ggplot(biomass_mean, aes(x=abs(Lat), y=log10(biomass.area+0.01))) + geom_point() + facet_wrap(~habitat,ncol=2) +
#   geom_smooth(method='lm')


# calculate average biomass by site and habitat
biomass_site <- ddply( biomass_total, .(Country,habitat), 
                       summarize, meanLat=mean(Lat), 
                       biomass.area=mean(biomass.area,na.rm=TRUE) )

ggplot(biomass_site, aes(x=abs(meanLat), y=log10(biomass.area+0.01))) + geom_point() + 
  facet_wrap(~habitat,ncol=2) #+ geom_smooth(se=F) + geom_smooth(method='lm')

# # graph biomass by habitat type
# ggplot(biomass_site, aes(x=habitat, y=log10(biomass.area+0.01), fill = habitat)) + 
#   geom_boxplot(width=0.5,notch=TRUE) +
#   geom_point(alpha=0.1) +
#   # facet_grid(.~Year) +
#   theme_classic() +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14),
#         legend.position = "none") +
#   labs(title="Fish biomass by habitat type",
#        x = "Type of Habitat", y="Biomass") 




# get family level biomass for each site
# first get total biomass within each family for each seine
fam.biomass <- consumer_biomass %>%
  dplyr::group_by( Country, Site.Name, Lat, Long, habitat, Time, Depth, Distance, width.transect, family ) %>%
  dplyr::summarise( biomass.area = sum(biomass.area,na.rm=T) )

# now get average biomass at the level of seining sites
fam.mass.mean <- fam.biomass %>%
  dplyr::group_by( Country, Site.Name, Lat, Long, habitat, family ) %>%
  dplyr::summarise( mean.biomass = mean(biomass.area,na.rm=T) ) %>%
  ungroup()

# join with rate.mean (note again that individual seines DO NOT match up with consumption assays)
fam.mass.rate <- left_join( fam.mass.mean, rate.mean )
#
# # look family by family at the relationship between abundance and consumption rate
# # windows(15,12)
# ggplot( fam.mass.rate, aes(x=log(mean.biomass),y=rate, col=habitat) ) + facet_wrap(~family, scales="free") +
#   geom_point() + geom_smooth(method='lm', aes(group=1))

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

# # look family by family at the relationship between abundance and consumption rate
# windows(10,7)
# seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
# unveg.color    <- "#d8b365"
# ggplot( fam.mass.rate.sel, aes(x=log(mean.biomass),y=rate, col=habitat) ) + facet_wrap(~family, scales="free_x") +
#   # geom_smooth(method='glm', aes(group=1), method.args=list(family=quasibinomial), col='black', se=T) +
#   stat_ellipse( aes(group=1), col='black' ) +
#   # geom_smooth( aes(group=1), se=F, col='black') +
#   # geom_smooth( method='lm', aes(group=1), se=F, col='black' ) +
#   geom_point(size=3) + ylim(c(0,1)) + 
#   xlab( expression(paste(log[10],"( mean biomass (g/",m^2,") )")) ) +
#   ylab( "consumption rate" ) +
#   scale_color_manual( values=c(seagrass.color,unveg.color)) +
#   theme( text = element_text(size=12),
#          axis.text=element_text(size=10) )


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
  do( fitFam = glm(rate ~ log10(mean.biomass+0.001), 
                           family="quasibinomial", data=.) )
# get coefficients and model summaries
library(broom)
fam.coef <- tidy(fam.glm, fitFam)
fam.coef %>% 
  filter( term=="log10(mean.biomass + 0.001)" ) %>% 
  select( -term ) %>% 
  arrange(estimate)
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

# two bins for temperature (choose a cutoff based on differences in fish composition or consumption rate...
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
