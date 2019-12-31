#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# summarise taxa in video data and get traits
# Code by Matt Whalen
# updated 2018.01.11
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################


###   REMOVE NSW1.1 AND 1.2 BECAUSE THESE WERE NOT CONDUCTED IN SEAGRASS



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
library( tidyverse )
library(plyr)
library( taxize )
library( reshape2 )

###########################################
## DATA FROM GOPROS                     ###
###########################################


## GOPRO DATA
go <- read.csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv", strip.white = TRUE, stringsAsFactors = FALSE )
# get unique taxa for each site
go$sciname <- with(go, paste(Genus,Species))
# Add dummy column to contain number of individuals seen (always 1)
go$Abundance <- 1

# get summary of taxa by country
gosumcountry <- ddply( go, .(Country,sciname), summarise, Abundance=sum(Abundance) )



#########################
# Fish Trophic Levels   #
#########################
library(rfishbase)
fishSpecies <- sort( unique( gosumcountry$sciname ) )

fishEcology <- ecology( fishSpecies, fields = c("SpecCode", "FeedingType", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph") )

# sometimes DietTroph (based on whole range of food items) unavailable but FoodTroph (based on individual food items) more often is available
fishEcology$Trophic.Level <- fishEcology$DietTroph
for( i in 1:nrow(fishEcology) ) {
  if( is.na(fishEcology$DietTroph[i]) )   fishEcology$Trophic.Level[i] <- fishEcology$FoodTroph[i]
}

# which taxa are missing trophic levels?
trophMiss <- fishEcology[ which(is.na(fishEcology$Trophic.Level)), ]


## join Bitemap fish data with trait data from Reef Life Survey
# read in fish trait data from (source = RLS?)
fishtraits <- read.csv("../Data/Fish Biomass + Traits//Traits_all-species_edit.csv")
# create column for joining with fishEcology
fishtraits$sciname <- fishtraits$CURRENT_TAXONOMIC_NAME


# first, make sure all bitemap species are represented
fishINFO <- full_join(fishEcology,data.frame(sciname=fishSpecies))
# now, left_join so we don't include taxa that aren't represented in Bitemap
fishINFO <- left_join( fishINFO, fishtraits ) # lots of gaps
# select relevant columns
fishINFO <- fishINFO %>%
  dplyr::select( sciname, FeedingType, Trophic.group, Water.column, Diel.Activity, Habitat, Trophic.Level )



#################
# Fish Taxonomy #
#################
# add taxonomy for fishes
library(taxize)
# strip sp from unknown species before running classification
scisplit <- strsplit( fishINFO$sciname, " ", fixed=TRUE )
which(lapply(scisplit,length)==1)
# scisplit[[1460]] <- c("Clupeidae","sp.")
genusSpecies <- lapply(scisplit,function(z) data.frame(genus=z[1],species=z[2]))
genusSpecies <- do.call(rbind,genusSpecies)
fishINFO <- cbind(fishINFO,genusSpecies)

# omit sp designations
for(i in 1:nrow(fishINFO)){
  if(fishINFO$species[i] %in% c("sp.","spp.","sp")) fishINFO$species[i] <- ""
}

# recombine genus and species without sp. designations (spaces should not matter here)
fishINFO$sciname2 <- with(fishINFO, paste(genus,species) )
fishINFO <- fishINFO[ order(fishINFO$sciname2), ]
# get rid of duplicate entries
fishINFO <- fishINFO[ !duplicated(fishINFO$sciname2), ]

# run classification on all fish taxa
fishtax <- classification(fishINFO$sciname2,db="ncbi")


# only keep non-NA ones
which( unlist(lapply( fishtax, function(z) all(!is.na(z)) )) )
cf <- fishtax[which( unlist(lapply( fishtax, function(z) all(!is.na(z)) )) )]

# combine the results
# cannot cbind with NA values, so rbind
cr <- do.call( rbind,cf )
cr$unique <- rownames(cr)
cr$sciname2 <- unlist(lapply( strsplit( cr$unique, split = ".", fixed=TRUE ), function(z) z[1] ))

# melt and recast
cm <- melt(cr)
# omit no rank
cm <- cm[cm$rank!="no rank",]
ccast <- dcast( cm, sciname2~rank, value.var = "name" )

# merge taxonomy with traits
fishes <- full_join( fishINFO, ccast, "sciname2" )
head(fishes)

# write this to disk
# write.csv( fishes, "Output Data/fish_traits+taxonomy.csv", row.names=FALSE )


# Join the go pro data with traits
goprotraits <- left_join( gosumcountry, fishes, by="sciname")

write.csv( goprotraits, "Output Data/GoPro_Traits.csv", row.names=FALSE )


# merge goprotrais with other traits from seining data
# read data from seines
seines <- read.csv( "Output Data/Bitemap_seine_taxa_bySite_traits.csv", stringsAsFactors = FALSE )
names(seines)[5] <- "Total"
names(seines)


# select columns from goPro
names(goprotraits)
goprotraitssel <- goprotraits %>%
  select( Country, sciname,Total=Abundance,family,FeedingType,Trophic.group,Water.column,Diel.Activity,Habitat)
names(goprotraitssel)

traits <- full_join( seines,goprotraitssel )

# write to disk
write.csv( traits, "Output Data/Bitemap_seine+video_bySite_traits.csv", row.names=FALSE )
