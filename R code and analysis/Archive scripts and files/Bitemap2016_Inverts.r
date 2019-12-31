#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Isolating invertebrate (esp. crabs) from seine and GoPro data
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2018.01.09
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################


###   REMOVE NSW1.1 AND 1.2 BECAUSE THESE WERE NOT CONDUCTED IN SEAGRASS

###  Sites where crabs should be important, but data do not exist
# Ireland (Nessa O'Connor) -- Nessa says she thinks Carcinus maenus was responsible for all predation event, although they lack the data
# Virginia (Eastern Shore) -- Callinectes sapidus should be an important player
# others to consider: Korea?



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
## DATA ON Predators                    ###
###########################################

## sources: seining, GoPros

## FISH SEINING DATA
seines <- read.csv( '../Data/Bitemap_Seine_ALL-DATA_20180109.csv', strip.white = TRUE, stringsAsFactors = FALSE )
names(seines)[4] <- "habitat"
# convert all vegetated and unvegetated sites to common categories
seines$habitat[ seines$habitat %in% c("Seagrass " )] <- "Seagrass"
seines$habitat[ seines$habitat %in% c("unveg","Unveg","Unvegetated" )] <- "Unveg"
seines <- droplevels(seines)
# for all organisms from France, UNC2, Wales, QLD2, QLD3, multiple lengths by 10 to convert from cm to mm
seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] <- 
  seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] * 10

# one genus listed as "Crab". Change this to Decapoda
seines$Genus[seines$Genus=="Crab"] <- "Decapoda"


# isolate inverts from the seines
seineinvert <- seines[  seines$Phylum!="Chordata" & !is.na(seines$Phylum), ]
invertgen <- unique( seineinvert$Genus )
sitax <- classification( invertgen, db="ncbi" )
sitax <- rbind(sitax)
# select which ranks to keep
taxkeepsi <- sitax[ sitax$rank %in% c("phylum","class","order","infraorder","family"), ]

# melt and recast taxonomy data
simelt <- melt(taxkeepsi)
sitaxtable <- dcast( simelt, query~rank, value.var="name" )
names(sitaxtable)[1] <- c("Genus")

# merge taxonomy with data
sinvert <- left_join( seineinvert, sitaxtable )


## GOPRO DATA
go <- read.csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv", strip.white = TRUE, stringsAsFactors = FALSE )
# get unique taxa for each site
go$sciname <- with(go, paste(Genus,Species))
# Add dummy column to contain number of individuals seen (always 1)
go$Abundance <- 1

# get summary of taxa by country and habitat
gosum <- go %>% 
  group_by( Country, habitat, Genus ) %>%
  summarise( sum(Abundance)  )

# get taxonomy for all unique genera
gen <- unique( go$Genus )
gen <- gen[ !is.na(gen) & gen!="" & gen!="NO ID"]
gotax <- classification( gen, db = "ncbi" )
gotax <- rbind(gotax)
# select which ranks to keep
taxkeep <- gotax[ gotax$rank %in% c("phylum","class","order","infraorder","family"), ]

# melt and recast taxonomy data
simelt <- melt(taxkeep)
gotaxtable <- dcast( simelt, query~rank, value.var="name" )
names(gotaxtable)[1] <- c("Genus")

# merge taxonomy with data
go <- left_join( go, gotaxtable )


###########################################
## Extract invertebrates
## Focus on crabs
###########################################

### Decapods from both data sets

## Seines
decaseine <- sinvert[ sinvert$order=="Decapoda", ]
decaseine <- decaseine %>%
  dplyr::select( Country,habitat,Genus,Species,class,family,order,phylum,Abundance )

## GOPRO
decago    <- go[ go$order=="Decapoda", ]
decago <- decago %>%
  dplyr::select( Country,habitat,Genus,Species,class,family,order,phylum,Abundance )

# merge these two together
deca <- full_join(decaseine,decago)

# summarize by site
decasum <- ddply( deca, .(Country,habitat,Genus,Species), summarise, Total=sum(Abundance) )
decasum <- ddply( deca, .(Country), summarise, Total=sum(Abundance) )

write.csv( data.frame(Country=decasum[!is.na(decasum$Total),] ), file = "Output Data/crustsum.csv",
           row.names=FALSE )

### Brachyurans from both data sets
## Seines
crabseine <- sinvert[ sinvert$infraorder=="Brachyura", ]
crabseine <- crabseine %>%
  dplyr::select( Country,habitat,Genus,Species,class,family,order,infraorder,phylum,Abundance )

## GOPRO
crabgo    <- go[ go$infraorder=="Brachyura", ]
crabgo <- crabgo %>%
  dplyr::select( Country,habitat,Genus,Species,class,family,order,infraorder,phylum,Abundance )

# merge these two together
crab <- full_join(crabseine,crabgo)

# summarize by site
crabsum <- ddply( crab, .(Country,habitat,Genus,Species), summarise, Total=sum(Abundance) )
ddply( crab, .(Country,habitat), summarise, Total=sum(Abundance) )
crabsum <- ddply( crab, .(Country), summarise, Total=sum(Abundance) )

crabsites <- sort( unique(crab$Country) )

write.csv( data.frame(Country=crabsites[!is.na(crabsites)]), file = "Output Data/crabsites.csv",
           row.names=FALSE )
write.csv( data.frame(Country=crabsum[!is.na(crabsum$Total),] ), file = "Output Data/crabsum.csv",
           row.names=FALSE )
