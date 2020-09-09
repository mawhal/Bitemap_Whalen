#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 



## Seine Data
# calculate different summaries of fish community that can be related to squidpop consumption intensity
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


### This script makes a variety of species list for the 2016 Bitemap project


library(tidyverse)


## read data

# seine abundance + biomass
seine <- read_csv( "R code and analysis/Output Data/Bitemap_seine_abundance_biomass.csv" )

# trait information that includes video data
video <- read_csv( "Data/Video Data/Bitemap_Video_Data_ALL.csv" )
# fix some names
video$Genus <- gsub( "Athrina", "Atherina", video$Genus )
video$Genus <- gsub( "Eucinostramus", "Eucinostomus", video$Genus )
video$Genus <- gsub( "Lactophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Lactrophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Alters", "Aluterus", video$Genus )
video$Genus <- gsub( "Sphyraema", "Sphyraena", video$Genus )
video$Genus <- gsub( "Chelio", "Cheilio", video$Genus )

phyla <- read_csv( "R code and analysis/Output Data/predator_families+phyla.csv" )



videosel <- video %>% 
  select( Country, habitat, Genus, Species ) %>% 
  mutate( habitat=tolower(habitat) ) %>% 
  filter( !is.na(Genus) ) %>% 
  distinct()
videosel$habitat[ videosel$habitat == "seagrass" ] <- "Seagrass"
videosel$habitat[ videosel$habitat == "unveg" ] <- "Unvegetated"
videosel$method <- "video"

genfam <- seine %>% select(family,Genus) %>% distinct()
video1 <- left_join( videosel, genfam )

# lookup taxa with taxize
taxa <- video1 %>% 
  filter( is.na(family) ) %>% 
  select( Genus ) %>% 
  distinct() %>%
  arrange( Genus ) %>% 
  unlist()
taxa <- na.omit(taxa)
taxa <- taxa[ !(taxa %in% c("NO ID","Not identified")) ]
# do the lookup
# library(taxize)
# options(ENTREZ_KEY = "a06da45de96c4b7d680015d4c5b46694d908")
# taxclass <- taxize::classification( taxa, db="ncbi", 
#                                     ENTREZ_KEY = "a06da45de96c4b7d680015d4c5b46694d908" )
# taxfam <- rbind(taxclass) %>% 
#   filter(rank=='family') %>% 
#   select( family=name, Genus=query )
# # write taxfam to disk so we don't have to look it up every time
# write_csv( taxfam, "Output Data/families_taxize_RDA.csv" )
taxfam <-  read_csv( "R code and analysis/Output Data/families_taxize_RDA.csv" )

famfam <- bind_rows( genfam, taxfam )
videofam <- left_join(videosel, famfam)

# vet taxa without family info
videofam %>% filter( Country %in% c("Mexico (ICML)","Mexico (ICML2)")) %>% 
  arrange( Country, family, Genus )
filter( videofam, is.na(family) )
videofam$family[ videofam$Genus %in% "Aluterus"] <- "Monacanthidae"
videofam$family[ videofam$Genus %in% "Caesio"] <- "Caesionidae"
videofam$family[ videofam$Genus %in% "Cheilio"] <- "Labridae"
videofam$phylum <- "Chordata"
videofam$phylum[ videofam$family %in% c("Cancridae","Palaemonidae")] <- "Arthropoda"

# summarize seine data to get abundance per site and habitat
seine.sum <- seine %>% 
  group_by( Country, Lat, Long, habitat, Date, Time, phylum, family, Genus, Species, area ) %>% 
  summarize( length=mean(Length,na.rm=T), catch.per.m2=sum(cpua,na.rm=T),
             biomass.per.m2=sum(biomass.area,na.rm=T) ) %>% 
  group_by( Country, habitat, phylum, family, Genus, Species ) %>% 
  summarize( area=sum(area),length=mean(length,na.rm=T), catch.per.m2=mean(catch.per.m2,na.rm=T),
             biomass.per.m2=mean(biomass.per.m2,na.rm=T) ) %>% 
  mutate( method="seine" )


# merge
sv <- full_join(seine.sum, videofam )


# combine names
sv <- sv %>% 
  mutate( taxon = paste(  Genus, Species ) )
sv$taxon <- gsub( " NA", "", sv$taxon )

# site info
sites <- read_csv( "R code and analysis/Output Data/Bitemap_BioORACLE_20190107.csv" )
sites <- sites %>% 
  select( Country, Site, meanLat, meanLong, basin, hemi, coast)
sv <- left_join( sv, sites )

svfinal <- sv %>% 
  ungroup() %>% 
  select( Country, Site, meanLat, meanLong, basin, hemi, coast, habitat, 
          phylum, family, taxon, method,
          area.seined.m=area, total.length.mean=length, catch.per.m2, biomass.per.m2, )

# write to disk
write_csv( svfinal, "R code and analysis/Output Data/Bitemap2016_consumer_seine_video_summary.csv" )

