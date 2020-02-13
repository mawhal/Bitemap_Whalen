#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Cleaning fish community data, calculating summaries
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 17 August 2019
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



# load libraries
# data.frame manipulation and joining
library(tidyverse)



## read data

# seine abundance + biomass
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )


# trait information that includes video data
video <- read_csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv" )
video <- video %>% 
  select( Country, habitat, Genus, Species ) %>% 
  mutate( habitat=tolower(habitat) )
video$habitat[ video$habitat == "seagrass" ] <- "Seagrass"
video$habitat[ video$habitat == "unveg" ] <- "Unvegetated"
# fix some names
video$Genus <- gsub( "Athrina", "Atherina", video$Genus )
video$Genus <- gsub( "Eucinostramus", "Eucinostomus", video$Genus )
video$Genus <- gsub( "Lactophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Lactrophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Alters", "Aluterus", video$Genus )
video$Genus <- gsub( "Sphyraema", "Sphyraena", video$Genus )
video <- video %>% 
  unite( SPECIES_NAME, Genus, Species, sep=" ", remove=F ) 

genfam <- seine %>% select(family,Genus) %>% distinct()
video1 <- left_join(video, genfam)

# lookup taxa with taxize
taxfam <-  read_csv( "Output Data/families_taxize_RDA.csv" )

famfam <- bind_rows( genfam, taxfam )
video <- left_join(video, famfam)

# vet taxa without family info
video %>% filter( Country %in% c("Mexico (ICML)","Mexico (ICML2)")) %>% 
  arrange( Country, family, Genus )
filter( video, is.na(family) )
video$family[ video$Genus %in% "Aluterus"] <- "Monacanthidae"
video$family[ video$Genus %in% "Caesio"] <- "Caesionidae"

video <- filter( video, !is.na(family) )


## summarize seine data
# sum of catch per seine
seine.totals <- seine %>% 
  dplyr::group_by(Country,Site.Name,Lat,Long,habitat,Date,Time,family,Genus,Species,SPECIES_NAME) %>% 
  dplyr::summarize( cpua=sum(cpua,na.rm=T) )
# because we want taxon specific catch we need to include zeros when the taxon was not caught
seine.all <- seine.totals %>% 
  dplyr::group_by(Country,Site.Name,Lat,Long,habitat,Date,Time,family,Genus,Species) %>% 
  tidyr::spread( SPECIES_NAME, cpua, fill=0 ) %>% 
  tidyr::gather( SPECIES_NAME, cpua, -Country,-Site.Name,-Lat,-Long,-habitat,-Date,-Time,-family,-Genus,-Species )
# avearage by site and habitat
mean.catch <- seine.all %>% 
  dplyr::group_by(Country,habitat,family,Genus,Species,SPECIES_NAME) %>% 
  dplyr::summarize( cpua=mean(cpua) )
  

# merge conusmer data with site information and taxonomy
mean.catch$family[ mean.catch$SPECIES_NAME == "Menticirrhus ophicephalus" ] <- "Sciaenidae"



## summarize data by family
catch.species <- mean.catch %>% 
  dplyr::group_by( Country, habitat, family,Genus,Species, SPECIES_NAME ) %>% 
  dplyr::summarize( P=sum(cpua) ) %>% 
  dplyr::mutate( P=ifelse(P>0,1,0) ) %>% 
  filter( P == 1 )


# merge with video data
video$P <- 1

catch.video <- bind_rows( catch.species, video )
catch.video <- catch.video %>% 
  distinct()


# remove unresolved taxa for certain families only (those for which other species are identified at the same site or nearby)
fam.sp <- catch.video %>% ungroup() %>% filter(Species %in% c("sp.","spp.")) %>% select(family)
catch.video %>% filter( family %in% fam.sp$family ) %>% arrange(family)
family c("Embiotocidae","Cottidae","Hyporhamphus","Labridae","Lethrinidae","Lutjanidae","Pholidae")
Genus c("Gobiosoma","Liza","Paralichthys")

# summarize to get family, number of taxa, number of sites, presence inside and outside of seagrass
taxon_summary <- catch.video %>% 
  ungroup() %>% 
  group_by( family ) %>% 
  summarize( taxa=length(unique(SPECIES_NAME)), sites=length(unique(Country)), habs=length(unique(habitat)), seagrass=ifelse(any(habitat=="Seagrass"),1,0) )


# read in the consumption rate correlations and add to the table
cap1 <- read_csv("Output Data/multivar_constr_taxa.csv") %>% 
  select(family,CAP1)

taxon_sum <- left_join( taxon_summary, cap1 )



# write to disk
write_csv( taxon_sum, "Output Data/taxon_list.csv")
