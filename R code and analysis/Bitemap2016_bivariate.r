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



########################################################################################
########################################################################################



# library(profvis)  - can look at places where the code is slow, bottlenecks, etc.
# profvis({

# load libraries
# data.frame manipulation and joining
library(tidyverse)
library(plyr)
library(ggrepel) # for plotting with text
# library(reshape2)
library(cowplot) # for arranging multiple plots
# geospatial data
# library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(vegan)
library(viridis)
library(here)



## read data

# seine abundance + biomass
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )

# environmental data from remote sensing and oceanographic expeditions
oracle <- read_csv( "Output Data/Bitemap_BioORACLE_20190107.csv")[,1:36]

# squidpop consumption rates
pop <- read_csv( "Output Data/Bitemap_rate.env.20190423.csv" )
pop <- pop %>% 
  dplyr::group_by( Country,Site,habitat ) %>% 
  dplyr::summarise( sstmean=mean(sstmean), temp=mean(temp), rate=mean(rate) )

# trait information that includes video data
video <- read_csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv" )
video <- video %>% 
  select( Country, habitat, Genus ) %>% 
  mutate( habitat=tolower(habitat) ) %>% 
  distinct()
video$habitat[ video$habitat == "seagrass" ] <- "Seagrass"
video$habitat[ video$habitat == "unveg" ] <- "Unvegetated"
video$Genus[ video$Genus=="Athrina" ] <- "Atherina"

genfam <- seine %>% select(family,Genus) %>% distinct()
video1 <- left_join(video, genfam)

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
taxfam <-  read_csv( "Output Data/families_taxize_RDA.csv" )

famfam <- bind_rows( genfam, taxfam )
video <- left_join(video, famfam)
video <- filter( video, !is.na(family) )


## summarize seine data
# sum of catch per seine
seine.totals <- seine %>% 
  dplyr::group_by(Country,Site.Name,Lat,Long,habitat,Date,Time,SPECIES_NAME) %>% 
  dplyr::summarize( cpua=sum(cpua,na.rm=T) )
# because we want taxon specific catch we need to include zeros when the taxon was not caught
seine.all <- seine.totals %>% 
  dplyr::group_by(Country,Site.Name,Lat,Long,habitat,Date,Time) %>% 
  tidyr::spread( SPECIES_NAME, cpua, fill=0 ) %>% 
  tidyr::gather( SPECIES_NAME, cpua, -Country,-Site.Name,-Lat,-Long,-habitat,-Date,-Time )
# avearage by site and habitat
mean.catch <- seine.all %>% 
  dplyr::group_by(Country,Site.Name,habitat,SPECIES_NAME) %>% 
  dplyr::summarize( cpua=mean(cpua) )
  
# spread out seine data
mean.catch %>% 
  dplyr::group_by(Country,Site.Name,habitat) %>% 
  tidyr::spread( SPECIES_NAME, cpua, fill=0 )

# merge conusmer data with site information and taxonomy
tax <- seine %>% 
  select( phylum, family, SPECIES_NAME ) %>% 
  distinct()
mean.catch.tax <- left_join( mean.catch, tax )

site <- select( oracle, Site, Country )

mean.catch.site <- right_join( site, mean.catch.tax )

## summarize data by family
catch.fam <- mean.catch.site %>% 
  dplyr::group_by( Site, habitat, family ) %>% 
  dplyr::summarize( P=sum(cpua) ) %>% 
  dplyr::mutate( P=ifelse(P>0,1,0) )


# merge with video data
video <- left_join(site,video)
video <- select(video, Site, habitat, family)
video <- video %>% filter( !is.na(family) )
video$P <- 1

catch.video <- bind_rows( catch.fam, video )
catch.video <- catch.video %>% 
  distinct()

# because we are merging seine and video data in some cases, we get discrepancies
# for example: Sparids were caught on video but not in the seine in Croatia
catch.video %>% 
  dplyr::group_by(Site,habitat,family) %>% 
  dplyr::summarise( count=length(P) ) %>% 
  filter( count>1 )

# since there are only ones and zeros, we can add up the presence values
catch.video <- catch.video %>% 
  dplyr::group_by(Site,habitat,family) %>% 
  dplyr::summarise( P=sum(P) ) 

# spread out again
fam.catch.wide <- catch.video %>% 
  group_by( Site, habitat ) %>% 
  spread( family, P, fill=0 )


# can average across habitat types, too

fam.data <- fam.catch.wide[,-c(1:2)]

fam.meta <- fam.catch.wide[,1:2]
fam.meta <- left_join( fam.meta, pop )

fam.meta <- fam.meta %>%  tidyr::unite( SH, Site, habitat, remove=FALSE )

rownames(fam.data) <- fam.meta$SH

# write to disk
fam.all <- bind_cols(fam.meta,fam.data)
write_csv( fam.all, "Output Data/consumer_presence_wide.csv")

## add ordination results and one trait, functional diversity
rdaunc <- read_csv( "Output Data/multivar_unconstr_sites.csv" )
rdaunc <- rdaunc %>% separate( SH, c("Site", "habitat") )
brda <- read_csv( "Output Data/biomass_RDAselected.csv" ) # 30 sites
brda <- brda %>%
  select( Site, habitat, biomass, cpua )
active <- read_csv( "Output Data/consumer_active_ratio2.csv" )
active <- active %>% 
  select( Site, habitat, act.ratio )
fds  <- read_csv( "Output Data/FunctionalDiversity_indices.csv" ) # only 21 sites
fds <- as.data.frame(fds)
fds[is.na(fds)] <- 0
cwm <- read_csv( "Output Data/FunctionalDiversity_CWM.csv" )

addon <- left_join(left_join(left_join(left_join(rdaunc,brda),active),fds),cwm) 
addon <- addon %>% 
  replace_na( replace=list(FRic=0,qual.FRic=0,FEve=0,FDis=0,RaoQ=0,FGR=0)  )

rate.mean <- left_join(fam.meta, addon)

rate.mean.pairs <- rate.mean %>% 
  mutate( log.cpua = log10(cpua+0.001),
          logit.rate = car::logit(rate) ) %>% 
  select( Site, habitat, `mean\nannual\nSST`=sstmean, MDS1, 
          activity=act.ratio.tax, `selected\nabundance`=log.cpua, 
          `consumption\nrate`=logit.rate, FRic, act )
# bivariate correlations
chart.Correlation <- source("chart.correlation.r")$value

chart.Correlation( rate.mean.pairs[,c("mean\nannual\nSST",
                                      "activity",
                                      "FRic",
                                      "MDS1",
                                      "selected\nabundance",
                                      "consumption\nrate") ],
                   histogram = F, method="spearman" )
