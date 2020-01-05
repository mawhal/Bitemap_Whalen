#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Isolate fish families that are identified in RDA analyses
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 20 September 2019
#####################################################################


## load libraries
library(tidyverse) 


## read seine data
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )

## read environmental and squidpop data -- also has site abbreviations
rate.env <- read_csv( "Output Data/Bitemap_DATA_analysis_20180825.csv" )
rate.mean <- rate.env %>%
  dplyr::group_by( Country, Site, hemi, basin, coast, habitat, sstmean, parmean, MDS1  ) %>%
  dplyr::summarize( rate=mean(rate,na.rm=T), temp=mean(temp,na.rm=T), sal=mean(sal,na.rm=T) )
rate.mean2 <- rate.env %>%
  dplyr::group_by( Country, Site, hemi, basin, coast, sstmean, parmean, MDS1  ) %>%
  dplyr::summarize( rate=mean(rate,na.rm=T), temp=mean(temp,na.rm=T), sal=mean(sal,na.rm=T) )


## merge these two
sr <- left_join(seine,rate.mean)



## read the families that were positively associated with predation rates in RDA
##----
# these include all families from both abundance weighted and presence-absence analysis
# fam <- read_csv( "Output Data/multivar_families_positive.csv" )
# DELETED THIS BECAUSE WE WANT JUST THE PRESENCE-ABBSENCE ANALYSIS
##----
#
fam <- read_csv( "Output Data/multivar_constr_taxa.csv" )
fam.pos <- fam %>% filter( CAP1 >0 )
top <- fam.pos$family
# top <- fam$taxon[-12] # gets rid of Serranids when more taxa are included (e.g., video data also included for presence-absence analysis)
# the abundance-weight families were a bit more conservative
# top <- c( 'Penaeidae','Tetraodontidae','Sillaginidae','Palaemonidae','Sphyraenidae','Terapontidae','Hemiramphidae',
#          'Monodactylidae','Gerreidae','Haemulidae','Tetrarogidae','Sparidae' )


## filter seine data to only include taxa identified in multivariate analysis
sf <- sr[ sr$family %in% top, ]




# which taxa don't have length-weight regression estimates
miss <- sf[ is.na(sf$A), ]
sort(unique(miss$SPECIES_NAME))
# misslw <- rfishbase::length_weight(sort(unique(miss$SPECIES_NAME)))
# lwadd <- misslw %>% filter( !is.na(Species) ) %>%
#   select( SPECIES_NAME=Species, A=a, B=b ) %>%
#   group_by( SPECIES_NAME ) %>%
#   summarize( A=mean(A,na.rm=T), B=mean(B,na.rm=T) )
# 
# for( i in 1:4 ){
#   pred_biom$A[ pred_biom$SPECIES_NAME == lwadd$SPECIES_NAME[i] ] <- lwadd$A[i]
#   pred_biom$B[ pred_biom$SPECIES_NAME == lwadd$SPECIES_NAME[i] ] <- lwadd$B[i]
# }
# pred_biom[ is.na(pred_biom$A), ]

## write pred_biom to disk, could be useful for downstream analysis (e.g., traits weighted by abundance or biomass)
# write_csv( pred_biom, "Output Data/consumer_biomass_seine.csv" )


## Summarize the seine data -- need to get totals per seine, then mean per habitat:site, then we can merge with rate data
# Calculate total fish biomass in each seine
biomass_total <- sf %>%
  dplyr::group_by(Site.Name,Country,habitat,Date,Time,Lat,Long,area) %>%
  dplyr::summarize( biomass=sum(biomass.area, na.rm=TRUE), cpua=sum(cpua,na.rm=T) ) %>%
  dplyr::group_by( Country,habitat ) %>%
  dplyr::summarize( biomass=mean(biomass), cpua=mean(cpua) )
summary(biomass_total)



## bring biomass back together with rate data
biomass.rate <- left_join( rate.mean, biomass_total )
# get rid of places that didn't seine
biomass.rate$Country[ !(biomass.rate$Country %in% seine$Country)]
biomass.rate <- biomass.rate[ biomass.rate$Country %in% seine$Country,]
biomass.rate$biomass[ is.na(biomass.rate$biomass) ] <- 0
biomass.rate$cpua[ is.na(biomass.rate$cpua) ] <- 0
# remove NC2 unvegetated - no seine data (nothing caught)
biomass.rate <- biomass.rate %>%filter( Site != "NC2" | habitat != "Unvegetated" )
plot( rate ~ log(biomass+0.001), biomass.rate )
plot( rate ~ log(cpua+0.003), biomass.rate )

biomass.rate %>%  filter( rate==max(rate) )
ggplot( biomass.rate, aes( x=log(biomass+0.001), y=rate)) + geom_point() + 
  geom_smooth(method='glm',method.args=list(family='quasibinomial'))



# outliers
biomass.rate[ !is.na(biomass.rate$biomass) & biomass.rate$biomass.area==max(biomass.rate$biomass.area, na.rm=T), ]
biomass.rate[ !is.na(biomass.rate$biomass) & biomass.rate$biomass.area < exp(-5), ]
biomass.rate[ !is.na(biomass.rate$biomass) & biomass.rate$biomass.area > exp(0), ]
out <- sf[sf$Site=="Baja" & sf$habitat=="Seagrass",]
out <- sf[sf$Site=="Wales" & sf$habitat=="Unvegetated",]
out <- sf[sf$Site=="CA" ,]
out <- sf[sf$Site=="CA3" & sf$habitat=="Seagrass",]
out <- sf[sf$Site=="QLD4" & sf$habitat=="Unvegetated",]
out <- sr[sr$Site=="QLD4" & sr$habitat=="Unvegetated",]
out <- sf[sf$Site=="WA2",]
out %>%
  group_by(SPECIES_NAME) %>%
  summarise(abundance=sum(Abundance))



# save to disk
write_csv( biomass.rate, "Output Data/biomass_RDAselected.csv" )

