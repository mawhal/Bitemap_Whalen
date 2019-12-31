#############################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Use traits of consumers in analyses of community stucture, consumption, etc.
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 30 August 2019
#############################################################################


## load libraries
library(tidyverse) 


## read trait data - traits imputed using package mice in other script
trait <- read_csv( "../Data/Fish Biomass + Traits/Bitemap_REQUEST_trait_siene+video_RESPONSES_SORTED_redux.csv" )
trait <- trait %>%
  # filter(eat.squid==1) %>% 
  select( Country, family, sciName, active.ambush, eat.squid ) %>%
  distinct()
trait$sciName <- gsub( " ", ".", trait$sciName  )
trait <- trait %>% 
  separate( sciName, into=c("Genus", "Species"), remove = FALSE )
trait.genus <- trait %>% 
  select( Genus, active.ambush, eat.squid ) %>% 
  distinct() %>% 
  arrange( Genus )


# read table with families and phyla 
phyla <-read_csv( "Output Data/predator_families+phyla.csv" )

# seine abundance + biomass
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )
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
# avearage by site and habitat, and convert to presence-absence
mean.catch <- seine.all %>% 
  dplyr::group_by(Country,habitat,SPECIES_NAME) %>% 
  dplyr::summarize( cpua=mean(cpua) ) %>% 
  dplyr::mutate( P=ifelse(cpua>0,1,0) ) %>% 
  filter( P==1 ) %>% 
  separate( SPECIES_NAME, into=c("Genus","Species"), sep=" ", remove=FALSE)

# trait information that includes video data
video <- read_csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv" )
video <- video %>% 
  select( Country, habitat, Genus, Species ) %>% 
  mutate( habitat=tolower(habitat) ) %>% 
  distinct()
video$habitat[ video$habitat == "seagrass" ] <- "Seagrass"
video$habitat[ video$habitat == "unveg" ] <- "Unvegetated"
video <- video %>% filter( !is.na(Genus) )
video$P <- 1

# merge consumer observations
consumers <- bind_rows(select(mean.catch, Country, habitat, Genus, Species, P),video)


# merge traits and consumers
trait.hab.pres <- left_join( consumers, trait.genus )



# read data on consumption rate, sites, environment
pop <- read_csv( "Output Data/Bitemap_rate.env.20190423.csv" )
pop <- pop %>% 
  dplyr::group_by( Country, Site,habitat,meanLat ) %>% 
  dplyr::summarise( rate=mean(rate) )





### summarize trait data by site
# filter out only fish?
suse <- trait.hab.pres
schord <- filter( s1, phylum=="Chordata" )
suse <- left_join( distinct(select(s1,Country,Site)), schord )
# ratio of active to ambush at each site by individual 
rat.tax <- suse %>%
  dplyr::group_by( Country, habitat ) %>%
  dplyr::summarize(act.ratio = length(active.ambush[active.ambush %in% c("active","ambush.active")]) / length(Genus)  ) 

# merge
rat.site <- left_join( rat.tax, pop )

# plot
ggplot( rat.site, aes(x=abs(meanLat), y=act.ratio)) + geom_point() + geom_smooth()
ggplot( rat.site, aes(x=act.ratio,  y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))
# there is some relationship here, but not very strong -- just one trait that is not very well identified across all taxa
hilo <- rat.site %>% 
  ungroup() %>% 
  filter( act.ratio==1, rate<0.5 ) %>% 
  arrange(rate) %>% 
  select( Site ) %>% 
  unlist()

# nice figure
windows(2.6,2.6)
ggplot( rat.site, aes(x=act.ratio,  y=rate)) + 
  geom_point(pch=16,alpha=0.5) +
  xlab("Proportion of active fish") + ylab("Consumption rate") +
  theme_classic()
summary(glm(rate~act.ratio,data=rat.site,family="quasibinomial"))


# write to disk
write_csv( rat.site, "Output Data/consumer_active_ratio.csv" )
