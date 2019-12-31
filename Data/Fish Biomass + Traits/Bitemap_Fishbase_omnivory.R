#####################################
### MarineGEO Bitemap
### Fish traits, revisted
### Matt Whalen
### created 06 March 2019
####################################

# goal -- take compiled trait data (especially trophic level)
#  and compare to full list of taxa from seines
#  fill in gaps where necessary
#  then compile a new family-level summary of trophic ecology 
#      based on those species sampled in Bitemap



## additional goal of this script (see Bitemap_Fishbase.R)
# assess herbivory at the site level, and integrate this into MDS analysis


# load libraries
library(tidyverse)


# read data
# existing traits
trait  <- read.csv( "../../R code and analysis/Output Data/fish_traits+taxonomy_20180122.csv", 
                    stringsAsFactors = FALSE )

# consider adding the species from individual sites and the range of traits from this analysis
partners <- read.csv( "Bitemap_REQUEST_trait_siene+video_RESPONSES_SORTED.csv",
                      stringsAsFactors = FALSE )


# # summarize Fishbase traits by family
# trait.sum <- trait %>%
#   group_by( family ) %>%
#   summarise( feedingtypes=list(unique(FeedingType)), trophicmean=mean(Trophic.Level,na.rm=T),
#              trophicsd = sd(Trophic.Level,na.rm=T), trophicn = length(unique(sciname)),
#              trophicmin = min(Trophic.Level, na.rm=T), trophicmax = max(Trophic.Level, na.rm=T) )

# summarize traits from Partners by Species
partner.sum <- partners %>%
  group_by( sciName ) %>%
  summarise( trophicgroups=paste(unique(tolower(trophicGroupRLS[!is.na(trophicGroupRLS)])),collapse="; "),
             watercolumn=paste(unique(tolower(waterColumnRLS)),collapse="; "),
             diel=paste(unique(tolower(dielActivityRLS)),collapse="; "),
             habitats=paste(unique(tolower(habitatRLS)),collapse="; ")
             )

# merge these two
traits <- full_join(trait,partner.sum, by=c("sciname"="sciName") )


# for each family, get instances of herbivory and omnivory (how is this defined??)
# does trophic group include omnivory?
traits$omnivory <- 0
traits$omnivory[ grep( "mnivor*", traits$trophicgroups ) ] <- 1
traits$omnivory[ traits$trophicgroups=="" ] <- NA
traits$omnivory[ grep( "variable*", traits$FeedingType ) ] <- 1

# does trophic group include herbivory?
traits$herbivory <- 0
traits$herbivory[ grep( "herbiv*", traits$trophicgroups ) ] <- 1
traits$herbivory[ traits$trophicgroups=="" ] <- NA
traits$herbivory[ grep( "plants", traits$FeedingType ) ] <- 1
traits$herbivory[ grep( "variable*", traits$FeedingType ) ] <- 1





# select relevant columns and get rid of duplicates
seine.spp <- partners %>%
  select( Country, sciname=sciName ) %>%
  distinct()

# merge back with seine data
seine <- left_join( seine.spp, traits )


# calculate site-level omnivory and herbivory totals
site.omni <- seine %>%
  group_by( Country ) %>%
  summarise( omnisum=sum(omnivory,na.rm=T),herbsum=sum(herbivory,na.rm=T),
             omnirich=length(omnivory[!is.na(omnivory)]), herbrich=length(herbivory[!is.na(herbivory)]) ) %>%
  mutate( omniprop=omnisum/omnirich, herbprop=herbsum/herbrich )




# write to disk
write.csv( site.omni, "../../R code and analysis/Output Data/omnivory_site.csv", row.names=FALSE )
