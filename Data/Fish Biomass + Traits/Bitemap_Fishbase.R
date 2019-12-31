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


# load libraries
library(tidyverse)


# read data
# existing traits
trait  <- read.csv( "../../R code and analysis/Output Data/fish_traits+taxonomy_20180122.csv", 
                    stringsAsFactors = FALSE )

# consider adding the species from individual sites and the range of traits from this analysis
partners <- read.csv( "Bitemap_REQUEST_trait_siene+video_RESPONSES_SORTED.csv",
                      stringsAsFactors = FALSE )


# summarize Fishbase traits by family
trait.sum <- trait %>%
  group_by( family ) %>%
  summarise( feedingtypes=list(unique(FeedingType)), trophicmean=mean(Trophic.Level,na.rm=T), 
             trophicsd = sd(Trophic.Level,na.rm=T), trophicn = length(unique(sciname)),
             trophicmin = min(Trophic.Level, na.rm=T), trophicmax = max(Trophic.Level, na.rm=T) )

# summarize traits from Partners by family
partner.sum <- partners %>%
  group_by( family ) %>%
  summarise( trophicgroups=list(unique(trophicGroupRLS)), watercolumn=list(unique(waterColumnRLS)),
             diel = list(unique(dielActivityRLS)), habitats=list(unique(habitatRLS)) )
             
# merge these two
traits <- full_join(trait.sum,partner.sum)
View(sort(unique( Hmisc::capitalize(tolower(unlist(traits$feedingtypes))) )))
View(sort(unique( Hmisc::capitalize(tolower(unlist(traits$trophicgroups))) )))
View(sort(unique( Hmisc::capitalize(tolower(unlist(traits$watercolumn))) )))
View(sort(unique( Hmisc::capitalize(tolower(unlist(traits$diel))) )))
View(sort(unique( Hmisc::capitalize(tolower(unlist(traits$habitats))) )))

# global mean trophic level
global <- trait %>% summarise( trophicmean = mean(Trophic.Level, na.rm=T), trophicsd=sd(Trophic.Level,na.rm=T),
                               trophicmin = min(Trophic.Level, na.rm=T), trophicmax = max(Trophic.Level, na.rm=T) )
global$family="ALL"
traits <- full_join( traits, global )

# color for ALL vs indidvual families
traits$color <- "2"
traits$color[ traits$family == "ALL" ] <- "1"



# select ones that were "strong" in NMDS
trait.choose <- traits %>% filter( family %in% c("Sparidae","Tetraodontidae","Labridae","Mullidae","Tetrarogidae",
                                 "Haemulidae","Terapontidae","Tetraodontidae","Sillaginidae",
                                 "Hemiramphidae", "ALL") )

# trait.choose %>% 
#   group_by(family)%>%
#   mutate( View(sort(unique( Hmisc::capitalize(tolower(unlist(feedingtypes))) ))) )

write.csv( trait.choose, "trait_summary_family_strong.csv", row.names=F )

windows(3,3 )
ggplot( trait.choose, aes(x=family,y=trophicmean,ymax=trophicmax,ymin=trophicmin, color=color) ) + 
  geom_point(size=2) + geom_errorbar(width=0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_text( aes(y=trophicmean,label=trophicn), nudge_x = -0.33 ) +
  scale_color_manual(values=c("firebrick","slateblue")) +
  ylab( "Mean trophic level ± range") + xlab("Family")



# try to plot them all
# need multiple panels
  
  
  
# look for omnivory and herbivory among all families
traits$troph2 <- apply( traits, 1, function(z) {
  tmp <- z$trophicgroups
  comb <- paste( sort(unique( tmp ) ), collapse="; " )
  return(comb)
} )

# does trophic group include omnivory?
traits$omnivory <- 0
traits$omnivory[ grep( "mnivor*", traits$troph2 ) ] <- 1
traits$omnivory[ is.na(traits$trophicgroups) ] <- NA


# does trophic group include herbivory?
traits$herbivory <- 0
traits$herbivory[ grep( "herbiv*", traits$troph2 ) ] <- 1
traits$herbivory[ is.na(traits$trophicgroups) ] <- NA

# select columns
traitsave <- traits %>%
  filter( !is.na(family) ) %>%
  select( family, trophicmean, trophicmin, trophicmax, omnivory, herbivory )

# write to disk
write.csv( traitsave, "../../R code and analysis/Output Data/trophic_family_strong_omni.csv", row.names = FALSE )
