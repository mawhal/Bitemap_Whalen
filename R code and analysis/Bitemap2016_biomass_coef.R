#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Prepare Biomass conversion table
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 20 September 2019
#####################################################################

# load libaries
library( tidyverse )

# read biomass coefficients
coef <- read.csv("../Data/Fish Biomass + Traits/20160710_RLS_biomass_coefs_20170921.csv")
# omit rows from coef for which we have no estimates
coef <- coef[ !is.na(coef$A), ]
# remove duplicates from coef
coef <- coef[ !duplicated(coef), ]

# FROM PREVIOUS SCRIPT
# 
# # identify species not in character list
# lookup <- sort(unique(fish.clean$SPECIES_NAME[!(fish.clean$SPECIES_NAME %in% coef$SPECIES_NAME)]))
# sort(unique(fish.clean$SPECIES_NAME[(fish.clean$SPECIES_NAME %in% coef$SPECIES_NAME)]))

# # use rfishbase to look up taxa not represented in the reef life survey biomass conversion data.frame (coef)
# fishbaseLW <- length_weight(lookup)
# sort(unique(fishbaseLW$sciname)) # added 67 taxa
# unrepresented <- sort(unique(fish.clean$SPECIES_NAME[!(fish.clean$SPECIES_NAME %in% c(fishbaseLW$sciname, as.character(coef$SPECIES_NAME)) )]))
# write.csv( fishbaseLW, "../Data/Fish Biomass/fishbaseLW.csv", row.names = FALSE)
# write.csv( data.frame(unrepresented), "../Data/Fish Biomass/fishbaseMISSING.csv", row.names = FALSE)

# read the length-weight relationships looked up in fishbase
fishbaseLW <- read.csv( "../Data/Fish Biomass + Traits/fishbaseLW_20170921.csv" )
# only keep rows that were chosen from the database (another option would be to average all estimates by taxon)
fishbaseLW <- fishbaseLW[ fishbaseLW$Keep == "Y", ]
# average by taxon
fishbaseLW <- fishbaseLW %>%
  group_by( sciname) %>%
  summarize( a=mean(a),aTL=mean(aTL),b=mean(b) )

# accept the estimates for a (intercept) that account for standard vs total length
for(i in 1:nrow(fishbaseLW) ){
  if( !is.na(fishbaseLW$aTL[i]) ) fishbaseLW$a[i] <- fishbaseLW$aTL[i]
}
# rename sciname to Species
names(fishbaseLW) <- c( "SPECIES_NAME", "A", "aTL", "B" )

# read in table for taxa that were looked up manually in fishbase
fishbaseManual <- read.csv( "../Data/Fish Biomass + Traits/Bitemap_FishbaseLW_manual.csv" )[,1:3]
names(fishbaseManual) <- c( "SPECIES_NAME", "A", "B" )

# combine all data from fishbase, including those provided by Reef Life Survey
fishbaseFull <- full_join( fishbaseLW, fishbaseManual )
biomassCoef <- full_join(  coef, fishbaseFull )
biomassCoef <- biomassCoef %>%
  dplyr::select( SPECIES_NAME, A, B )
length(unique(biomassCoef$SPECIES_NAME))
# omit duplicated rows (separate estimates for a taxon), default to use Reef Life Survey estimates for consistency
biomassCoef <- biomassCoef[ !duplicated(biomassCoef$SPECIES_NAME), ]


## Write to disk
write_csv( biomassCoef, "Output Data/Bitemap_biomass_coef.csv" )
