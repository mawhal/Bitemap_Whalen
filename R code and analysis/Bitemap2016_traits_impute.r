#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Impute traits of consumers based on values from relatives
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 30 August 2019
#####################################################################


## load libraries
library(tidyverse) 
library(mice)
library(rfishbase)

## read trait data
draw <- read_csv( "Data/Fish Biomass + Traits/Bitemap_trait_siene+video_RESPONSES_SORTED_redux.csv" )
traits <- read_csv( "Data/processed/consumer_eat_squid.csv" )
# merge eat squidpop trait
dmerge <- left_join( select(draw,-eat.squid), select(traits, Country, sciName, eat.squid) )

## mice package vignette
dmerge %>%
  select( feedingTypeFishbase, trophicGroupRLS, waterColumnRLS, eat.squid ) %>%
  md.pattern()

d <- dmerge

## for species that appear multiple times in a region, fill in the trait data
ddup <- d %>%
  group_by( Country, sciName ) %>%
  summarize( n=length(totalCount) ) %>%
  filter( n>1 )
ddup <- left_join( ddup, d )
ddup <- ddup %>%
  filter( is.na(feedingTypeFishbase) | is.na(trophicGroupRLS) | is.na(waterColumnRLS) ) %>%
  arrange( family, sciName, Country )
# note: Whalen used this to go back to the data loaded in this script and change it.





## add more traits from fishbase
fish <-  c( gsub("[.]"," ",tail(trait.keep$sciName)), "Lagodon rhomboides" )
fish <- sort(unique(d$sciName))
inverts <- sort( unique(d$sciName[d$family %in% 
                                    c("Portunidae","Penaeidae","Palaemonidae", 'Cancridae')]) )
eco <- ecology( fish, fields = c('Species','Herbivory2', 'DietTroph', 'DietSeTroph', 'FoodTroph', 
                                 'FoodSeTroph', 'SoftBottom', 'Sand', 'Mud', 'HardBottom', 
                                 'SeaGrassBeds', 'CoralReefs', 'Estuaries', 'Mangroves', 
                                 'Intertidal', 'Saltmarshes') )  # Herbivory2, DietTroph, DietSeTroph, FoodTroph, FoodSeTroph, SoftBottom, Sand, Mud, HardBottom, SeaGrassBeds, CoralReefs, Estuaries, Mangroves, Intertidal, Saltmarshes
eco2 <- ecosystem( fish )  # Climate
eco2 %>% select(Species,Climate)
food <- fooditems(  fish, fields=c('Species','FoodI','FoodII'))  # FoodI, FoodII
swim <- swimming(fish, fields = c('Species','AdultType','AdultMode')) # AdultType, AdultMode
morph <- morphology(fish,fields = c('Species','BodyShapeI', 'BodyShapeII', 'Forehead', 'TypeofMouth', 'PosofMouth', 'CShape', 'Attributes') ) #BodyShapeI, BodyShapeII, Forehead, TypeofMouth, PosofMouth, Cshape, Atributes, 
morph.met <- morphometrics(fish, fields = c( 'Species','TL', 'HL', 'BD', 'SnoutTipX', 'SnoutTipY', 'AspectRatio') ) #Total Length, Head Length, Body Depth, SnoutTipX, SnoutTipY, AspectRatio
brain <- brains( fish, fields = c('Species','EncIndex') ) # EncIndex
# repeat for inverts
ieco <- ecology( inverts, fields = c('Species','Herbivory2', 'DietTroph', 'DietSeTroph', 'FoodTroph', 
                                     'FoodSeTroph', 'SoftBottom', 'Sand', 'Mud', 'HardBottom', 
                                     'SeaGrassBeds', 'CoralReefs', 'Estuaries', 'Mangroves', 
                                     'Intertidal'),
                 server='sealifebase' )  # Herbivory2, DietTroph, DietSeTroph, FoodTroph, FoodSeTroph, SoftBottom, Sand, Mud, HardBottom, SeaGrassBeds, CoralReefs, Estuaries, Mangroves, Intertidal, Saltmarshes
ieco2 <- ecosystem( inverts,
                    server='sealifebase' )  # Climate
ieco2 %>% select(Species,Climate)
ifood <- fooditems(  inverts, fields=c('Species','FoodI','FoodII'),
                     server='sealifebase')  # FoodI, FoodII

# summarize the ones that differ from number of species
Mode = function(x){ 
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = names(ta)[ta == tam]
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(paste(mod,collapse = ";"))
}

# fish
eco[which(table(eco$Species)>1),]
ecosum <- eco %>%
  select(-Herbivory2) %>%
  group_by( Species ) %>%
  summarize_all( mean, na.rm=T )
herbsum <- eco %>%
  select( Species, Herbivory2 ) %>%
  group_by( Species) %>%
  summarize( Herbivory2=Mode(Herbivory2))
ecosum <- full_join(ecosum,herbsum)
ecosum$Intertidal[ !is.na(ecosum$Intertidal) & ecosum$Intertidal == -0.5 ] <- -1
ecosum$Estuaries[ !is.na(ecosum$Estuaries) & ecosum$Estuaries == -0.5 ] <- -1

eco2sum <- eco2 %>%
  mutate( Climate = tolower(Climate) ) %>%
  group_by( Species ) %>%
  summarize( Climate=Mode(Climate) )

foodsum <- food %>%
  group_by( Species ) %>%
  summarise( FoodI=paste(sort(unique(FoodI)),collapse=';'),
             FoodII=paste(sort(unique(FoodII)),collapse=';') )

rem <- filter(morph,Species=="Clupea harengus")[2,]
morphsum <- anti_join(morph,rem)

morph.met[which(table(morph.met$Species)>1),]
morphmetsum <- morph.met %>%
  group_by(Species) %>%
  summarise_all( mean, na.rm=T )
  
brainsum <- brain %>%
  group_by( Species ) %>%
  summarize( EncIndex=mean(EncIndex,na.rm=T) )

# inverts

ieco[which(table(ieco$Species)>1),]
iecosum <- ieco %>%
  select(-Herbivory2) %>%
  group_by( Species ) %>%
  summarize_all( mean, na.rm=T )
iherbsum <- ieco %>%
  select( Species, Herbivory2 ) %>%
  group_by( Species) %>%
  summarize( Herbivory2=Mode(Herbivory2))
iecosum <- full_join(iecosum,iherbsum)
iecosum$Intertidal[ !is.na(iecosum$Intertidal) & iecosum$Intertidal == -0.5 ] <- -1
iecosum$Estuaries[ !is.na(iecosum$Estuaries) & iecosum$Estuaries == -0.5 ] <- -1

ieco2sum <- ieco2 %>%
  mutate( Climate = tolower(Climate) ) %>%
  group_by( Species ) %>%
  summarize( Climate=Mode(Climate) )

ifoodsum <- ifood %>%
  group_by( Species ) %>%
  summarise( FoodI=paste(sort(unique(FoodI)),collapse=';'),
             FoodII=paste(sort(unique(FoodII)),collapse=';') )






# merge new traits
trait.new <- full_join( full_join( full_join( full_join( full_join( morphsum, morphmetsum ), brainsum), ecosum), eco2sum), foodsum)
itrait.new <- full_join( full_join( iecosum, ieco2sum), ifoodsum)
# join fish and invert traits
trait.full <- full_join( trait.new, itrait.new )
# merge with data
dplus <- left_join(d,trait.full, by=c('sciName'='Species'))

# # convert to numeric
# dnum <- do.call( cbind, lapply( dtrait, function(z) as.numeric(z) ) )
# 

## IMPUTATION
# using mean of values
# dimp <- mice( dnum, method = "mean", m = 1, maxit = 1)
# complete( dimp )

# using the default method for factors with >2 unordered levels
# impute feeding and trophic things together, water column use separately
dtoimp <- dplus %>%
  select(sciName, family, active.ambush, feedingTypeFishbase, trophicGroupRLS, waterColumnRLS, eat.squid ) 
dtoimp[,4:6] <- apply( dtoimp[,4:6], 2, tolower )
dtoimp[,4:6]  # Need to make columns to be imputed factors!
dtoimpfct <- dtoimp %>% 
  mutate_all( as.factor )

md.pattern( dtoimpfct )
  
dimp2 <- mice( dtoimpfct, blocks=list(c("feedingTypeFishbase","trophicGroupRLS",
                                        "waterColumnRLS", "eat.squid")) )

complete( dimp2 )
dimp2$predictorMatrix
# if family column is included, it will be used to predict the imputed value?


## clean data and write to disk
dnew <- as_tibble(data.frame(complete( dimp2 )))
names(dnew) <- c("scieName", "family", "active.ambush",
                 "feedingImpute", "trophicImpute", "waterColumnImpute","eat.squidImpute")

dcomb <- data.frame( d[,1:8], dnew[,4], d[,9], dnew[,5], d[,10], dnew[,6], d[,14], dnew[,7], 
                     dplus[15:46] )

write_csv( dcomb, "../Data/Fish Biomass + Traits/Bitemap_trait_impute.csv" )
