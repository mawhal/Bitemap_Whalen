#############################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Use traits of consumers in analyses of community stucture, consumption, etc.
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 09 January 2020
#############################################################################


## load libraries
library(tidyverse) 

## TRAITS
## read trait data - traits imputed using package mice in other script
trait <- read_csv( "Data/Fish Biomass + Traits/Bitemap_trait_impute.csv" )
# clean up this dataset a bit
trait$family[ trait$family=="Channidae"] <- "Sciaenidae"
trait$family[ trait$family=="Sparisoma viride"] <- "Scaridae"
# filter taxa
trait <- trait %>% filter( sciName != "Aphanius fasciatus REDUNDANT" )
#  convert benthopelagic to demersal
trait$waterColumnRLS <- tolower(trait$waterColumnRLS)
trait$waterColumnRLS <- gsub( ".*benthopelagic.*","demersal",trait$waterColumnRLS)
trait$waterColumnRLS[ trait$waterColumnRLS == "pelagic"] <- 'pelagic non-site attached'
trait$waterColumnRLS[ trait$family == "Cottidae"] <- 'benthic'
# name change
trait$sciName[ trait$sciName=="Romaleon polyodon" ] <- "Romaleon setosum"
trait[ trait$feedingTypeFishbase %in% "filtering plankton", ]
trait[ trait$feedingTypeFishbase %in% "filtering plankton", ]

# summary functions to get a list of unique traits
num <- function(x) length(sort(unique(tolower(x))))
numlist <- function(x) paste(sort(unique(tolower(x))),collapse=";")

trait.unique.num <- trait %>% 
  group_by( sciName, family ) %>% 
  summarize( act=num(active.ambush), feed=num(feedingTypeFishbase),
             troph=num(trophicGroupRLS), watercol=num(waterColumnRLS),
             clim=num(Climate), food=num(FoodTroph), asp=num(AspectRatio),
             tl = num(TL), body=num(BodyShapeI), snoutx=num(SnoutTipX), snouty=num(SnoutTipY) )
trait.unique <- trait %>% 
  group_by( sciName, family ) %>% 
  summarize( act=numlist(active.ambush), feed=numlist(feedingTypeFishbase),
             troph=numlist(trophicGroupRLS), watercol=numlist(waterColumnRLS),
             clim=numlist(Climate), food=numlist(FoodTroph), asp=numlist(AspectRatio),
             tl = numlist(TL), body=numlist(BodyShapeI), snoutx=numlist(SnoutTipX), snouty=numlist(SnoutTipY) )
multis.row <- apply( trait.unique.num[,-c(1,2)], 1, function(x) any(x>1) )
multis.col <- apply( trait.unique.num[,-c(1,2)], 2, function(x) any(x>1) )
# manually adjust some categories
trait.unique$feed[ trait.unique$feed == "browsing on substrate;variable" ] <- "browsing on substrate"
trait.unique$feed[ trait.unique$feed == "hunting macrofauna (predator);selective plankton feeding" ] <- "selective plankton feeding"

trait.unique[ multis.row, c(1, which(multis.col)+2) ]

trait.unique.family <- trait %>% 
  group_by( family ) %>% 
  summarize( act=numlist(active.ambush), feed=numlist(feedingTypeFishbase),
             troph=numlist(trophicGroupRLS), watercol=numlist(waterColumnRLS),
             clim=numlist(Climate), food=numlist(FoodTroph), asp=numlist(AspectRatio),
             tl = numlist(TL), body=numlist(BodyShapeI), snoutx=num(SnoutTipX), snouty=num(SnoutTipY) )
tp <- trait.unique.family %>% 
  select( family, act, feed, troph, watercol, clim, body )
write_csv( tp, "Data/processed/traits_family_raw_summary.csv" )

# edited summary
tp_edit <- read_csv( "Data/processed/traits_family_edit_summary.csv" )



# fill in missing cells with family-level unique traits
tp_merge <- left_join( trait.unique, tp_edit, by=c("family","act") )
tp_fill <- tp_merge %>% 
  mutate( feed = replace(feed.x, feed.x=="", feed.y),
          troph= replace(troph.x, troph.x=="", troph.y),
          watercol= replace(watercol.x, watercol.x=="", watercol.y),
          clim= replace(clim.x, clim.x=="", clim.y),
          body= replace(body.x, body.x=="", body.y) )
trait.final <- tp_fill %>% 
  select(sciName,family,act,food,asp,tl,feed,troph,watercol,clim,body,snoutx,snouty) %>% 
  mutate( food=replace(food,food=="",NA),
          asp=replace(asp,asp=="",NA),
          tl=replace(tl,tl=="",NA),
          snoutx=replace(snoutx,snoutx=="",NA), 
          snouty=replace(snouty,snouty=="",NA) )


trait.final <- ungroup(trait.final)


# simplify scinames
trait.final <- trait.final %>% 
  mutate( name=vegan::make.cepnames( trait.final$sciName ) )
# clean up duplicate entries
options( stringsAsFactors = FALSE )
clean <- function(column){
  unlist( lapply(strsplit( column, split=';' ),function(z) unlist(paste(sort(unique(z)),collapse=";")) ) )
}
trait.final <- data.frame(apply(trait.final,2,clean))


# write to disk and potentially update manuall
write_csv( trait.final, "Data/processed/traits_clean_intial.csv")



### look at range of traits - can we break these down at all
# for each categorical trait of interest, which ones have multiple entries (marked by semicolon)
# traits of interest
columns <- c("feed","troph","watercol","body")

# Multiple trait listings
trait.final %>% as_tibble() %>%  select( sciName, family, columns ) %>% 
  filter( grepl(";",feed) |  grepl(";",troph)  | grepl(";",watercol) | grepl(";",body) ) %>% 
  arrange( family,sciName ) %>% 
  View()


### Water Columns
# Multiple trait listings
trait.final %>% as_tibble() %>%  select( sciName, family, columns ) %>% 
  filter( grepl(";",watercol) ) %>% 
  arrange( watercol,family,sciName ) %>% 
  View()
# get rid of benthopelagic
trait.final$watercol <- gsub( ".*benthopelagic.*","demersal",trait.final$watercol)
trait.final$watercol[ trait.final$family == "Palaemonidae" ] <- "benthic"
# benthic;demersal <- demersal
trait.final$watercol <- gsub( "benthic;demersal","demersal",trait.final$watercol)
# triggerfish
trait.final$watercol[ trait.final$sciName == "Balistes vetula"] <- "demersal"
trait.final$troph[ trait.final$sciName == "Balistes vetula"] <- "benthic invertivore"
trait.final$watercol[ trait.final$sciName == "Canthidermis sufflamen"] <- "pelagic site attached"
trait.final$troph[ trait.final$sciName == "Canthidermis sufflamen"] <- "benthic invertivore;planktivore"
# demersal
trait.final$watercol[trait.final$sciName %in% c( "Elops saurus","Sphoeroides testudineus","Arothron sp.") ] <- "demersal"
# Embiotocidae
trait.final$watercol[trait.final$family=="Embiotocidae"] <- "demersal;pelagic site attached"
# monodactylidae
trait.final$watercol[trait.final$family=="Monodactylidae"] <- "demersal"


### Feeding
# any mention of omnivore make it omnivore
trait.final$troph <- gsub( ".*omnivore.*","omnivore", trait.final$troph )

# anything featuring variable make variable
trait.final$feed <- gsub( ".*variable.*", "variable", trait.final$feed )

# find mentions of herbivory along with others
t.herb <- trait.final %>% as_tibble() %>%  select( sciName, family, columns ) %>% 
  filter( grepl( "herbivore", trait.final$troph ) & grepl(";",troph)  ) 




# multiple entries under feed <- variable?
# remove omnivory + cleaner
remove_omni_cleaner <- c("Cheilio inermis","Coris julis","Coris caudimacula","Coris pictoides","Coris sp.",
                 "Halichoeres bivittatus", "Halichoeres chloropterus","Halichoeres garnoti",
                 "Halichoeres trimaculatus","Thalassoma sp.")
trait.final$troph[ trait.final$sciName %in% remove_omni_cleaner ] <- gsub(";scraping herbivore","",trait.final$troph[ trait.final$sciName %in% remove_omni_cleaner ] )
trait.final$troph[ trait.final$sciName %in% remove_omni_cleaner ] <- gsub(";cleaner","",trait.final$troph[ trait.final$sciName %in% remove_omni_cleaner ] )

# selective plankton feeder
select_plankt <- "Leptojulis sp."
trait.final$feed[ trait.final$sciName %in% select_plankt ] <- "selective plankton feeding"
trait.final$troph[ trait.final$sciName %in% select_plankt ] <- "planktivore"


# benthic invertivores
benthic_invert <- c("Halichoeres chloropterus","Halichoeres garnoti","Halichoeres trimaculatus",
                    "Halichoeres sp.","Labrus bergylta","Symphodus cinereus","Symphodus melops",
                    "Himanthalia elongata")

trait.final$troph[ trait.final$sciName %in% benthic_invert ] <- "benthic invertivore"

# browsing hervivore
brows_herb <- c("Sparisoma viride") # browsing on substrate, herbivore
trait.final$feed[ trait.final$sciName %in% brows_herb ] <- "browsing on substrate"
trait.final$troph[ trait.final$sciName %in% brows_herb ] <- "scraping herbivore"


# Carangidae
# benthic invertivore;higher carnivore
trait.final$troph[ grep( "Caran*", trait.final$sciName )] <- "benthic invertivore;higher carnivore"
# higher carnivore 
trait.final$troph[ trait.final$sciName=="Scomberoides tala"] <- "higher carnivore"
# benthic invertivore
trait.final$troph[ grep( "Trachinotus*", trait.final$sciName )] <- "benthic invertivore"



# Clupeidae
trait.final$watercol[ trait.final$family=="Clupeidae" & 
                        trait.final$watercol=="pelagic;pelagic non-site attached;pelagic site attached"] <- "pelagic non-site attached"
# planktivore - variable (can switch filtering and seletive)
# Clupea harengus
trait.final$feed[ trait.final$sciName=="Clupea harengus"] <- "variable"
trait.final$troph[ trait.final$sciName=="Clupea harengus"] <- "planktivore"
# selective plankton
# Harengula jaguana
# Sprattus sprattus
trait.final$feed[ trait.final$sciName %in% c("Sprattus sprattus","Harengula jaguana")] <- "selective plankton feeding"
trait.final$troph[ trait.final$sciName %in% c("Sprattus sprattus","Harengula jaguana")] <- "planktivore"


# Gadidae - check length of gadids and guess feeding (plankton vs inverts vs higher carnivore)
trait.final$troph[trait.final$family=="Gadidae"] <- "benthic invertivore;higher carnivore"



### Body shape
trait.final %>% as_tibble() %>%  select( sciName, family, columns ) %>% 
  filter( grepl(";",body) ) %>% 
  arrange( body,family,sciName ) %>% 
  View()
# Cottidae
# trait.final$body[trait.final$family=="Cottidae" & ] <- 




### write cleaned version to disk
write_csv( trait.final, "Data/processed/Bitemap_traits_clean_final.csv")

