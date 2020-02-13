#############################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Use traits of consumers in analyses of community stucture, consumption, etc.
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 30 August 2019
#############################################################################


## load libraries
library(tidyverse) 
library(FD)


## TRAITS
## read trait data - see traits_imput.R and traits_clean.r for details
trait.final <- read_csv( "Output Data/traits_clean_final.csv" )


### switch which traits to use
trait.use <- trait.final %>% 
  select( sciName, name, act, feed, troph, watercol, body  )



## mice package 
md <- trait.use %>%
  mice::md.pattern()




## read seine data
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )
seine <- seine %>% 
  select( Country, habitat, Lat, Long, Date, Time, 
          sciName=SPECIES_NAME, family, Genus,
          length=Length, biomass, abun=Abundance )
# change some taxon names
seine$sciName[ seine$sciName=="Pelates sexlineatus"] <- "Helotes sexlineatus"
seine$family[ seine$sciName %in% "Carcinus maenas" ] <- "Carcinidae"
unique(seine$sciName[ !(seine$sciName %in% trait.final$sciName) ]) # still a few missing taxa
seine$family[ seine$sciName %in% "Menticirrhus ophicephalus" ] <- "Sciaenidae" 
seine$family[ seine$sciName %in% "Mallotus villosus" ] <- "Osmeridae" 
# summarize median per capita length and biomass per seine
seine_per <- seine %>% 
  group_by( Country, habitat, Lat, Long, Date, Time, family, sciName ) %>% 
  summarize( length=median(length,na.rm=T), biomass=median(biomass,na.rm=T), abun=sum(abun,na.rm=T) )
# summarize mean across multiple seines within habitats
seine_hab <- seine_per %>% 
  group_by( Country, habitat, family, sciName ) %>% 
  summarize( length=mean(length,na.rm=T), biomass=mean(biomass,na.rm=T), abun=mean(abun,na.rm=T) )
# summarize site info for merging below
seine_merge <- seine_per %>% 
  group_by( Country, habitat, family, sciName ) %>% 
  summarize( present=1 )

# trait information that includes video data
video <- read_csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv" )
# some sites list the full species name under Species instead of the epithet
video$Species[ grep( "[A-Z]", video$Species ) ] <- 
  unlist( lapply( strsplit( video$Species[ grep( "[A-Z]", video$Species ) ], " " ), 
                  function(z) paste(z[-1],collapse=" ") ) )
# fix some names
video$Genus <- gsub( "Athrina", "Atherina", video$Genus )
video$Genus <- gsub( "Eucinostramus", "Eucinostomus", video$Genus )
video$Genus <- gsub( "Lactophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Lactrophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Alters", "Aluterus", video$Genus )
video$Genus <- gsub( "Sphyraema", "Sphyraena", video$Genus )

# 
video <- video %>% 
  unite( sciName, Genus, Species, sep=" ", remove = F ) %>% 
  select( Country, habitat, sciName, Genus ) %>% 
  mutate( habitat=tolower(habitat) ) %>% 
  distinct()
video$habitat[ video$habitat == "seagrass" ] <- "Seagrass"
video$habitat[ video$habitat == "unveg" ] <- "Unvegetated"
# fix more names
video$sciName <- gsub( " sp", " sp.", video$sciName )
video$sciName <- gsub( " sp..", " sp.", video$sciName )
video$sciName <- gsub( " NA", "", video$sciName )

genfam <- seine %>% select(family,Genus) %>% distinct()
video1 <- left_join(video, genfam)


# lookup taxa with taxize
taxa <- video1 %>% 
  filter( is.na(family) ) %>% 
  separate( sciName, c("Genus","species")) %>% 
  select( Genus ) %>% 
  distinct() %>%
  arrange( Genus ) %>% 
  unlist()
taxa <- na.omit(taxa)
taxa <- taxa[ !(taxa %in% c("NO","Not", "NA")) ]
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
video_fam <- left_join(video, famfam, by="Genus")
video_fam <- video_fam %>% filter( sciName != "NA" ) %>% 
  filter( sciName != "NO ID" )

# name change Romaleon polyodon to Romaleon  
video_fam$sciName[ video_fam$sciName=="Romaleon polyodon" ] <- "Romaleon setosum"
video_fam$sciName[ video_fam$sciName=="Sparisoma viridae" ] <- "Sparisoma viride"
video_fam$sciName[ video_fam$sciName=="Diodon histrix" ] <- "Diodon hystrix"
video_fam$sciName[ video_fam$sciName=="Zosterisessor offiocephalus" ] <- "Zosterisessor ophiocephalus"
# other names to fix
# z ophiocephalus appears twice
video_fam <- video_fam %>% distinct()
# add fmaily info for Aluterus scripus
video_fam$family[ video_fam$sciName %in% "Aluterus scriptus" ] <- "Monacanthidae"
# cleanup
# taxa with traits needed
video_fam %>% filter(is.na(family))


## read site data
sites <- read_csv( "Output Data/Bitemap_BioORACLE_20190107.csv" )

## select columns
sites <- sites %>%
  select( Country, Site, basin, coast, meanLat, meanLong, chlomax, chlomean, chlomin, chlorange, cloudmax, cloudmean, cloudmin,
          parmax, parmean, salinity, sstmax, sstmean, sstmin, sstrange )

## consumption rate data
rate.env <- read_csv( "Output Data/Bitemap_rate.env.20190423.csv" ) %>% 
  select( Country, Site, habitat, rate )
rate.hab <- rate.env %>% 
  group_by( Country, Site, habitat ) %>% 
  summarize( rate=mean(rate,na.rm=T) )

## merge
ss <- left_join( rate.hab, sites )
all_taxa <- full_join( video_fam, seine_merge )
all_taxa %>% filter(is.na(family))
s1 <- left_join( all_taxa, ss )
s1 <- left_join( s1, trait.use )

# manually fill in missing taxa
# library(rfishbase)
# eco <- ecology(unlist(missingtaxa$sciName), fields = "FeedingType")
# morph <- morphology( unlist(missingtaxa$sciName),fields = c('Species','BodyShapeI') )
s1[ s1$sciName == "Aluterus scriptus", c("act","feed","troph","watercol","body")] <- 
  c("active","variable","omnivore","demersal","short and / or deep")
s1[ s1$sciName == "Arripis georgianus", c("act","feed","troph","watercol","body")] <- 
  c("active","selective plankton feeding","planktivore;higher carnivore",
    "pelagic site attached", "fusiform / normal" )
s1[ s1$sciName == "Mallotus villosus", c("act")] <- "active"
s1[ s1$sciName == "Mallotus villosus", c("feed")] <- "selective plankton feeding"
s1[ s1$sciName == "Mallotus villosus", c("troph")] <- "planktivore"
s1[ s1$sciName == "Mallotus villosus", c("watercol")] <- "pelagic non-site attached"
s1[ s1$sciName == "Mallotus villosus", c("body")] <- "elongated"
s1[ s1$sciName == "Hyporhamphus sajori", "act" ] <- "active"
s1[ s1$sciName == "Hyporhamphus sajori", "feed" ] <- "selective plankton feeding"
s1[ s1$sciName == "Hyporhamphus sajori", "troph" ] <- "planktivore"
s1[ s1$sciName == "Hyporhamphus sajori", "watercol" ] <- "pelagic site attached"
s1[ s1$sciName == "Hyporhamphus sajori", "body" ] <- "elongated"
# fill missing traits for video only taxa
missingtaxa <- s1 %>% filter(is.na(act)) %>% 
  select(sciName, Genus, family, act, feed, troph, watercol, body ) %>%  
  distinct()


# remove rows with NA for activity - omit things like barracuda and green sea turtles
s1 <- s1 %>% filter( !is.na(act) )


# duplicate entries?
s1 %>%
  arrange( Site, habitat, sciName ) %>%
  group_by( Site, habitat, sciName) %>%
  summarise( entries=length(sciName) ) %>%
  filter( entries>1 )


trait.use <- ungroup(trait.use)


# filter out taxa not in trait.final and vice versa
s1 <- s1 %>% 
  filter( sciName %in% trait.use$sciName )
trait.use <- trait.use %>% 
  filter( sciName %in% s1$sciName )
s1$sciName %in% trait.use$sciName
trait.use$sciName %in% s1$sciName

# simplify scinames
# trait.use <- trait.use %>% 
#   mutate( name=vegan::make.cepnames( trait.use$sciName ) )
# s1 <- left_join( s1, select(trait.use,sciName,name))
  
unique(s1$name) %in% unique(trait.use$name)

## write the seine-trait data to disk as well as the trait data themselvees
write_csv( s1, "Output Data/Bitemap_trait_presence+absence.csv" )
write_csv( trait.use, "Output Data/Bitemap_trait_use_presence+absence.csv" )




### summarize trait data by site
# filter out only fish?
suse <- s1 #filter( s1, phylum=="Chordata" )

# indvidual and taxon based ratios below are currently the same because only one individual per site*habitat combination
# ratio of active to ambush at each site by individuals 
rat.ind <- suse %>%
  dplyr::group_by( Site, habitat ) %>%
  dplyr::summarize(act.ratio.ind = length(act[act %in% c("active","ambush.active")]) / length(sciName)  ) 
# and by taxon
rat.tax <- suse %>%
  dplyr::select( Site, habitat, sciName, act ) %>%
  dplyr::distinct() %>%
  dplyr::group_by( Site, habitat ) %>%
  dplyr::summarize( act.ratio.tax = length(act[act %in% c("active","ambush.active")]) / length(sciName),
                    richness=length(sciName)) 
# merge
rat <- left_join( rat.ind, rat.tax )
rat.site <- left_join( rat, ss )

# write to disk
write_csv( rat, "Output Data/consumer_active_ratio_PA.csv" )


# plot
ggplot( rat.site, aes(x=abs(meanLat), y=act.ratio.ind)) + geom_point() + geom_smooth()
# ggplot( rat.site, aes(x=abs(meanLat), y=act.ratio.tax)) + geom_point() + geom_smooth()
# ggplot( rat.site, aes(x=act.ratio.tax,  y=act.ratio.ind)) + geom_point() + geom_smooth()
ggplot( rat.site, aes(x=act.ratio.ind,  y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))

# nice figure
# windows(3,3)
ggplot( rat.site, aes(x=act.ratio.tax,  y=rate)) + 
  geom_point(pch=16,alpha=0.5) +
  xlab("Ratio of active to restful taxa") + ylab("Consumption rate") +
  theme_classic()
summary(glm(rate~act.ratio.tax,data=rat.site,family="quasibinomial"))






## community weighted mean traits
# FD package
trait.xuse <- data.frame( trait.use[ trait.use$sciName %in% s1$sciName,] )
xuse <- trait.xuse %>% 
  select( act, feed, troph, watercol, body )
xuse[xuse==""] <- NA
xuse$act[ xuse$act == "active;ambush" ] <- "ambush"
# to numeric
# xuse[,c(2)] <- apply(xuse[,c(2)],2,as.numeric)
xuse$act <- as.numeric( as.factor(xuse$act))-1
summary(xuse)
row.names(xuse) <- trait.xuse$name
xuse <- xuse[ order(rownames(xuse)), ]
# xuse$squidpop <- as.numeric(as.character(xuse$squidpop))
mice::md.pattern( xuse )
# seine.keep <- s1[ s1$sciName %in% trait.keep$sciName, ]
  


### Ignore abundance patterns for Presence-Absence analysis,
# but can fill a matrix with ones and zeros
# abundance patterns
scomm <- s1 %>% mutate( abun=1 )
comm <- scomm %>%
  ungroup() %>%
  unite( "SH", Site, habitat  ) %>%
  select( SH, name, abun ) %>%
  spread( name,abun,fill=0  )
ause <- as.matrix( comm[,-1] )
rownames(ause) <- comm$SH

colnames(ause) == rownames(xuse)





gd <- gowdis( xuse, asym.bin = 1 )
m1 <- dbFD( gd, ause, w.abun = T, 
            corr="lingoes", print.pco=T, calc.FGR = T, clust.type="kmeans",
            km.inf.gr=3,km.sup.gr=20,
            calc.CWM = FALSE, m="min", calc.FDiv = FALSE, stand.FRic = TRUE )

# Community weight means
cwm <- functcomp( xuse, ause )
cwm$SH <- comm$SH
cwm <- cwm %>% 
  separate( SH, c("Site","habitat"))
cwm.ss <- left_join( cwm, ss )

write_csv( cwm, "Output Data/FunctionalDiversity_CWM_PA.csv")

# quantitative traits
ggplot( cwm.ss, aes( x=abs(meanLat), y=as.numeric(as.character(act)) )) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))
ggplot( cwm.ss, aes( x=as.numeric(as.character(act)), y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))

# qualitative traits
ggplot( cwm.ss, aes( x=abs(meanLat), y=feed)) + geom_point() 
ggplot( cwm.ss, aes( x=abs(meanLat), y=troph)) + geom_point() 
ggplot( cwm.ss, aes( x=abs(meanLat), y=watercol)) + geom_point() 

ggplot( cwm.ss, aes( x=feed, y=rate)) + geom_point() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
ggplot( cwm.ss, aes( x=troph, y=rate)) + geom_boxplot()+ geom_point() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
ggplot( cwm.ss, aes( x=watercol, y=rate)) + geom_point() +
  theme(axis.text.x = element_text(angle=45,hjust=1))




## Other dimensions of functional diversity
names(m1)
# number of functional groups
fd <- with(m1, data.frame(FRic,qual.FRic,FEve,FDis, RaoQ,FGR) )
fd$SH <- rownames(fd)
fd <- fd %>% 
  separate( SH, c("Site", "habitat"))

fds <- left_join(fd,ss)
psych::pairs.panels( fds[,c("FRic","FEve","FDis","RaoQ","FGR","rate")])
#RaoQ, FGR, Fric
# from the help page: Rao's quadratic entropy (Q) is computed from the uncorrected species-species distance matrix via divc. See Botta-Dukát (2005) for details. Rao's Q is conceptually similar to FDis, and simulations (via simul.dbFD) have shown high positive correlations between the two indices (Laliberté and Legendre 2010). Still, one potential advantage of FDis over Rao's Q is that in the unweighted case (i.e. with presence-absence data), it opens possibilities for formal statistical tests for differences in FD between two or more communities through a distance-based test for homogeneity of multivariate dispersions (Anderson 2006); see betadisper for more details.
psych::pairs.panels( fds[,c("FRic","RaoQ","FGR","rate")])
write_csv( fd, "Output Data/FunctionalDiversity_indices_PA.csv" )

with(fds,plot(FRic~cwm$act))
with(fds,plot(FRic~cwm$feed))
with(fds,plot(FRic~cwm$troph))
with(fds,plot(FRic~cwm$watercol))
with(fds,plot(FRic~cwm$body))
