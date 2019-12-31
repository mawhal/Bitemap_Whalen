#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Cleaning fish community data, calculating summaries
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2017.08.11
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################

## UPDATES FOR THIS SCRIPT
# 2017.09.21: incorporate fish biomass coefficients to calculate mass from length
# 2017.09.26: finished first draft of code to compile length-weight relationships and calculate biomass
# 2017.09.27: construct fish size distribution figures
# 2017.12.27: collect taxonomic information for invertebrates
# 2018.03.08: include information from partners about traits related to feeding and likelihood of interaction with squidpops
# 2018.03.09: calculate abundance totals and diversity (Effective Number of Species PIE) using all taxa and fish only
# 2018.04.10: look at family level diversity (maybe genus also?) in terms of presence absence of families
# 2018.04.13: remove more fish families from final list (Pholidae, Agonidae, )
# 2018.04.18: save warm water species richness for regression with predation rate
# 2018.05.06: investigate relationship bw fish diversity and biomass at different temperatures
# 2018.08.11: fish family diversity. Deal with fact that there could be over and underestimates of richness at a site given taxonomic uncertainty
# 2018.08.11: More on fish family NMDS, 
              # rotating to environmental predictors -- rotating to sstmean works well
              # including looking at numbers of axes, 
# 2018.08.25: Fish NMDS: look at number of axes. 4 axes are strongly supported and model works well
              # add all four of these axes onto the environmental data.frame -- to be used in glmer
# 2018.08.26: Different distance indices. Default with metaMDS is Bray-Curtis. Try Jaccard, Raup-Crick
              # Jaccard results basically the same as Bray. Raup-Crick more interesting, better theoretical foundations (Chase et al 2011)
# 2018.08.27: add predation rate data (summary for each site + habitat)

## Seine Data
# calculate different summaries of fish community that can be related to predation intensity
# summaries: total abundance, functional group abundance, richness, diversity metrics
#            fish biomass (need to calculate based on size and taxonomy)













########################################################################################
###  Seines were pulled for different distances. Need to take this into account!!
###  We also lack information on how wide seines were!
########################################################################################
########################################################################################
########################################################################################















# load libraries
# data.frame manipulation and joining
library(tidyverse)
library(plyr)
library(ggrepel) # for plotting with text
library(reshape2)
# geospatial data
# library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(rfishbase)
library(taxize)



###########################################
## DATA ON PREDATION AND FISH COMMUNITIES #
###########################################

# FISH SEINING DATA
seines <- read.csv( '../Data/Bitemap_Seine_ALL-DATA_20180322.csv', stringsAsFactors = FALSE, strip.white = TRUE)

# SITE LOCATIONS
siteGPS <- read.csv( '../Data/Bitemap_sites.csv')
# calculate distances between all points
siteGPS[,3:4]
# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)
rads    <- deg2rad(siteGPS[,3:4])
# calculate pairwise distances for all points
source("VincentryInverseFunction.R")
# empty matrix to hold data

dist.vif <- function( data ) {
  vif.mat <- matrix( NA, nrow=nrow(data), ncol=nrow(data) )
  for( i in 1:nrow(data) ){
    for( j in 1:nrow(data) ){
      vif.mat[i,j] = gcd.vif( data[i,2], data[i,1], data[j,2], data[j,1] )
    }
  }
}

# Predation Rates from other script
rate.env <- read.csv( "Output Data/Bitemap_DATA_analysis_20180825.csv", stringsAsFactors = FALSE )
# below, we deal with site and habitat level data, not at the same scale
# summarize to average predation rates at site and habitat level
rate.mean <- ddply( rate.env, .(Country), summarise, rate=mean( rate,na.rm=TRUE ) )
ggplot( rate.mean, aes(x=Country,y=rate)) + geom_point()#+ facet_wrap(~Country)

###########################################################################################################
######                                    Cleaning the data                                           #####
###########################################################################################################

# clean up  data
names(seines)[4] <- "habitat"
# convert all vegetated and unvegetated sites to common categories
seines$habitat[ seines$habitat %in% c("Seagrass ","seagrass" )] <- "Seagrass"
seines$habitat[ seines$habitat %in% c("unveg","Unveg","Unvegetated" )] <- "Unveg"
seines <- droplevels(seines)
# for all organisms from France, UNC2, Wales, QLD2, QLD3, multiple lengths by 10 to convert from cm to mm
seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] <- 
  seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] * 10


# select relevant columns
seines <- seines %>%
  dplyr::select( Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, Q.PA, Phylum,
                 Genus, Species, Length, Abundance, Country )


# 
# # DON'T REMOVE INVERTS AT THIS STAGE
# # which sites sampled crabs?
# pods <- seines[ seines$Phylum=="Arthropoda" & !is.na(seines$Phylum), ]
# # summarize by site
# ddply( pods, .(Country,Genus,Species), summarise, Total=sum(Abundance) )
# 
# # remove inverts
# fish <- seines[ seines$Phylum=="Chordata" & !is.na(seines$Phylum), ]


### Fish length info
# What to do with rows where we have fish abundance but no length info?
# options: ignore these rows (underestimate biomass, potentially substantially)
#  use an average size based on other length estimates for biomass (unclear if under- or overestimate)
#  use largest or smallest lengths
# for now, use average (median) size for all other length estimates

# find rows with NA for size, but omit rows where nothing was found
seinesNA <- seines[ is.na(seines$Length) & seines$Abundance>0, ]

# calculate median sizes for all other seines species from each seine
seinesL  <- seines[ !is.na(seines$Length) & seines$Abundance>0, ]
aveLengths <- ddply( seinesL, .(Site.Name,habitat,Country,Date,Genus,Species), 
                     summarize, medLength=median(Length), meanLength=mean(Length) )
# remove NA values for unknown Genera
aveLengths[is.na(aveLengths$Genus),]
aveLengths <- aveLengths[!is.na(aveLengths$Genus),]

# replace NA length values with calculated median
# rename seinesNA as seinesEST for adding length ESTimates
# Note this will only work for fishes
seinesEST <- seinesNA
# loop over all rows and get median lengths
for( i in 1:nrow(seinesEST) ){
  if(seinesEST$Phylum[i] == "Chordata" ) seinesEST$Length[i] <- aveLengths$meanLength[ aveLengths$Country==seinesEST$Country[i] &
                                                aveLengths$habitat==seinesEST$habitat[i] &  
                                                aveLengths$Site.Name==seinesEST$Site.Name[i] & 
                                                aveLengths$Date==seinesEST$Date[i] &
                                                aveLengths$Genus==seinesEST$Genus[i] &
                                                aveLengths$Species==seinesEST$Species[i] ]
}
### errors from this step that were 'corrected'
## Syngnathus from NC, no other pipeseines were caught but this one not measured, 
# so use median of all syngnathus captured from that site
## Ambassis jacksoniensis from NSW2 20170223 + 201702024 not measured to maximize survival. 
# Notes state that all were between 15 and 60 mm, so using mean of these values (37.5mm)

# merge size estimates back with original seines data
seines.clean <- full_join( seinesL, seinesEST )

# combine genus and species with name to match length-weight coefficients (see below)
seines.clean$SPECIES_NAME <- with(seines.clean, paste(Genus, Species) )


## include data from partners (using trait information) to assess likelihood of interaction with squidpops
# Note this data is qualitative in nature. It is based on expert opinion, known feeding traits, but not yet to direct observation
traits <- read.csv( '../Data/Bitemap_REQUEST_trait_siene+video_EAT_SQUID.csv', 
                    stringsAsFactors = FALSE, strip.white = TRUE )
str(traits)
# select columns
traits <- traits %>%
  dplyr::select( Country, siteName, sciName, totalCount, family, feeding=feedingTypeFishbase, trophic=trophicGroupRLS,
          waterColumn=waterColumnRLS, diel=dielActivityRLS, habitat=habitatRLS, eat.squid )
# which families are missing squid eating info?
# add phylum
# tax <- classification( sort(unique(traits$family)), db="worms" )
# taxbind <- rbind(tax)
# phyla <- taxbind %>% filter( rank=="Phylum" ) %>% 
#   dplyr::select( phylum=name, family=query )
# write.csv( phyla, "Output Data/predator_families+phyla.csv", row.names = FALSE )
phyla <- read.csv( "Output Data/predator_families+phyla.csv", stringsAsFactors = FALSE )
traits <- left_join( traits, phyla )

# Note some were left ambiguous because a clear decision could not be made
sort( unique( traits$family[is.na(traits$eat.squid)] ) )

# decision for different taxa
# not likely to eat squid
no.squid <- c( "Pleuronectidae", "Paralichthyidae", "Rhombosoleidae", "Hippolytidae", 
               "Gasterosteidae", "Embiotocidae","Pholidae","Agonidae",
               "Blenniidae", "Gobiidae" )
# likely to eat squid
yes.squid <- c( "Labridae", "Penaeidae" ) # not some labrids already marked as not eating squid because they are herbivorous

# impose decisions on data
traits$eat.squid[ is.na(traits$eat.squid) & traits$family %in% no.squid ] <- 0
traits$eat.squid[ traits$family %in% no.squid ] <- 0
traits$eat.squid[ is.na(traits$eat.squid) & traits$family %in% yes.squid ] <- 1

traits[ traits$family %in% c("Pholidae","Agonidae"), ]
traits[ traits$family %in% c("Sparidae"), ]
traits[ traits$family %in% c("Labridae") & traits$Country == "Australia (QLD3)", ]


# add site information
sites <- read.csv( "Output Data/Bitemap_BioORACLE_2018.03.13.csv", stringsAsFactors = FALSE )
sites <- left_join( sites,siteGPS[,1:2], by="Country" )

# add ocean basin information here
sites$Country
sites$basin <- c( 1,1,1,1,1,1,1,1,1,2,2,2,2,1,3,3,4,2,3,1,1,2,2,2,2,1,1,1,1,2,2,2,2,2,1,2,2,2,1,1,2)
sites$basin <- factor( sites$basin, levels=1:4, labels=c("Pacific","Atlantic", "Mediterranean", "Indian") )
sites$coast <- c( 1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,2,2,2,2,1,2,1,1,2,1,2,2,2,2,1,1,1,1,1,2,1,1,1,2,2,2)
sites$coast <- factor( sites$coast, levels=1:2, labels=c("West","East") )

# merge sites with traits
traits <- left_join( traits, sites, by=c("Country") )

# select columns
traits.clean <- traits %>%
  dplyr::select( SPECIES_NAME=sciName, phylum, family, eat.squid ) 

# remove duplicates
traits.unique <- traits.clean[!duplicated(traits.clean[,1]),]

# merge trait information into seine data
seine <- left_join( seines.clean, traits.unique, by=c("SPECIES_NAME") )


#### filter predators by likelihod of interacting with squidpops
predators <- seine[ seine$eat.squid==1, ]
# remove NA values
predators <-  predators[ !is.na(predators$habitat), ]

# How certain are we that we have a good estimate of species richness at each site?
# If a species is labeled with a genus or family designation, do we count one species? YES
# If a species is labelled as one potential member of a genus, but other known taxa
#      are present, should we count all species??
# aggregate species list for each site, then interrogate each list
by( predators, predators$Site.Name, function(z) sort(unique( z$SPECIES_NAME )) )
by( predators, predators$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# these all look good


# Build family level presence-absence matrices for those thought to eat and not thought to eat squid

# select Country, family, totalCount, eat.squid
trait.fam <- traits %>%
  dplyr::select( Country, Site,hemi,basin,coast, phylum, family, meanLat, sstmean, eat.squid, totalCount )
trait.pred <- left_join( trait.fam, rate.mean )

# aggregate to get total counts by family
family.all <- plyr::ddply( trait.fam, .(Site,family,eat.squid), summarize,
                     totalCount=sum(totalCount,na.rm=TRUE) )

family.squid <- ddply( trait.fam[ trait.fam$eat.squid==1, ], .(Site,family), summarize,
                       totalCount=sum(totalCount,na.rm=TRUE) )

# calculate family richness
famrich.all   <- ddply( family.all, .(Site), summarize, famrich.all = length(unique(family)) )
famrich.squid <- ddply( family.squid, .(Site), summarize, famrich.squid = length(unique(family)) )

famrich <- full_join( famrich.all, famrich.squid )
famrich.geo <- left_join( famrich, sites )

# number of predator families (likely to eat squidpops) as a function of distance from equator
ggplot( data=famrich.geo, aes(x=abs(meanLat),y=famrich.squid)) + geom_point() # India is low, but based on 

# presence-absence matrix for families (all and thought to eat squid)
# melt
fam.melt <- melt( trait.pred, id.vars = c(1:10,12), measure.vars=11 )
# remove Taxon NA
fam.melt <- fam.melt[ !is.na(fam.melt$family), ]
# cast
fam.pa   <- dcast( fam.melt, Site+hemi+basin+coast+meanLat+sstmean+rate~family )

# can select just those families included in biomass analysis
fam.melt2 <- fam.melt[ fam.melt$eat.squid == 1,]# & fam.melt$phylum=="Chordata", ]
# remove Taxon NA
fam.melt2 <- fam.melt2[ !is.na(fam.melt2$family), ]
# cast
fam.pa2   <- dcast( fam.melt2, Site+hemi+basin+coast+meanLat+sstmean+rate~family )

# separate community data from site data
group <- fam.pa2
fam.meta <- group[,1:7]
fam.data <- group[,-c(1:7)]
# make data presence-absence
fam.data[ fam.data>0 ] <- 1
colSums(fam.data)
rowSums(fam.data)

# which sites only had 2 families
group[ rowSums(fam.data) <= 2, ]  


# NMDS on families
library(vegan)
library(viridis)

# test a number of dimensions
mds2 <- metaMDS( fam.data, "raup", k = 2, binary=TRUE, trymax = 100 )
mds3 <- metaMDS( fam.data, "raup", k = 3, binary=TRUE, trymax = 100 )
mds4 <- metaMDS( fam.data, "raup", k = 4, binary=TRUE, trymax = 100 )
mds5 <- metaMDS( fam.data, "raup", k = 5, binary=TRUE, trymax = 100 )
mds6 <- metaMDS( fam.data, "raup", k = 6, binary=TRUE, trymax = 100 )
mds7 <- metaMDS( fam.data, "raup", k = 7, binary=TRUE, trymax = 100 )
mds8 <- metaMDS( fam.data, "raup", k = 8, binary=TRUE, trymax = 100 )
mds9 <- metaMDS( fam.data, "raup", k = 9, binary=TRUE, trymax = 100 )
mds10 <- metaMDS( fam.data, "raup", k = 10, binary=TRUE, trymax = 100 )
mds11 <- metaMDS( fam.data, "raup", k = 11, binary=TRUE, trymax = 100 )
mds12 <- metaMDS( fam.data, "raup", k = 12, binary=TRUE, trymax = 100 )
mds13 <- metaMDS( fam.data, "raup", k = 13, binary=TRUE, trymax = 100 )
mds.comp <- list( mds2, mds3, mds4, mds5, mds6, mds7, mds8, mds9, mds10, mds11, mds12, mds13 )
stress <- unlist( lapply( mds.comp, function(z) z$stress ) )
plot( stress ) # no break point

# windows(7,4)
mds <- metaMDS( fam.data, "raup", k = 3, binary=TRUE, trymax = 100 )
mds
stressplot( mds )


# rotate the NMDS to mean annual SST
mds.rot <- MDSrotate( mds, group$sstmean, na.rm=TRUE )
mds.use <- mds.rot

plot( mds.use, "sites", cex=(fam.meta$rate*3)+1 )
vec.fam <- envfit( mds.use$points, fam.data, perm=1000 )
plot( vec.fam, p.max=0.01, col="blue" )
vec.fam.df <- as.data.frame(vec.fam$vectors$arrows*sqrt(vec.fam$vectors$r))
vec.fam.df$family <- rownames(vec.fam.df)
vec.fam.df$P      <- vec.fam$vectors$pvals

fam.sel <- rownames(vec.fam$vectors[[1]])[ vec.fam$vectors$pvals < 0.05 ]

# plot species richness within each families
by( predators, predators$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# only select relevant columns of predator data.frame
pred.sel <- predators %>%
  select( Country, Site.Name, Lat, Long, habitat, SPECIES_NAME, family ) %>%
  distinct( )
# aggregate to count the number of rows in each family at each site
pred.richness <- ddply( pred.sel, .(Country, Lat, Long, habitat, family), 
                        summarize, S = length(SPECIES_NAME) )
summary(pred.richness)
# tend to be very few species within families
# look at species list for sites with high species richness within particular families
pred.richness[ pred.richness$S >2, ]
predators[ predators$Country=="Canada (BC)" & predators$habitat=="Seagrass" & predators$family=="Cottidae",]
predators[ predators$Country=="USA (CA)" & predators$habitat=="Seagrass" & predators$family=="Serranidae",]
predators[ predators$Country=="USA (TX)" & predators$habitat=="Seagrass" & predators$family=="Sciaenidae",]
# merge with the environmental data
pred.rich.env <- full_join( pred.richness, sites )
# select influential families
pre.sel <- pred.rich.env[ pred.rich.env$family %in% fam.sel, ]


# pad with zeros
S <- 0
zeros <- expand.grid( Site = fam.meta$Site, family=unique(phyla$family), S=S )
pre.pad <- full_join( pre.sel,zeros, by=c("Site","family","S") )

# this is actually abundance within families, I think
group.sel <- cbind( fam.meta, group[,fam.sel] )
group.sel.melt <- melt( group.sel, id.vars=1:6, value.name = "S", variable.name="family"  ) 
# add phylum back in
group.sel.melt <- left_join( group.sel.melt, phyla )

# 
# # windows( 6,3 )
# ggplot( group.sel.melt, aes(x=sstmean, y=S ) ) + facet_wrap(phylum~family, ncol=5) +
#   geom_smooth(se=F) + 
#   geom_point( pch=1 ) +
#   ylab( "Abundance of each family" ) + xlab("Mean Annual SST (°C)") 
# 
# ggplot( pre.pad, aes(x=sstmean, y=S ) ) + facet_wrap(~family, ncol=10) +
#   geom_smooth(se=F) + 
#   geom_point( pch=1 ) +
#   ylab( "Species richness" ) + xlab("Mean Annual SST (°C)")

# dev.off()


# summarize species richness of families associated with colder or warmer temps
vec.pos <- vec.fam.df %>% filter( MDS1>0 )
vec.neg <- vec.fam.df %>% filter( MDS1<0 )

# aggregate and combine
pos.richness <- plyr::ddply( group.sel.melt[ group.sel.melt$family %in% vec.pos$family, ], .(sstmean), 
                             summarize, S=sum(S) )
pos.richness$thermal <- "S.warm"
neg.richness <- plyr::ddply( group.sel.melt[ group.sel.melt$family %in% vec.neg$family, ], .(sstmean), 
                             summarize, S=sum(S) )
neg.richness$thermal <- "S.cold"
total.richness <- plyr::ddply( group.sel.melt, .(sstmean), 
                               summarize, S=sum(S) )
total.richness$thermal <- "S.all"
richness.groups <- rbind(pos.richness,neg.richness,total.richness)

# merge richness and sites information
# first melt and cast richness.groups
richness.melt <- melt( richness.groups, id=c(1,3) )
richness.cast <- dcast( richness.melt, sstmean~thermal )
sites.rich <- full_join( sites, richness.cast )



# plot them
# windows( 3.5,3 )
# ggplot( richness.groups, aes(x=sstmean,y=S) ) + facet_wrap( ~thermal ) + 
#   geom_smooth( se=F ) + geom_point( pch=1 ) +
#   ylab( "Species Richness" ) + xlab("Mean Annual SST (°C)")

adonis2( fam.data ~ coast:(hemi+basin), data=group )
adonis2( fam.data ~ rate, data=group )
adonis2( fam.data ~ coast:(hemi+basin)+rate, data=group )

# create interaction dummies
group$hemicoast  <- with( group, paste(hemi,coast) )
group$coastbasin <- with( group, paste(coast,basin) )

# combine data and NMDS results
mdsfd <- cbind( group, mds.use$points )

# now merge in the rotated nmds results
sites.nmds <- full_join( sites.rich, mdsfd )
# strong correlation with mean annual SST
plot( MDS1 ~ sstmean, sites.nmds )



### Ellipses for groups
plot( mds, "sites" )
ord<-ordiellipse( mds, fam.pa$coastbasin, display = "sites", 
                  kind = "se", conf = 0.95, label = T )
## new data.frame to show ellipses
# https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
# define function for ellipse
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# put everything together into a data.frame
df_ell <- data.frame()
for(g in unique(fam.pa$coastbasin)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(mdsfd[mdsfd$coastbasin==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
} # okay if this fails (it does for Indian ocean)


## unwrap ordisurf object
# https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/
#ordisurf:
ordi<-ordisurf( mds.use, group$sstmean, display = "sites", 
               kind = "se", conf = 0.95, label = T )
ordi.grid <- ordi$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.pred <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.pred$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.na <- data.frame(na.omit(ordi.pred)) #gets rid of the nas
# ordi.na #looks ready for plotting!


# the plots
# windows(6,4)
ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
  stat_contour( data = ordi.na, aes(x = x, y = y, z = z, colour = ..level..),
               binwidth = 1.5, size=0.8 ) + #can change the binwidth depending on how many contours you want
  # geom_text_repel( aes(label=Site) ) + 
  geom_point( aes(size=rate) ) +
  scale_color_viridis( name = "Mean\nAnnual SST" ) +
  xlim(c(-1.2,0.9)) + ylim(c(-0.8,0.85)) 



with( mdsfd, cor.test( MDS1, sstmean ) )
cor.test( mds$points[,1], mdsfd$sstmean )
cor.test( mds.rot$points[,1], mdsfd$sstmean )

# 
vec.fam.sel <- vec.fam.df[vec.fam.df$P < 0.02,]
ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
  geom_point(  col='slateblue' ) +
  geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black") + 
  geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family),
                  point.padding=1, box.padding = 0.1, 
                  size=4, col="black") +
  xlim(c(-1.2,0.9)) + ylim(c(-0.8,0.85)) 
# 
# ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
#   geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, col=group), size=0.5) +
#   geom_text_repel( aes(label=Site) ) + geom_point() +
#   scale_color_discrete( name = "Coast" )  +
#   xlim(c(-1.5,1.5)) + ylim(c(-1.2,0.75))


# what's driving communities?
# SIMPER? Bubble Plots

# try PCA, best fits??
# RDA? Correspondence analysis


## Redundant Descriminant analysis
rda1 <- rda(  fam.data ~ sstmean + meanLat, group, na.action=na.exclude )
plot(rda1)

## Canonical Correspondence analysis
cca1 <- cca(  fam.data ~ sstmean, group, na.action=na.exclude )
cca1 <- cca(  fam.data ~ sstmean, group, na.action=na.exclude )
plot(cca1)

# get rid of NA
fam.na <- na.omit( fam.pa )
# separate community data from site data
fam.meta <- fam.na[,1:6]
fam.data <- fam.na[,-c(1:6)]
## environmental data matrix
vare.cca <- cca( fam.data ~ sstmean + basin, fam.meta )
vare.cca
plot(vare.cca)

## RDA
rda1 <- rda( fam.data ~ sstmean + hemi*basin*coast, fam.meta )
anova(rda1)
summary(rda1)
plot(rda1) 





### add a switch to only allow taxa associated with warm waters in the calculation
pred.all  <- predators
pred.warm <- predators[ predators$family %in% vec.pos$family, ]

# decide which predators set to use
pred.use <- pred.all



##############################
# Predator size distribution #
##############################


# prob first need to omit rows lacking a size, then use rep() to make a vector of lengths for for each count


# only include rows with a length estimate (pred biomass assumed average lengths for preds not measured...see below)
pred.size <- pred.use[ !is.na(pred.use$Length), ]


# use ggplot to make violin plots for each site x habitat combination
ggplot( pred.size, aes(x=habitat, y=log10(Length) )) + geom_boxplot(width=0.5) + facet_wrap(~Country)
ggplot( pred.size, aes(x=habitat, y=log10(Length) )) + geom_violin() + facet_wrap(~Country)


# pred length by latitude
ggplot( pred.size, aes(x=abs(Lat), y=log10(Length) )) + geom_point() + facet_wrap(~habitat) + geom_smooth()
pred.size.median <- ddply( pred.size, .(Country,Lat,habitat), summarize, Length=median(Length,na.rm=T) )
ggplot( pred.size.median, aes(x=abs(Lat), y=Length)) + 
  geom_smooth(se=F) + geom_point()


# calculate median size for each seine
pred.med.len <- ddply( pred.size, .(Country,Site.Name,Lat,Date,Time),
                       summarize, medLength=median(Length,na.rm=T) )

ggplot( pred.med.len, aes(x=abs(Lat), y=medLength )) + geom_point() + geom_smooth() +
  xlab('degrees latitude from equator') + ylab('median length of fish per seine')


# #########################
# # Fish Trophic Levels   #
# #########################
# library(rfishbase)
# fishSpecies <- sort( unique( fish.clean$SPECIES_NAME ) )
# 
# fishEcology <- ecology( fishSpecies, fields = c("SpecCode", "FeedingType", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph") )
# 
# # sometimes DietTroph (based on whole range of food items) unavailable but FoodTroph (based on individual food items) more often is available
# fishEcology$Trophic.Level <- fishEcology$DietTroph
# for( i in 1:nrow(fishEcology) ) {
#   if( is.na(fishEcology$DietTroph[i]) )   fishEcology$Trophic.Level[i] <- fishEcology$FoodTroph[i]
# }
# 
# # which taxa are missing trophic levels?
# trophMiss <- fishEcology[ which(is.na(fishEcology$Trophic.Level)), ]
# write.csv( trophMiss, "../Data/Fish Biomass + Traits/Bitemap_2017_fishBase_missingTrophic.csv", row.names = FALSE )
# # read in data.frame with info added from fishbase
# trophFill <- read.csv( "../Data/Fish Biomass + Traits/Bitemap_2017_fishBase_missingTrophic.csv" )
# trophFill$Trophic.Level <- trophFill$Trophic
# trophFill$fishbaseTrophEst <- c(2.9,3.2,3.3,3.3,3.4,3.4,3.3,3.2,3.2,3.4,3.7)
# fishEcology <- full_join( fishEcology, trophFill[,c('sciname','fishbaseTrophEst')] )
# for( i in 1:nrow(fishEcology) ) {
#   if( is.na(fishEcology$Trophic[i]) )   fishEcology$Trophic.Level[i] <- fishEcology$fishbaseTrophEst[i]
# }
# 
# 
# # 
# # some trophic level info added manually
# trophManual <- read.csv( "../Data/Fish Biomass + Traits/Bitemap_2017_fishBase_noEcology_ecologyAdded.csv" )
# names(trophManual) <- c("sciname","Trophic.Level")
# #
# fishINFO <- full_join( fishEcology, trophManual )
# fishSpecies[ !(fishSpecies %in% fishINFO$sciname)] # several species not represented
# 
# ## join Bitemap fish data with trait data from Reef Life Survey
# # read in fish trait data from (source = RLS?)
# fishtraits <- read.csv("../Data/Fish Biomass + Traits/Traits_all-species_edit.csv")
# # create column for joining with fishEcology
# fishtraits$sciname <- fishtraits$CURRENT_TAXONOMIC_NAME
# 
# 
# # first, make sure all bitemap species are represented
# fishINFO <- full_join(fishINFO,data.frame(sciname=fishSpecies))
# # now, left_join so we don't include taxa that aren't represented in Bitemap
# fishINFO <- left_join( fishINFO, fishtraits ) # lots of gaps
# # select relevant columns
# fishINFO <- fishINFO %>%
#   dplyr::select( sciname, FeedingType, Trophic.group, Water.column, Diel.Activity, Habitat, Trophic.Level )
# 
# 
# # (fishbaseMiss <- fishSpecies[ !(fishSpecies %in% fishINFO$sciname)])
# #  write.csv( data.frame(fishbaseMiss), "../Data/Fish Biomass/Bitemap_2017_fishBase_noEcology.csv", row.names=FALSE )
# 
# 
# 
# 
# 
# 
# #################
# # Fish Taxonomy #
# #################
# # add taxonomy for fishes
# library(taxize)
# # strip sp from unknown species before running classification
# scisplit <- strsplit( fishINFO$sciname, " ", fixed=TRUE )
# which(lapply(scisplit,length)==1)
# # scisplit[[1460]] <- c("Clupeidae","sp.")
# genusSpecies <- lapply(scisplit,function(z) data.frame(genus=z[1],species=z[2]))
# genusSpecies <- do.call(rbind,genusSpecies)
# fishINFO <- cbind(fishINFO,genusSpecies)
# 
# # omit sp designations
# for(i in 1:nrow(fishINFO)){
#   if(fishINFO$species[i] %in% c("sp.","spp.","sp. [Belize]","sp. [bottomei]")) fishINFO$species[i] <- ""
# }
# 
# # recombine genus and species without sp. designations (spaces should not matter here)
# fishINFO$sciname2 <- with(fishINFO, paste(genus,species) )
# fishINFO <- fishINFO[ order(fishINFO$sciname2), ]
# # get rid of duplicate entries
# fishINFO <- fishINFO[ !duplicated(fishINFO$sciname2), ]
# 
# # run classification on all fish taxa
# fishtax <- classification(fishINFO$sciname2,db="ncbi")
# 
# 
# # only keep non-NA ones
# which( unlist(lapply( fishtax, function(z) all(!is.na(z)) )) )
# cf <- fishtax[which( unlist(lapply( fishtax, function(z) all(!is.na(z)) )) )]
# 
# # combine the results
# # cannot cbind with NA values, so rbind
# cr <- do.call( rbind,cf )
# cr$unique <- rownames(cr)
# cr$sciname2 <- unlist(lapply( strsplit( cr$unique, split = ".", fixed=TRUE ), function(z) z[1] ))
# 
# # melt and recast
# cm <- melt(cr)
# # omit no rank
# cm <- cm[cm$rank!="no rank",]
# ccast <- dcast( cm, sciname2~rank, value.var = "name" )
# 
# # merge taxonomy with traits
# fishes <- full_join( fishINFO, ccast, "sciname2" )
# head(fishes)

# write this to disk
# write.csv( fishes, "Output Data/fish_traits+taxonomy_20180122.csv", row.names=FALSE )




##################
# Fish diversity #
##################


# combine seines and sites
pred.site <- left_join( pred.use, siteGPS )


# summarize seine data based on counts of each taxon
pred.abund <- ddply( pred.site, .(Country,meanLat,meanLong,habitat,phylum,Genus,Species), 
                      summarise, Abundance=sum(Abundance) )
# summarize seine data based on each indivdual seine (use site.name, date, time)
pred.abund.each <- ddply( pred.site, .(Country,Site.Name,Distance,Date,Time,Lat,Long,habitat,phylum,Genus,Species), 
                      summarise, Abundance=sum(Abundance) )



# for species diversity metrics, we need counts of each species in each habitat
# # omit invertebrates
# pred.abund <- seine.abund[ seine.abund$Phylum=="Chordata", ]
# pred.abund.each <- seine.abund.each[ seine.abund.each$Phylum=="Chordata", ]

# consider using different diversity metrics (e.g. ENSpie as in Chase & Knight 2013)
ENSPIE <- function(prop){
  ifelse( sum(prop,na.rm=T)>0, 1 / sum(prop^2, na.rm=T), NA ) 
}   
prop <- c(0.5,0.4,0.1,NA)
ENSPIE(prop)
ENSPIE(0)
# ENSPIE is less scale dependent than many other diversity metrics (arguable important here), 
# but tends to be sensitive to aggregation (which can be affected by spatial scale, 
# and likely an issue for schooling predators and seining)

# convert abundance to proportion
# calculate ENSPIE for each seine
# calculate total abundance for each seine
pred.totals.each <- ddply( pred.abund.each, .(Country,Site.Name,Distance,Date,Time,Lat,Long,habitat),
                      summarize, Total = sum(Abundance), TotalFish = sum(Abundance[phylum=="Chordata"]) )
head(pred.comm.each[ with(pred.comm.each,order(Country,Site.Name,Distance,Date,Time,habitat)), ])

pred.comm.each <- left_join(pred.abund.each,pred.totals.each)
# remove fish total information for taxa that are not fish
pred.comm.each$TotalFish[ pred.comm.each$phylum != "Chordata" ] <- NA
# calculate the proportion of each species in each site
pred.comm.each$prop <- pred.comm.each$Abundance/pred.comm.each$Total
pred.comm.each$propFish <- pred.comm.each$Abundance/pred.comm.each$TotalFish
# calculate ENSPIE for each seine
ENSPIE.seine <- ddply( pred.comm.each, .(Country,Site.Name,Distance,Date,Time,Lat,Long,habitat), 
                      summarize, ENSPIE = ENSPIE(prop), ENSPIEfish = ENSPIE(propFish) )


# calculate total abundance for each site (lump all replicate seines together)
pred.totals <- ddply( pred.abund, .(Country,meanLat,meanLong,habitat),
       summarize, Total = sum(Abundance), TotalFish = sum(Abundance[phylum=="Chordata"]))
# get rid of NA rows
pred.totals <- pred.totals[!is.na(pred.totals$Country),]
# include totals as a column in pred.abund
pred.comm <- left_join(pred.abund.each,pred.totals)
# calculate the proportion of each species in each site
pred.comm$prop <- pred.comm$Abundance/pred.comm$Total

# calculate ENSPIE for each site
ENSPIE.site <- ddply( pred.comm, .(Country,meanLat,meanLong,habitat), 
       summarize, ENSPIE = ENSPIE(prop))
# get rid of NA rows
ENSPIE.site <- ENSPIE.site[!is.na(ENSPIE.site$Country),]

# calculate effect size as difference between ENSPIE estimates
ENSPIE.site <- ENSPIE.site[ with(ENSPIE.site,order(Country,habitat)), ]

# make all Unveg sites negative
ENSPIE.site$ENSPIEdiff <- ENSPIE.site$ENSPIE
ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] <- -(ENSPIE.site$ENSPIEdiff[ENSPIE.site$habitat=="Unveg"] )
ENSPIE.diff <- ddply( ENSPIE.site, .(Country,meanLat,meanLong), summarise, ENSPIEdiff=sum(ENSPIEdiff) )


ggplot( ENSPIE.site, aes(x=habitat,y=ENSPIE)) + geom_point() + 
  geom_line(aes(group=Country)) + facet_wrap(~Country)
ggplot( ENSPIE.site, aes(x=abs(meanLat),y=ENSPIE)) + geom_point() + facet_wrap(~habitat, ncol=2) + geom_smooth()
ggplot( ENSPIE.diff, aes(x=abs(meanLat),y=ENSPIE)) + geom_point() + geom_smooth()

ggplot(ENSPIE.site, aes(x=habitat, y=ENSPIE, fill = habitat)) + geom_boxplot(width=0.5, notch=FALSE) +
  geom_point(alpha=0.1) +
  # facet_grid(.~Year) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  labs(title="Diveristy of Select Predators By Habitat Type",
       x = "Type of Habitat", y="Diversity\n(effective number of species)") 


# average by site
ggplot( pred.totals.each, aes(x=abs(Lat),y=log10(Total)) ) + geom_point() + geom_smooth()
ggplot( pred.totals, aes(x=abs(meanLat),y=log10(Total)) ) + geom_point() + 
  geom_smooth(se=F) + geom_smooth(method='lm')





###################################################################################
# PREDATOR BIOMASS                                                                #
###################################################################################

# read biomass coefficients
coef <- read.csv("../Data/Fish Biomass + Traits/20160710_RLS_biomass_coefs_20170921.csv")
# omit rows from coef for which we have no estimates
coef <- coef[ !is.na(coef$A), ]
# remove duplicates from coef
coef <- coef[ !duplicated(coef), ]



# identify species not in character list
lookup <- sort(unique(fish.clean$SPECIES_NAME[!(fish.clean$SPECIES_NAME %in% coef$SPECIES_NAME)]))
sort(unique(fish.clean$SPECIES_NAME[(fish.clean$SPECIES_NAME %in% coef$SPECIES_NAME)]))

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
fishbaseLW <- ddply( fishbaseLW, .(sciname), 
                     summarize, a=mean(a),aTL=mean(aTL),b=mean(b) )
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

# match seine data with biomass conversion coefficients
pred_biom <- left_join(pred.use, biomassCoef, by="SPECIES_NAME" )

# remove extraneous columns
pred_biom <- pred_biom %>%
  dplyr::select(Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, Q.PA, Phylum, Genus, Species, 
         Length, Abundance, Country, SPECIES_NAME, A, B )


# which taxa don't have length-weight regression estimates
miss <- pred_biom[ is.na(pred_biom$A), ]
sort(unique(miss$SPECIES_NAME))

# create biomass column -- note that I've divided all lengths by 10 to convert mm to cm, which is the unit 
pred_biom$biomass <- with(pred_biom, Abundance*(A*(Length/10)^B))

# which ones have no estimate
pred_biom[ is.na(pred_biom$biomass), ]

# # make biomass into kg
# pred_biom$biomass <- pred_biom$biomass*0.001

# Calculate total fish biomass in each seine
biomass_total <- ddply( pred_biom, .(Site.Name,Country,Date,Time,Lat,Long,Distance,habitat), 
                        summarize, biomass=sum(biomass, na.rm=TRUE) )
sort(unique(biomass_total$Country))
sort(unique(predators$Country))

# standardize biomass by seine distance
biomass_total$biomass.std <- with(biomass_total, biomass/Distance)

# graph biomass by habitat type
ggplot(biomass_total, aes(x=habitat, y=log10(biomass.std*1000), fill = habitat)) + geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  # facet_grid(.~Year) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  labs(title="Fish biomass by habitat type",
       x = "Type of Habitat", y="Biomass") 
  

# total biomass by Latitude
ggplot(biomass_total, aes(x=Lat, y=log10(biomass.std*1000))) + geom_point() + facet_wrap(~habitat,ncol=1)
ggplot(biomass_total, aes(x=abs(Lat), y=log10(biomass.std*1000))) + geom_point() + facet_wrap(~habitat,ncol=1) +
  geom_smooth(method='lm')

ggplot( biomass_total, aes(x=abs(Lat),y=log10(biomass))) + geom_point() + facet_wrap(~habitat, ncol=2) + geom_smooth()

# calculate average biomass by seine and habitat
biomass_mean <- ddply( biomass_total, .(Site.Name,Country,Date,Lat,Long,Distance,habitat), 
                       summarize, biomass=mean(biomass, na.rm=TRUE), 
                       biomass.std=mean(biomass.std,na.rm=TRUE) )
ggplot(biomass_mean, aes(x=abs(Lat), y=log10(biomass.std*1000))) + geom_point() + facet_wrap(~habitat,ncol=2) +
  geom_smooth(method='lm')


# calculate average biomass by site and habitat
biomass_site <- ddply( biomass_total, .(Country,habitat), 
                       summarize, meanLat=mean(Lat), biomass=mean(biomass, na.rm=TRUE), 
                       biomass.std=mean(biomass.std,na.rm=TRUE) )

ggplot(biomass_site, aes(x=abs(meanLat), y=log10(biomass.std*1000))) + geom_point() + facet_wrap(~habitat,ncol=2) +
  geom_smooth(se=F) + geom_smooth(method='lm')

# graph biomass by habitat type
ggplot(biomass_site, aes(x=habitat, y=log10(biomass.std*1000), fill = habitat)) + 
  geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  # facet_grid(.~Year) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  labs(title="Fish biomass by habitat type",
       x = "Type of Habitat", y="Biomass") 





################################################
## Combine data summaries and save to disk    ##
################################################

# summarize these at the highest level of detail (smalles scale) possible
# by country, habitat, replicate seine (lat, long, date): not resolvable for squidpops at every site
# only means by habitat and country will be directly comparable to squidpops at every site


# Abundance ('Total' above, but rename to abundance here)
abundance <- pred.totals.each 
abundance <- abundance %>%
  mutate( Total.std=Total/Distance, TotalFish.std=TotalFish/Distance )

# diversity (ENSPIE)
diversity <- ENSPIE.seine 

# Biomass
biomass <- biomass_total 


## merge
dim(abundance); dim(diversity); dim(biomass)
predator1 <- left_join(abundance,diversity)
predator  <- left_join(predator1,biomass)
# when different ways of counting species are used, some sites are dropped
# merge with a bigger data set, the select relevant column
pred.all <- left_join( predators,predator )

head(pred.all)
unique( pred.all$Country[ is.na(pred.all$biomass) ] )
# convert NA to 0 because there is no abundance, diversity, or biomass if enough species filtered out
pred.all[ is.na(pred.all$biomass), 
          c("Total","TotalFish","Total.std","TotalFish.std","ENSPIE","ENSPIEfish","biomass") ] <- 0


# relationships between predator summaries
ggplot( pred.all, aes(x=log10(TotalFish),y=log10(biomass*1000))) + geom_point() +
  geom_smooth()
ggplot( pred.all[pred.all$biomass>0,], aes(x=ENSPIEfish,y=(biomass))) + geom_point() +
  geom_smooth( ) 
ggplot( pred.all, aes(x=log10(TotalFish),y=ENSPIE)) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=log10(Total),y=ENSPIE)) + geom_point() +geom_smooth(method='lm')

# latitudinal gradients
ggplot( pred.all, aes(x=abs(Lat),y=ENSPIEfish)) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=(TotalFish))) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=log10(Total/Distance))) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=log10(biomass))) + geom_point() +geom_smooth(method='lm')
ggplot( pred.all, aes(x=abs(Lat),y=log10(biomass/Distance))) + geom_point() +geom_smooth(method='lm')

# boxplots of abundance, diversity, and biomass by habitat
# point colors
seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
unveg.color    <- "#d8b365"
# seagrass.color <- "#4fc179"    # see http://www.color-hex.com/color-palette/31668
# unveg.color    <- "blanchedalmond"

# graph predator abundance by habitat type
ggplot(pred.all, aes(x=habitat, y=Total.std*25, fill = habitat)) + 
  geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  scale_y_log10(breaks=c(1,10,100,1000)) +
  scale_x_discrete( labels=c("Seagrass", "Unvegetated")) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) + 
  labs( title="Predator abundance\nby habitat type",
       x = "Type of Habitat", y="Abundance (CPUE)" ) 

# graph fish abundance by habitat type
ggplot(pred.all, aes(x=habitat, y=TotalFish.std*25, fill = habitat)) + 
  geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  scale_y_log10(breaks=c(1,10,100,1000)) +
  scale_x_discrete( labels=c("Seagrass", "Unvegetated")) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) + 
  labs(title="Fish abundance\nby habitat type",
       x = "Type of Habitat", y="Abundance (CPUE)") 

# graph predator diversity by habitat type
ggplot(pred.all, aes(x=habitat, y=ENSPIE, fill = habitat)) + 
  geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  scale_x_discrete( labels=c("Seagrass", "Unvegetated")) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) + 
  labs(title="Predator diversity\nby habitat type",
       x = "Type of Habitat", y="Effective Number\nof Species") 

# graph fish diversity by habitat type
ggplot(pred.all, aes(x=habitat, y=ENSPIEfish, fill = habitat)) + 
  geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  scale_x_discrete( labels=c("Seagrass", "Unvegetated")) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) + 
  labs(title="Fish diversity\nby habitat type",
       x = "Type of Habitat", y="Effective Number of Species") 

# graph fish diversity by habitat type
# windows(3,3)
ggplot(pred.all, aes(x=habitat, y=biomass.std*1000+1, fill = habitat)) + 
  geom_boxplot(width=0.5,notch=TRUE) +
  geom_point(alpha=0.1) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none") +
  scale_y_log10(breaks=c(1,10,100), labels=c(0,10,100)) +
  scale_x_discrete( labels=c("Seagrass", "Unvegetated")) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) + 
  labs( x = "Type of Habitat", y="Biomass\n(grams per meter)") 




############################################################################################
# relationship between biomass and diversity and temperature
# merge predator summaries with environmental data
pred.env <- left_join( pred.all, sites, by=c("Country") )
# select relevant columns and remove duplicate entries
pred.env <- pred.env %>%
  select( Site, Site.Name, Country, meanLat, meanLong, habitat, 
          TotalFish.std, ENSPIEfish, biomass.std, sstmean, hemi ) %>%
  distinct()
# get rid of rows without sstmean estimates
pred.env <- pred.env[ !is.na(pred.env$sstmean), ]
  

# plot biomass ~ diversity
ggplot( data=pred.env, aes( x=ENSPIEfish, y=log10(biomass.std) )) +
  geom_point() + geom_smooth()

# two bins for temperature (choose a cutoff based on differences in fish composition or predation rate...
#    current cutoff is 15.6 C)
ggplot( data=pred.env, aes( x=ENSPIEfish, y=log10(biomass.std+1), 
                            col=cut(sstmean,breaks=2), group=cut(sstmean,breaks=2) )) +
  geom_point( size=3, alpha=0.5 ) + geom_smooth( method="lm", se=F) +
  ylab( expression(paste(log[10],"(biomass)") )) + xlab( "Effective Number of Species (PIE)" ) +
  scale_color_manual( name="SST bins (°C)", values=c("blue","red") ) +
  theme(legend.justification=c(1,0), legend.position=c(1,0.77))

# maybe not surprising, since lots of cold water species are removed...what if we look at all species?
library(lmerTest)
pred.env$tempbin <- cut(pred.env$sstmean,breaks=2)
pred.env$tempbin <- factor( pred.env$tempbin, levels=c( "(15.6,26]","(5.3,15.6]" ) )
summary( lmer( log10(biomass.std+1) ~ ENSPIEfish * tempbin + (1|Country), data=pred.env ))




# sites with high abundance
pred.all[pred.all$Total>300,] 

## write predator summaries to disk
write.csv(  pred.all, "Bitemap_SEINE_summaries_20180825.csv", row.names=F )
write.csv(  sites.nmds, "Output Data/Bitemap_BioORACLE_20180825.csv", row.names = FALSE )
