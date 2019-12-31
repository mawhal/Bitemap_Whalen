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



## Seine Data
# calculate different summaries of fish community that can be related to squidpop consumption intensity
# summaries: total abundance, functional group abundance, richness, diversity metrics
#            fish biomass (need to calculate based on size and taxonomy)

########################################################################################
###  Seines were pulled for different distances. We take this into account,
###  but we still lack information on how wide seines were!
########################################################################################


## Video Data
# only use video data to identify which fish families were present


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

# #### SEINE and TRAIT data
# here::here()
# seine <- read_csv( "Output Data/Bitemap_SEINE_clean.csv" )
# 
# ### trait and environment data
# traits <- read_csv("Output Data/Bitemap_seine_trait.csv")
# 
# ### squidpop consumption rates
# rate.env <- read_csv( "Output Data/Bitemap_DATA_analysis_20180825.csv")
# # below, we deal with site and habitat level data, not at the same scale
# # summarize to average predation rates at site and habitat level -- need to be careful with habitat because this has problems with traits
# rate.mean <- rate.env %>%
#   dplyr::group_by(Country,  habitat) %>%
#   dplyr::summarise( rate=mean( rate,na.rm=TRUE ) )
# rate.mean2 <- rate.env %>%
#   dplyr::group_by(Country, Site) %>%
#   dplyr::summarise( rate=mean( rate,na.rm=TRUE ),
#                     sstmean=mean( sstmean, na.rm=TRUE ),
#                     salinity=mean( salinity, na.rm=TRUE ),
#                     temp=mean( temp, na.rm=TRUE ),
#                     sal=mean( sal, na.rm=TRUE ),
#                     meanLat=mean( meanLat, na.rm=TRUE ),
#                     meanLong=mean( meanLong, na.rm=TRUE ),
#                     chlomean=mean( chlomean, na.rm=TRUE ),
#                     hemi=unique( hemi),
#                     hemicoast=unique( hemicoast),
#                     coastbasin=unique( coastbasin),
#                     basin=unique( basin ) )
# 
# ### site information
# sites <- read_csv( "Output Data/Bitemap_BioORACLE_20190107.csv")
# 
# 
# 
# 
# #### filter consumers by likelihod of interacting with squidpops
# seine[is.na(seine$eat.squid),]
# 
# ### SWITCH next two rows - all species or just the ones we think will eat squidpops
# # consumers <- seine[ seine$eat.squid==1, ]
# consumers <- seine
# 
# 
# # remove NA values
# consumers <-  consumers[ !is.na(consumers$habitat), ]
# 
# # How certain are we that we have a good estimate of species richness at each site?
# # If a species is labeled with a genus or family designation, do we count one species? YES
# # If a species is labelled as one potential member of a genus, but other known taxa
# #      are present, should we count all species??
# # aggregate species list for each site, then interrogate each list
# # by( consumers, consumers$Site.Name, function(z) sort(unique( z$SPECIES_NAME )) )
# by( consumers, list(consumers$Country,consumers$habitat), function(z) sort(unique( z$SPECIES_NAME )) )
# # these all look good
# # Sparidae
# by( consumers[consumers$family=="Sparidae",], 
#     consumers[consumers$family=="Sparidae",]$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# # Cottidae
# by( consumers[consumers$family=="Cottidae",], 
#     consumers[consumers$family=="Cottidae",]$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# 
# # standardize predator abundance by seined area and volume
# consumers <- consumers %>%
#   # calculate area and approximate volume seined
#   # assuming the seined volume is a wedge-shaped object (triangular in cross-section) because seines pulled from depth to shore, hence the division by 2
#   dplyr::mutate( area = Distance*width.transect, volume=(Distance*Depth)/2 * width.transect  ) %>%
#   # standardize abundance by area and volume (Catch per unit effort) units are meters squared and cubed
#   dplyr::mutate( cpua=Abundance/area, cpuv=Abundance/volume )
# 
# # get total abundance by Site and seine 
# abundance.all <- consumers %>%
#   dplyr::group_by( Site.Name, Country, Date, Time, Lat, Long, Distance, area, volume, habitat) %>%
#   dplyr::summarise( total.cpua=sum(cpua, na.rm=T), total.cpuv=sum(cpuv, na.rm=T))
# 
# # NEED TO DO THIS FOR FISH SEPARATELY?
# abundance.fish <- consumers %>%
#   filter( phylum=="Chordata" ) %>%
#   dplyr::group_by( Site.Name, Country, Date, Time, Lat, Long, Distance, area, volume, habitat) %>%
#   dplyr::summarise( total.cpua.fish=sum(cpua, na.rm=T), total.cpuv.fish=sum(cpuv, na.rm=T))
# 
# # join these
# abundance <- left_join( abundance.all, abundance.fish )
# abundance$total.cpua.fish[ is.na(abundance$total.cpua.fish) ] <- 0
# abundance$total.cpuv.fish[ is.na(abundance$total.cpuv.fish) ] <- 0
# 
# # Build family level abundance (+presence-absence) matrices for families thought to eat squid
# # For each seine
# # Note: totalCount from data.frame traits is not valid and must be standardized
# # use consumers
# families <- consumers %>%
#   dplyr::group_by( Country, Site.Name, Lat, Long, habitat, Date, Time, Depth, Distance, width.transect,
#             Phylum, family ) %>%
#   dplyr::summarise( cpua = sum(cpua), cpuv = sum(cpuv) ) 
#   
# # add some trait information back in
# traits.clean2 <- traits %>%
#   dplyr::select( Country, meanLat, sstmean, Site, basin, coast, family, eat.squid ) %>%
#   dplyr::distinct() %>%
#   dplyr::arrange( Country )
#   
# fam.trait <- left_join( families, traits.clean2 )
# # NA values associated with previously unknown families
# fam.trait[ is.na(fam.trait$Site), c("meanLat","sstmean","Site","basin","coast")] <- fam.trait[ which(is.na(fam.trait$Site))+1, c("meanLat","sstmean","Site","basin","coast")]
# 
# # a few families to consider
# fam.trait[fam.trait$family=="Fundulidae",]
# fam.trait[fam.trait$family=="Palaemonidae",]
# fam.trait[fam.trait$family=="Portunidae",]
# fam.trait[fam.trait$family=="Cancridae",]
# fam.trait[fam.trait$family=="Cottidae",]
# fam.trait[fam.trait$family=="Haemulidae",]
# fam.trait[fam.trait$family=="Lutjanidae",]
# 
# # use gather and spread instead of melt and cast
# # first get average abundance across all seines within a particular habitat at a site
# fam.cpua <- fam.trait %>%
#   dplyr::group_by( Country, Site, family ) %>%
#   dplyr::summarise( mean.cpua = mean(cpua,na.rm=T) ) %>%
#   ungroup()
# 
# fam.cpua.spread <- fam.trait %>%
#   dplyr::group_by( Country, Site, family ) %>%
#   dplyr::summarise( mean.cpua = mean(cpua,na.rm=T) ) %>%
#   ungroup() %>%
#   # then spread out the families
#   spread( family, mean.cpua, fill = 0 )
# 
# 
# 
# 
# # join with rate.mean (note again that individual seines DO NOT match up with squidpop consumption assays)
# fam.cpua.rate <- left_join( fam.cpua, rate.mean2 )
# #
# # look family by family at the relationship between abundance and squidpop consumption rate
# # windows()
# # ggplot( fam.cpua.rate, aes(x=mean.cpua,y=rate, col=habitat) ) + facet_wrap(~family, scales="free") +
# #   geom_point() + geom_smooth(method='lm',se=F)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # separate community data from site data
# group <- fam.cpua.spread
# fam.meta <- group[,1:2]
# # fam.meta$ID <- with( fam.meta, paste( substr(Site,1,4), substr(habitat,1,1),sep="_") )
# fam.data <- as.matrix(group[,-c(1:2)])
# row.names(fam.data) <- fam.meta$Site
# # join with rate
# fam.meta <- left_join(fam.meta, rate.mean2)
# row.names(fam.data)
# 
# 
# # site and family summaries
# sort( colSums(fam.data) )
# rowSums(fam.data)
# 
# 
# # add omnivory and herbivory?
# 
# # add summarized trophic trait information from FishBase script
# trophic <- read.csv( "Output Data/omnivory_site.csv", stringsAsFactors = FALSE )
# trophic$meanprop <- apply( trophic[,c("omniprop","herbprop")], 1, mean )
# trophsite <- left_join(sites,trophic)
# fam.meta2 <- left_join( fam.meta, trophsite )
# 


## read data

# seine abundance + biomass
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )

# environmental data from remote sensing and oceanographic expeditions
oracle <- read_csv( "Output Data/Bitemap_BioORACLE_20190107.csv")[,1:36]

# in situ measurements
directenv <- read_csv( '../Data/Environmental Data/BiteMap_EnvironmentalData_Compilation_20190228.csv' )
# combine in situ temperature and aqua
directenv$Temp[ is.na(directenv$Temp) ] <- directenv$aqua[ is.na(directenv$Temp) ]
names(directenv)[names(directenv)=="Seagrass.Unveg"] <- "habitat"
# directenv$habitat[ directenv$habitat=="seagrass"] <- "Seagrass"
directenv$habitat[ directenv$habitat=="Unveg"] <- "Unvegetated"

# get average temp and salinity for each site on each day
directenvmean <- directenv %>% 
  dplyr::group_by(Site.Name,Lat,Long,habitat,Date) %>% 
  dplyr::summarize( temp=mean(Temp, na.rm=T), sal=mean(Sal, na.rm=T) )

# get average temp and salinity for each Country and habitat
directenvsite <- ddply( directenv, .(Country,habitat), summarize,
                        temp=mean(Temp), sal=mean(Sal) )

# merge Bio-ORACLE data and in situ measurements
env <- full_join( oracle, directenvsite )


# also need population size and fishind pressure index



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
mean.catch.fam <- mean.catch.site %>% 
  dplyr::group_by( Site, habitat, family ) %>% 
  dplyr::summarize( cpua=sum(cpua) )

# spread out again
fam.catch.wide <- mean.catch.fam %>% 
  group_by( Site, habitat ) %>% 
  spread( family, cpua )


# can average across habitat types, too

fam.data <- fam.catch.wide[,-c(1:2)]

fam.meta <- fam.catch.wide[,1:2]
fam.meta <- left_join( fam.meta, env )

fam.meta <- fam.meta %>%  tidyr::unite( SH, Site, habitat, remove=FALSE )

rownames(fam.data) <- fam.meta$SH

##### MULTIVARIATE SECTION

# NMDS on families


# windows(7,4)
mds <- metaMDS( fam.data, "bray", k = 3, binary=TRUE, trymax = 100 )
mds
 
stressplot( mds )
plot(mds, display = "sites", cex=fam.meta$rate*10+1 )
points(mds, display = "species", cex=1, col='red', pch=3)

mdsrate <- data.frame( scores(mds), fam.meta )
ggplot( data=mdsrate, aes(x=NMDS1,y=rate) ) + geom_point(size=5) +
  geom_smooth( aes(group=1), method='glm', method.args=list(family='quasibinomial')) 

# capscale
# row.names(fam.data) <- fam.meta$Site
# reduce number of families to most abundant
# also get rid of Embiotocids because they are super abundant at a few sites and skew the pattern
# Y  <- fam.data[, which( colSums(fam.data) > 0.02 & colSums(fam.data) < 4 ) ] 
Y <- fam.data
Y2 <- ifelse( Y > 0, 1, 0 ) 
# # also allow for more sites and taxa with presence-absence data by allowing for video data
# traits <- read_csv( '../Data/Bitemap_REQUEST_trait_siene+video_EAT_SQUID.csv' )
# famlong <- traits %>% 
#   dplyr::group_by( Country, family ) %>%
#   dplyr::summarize( n=sum(totalCount) ) %>%
#   dplyr::mutate( n=ifelse(n>0,1,0) )
# famwide <- famlong %>%
#   dplyr::group_by( Country ) %>%
#   spread( key=family,value=n, fill=0 ) %>%
#   select( -'<NA>', -Aplysiidae) %>%
#   filter( Country != "USA (DE)")
# # merge with sites
# sites.site <- sites %>%
#   select(Country,Site,sstmean)
# famsite <- right_join( sites.site, famwide )
# colSums(famsite[,-c(1,2, 3)])
# Y2 <- as.matrix(famsite[,-c(1,2, 3)])
# rownames(Y2) <- famsite$Site
# # merge with rates
# rate2 <- left_join(left_join(sites.site, famsite[,c(1,2,3)]), rate.mean2)
# fam.meta2 <- left_join(famsite[,c(1,2)], rate2)
# rownames(fam.meta2) <- fam.meta2$Site
# 

# Y should include sites where we have video data when presence absence is being used




## PCA on consumer densities
p  <- prcomp( Y,scale. = T )
summary(p)
biplot(p)


# RDA, CAP
# abundance
r  <- rda( Y~rate,fam.meta, scale=T )
r2 <- capscale( Y~rate,fam.meta, distance="bray", scale=T )
# presence-absence
j  <- rda( Y2~rate,fam.meta, scale=T )
j2 <- capscale( Y2~scale(rate),fam.meta, distance="raup", scale=T )
# r  <- rda( y, scale=T )
# summaries
custum <- function(model){
  list( coef(model),
        RsquareAdj(model),
        summary(model)$cont$importance[,1:6] )
}
custum(j2)

# plots
plot(r, scaling=3)
plot(r2, scaling=3)
plot(j, scaling=3)
plot(j2, scaling=3)
par(mar=c(5,4,2,2))
plot(j2,choices=c(1,2), scaling=2)
# extract axes
nax <- 1:2
s1 <- scores(r,choices = nax)$sites
s2 <- scores(r2,choices = nax)$sites
s3 <- scores(j,choices = nax)$sites
s4 <- scores(j2,choices = nax)$sites
# colnames(s2)[1] <- paste( colnames(s2)[1],"2",sep="_" )
# extract taxa
t1 <- scores(r,choices = nax)$species
t2 <- scores(r2,choices = nax)$species
t3 <- scores(j,choices = nax)$species
t4 <- scores(j2,choices = nax)$species


## Merge and save to disk
# make site a column in each site scores matrix
s1.df <- data.frame( site=row.names(s1), RDA1.s1=s1[,1], PC1.s1=s1[,2] )
s2.df <- data.frame( site=row.names(s2), CAP1.s2=s2[,1], MDS1.s2=s2[,2] )
s3.df <- data.frame( site=row.names(s3), RDA1.s3=s3[,1], PC1.s3=s3[,2] )
s4.df <- data.frame( site=row.names(s4), CAP1.s4=s4[,1], MDS1.s4=s4[,2] )
s.df <- full_join( full_join( full_join( s1.df, s2.df ), s3.df), s4.df ) 
# make taxon a column in each species loadings matrix
t1.df <- data.frame( taxon=row.names(t1), RDA1.t1=t1[,1], PC1.t1=t1[,2] )
t2.df <- data.frame( taxon=row.names(t2), CAP1.t2=t2[,1], MDS1.t2=t2[,2] )
t3.df <- data.frame( taxon=row.names(t3), RDA1.t3=t3[,1], PC1.t3=t3[,2] )
t4.df <- data.frame( taxon=row.names(t4), CAP1.t4=t4[,1], MDS1.t4=t4[,2] )
t.df <- full_join( full_join( full_join( t1.df, t2.df ), t3.df), t4.df ) 

write_csv( s.df, "Output Data/multivar_compare_sites.csv" )
write_csv( t.df, "Output Data/multivar_compare_taxa.csv" )

# also write to disk the vector for the constrained axes
vec.dir <- data.frame( model=c("r","r2","j","j2"), direction= c( scores(r, display="bp")[1], scores(r2, display="bp")[1],
   scores(j, display="bp")[1], scores(j2, display="bp")[1] ) )
write_csv( vec.dir, "Output Data/multivar_compare_direction.csv" )


# join with rate.mean
sr <- data.frame( fam.meta, s4 )
# sr$rate2 <- round(sr$rate*25)
# sr$rate3 <- cbind(sr$rate2,25-sr$rate2)
sr$CAP1 <- sr$CAP1 * vec.dir[4,2]

ggplot( data=sr, aes(x=CAP1,y=rate) ) + geom_point(size=5) +
  geom_smooth( aes(group=1), method='glm', method.args=list(family='quasibinomial')) 


# customize an RDA figure
# points make sure rate points to the right
capraup <- s4
capraup[,1] <- capraup[,1] * vec.dir[4,2]
capspec <- t4
capspec[,1] <- capspec[,1] * vec.dir[4,2]
# rates
capraup <- data.frame(capraup,fam.meta)
coef(j2) # Function coef will give the regression coefficients from centred environmental variables (constraints and conditions) to linear combination scores. The coefficients are for unstandardized environmental variables. The coefficients will be NA for aliased effects.
# proportion explained
RsquareAdj(j2)  # explains 20-25% of variation in consumption rate?
R2 <- eigenvals(j2)/sum(eigenvals(j2))
R2
summary(j2)
cummr2 <- R2
for(i in 2:length(eigenvals(j2))){
  cummr2[i] <- cummr2[i]+cummr2[i-1]  
}

asite <- ggplot( capraup, aes(x=CAP1,y=MDS1,col=habitat)) + 
  geom_hline(yintercept = 0, col='gray', lty=2 ) +
  geom_vline(xintercept = 0, col='gray', lty=2 ) +
  geom_point(aes(size=rate),  alpha=0.5) +
  xlab(paste0('CAP1 (',round(R2[1],3)*100, '%)')) +
  ylab(paste0('MDS1 (',round(R2[2],3)*100, '%)')) +
  scale_color_manual(values=c("green","gray25"))#+ 
  # geom_text_repel(aes(label=Site), point.padding = 0.5)
# add species
capspec <- data.frame( capspec, family=rownames(capspec) )
capspec2 <- capspec[ abs(capspec$CAP1) > sd(capspec$CAP1)*2 , ]  #abs(capspec$MDS1) > sd(capspec$MDS1)*2

# custom
capspec2$yadj <- c(0.1, 0.1, 0.05, 0, -0.025, -0.025, 0, 0 )
capspec2 <- capspec2 %>% 
  mutate( ynew = MDS1+yadj )

aspec <- asite + 
  # geom_text_repel( data=capspec2, aes(x=CAP1,y=MDS1,label=family), col='slateblue',
  #                  box.padding = 0.1 ) +
  geom_text( data=capspec2, aes(x=CAP1,y=ynew,label=family), col='slateblue',
             hjust = c(0,0,0,1,0,0,0,0), nudge_x = c(0.1,0.1,0.1,-0.1,0.1,0.1,0.1,0.1 ) ) +
  geom_segment( data=capspec2, aes(x=0,y=0,xend=CAP1,yend=MDS1), col='slateblue',
                arrow = arrow(length = unit(0.2, "cm")) ) 
b <- aspec + theme_minimal() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) 

# something like a Carangid might go after squid in the water column, but maybe not for small fish individuals
# Lutjanids are also quite active but they might only go after live/lively prey


## which families are responsible for high rates?  
# just get all taxa that fall out in positive association with the rda axes
# first make constrained RDA axes positive
t.rda <- t.df %>%
  select( RDA1.t1, CAP1.t2, RDA1.t3, CAP1.t4 )
vec.dir$direction
t.pos <- mapply( function(x,y) x*y, t.rda, vec.dir$direction)
rownames(t.pos) <- t.df$taxon
pos.fam <- lapply( apply( t.pos, 2, function(z) which(z>0) ), names )  #sd(z,na.rm=T)
pos.fam[[1]] == pos.fam[[2]]
pos.fam[[3]] == pos.fam[[4]]
pos.fam[[1]][!(pos.fam[[1]] %in% pos.fam[[4]])]
pos.fam[[1]] %in% pos.fam[[4]]

pf.un <- sort( unique( unlist( pos.fam ) ) )

# write to disk
write_csv( data.frame(taxon=pf.un), "Output Data/multivar_families_positive.csv" )




# another option to highlight the very strongest associated taxa
x <- colnames(s4)[1]
y <- "MDS1"
sort( t2[ t2[,x] > 0.09, x ] )
tu <- t2[ t2[,x] < -0.01, ]
plot( formula(paste(y,x,sep="~")), sr, 
      cex=rate*12+1, pch=21 )
scale = 2
text( tu[,x]*scale, tu[,y]*scale, labels=row.names(tu), col='blue' )
text( formula(paste(y,x,sep="~")), data=sr, labels=sr$Site, col='black' )

# which taxa to keep - will use this as an explanatory variable in models of consumption rate
fam.keep <- row.names(tu)
# also use RDA1, MDS axes, best or top PC axes (they are orthogonal), density of all taxa, 
# density of taxa selected in RDA, count of selected taxa (family richness), 
# then look at composition of sites that are not well-predicted by RDA



# make it long, then ggplot
sl <- data.frame(s1) %>%
  gather( key='axis', value='x')
sl$rate <- rep( fam.meta$rate, 17 )

ggplot( sl[sl$axis %in% paste0("PC",1:7),], aes(x=x, y=rate)) + facet_wrap(~axis) + geom_point() +
  geom_smooth(method='glm',se=T, method.args=list(family='quasibinomial'))








### for an analysis of all taxa, which predictor variables best explain dissimilarity?
e  <- rda( Y~sstmean,fam.meta, scale=T, na.action = na.omit )
e2 <- capscale( Y~rate,fam.meta, distance="bray", scale=T )
# presence-absence
h <- capscale( Y2~1,fam.meta, 
                distance="raup", scale=T, na.action = na.omit )
h2 <- capscale( Y2~sstmean+basin+hemi+coast+habitat,fam.meta, 
                distance="raup", scale=T, na.action = na.omit )
h3 <- capscale( Y2~chlomean+temp+sal,fam.meta, 
                distance="raup", scale=T, na.action = na.omit )
h4 <- capscale( Y2~scale(sstmean),fam.meta, 
                distance="raup", scale=T, na.action = na.omit )
h5 <- capscale( Y2~scale(temp),fam.meta, 
                distance="raup", scale=T, na.action = na.omit )
h6 <- capscale( Y2~scale(rate)+scale(temp),fam.meta, 
                distance="raup", scale=T, na.action = na.omit )
custum(h3)
custum(h4)
custum(h5)
custum(h6)
plot(h5)
plot(h6)
# 

# get scores
model <- h
ts <- scores(model,choices = nax)$sites
# colnames(s2)[1] <- paste( colnames(s2)[1],"2",sep="_" )
# extract taxa
tt <- scores(model,choices = nax)$species
tv <- scores(model, display="bp")[1]
# rsware
R2 <- eigenvals(model)/sum(eigenvals(h5))
R2

# write unconstrained ordination to disk
tsdf <- data.frame( SH=fam.meta$SH, ts )
write_csv( tsdf, "Output Data/multivar_unconstr_sites.csv" )
  
# make a figure
ts[,1] <- ts[,1] * tv
tt[,1] <- tt[,1] * tv

# rates
capunc <- data.frame(ts,fam.meta)
a <- ggplot( capunc, aes(x=MDS1,y=MDS2,col=temp)) + 
  geom_hline(yintercept = 0, col='gray', lty=2 ) +
  geom_vline(xintercept = 0, col='gray', lty=2 ) +
  geom_point(aes(size=rate),  alpha=0.75) +
  xlab(paste0('MDS1 (',round(R2[1],3)*100, '%)')) +
  ylab(paste0('MDS2 (',round(R2[2],3)*100, '%)')) +
  scale_color_viridis(option = "C") +
  theme_minimal() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) 

plot_grid( a, b, ncol=2, labels = "AUTO", align = 'hv' )


captemp <- data.frame(ts,fam.meta[ !is.na(fam.meta$temp), ])
ggplot( captemp, aes(x=CAP1,y=MDS1,col=habitat)) + 
  geom_hline(yintercept = 0, col='gray', lty=2 ) +
  geom_vline(xintercept = 0, col='gray', lty=2 ) +
  geom_point(aes(size=rate),  alpha=0.5) +
  xlab(paste0('CAP1 (',round(R2[1],3)*100, '%)')) +
  ylab(paste0('MDS1 (',round(R2[2],3)*100, '%)')) +
  scale_color_manual(values=c("green","gray25")) +
  theme_minimal() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) 
  


ggplot( data=capunc, aes(x=MDS1,y=rate) ) + geom_point(size=5) +
  geom_smooth( aes(group=1), method='glm', method.args=list(family='quasibinomial')) 

