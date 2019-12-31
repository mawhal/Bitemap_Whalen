#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Preliminary data analysis and exploration
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2018.03.05
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################

## UPDATES
# 2017.09.21: incorporate fish biomass coefficients to calculate mass from length
# 2017.09.26: finished first draft of code to compile length-weight relationships and calculate biomass
# 2017.09.27: construct fish size distribution figures
# 2017.10.23: incorporate environmental data from the dates of squidpop assays
# 2017.10.23: quantitative comparison of 1hr vs 24hr predation data
# 2017.11.07: Finish collating fish trophic and trait information  (see "Bitemap2016_Fish_Clean+Calculations.r")
# 2017.12.07: Added analysis of replicate level data for predation intentsity (squidpops)
# 2017.12.10: More in situ environmental data added (see directenv)
# 2017.12.14: AIC comparison of models with in situ and satellite temperature using same dataset 
# 2018.02.03: add a new site (recalculate satellite data); new predation rate calculation that accounts for two time points (1+24hours)
# 2018.02.14: model squidpops as zeros and ones rather than proportions
# 2018.02.16: incorporate seagrass bed characteristics into analysis
# 2018.02.18: model number of squidpops remaining as exponential decay
# 2018.03.06: use exponential rates only, remove old code
# 2018.03.06: added satellite data for Italy on days assays conducted (use for in situ temperature)

###   REMOVE NSW1.1 AND 1.2 AND USA(DE) BECAUSE THESE WERE NOT CONDUCTED IN SEAGRASS

##################################
## Quetions, hypotheses, theory ##
##################################

# Big questions:
# How does predation in marine systems vary on a global scale?
# What controls (or at least is correlated with) predation rates across the world ocean?
# Does the presence (and quality??) of structured habitat influence the predator community and predation rates?

# Hypotheses for patterns:
# stronger predation at lower latitudes because of metabolism, patterns and history of diversification, ...
    # perhaps expectation is a normal-ish distribution of predation intensity across latitude, regardless of hemisphere?

# Latitudinal gradient in fish diversity known, but we expect to see this, too. 
# Does latitudinal diversity gradient hold at size classes and relevant to mesograzer control in seagrass communities?
# Similarly, is there a latitudinal gradient when traits are considered (e.g. only consider mesopredators)?

# Hypotheses for controls:
# environment: biogeography (includes evolutionary history, climate history), modern climate, biodiversity

# Interactions between predation rate and fish diversity. Possible that these are correlated, but can we know
# whether this potential relationship is causal or that both are responses to a different, shared driver?

########################################################################


###############################
## Analysis goals/strategies ##
###############################

# consider two datasets, one for predation assay (squidpops) and other for fish community observations (seines)
# some analysis on each independently becuase conducted at different scales, times, reps
# obtain summaries for each habitat type and place and merge the two datasets together 
  # look for statistical relationships between predation and fish community (and environment)

## Environmental data: for each site (mean of all latitudes/longitutde from each site)
# Source: Bio-ORACLE, 
# temperature: mean, max, min, range
# salinity, precipitation
# nutrients, productivity, sunlight

## Seine Data
# calculate different summaries of fish community that can be related to predation intensity
# summaries: total abundance, functional group abundance, richness, diversity metrics
#            fish biomass (need to calculate based on size and taxonomy)




#############################################################
#################################################
#######################################
###############################
###############################
#######################################
#################################################
#############################################################
########################################################################
#############################################################
#################################################
#######################################
###############################



# load libraries
# data.frame manipulation and joining
library(tidyverse)
library(plyr)
library(reshape2)
# geospatial data
library(raster) # note that select() also occurs in dplyr
# statistics
library(bbmle)
library(lme4)
library(psych)


# define function vif
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


###########################################
## DATA ON PREDATION AND FISH COMMUNITIES #
###########################################

# read in data
# SQUID POP PREDATION ASSAYS
pops   <- read.csv( '../Data/OceanBitemap_Squidpop_Data_20180203.csv', stringsAsFactors = FALSE )
# rename pops columns
names(pops) <- c("timeStamp","name","email","Lat","Long","dateRetrieved","timeRetrieved",
                 "dateDeployed","timeDeployed","habitat","habitatDescription","numberDeployed","missing1",
                 "numberRetrieved","missing24","Country","Institution","prop1","prop24","notes")
# convert vegetated and unvegetated sites to common categories
pops$habitat[ pops$habitat %in% c("Seagrass","Seagrass Meadow","seagrass")] <- "Seagrass"
pops$habitat[ pops$habitat %in% c("Muddy Bottom","Sandy Bottom", "unvegetated ","unvegetated", "Unveg")] <- "Unveg"
pops$habitat[ pops$habitat %in% c("Artificial Habitat (dock, breakwater, weir, etc.)","Rocky Reef")] <- NA
pops <- pops[ !is.na(pops$habitat), ]
with(pops, table(Country,habitat))
pops$habitat <- factor(pops$habitat)

## reduce the number of columns
pops <- pops %>%
  dplyr::select( Country, Lat, Long, dateDeployed, habitat, N1=numberDeployed, N24=numberRetrieved, prop1, prop24 )


#####---------------------------------------------------------------------------
##### Consider 1 hour and 24 hour sampling points together in the same analysis

## melt and recast data so that prop1 and prop24 are a single column
popmelt <- melt( pops, id.vars = 1:7, measure.vars = 8:9 )

dcast( popmelt, Country~variable )

pop <- popmelt %>%
  dplyr::select( Country, Lat, Long, dateDeployed, habitat, N1,N24, hour=variable, prop=value ) %>%
  mutate( hour = as.numeric(as.character(dplyr::recode(hour, "prop1"=1, "prop24"=24))) )

# add a count column and fill based on whether hour is 1 or 24
pop$N <- ifelse( pop$hour==1, pop$N1, pop$N2 )

# get rid of extra columns for total numbers of squidpops
pop <- pop %>% dplyr::select( Country, Lat, Long, dateDeployed, habitat, hour, N, prop )

# add columns for successes (eaten) and failures (not eaten)
pop$eaten    <- round( with(pop, N*prop) )
pop$noteaten <- with(pop, N-eaten)


## Seagrass bed characteristics
seagrass <- read.csv( "../Data/Environmental Data/BiteMap_Seagrass_bed_character.csv", 
                      stringsAsFactors = FALSE, strip.white = TRUE )
# get all species
species <- strsplit( seagrass$seagrassSpecies, ", ", fixed=T )
seagrass$richness <- unlist(lapply( species, length ))
splist <- sort(unique(do.call(c,species)))
# species as rows
spmat <- matrix(NA,ncol=length(splist),nrow=nrow(seagrass))
for( i in 1:nrow(seagrass) ){
  for(j in 1:length(splist)){
    spmat[i,j] <- ifelse( splist[j] %in% species[[i]], 1, 0)
  }
}
splist <- gsub(pattern = " ",".",splist)
colnames(spmat) <- splist
colSums(spmat)
rowSums(spmat)
sea <- cbind( seagrass, spmat )
# # seagrass species occurence with latitude
# windows()
# par(mfrow=c(4,4))
# for( i in 14:29 ){
#   plot( sea[,i] ~ sea$latDecimal )
# }



## FISH SEINING DATA
# seines <- read.csv( '../Data/Bitemap_Seine_ALL-DATA_20180122.csv', strip.white = TRUE)
# names(seines)[4] <- "habitat"
# # convert all vegetated and unvegetated sites to common categories
# seines$habitat[ seines$habitat %in% c("Seagrass " )] <- "Seagrass"
# seines$habitat[ seines$habitat %in% c("unveg","Unveg","Unvegetated" )] <- "Unveg"
# seines <- droplevels(seines)
# # for all organisms from France, UNC2, Wales, QLD2, QLD3, multiple lengths by 10 to convert from cm to mm
# seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] <- 
#   seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] * 10


# # merging data
# # this is difficult because the replicate seines do not pair up with squidpop reps
# # so, in order to merge the two we need to use summaries at the site level
# # try to match things based on Lat/Long info
# popGPS   <- unique( pops[,c('Lat','Long','Institution','Country')])
# seineGPS <- unique( seines[,c('Lat','Long','Country','Site.Name','habitat')])
# # compare sites used in seining and predation assays
# GPSjoin <- full_join( popGPS,seineGPS )
# GPSjoin <- GPSjoin[ with(GPSjoin, order(Lat,Long)), ]
# # write.csv( GPSjoin, "Output Data/pops_sienes_GPSmatch.csv", row.names=F )
# # to complete merge, we need to average predation and fish metrics at site level
# 
# # combine all squidpop and seine replicates from a given location (e.g. NSW2)
# siteGPS <- ddply( GPSjoin, .(Country),
#                   summarise, meanLat=mean(Lat), meanLong=mean(Long) )
# # for each site within a country
# siteGPS.site <- ddply( GPSjoin, .(Country,Site.Name),
#                   summarise, meanLat=mean(Lat), meanLong=mean(Long) )
# write.csv( siteGPS, "../Data/Bitemap_sites.csv", row.names=FALSE )
siteGPS <- read.csv( "../Data/Bitemap_sites.csv" )

# quick map
library(ggalt)
library(ggmap)
library(rgdal)
world <- fortify(map_data("world"))

WorldData <- map_data('world')
WorldData %>% filter(region != "Antarctica") -> WorldData
WorldData <- fortify(WorldData)

p <- ggplot()
p <- p + geom_map(data=WorldData, map=WorldData,
                  aes(x=long, y=lat, group=group, map_id=region),
                  fill="white", colour="#7f7f7f", size=0.5)
p + coord_proj("+proj=wintri")


mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mp <- ggplot() +   mapWorld

gg <- ggplot( world, aes(map_id=region) ) +
  geom_map( map=world,
                    ,
                    color="black", fill="white", size=0.25)
gg <- gg + geom_point(data=(sites), 
                      aes(x=coords.x1, y=coords.x2), 
                      color="#b2182b")
gg <- gg + coord_proj("+proj=wintri")
gg <- gg + theme_map()
gg





###############################
#######################################
#################################################
#############################################################
########################################################################
#############################################################
#################################################
#######################################
###############################


###################################################################################
# EXTRACT BIO-ORACLE DATA                                                         #
###################################################################################

# This part of the script prepares Bio-ORACLE and WorldClim data for input sites,
# then merges these data with data from the input
# http://www.oracle.ugent.be/download.html (Shift-click to follow link)
# downloaded .rar data were extracted in Ubuntu using unrar function
# Note that BioOracle offers many useful oceanographic predictors

# # function to extract raster data within a radius defined by buffer
# buff <- function(rast,buffer=10000) {
#   unlist(lapply(extract( rast, input, buffer=buffer ),mean,na.rm=T))
# }


########### Bio-ORACLE ##################
# get a list of files to work with
# files <- list.files( "../Data/BioOracle Data/",pattern = ".asc")
# # read in the raster data files (will assume lat,long for projection)
# r <- lapply( paste0("../Data/BioOracle Data/",files), raster )
# # crop all rasters to same extent, keeping southern hemisphere
# e <- extent(-180,180,-70,70)
# r2 <- lapply( r, function(rast) crop(rast,e) )
# 
# ########### WorldClim Precipitation ###############
# # p <- raster( "C:/Users/mawha/Dropbox/Global Databases/WorldClim/WorldClim_precip_2-5.tif" )
# 
# ################## input data ######################
# # extract lat and long from input
# input <- siteGPS[,c("meanLong","meanLat")]
# 
# 
# ############ extract environmental data #############
# # average values in all raster cells within a given radius from the GPS pointS
# # Use a "buffer" or radius over which to look for raster cells surrounding each GPS point
# # that has data (note there is no data on land or freshwater for Bio-ORACLE, and no ocean data for WorldClim
# # This will take a while for many variables
# buffer <- 10000  # this is in meters if the map is projected correctly
# # if raster coordinate reference system (CRS) is undefined, it will assume lat/long, which is correct in this case
# 
# # use the buff function (defined above) to all rasters
# oracle <- lapply( r2, buff ) # this takes a very long time, average Lat and Long by site first
# # precip <- buff( p )
# 
# # combine all of these into a data.frame and give them names
# Environmentals <- data.frame( do.call( cbind, oracle ) )
# names(Environmentals) <- unlist( strsplit( files, ".asc") )
# # Environmentals$precip <- precip
# 
# sites <- cbind(siteGPS,Environmentals)
# # factor for hemisphere
# sites$hemi <- factor(ifelse( sites$meanLat>0,"North","South"))
# 
# # write Environmentals to disk
# write.csv( sites, "Output Data/Bitemap_BioORACLE.csv", row.names = FALSE )
sites <- read.csv( "Output Data/Bitemap_BioORACLE_2018.04.02.csv", stringsAsFactors = FALSE )

# add ocean basin information here
sites$Country
sites$basin <- c( 1,1,1,1,1,1,1,1,2,2,2,2,1,3,3,4,2,3,1,1,2,2,2,2,1,1,1,1,2,2,2,2,2,1,2,2,2,1,1,2)
sites$basin <- factor( sites$basin, levels=1:4, labels=c("Pacific","Atlantic", "Mediterranean", "Indian") )

# MISSING VALUES FROM Bio-ORACLE
unique( sites$Country )
sites$Country[ is.na(sites$sstmean) ]
# no data for
# New South Wales, NSW1.1
# University of Delaware





par(mfrow=c(2,2), mar=c(5,4,1,2)+0.1 )
with(sites, plot( sstrange~sstmean) )
with(sites, plot( sstrange~sstmin) )
with(sites, plot( sstrange~sstmax) )
with(sites, plot( sstrange~parmean) )
dev.off()
with(sites, plot( sstmean~parmean) )

m1 <- lm(sstrange~sstmean,sites)
m2 <- lm(sstrange~sstmax,sites)
m3 <- lm(sstrange~sstmin,sites)
AICctab( m1,m2,m3,nobs=nrow(na.omit(sites)) )
# temperature range is much better predicted by minimum temperature that max or mean



par(mar=c(5,4,2,2)+0.1)
plot( sstmin~sstmean, data=sites, type='n', ylab="Min or max SST", xlab="Mean annual SST",
      ylim=c(0,35))
  points( sstmax~sstmean, data=sites, col='darkorange' )
  points( sstmin~sstmean, data=sites, col='blue' )
  abline( lm(sstmax~sstmean, data=sites), col='darkorange' )
  abline( lm(sstmin~sstmean, data=sites), col='blue' )
  
ggplot( data=sites, aes(x=abs(meanLat),y=sstmean)) + facet_wrap(~hemi) + geom_point()
ggplot( data=sites, aes(x=abs(meanLat),y=sstmean,col=hemi)) + geom_point() + geom_smooth(se=F)
ggplot( data=sites, aes(x=sstmean,y=sstrange)) + facet_wrap(~hemi) + geom_point()


# read temperature and salinity data collected during predation assays
directenv <- read.csv( '../Data/Environmental Data/BiteMap_EnvironmentalData_Compilation_20180306.csv' )
names(directenv)[names(directenv)=="Seagrass.Unveg"] <- "habitat"
directenv$habitat[ directenv$habitat=="seagrass"] <- "Seagrass"
directenv$habitat[ directenv$habitat=="unveg"] <- "Unveg"
directenv <- droplevels(directenv)

# get average temp and salinity for each site on each day
directenvmean <- ddply( directenv, .(Site.Name,Lat,Long,habitat,Date), summarize,
                        temp=mean(Temp), sal=mean(Sal) )

# get average temp and salinity for each Country
directenvsite <- ddply( directenv, .(Country,habitat), summarize,
                        temp=mean(Temp), sal=mean(Sal) )

# combine satellite and in situ measurements
popenv <- left_join( pops, directenvsite )
allenv <- left_join( popenv, sites )



### ----------------------------------------------------------------------------------
### Figures of environmental data (latitude, temperature, etc.)

ggplot( allenv, aes(x=meanLat,y=sstmax) ) + geom_point(size=2.5) + 
  # geom_smooth( method='lm', formula = y~x^2, se=T ) +
  geom_smooth( aes(group=1), col='black', se=T ) +
  xlab("Degrees latitude") + 
  ylab(expression(paste("Maximum annual sea surface temperature (",degree,"C)",sep="")))

# latitude on y axis for nice map comparison
windows(1.75,3)
ggplot( allenv, aes(x=meanLat,y=sstmean) ) +  
  geom_smooth( aes(group=1),col='black', se=F, lwd=0.5 ) +
  geom_point(col="slateblue", alpha=0.8, size=2) +
  xlab("Degrees latitude") + 
  ylab("Mean annual\nSST (°C)") +
  theme_classic() +
  # Flip axes to let match up with a map
  coord_flip()


# range
ggplot( allenv, aes(x=abs(meanLat),y=sstrange,fill=hemi) ) + geom_point(size=2,pch=21) + 
  # geom_smooth( method='lm', formula = y~x^2, se=T ) +
  xlab("Degrees from equator") + 
  ylab(expression(paste("Range annual sea surface temperature (",degree,"C)",sep=""))) +
  scale_fill_manual(values=c("black","white"))

ggplot( allenv, aes(x=sstmean,y=sstrange,col=basin) ) + geom_point(size=3) + 
  geom_line(alpha=0.5,size=2)


# instantaneous temperature as a function of SST
windows(3,3)
ggplot( allenv, aes(x=sstmean,y=temp) ) + geom_smooth(se=F,col='black',lwd=0.5) + 
  geom_point(size=2,col='slateblue') + 
  ylab("Instantaneous water\ntemperature (°C)") + xlab("Mean annual SST (°C)") + theme_classic()
ggplot( rate.env, aes(x=sstmean,y=temp) ) + geom_point() +geom_smooth()

###







####-------------------------------------------------------------------------
## Model predation for each site as exponential decay of remaining pops (noteaten)
#

# add zeros
zeros <- pop[pop$hour==1,]
zeros$hour  <- 0
zeros$eaten <- 0
zeros$prop  <- 0
zeros$noteaten <- zeros$N

pop0 <- rbind( pop,zeros,zeros,zeros )

# make factor for Country, date, and habitat
pop0$ID <- with( pop0, paste( Country,dateDeployed,habitat, sep="." ))

windows(12,7)
ggplot( popall0, aes(x=hour,y=noteaten,col=habitat, group=ID ) )+ facet_wrap(~Country,ncol=8) + geom_point() +
  geom_smooth( method='glm', method.args=list(family=poisson), se=F, lwd=0.5 ) +
  xlab("Hour of predation assay") + ylab("Number of prey remaining")
dev.off()

### Get estimates for individual Poisson models
mods <- dlply( pop0, .(ID), glm, formula=noteaten~hour, family=poisson )
coefs <- ldply( mods, coef )

## merge these back with environment data
popmerge <- pop0[,c(1,4,5,11)]
rates <- left_join( coefs, popmerge, by="ID" )
# for some reason I get a lot of duplicates this way
rates <- rates[ !duplicated( rates ), ] 
# merge environment
rate.env  <- left_join( rates, allenv )

# convert rate to 1-hour so that rate makes more sense (bigger equals faster)
rate.env$rate <- 1-exp(rate.env$hour)
##


# Add in fish data
fish <- read.csv("Bitemap_predation+fish+env_summary.csv", stringsAsFactors = FALSE)
fish <- fish %>%
  dplyr::select( Country, habitat, ENSPIE,ENSPIEdiff, Total, biomass )
rate.fish <- left_join(rate.env,fish)

## For modeling,
# restrict dataset to rows with both in situ and satellite measurements
rateclean <- rate.env[ !is.na(rate.env$temp) & !is.na(rate.env$sstmean), ]
rateclean$abLat <- abs(rateclean$Lat)
# which sites do not have in situ or satellite measurements
unique( rate.env$Country[ is.na(rate.env$temp) ] )
unique( rate.env$Country[ is.na(rate.env$sstmean) ] )
rate.env[ rate.env$Country=="Italy", ]


## Compare rates inside and outside of seagrass across all sites
# boxplot of rates by country and habitat type
windows(12,7)
ggplot( rate.env, aes(x=habitat,y=rate, col=habitat)) + facet_wrap(~Country, ncol=8) + geom_boxplot() + geom_point()
dev.off()

# Sites where rates noticeably higher in seagrass
c( "Australia (NSW2)", "Australia (QLD3)", "Mexico (BN)", "USA (CA)", "USA (CA3)", "Wales" )

# Sites where rates outside of seagrass higher
c( "Belize", "USA (FL)", "Korea" )

# Sites with high predation across the board
c( "Australia (NSW)", "Italy", "USA (NC)", "USA (NC2)" )

# Sites with low predation across the board
c( "Australia (QLD)", "Australia (VIC)", "Brazil", "Canada (BC)", "Canada (QC)", "Chile", "Croatia", 
   "France", "India", "Ireland", "Mexico (ICML)", "Mexico (ICML2)", "Norway", "Panama", "USA (AK)", 
   "USA (CA2)", "USA (MA)", "USA (OR)", "USA (VA)", "USA (VA2)", "USA (WA)", "USA (WA2)" )







#### ---------------------------------------------------------------------------------------
#### look at predation data without averaging at site level


# predation rate data
ggplot(rate.env, aes(x=meanLat, y=rate )) + geom_point() + facet_grid(hemi~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial))  

ggplot(rate.env, aes(x=meanLat, y=rate, group=hemi )) + geom_point() + facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial))  

ggplot(rate.env, aes(x=abs(meanLat), y=rate, group=hemi )) + geom_point() + 
  geom_smooth(method='glm',method.args=list(family=quasibinomial))  

ggplot(rate.env, aes(x=abs(meanLat), y=rate, group=hemi )) + geom_point() + 
  geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2))  

ggplot(rate.env, aes(x=abs(meanLat), y=rate )) + geom_point() + 
  geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2))  

# evidence for peak at warm sites, but not the warmest sites
# peak is more extreme in southern hemisphere
ggplot(rate.env, aes(x=sstmean, y=rate, group=hemi )) + geom_point() + 
  geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2))  
ggplot(rate.env, aes(x=sstmean, y=rate, group=habitat )) + geom_point() + facet_wrap(~hemi) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2))

ggplot(rate.env, aes(x=sstrange, y=rate )) + geom_point() +geom_smooth()

ggplot(rate.env, aes(x=temp, y=rate )) + geom_point() + facet_wrap(~hemi) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial))

# no difference between habitat types at this scale of analysis
ggplot(rate.env, aes(x=temp, y=rate, group=habitat )) + geom_point() + 
  geom_smooth(method='glm',method.args=list(family=quasibinomial))
ggplot(rate.env, aes(x=temp, y=rate, group=habitat )) + geom_point() + facet_wrap(~hemi) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial))


# temperature range
windows(4,3)
ggplot( rate.env, aes(x=abs(meanLat),y=sstrange,size=rate) ) + geom_point(pch=21) + 
   # geom_smooth( data=rate.env[rate.env$rate<0.1 & abs(rate.env$Lat) < 40,], aes(group=1), method='lm', se=F ) +
   # geom_smooth( data=rate.env[rate.env$rate>0.5 & abs(rate.env$Lat) < 40,], aes(group=1), method='lm', se=F ) +
  xlab("Degrees from equator") + 
  ylab( "Range annual sea surface\ntemperature (°C)" ) +
  theme_classic()


# use all fish data available from seines to look at effect on predation
ggplot(rate.fish, aes(x=ENSPIE, y=rate, group=habitat )) + geom_point() + facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic() +
  xlab("Fish diversity (ENSPIE)") + ylab("Predation rate")
ggplot(rate.fish, aes(x=Total, y=rate, group=habitat )) + geom_point() + facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial))
ggplot(rate.fish, aes(x=biomass, y=rate )) + geom_point() +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()





###--------------------------------------------------------------------------------------------------
### MODELS

## model temp replicates using mixed effects model
tempm1 <- glmer( rate~temp + (1|Country), rateclean, family='binomial' )
summary(tempm1)
plot( y=resid(tempm1), x=abs(rateclean$Lat[!is.na(rateclean$temp)]) )

tempm2 <- glmer( rate~temp+habitat + (1|Country), rateclean, family='binomial' )
summary(tempm2)
vif.mer(tempm2)

tempm3 <- glmer( rate~poly(temp,2)*habitat + (1|Country), rateclean, family='binomial' )
summary(tempm3)
vif.mer(tempm3)

tempm4 <- glmer( rate~poly(temp,2) + (1|Country), rateclean, family='binomial' )
summary(tempm4)
vif.mer(tempm4)

tempm5 <- glmer( rate~temp*habitat + (1|Country), rateclean, family='binomial' )
summary(tempm5)
vif.mer(tempm5)

# model comparison
AICctab( tempm1,tempm2,tempm5,tempm3, nobs=nrow(rateclean[!is.na(rateclean$temp),]) )



## responses to annual temperature


ggplot(rateclean, aes(x=sstmean, y=rate )) + geom_point() + facet_wrap(~habitat,ncol=2) +
  geom_smooth(method='glm',  method.args=list(family=quasibinomial) ) 


## model sst replicates using mixed effects model
sstm2 <- glmer( rate~sstmean+habitat + (1|Country), rateclean, family='binomial' )
summary(sstm2)
summary( glm( rate~temp+habitat, rateclean, family='binomial' ) )

sstm3 <- glmer( rate~poly(sstmean,2)*habitat + (1|Country), rateclean, family='binomial' )
summary(sstm3)
vif.mer(sstm3)

sstm4 <- glmer( rate~sstmean + (1|Country), rateclean, family='binomial' )
summary(sstm4)
vif.mer(sstm4)

sstm5 <- glmer( rate~sstmean+sstrange + (1|Country), rateclean, family='binomial' )
summary(sstm5)
vif.mer(sstm5)

sstm6 <- glmer( rate~sstmin*sstrange + (1|Country), rateclean, family='binomial' )
summary(sstm6)
vif.mer(sstm6)

sstm7 <- glmer( rate~sstmin:sstrange + (1|Country), rateclean, family='binomial' )
summary(sstm7)
vif.mer(sstm7)

sstm8 <- glmer( rate~poly(sstmean,2) + (1|Country), rateclean, family='binomial' )
summary(sstm8)
vif.mer(sstm8)

m0 <- glmer( rate~1 + (1|Country), rateclean, family='binomial' )
summary(m0)

AICctab( m0,sstm2,sstm3,sstm4,sstm5,sstm6,sstm7,sstm8, nobs=nrow(rateclean[!is.na(rateclean$sstmean),]))
 

# effects of latitude and temp
latm1 <- glmer( rate~abLat + (1|Country), rateclean, family='binomial' )
summary(latm1)

latm2 <- glmer( rate~abLat+abLat:temp + (1|Country), rateclean, family='binomial' )
summary(latm2)
vif.mer(latm2)

latm3 <- glmer( rate~abLat:sstmean + (1|Country), rateclean, family='binomial' )
summary(latm3)
vif.mer(latm3)

latm4 <- glmer( rate~abLat+abLat:sstmin + (1|Country), rateclean, family='binomial' )
summary(latm4)
vif.mer(latm4)

latm5 <- glmer( rate~abLat+abLat:sstrange + (1|Country), rateclean, family='binomial' )
summary(latm5)
vif.mer(latm5)

latm6 <- glmer( rate~poly(abLat,2) + (1|Country), rateclean, family='binomial' )
summary(latm6)
vif.mer(latm6)

latm7 <- glmer( rate~abLat+abLat:sstmax + (1|Country), rateclean, family='binomial' )
summary(latm7)
vif.mer(latm7)

latm8 <- glmer( rate~abLat+hemi + (1|Country), rateclean, family='binomial' )
summary(latm8)
vif.mer(latm8)

# compare latitude models
AICctab( latm1,latm2,latm3,latm4,latm5,latm6,latm7,latm8, nobs=nrow(rateclean[!is.na(rateclean$temp),]))


# compare all models
AICctab( m0,latm3,latm4,sstm3,sstm8,tempm2,tempm1, nobs=nrow(rateclean) )

# compare all models WITHOUT quadratic terms
AICctab( m0,latm3,latm4,sstm2,sstm4,tempm2,tempm1, nobs=nrow(rateclean) )





###--------------------------------------------------------------------------------------------------------
### model predictions

# Annual mean sea surface temperature with quadratic AND habitat
pred.frame <- with(rateclean, expand.grid( sstmean=seq(min(sstmean,na.rm=T),max(sstmean,na.rm=T),len=100),
                                           habitat=unique(habitat)) )
X <- model.matrix(~poly(sstmean,2)*habitat,data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
pred <- data.frame(pred.frame,mean=(X%*%fixef(sstm3)))  # these are the point predictions
# Calculate variance between observations within each site (for each combination of fixed effects)
V <- vcov(sstm3)
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T

# Attach to dataframe, calculate standard errors and confidence intervals (using 1.96*sigma may be anti-conservative...)
predictions <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictions <- with(predictions, data.frame(predictions,
                                            p1.lo = mean-1.96*pred.se1,
                                            p1.hi = mean+1.96*pred.se1))
windows(5,3)
ggplot(rateclean, aes(x=sstmean, y=rate )) + geom_point(alpha=0.2,col='slateblue') + facet_wrap(~habitat,ncol=2) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(mean) )) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.lo)),lty=2) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.hi)),lty=2) +
  # geom_smooth( method='glm', method.args=list(family=binomial), se=F ) +
  ylab(expression(paste('Predation rate (',hr^-1,')'))) +
  xlab(expression(paste('Annual mean temperature (',degree,'C)'))) +
  theme_bw() + theme( panel.grid = element_blank()  )



# Latitude with sstmean
pred.frame <- with(rateclean, expand.grid( abLat=seq(min(abLat),max(abLat),len=20), 
                                          # sstmean=seq(min(sstmean),max(sstmean),len=5) ) )
                                          sstmean=c(min(sstmean),mean(sstmean),max(sstmean)) ) )
X <- model.matrix(~abLat*sstmean,data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
pred <- data.frame(pred.frame,mean=(X%*%fixef(latm3)))  # these are the point predictions
# Calculate variance between observations within each site (for each combination of fixed effects)
V <- vcov(latm3)
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T

# Attach to dataframe, calculate standard errors and confidence intervals (using 1.96*sigma may be anti-conservative...)
predictions <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictions <- with(predictions, data.frame(predictions,
                                            p1.lo = mean-1.96*pred.se1,
                                            p1.hi = mean+1.96*pred.se1))

ggplot(rateclean, aes(x=abLat, y=rate )) + geom_point(alpha=0.2,col='slateblue') +  
  geom_line( data=predictions, aes(x=abLat,y=logistic(mean), group=sstmean )) +
  # geom_line( data=predictions, aes(x=abLat,y=logistic(p1.lo),group=sstmean),lty=2) +
  # geom_line( data=predictions, aes(x=abLat,y=logistic(p1.hi),group=sstmean),lty=2) +
  ylab('Predation intensity\n(rate. squid missing after 24hr)') +
  xlab("Degrees from equator") + 
  theme_bw() + theme( panel.grid = element_blank()  )


# Latitude alone
pred.frame <- with(rateclean, expand.grid( abLat=seq(min(abLat,na.rm=T),max(abLat,na.rm=T),len=20) ))
X <- model.matrix(~abLat,data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
pred <- data.frame(pred.frame,mean=(X%*%fixef(latm1)))  # these are the point predictions
# Calculate variance between observations within each site (for each combination of fixed effects)
V <- vcov(latm1)
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T

# Attach to dataframe, calculate standard errors and confidence intervals (using 1.96*sigma may be anti-conservative...)
predictions <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictions <- with(predictions, data.frame(predictions,
                                            p1.lo = mean-1.96*pred.se1,
                                            p1.hi = mean+1.96*pred.se1))

ggplot(rateclean, aes(x=abLat, y=rate )) + geom_point(alpha=0.2,col='slateblue') + 
  geom_line( data=predictions, aes(x=abLat,y=logistic(mean))) +
  # geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.lo)),lty=2) +
  # geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.hi)),lty=2) +
  # geom_smooth( method='glm', method.args=list(family=binomial), se=F ) +
  ylab('Predation intensity\n(rate. squid missing after 24hr)') +
  xlab("Degrees from equator") + 
  theme_bw() + theme( panel.grid = element_blank()  )



# in situ temperature 
pred.frame <- with(rateclean, expand.grid( temp=seq(min(temp,na.rm=T),max(temp,na.rm=T),len=20) ) )
X <- model.matrix(~temp,data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
pred <- data.frame(pred.frame,mean=(X%*%fixef(tempm1)))  # these are the point predictions
# Calculate variance between observations within each site (for each combination of fixed effects)
V <- vcov(tempm1)
# # model predictions
# pred.frame <- with(rateclean, expand.grid( temp=seq(min(temp,na.rm=T),max(temp,na.rm=T),len=20), habitat=unique(habitat) ) )
# X <- model.matrix(~temp*habitat,data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
# pred <- data.frame(pred.frame,mean=(X%*%fixef(tempm1)))  # these are the point predictions
# # Calculate variance between observations within each site (for each combination of fixed effects)
# V <- vcov(tempm1)
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T

# Attach to dataframe, calculate standard errors and confidence intervals (using 1.96*sigma may be anti-conservative...)
predictions <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictions <- with(predictions, data.frame(predictions,
                                            p1.lo = mean-1.96*pred.se1,
                                            p1.hi = mean+1.96*pred.se1))
windows(3,3)
ggplot(rateclean, aes(x=temp, y=rate )) + geom_point(alpha=0.2,col='slateblue') + #facet_wrap(~habitat,ncol=2) +
  geom_line( data=predictions, aes(x=temp,y=logistic(mean) )) +
  geom_line( data=predictions, aes(x=temp,y=logistic(p1.lo)),lty=2) +
  geom_line( data=predictions, aes(x=temp,y=logistic(p1.hi)),lty=2) +
  # geom_smooth( method='glm', method.args=list(family=binomial), se=F ) +
  ylab(expression(paste('Predation rate (',hr^-1,')'))) +
  xlab(expression(paste('in situ temperature (',degree,'C)'))) + 
  theme_bw() + theme( panel.grid = element_blank()  )





####-----------------------------------------------------------------------------------------------------------------
## Does in situ temp explain residual variation in relationship with mean annual SST?
