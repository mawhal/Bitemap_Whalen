#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Preliminary data analysis and exploration
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2018.09.15
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
# 2018.03.12: add updated predator data summaries after filtering taxa that are unlikely to interact with squidpops
# 2018.03.22: updated GPS points for France
# 2018.04.08: Add more latitude models and figures. Test for quadratic latitude effect for contiguous ocean transects
# 2018.04.08: Add site abbreviations 
# 2018.04.18: Add species richness of warm families (added onto object sites in fish script)
# 2018.04.18: Matt realizes all of the above changes should have just been tracked with git
# 2018.08.25: Add new models that include fish community data (using family presence-absence)
# 2018.09.13: Finalize model set for AIC comparison
# 2018.09.15: Simple mediation model to test for effects of temperature and predator composition
# 2018.10.15: Matt adds Bitemap project to GitHub. Futher changes tracked there

###   REMOVE NSW1.1 AND 1.2 AND USA(DE) BECAUSE THESE WERE NOT CONDUCTED IN SEAGRASS?

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
library(ggrepel) # for plotting text on scatterplots

# statistics
library(bbmle)
library(lme4)
library(psych) # for function logistic
library(MuMIn)

# plotting
library(lattice)
library(viridis)


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


###################
## SQUIDPOP DATA ##
###################

# read in data
# SQUID POP PREDATION ASSAYS
pops   <- read.csv( '../Data/OceanBitemap_Squidpop_Data_20180313.csv', stringsAsFactors = FALSE )
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

# reduce number of columns
sea <- sea %>%
  dplyr::select( Country, siteName, contact=contactName, density=shootDensity.m2, 
                 cover=percentCover, canopy=canopyHeight, richness,
                 14:29)
head(sea)  
  
  

## PREDATOR DATA
predator <- read.csv( "Bitemap_SEINE_summaries_20180915.csv", stringsAsFactors = FALSE )



# # merging data
# # this is difficult because the replicate seines do not pair up with squidpop reps
# # so, in order to merge the two we need to use summaries at the site level
# # try to match things based on Lat/Long info
# popGPS   <- unique( pops[,c('Lat','Long','Country')])
# seines   <- read.csv( "../Data/Bitemap_Seine_ALL-DATA_20180322.csv", stringsAsFactors = FALSE )
#   names(seines)[4] <- "habitat"
#   # convert all vegetated and unvegetated sites to common categories
#   seines$habitat[ seines$habitat %in% c("Seagrass ","seagrass" )] <- "Seagrass"
#   seines$habitat[ seines$habitat %in% c("unveg","Unveg","Unvegetated" )] <- "Unveg"
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


# merge predator and squidpops
head(predator)

# standardize abundance based on seine length
predator$Total.std <- with( predator, Total/Distance )
predator$TotalFish.std <- with( predator, TotalFish/Distance )
predator$biomass.std <- with( predator, biomass/Distance )

# summarize predator data to site and habitat
pred.site <- ddply( predator, .(Country, habitat), summarize, Total=mean(Total), TotalFish=mean(TotalFish),
       Total.std=mean(Total.std), TotalFish.std=mean(TotalFish.std),
       ENSPIE=mean(ENSPIE), ENSPIEfish=mean(ENSPIEfish), 
       biomass=mean(biomass), biomass.std=mean(biomass.std) )

# 
# merge predators
pops.pred <- left_join( pops, pred.site, by=c("Country","habitat") )
# merge seagrass bed characteristics
pops.pred.hab <- left_join( pops.pred, sea, by="Country" )


##


# # Map of Sites
# # Winkel-triple projection
# library(maps)
# library(rnaturalearth) # for countries...see https://bhaskarvk.github.io/user2017.geodataviz/notebooks/02-Static-Maps.nb.html
# library(sp)
# 
# sites <- siteGPS
# coordinates(sites) <- ~meanLong+meanLat 
# proj4string(sites) <- '+init=epsg:4326'
# 
# world <- rnaturalearth::countries110
# world <- world[world$name != 'Antarctica',]
# 
# grid.lines.mj <- sp::gridlines(world,easts = seq(-180,180,by=30), norths = seq(-90,90,by=30))
# grid.lines.mi <- sp::gridlines(world,easts = seq(-165,195,by=15), norths = seq(-90,90,by=15))
# sites <- spTransform(sites, CRS("+proj=wintri"))
# world <- spTransform(world, CRS("+proj=wintri"))
# grid.lines.mj <- spTransform(grid.lines.mj,CRS("+proj=wintri"))
# grid.lines.mi <- spTransform(grid.lines.mi,CRS("+proj=wintri"))
# 
# # make the figure
# windows(12,7)
# par(mar = c(1, 1,0,0) + 0.1 )
# plot(methods::as(world, 'Spatial'), expandBB=c(0,0,0.05,0.05))
# 
# # plot gridlines
# plot(grid.lines.mi, col=grey(0.95), add=T)
# plot(grid.lines.mj, col=grey(0.9), add=T)
# # text(labels(grid.lines.mj, side=c(1,4), labelCRS = CRS("+init=epsg:4326")), 
#      # col = grey(.6), offset=0.3)
# 
# # plot land masses and points
# plot(world, add=TRUE, border="white", col=grey(0.8))   # previous border grey(0.2)
# points(sites, pch=21, bg=scales::alpha("slateblue",0.75), cex=2)
#  
# # add labels
# text(sites, labels=sites$Site, cex=1)
# 
# dev.off()

###############################
#######################################
#################################################
#############################################################
########################################################################
#############################################################
#################################################
#######################################
###############################


# ###################################################################################
# # EXTRACT BIO-ORACLE DATA                                                         #
# ###################################################################################
# 
# # This part of the script prepares Bio-ORACLE and WorldClim data for input sites,
# # then merges these data with data from the input
# # http://www.oracle.ugent.be/download.html (Shift-click to follow link)
# # downloaded .rar data were extracted in Ubuntu using unrar function
# # Note that BioOracle offers many useful oceanographic predictors
# library(raster)
# # function to extract raster data within a radius defined by buffer
# buff <- function(rast,buffer=10000) {
#   unlist(lapply(extract( rast, input, buffer=buffer ),mean,na.rm=T))
# }
# 
# ########### Bio-ORACLE ##################
# # get a list of files to work with
# files <- list.files( "../../Global Databases/BioOracle Data/",pattern = ".asc")
# # read in the raster data files (will assume lat,long for projection)
# r <- lapply( paste0("../../Global Databases/BioOracle Data/",files), raster )
# # crop all rasters to same extent, keeping southern hemisphere
# e <- extent(-180,180,-70,70)
# r2 <- lapply( r, function(rast) crop(rast,e) )
# 
# ########### WorldClim Precipitation ###############
# p <- raster( "C:/Users/mawha/Dropbox/Global Databases/WorldClim/WorldClim_precip_2-5.tif" )
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
# precip <- buff( p )
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
# write.csv( sites, "Output Data/Bitemap_BioORACLE_2018.03.13.csv", row.names = FALSE )
sites <- read.csv( "Output Data/Bitemap_BioORACLE_20180915.csv", stringsAsFactors = FALSE )
sites <- left_join( sites,siteGPS[,1:2], by=c("Country","Site") )

# add ocean basin information here
sites$Country
sites$basin <- c( 1,1,1,1,1,1,1,1,1,2,2,2,2,1,3,3,4,2,3,1,1,2,2,2,2,1,1,1,1,2,2,2,2,2,1,2,2,2,1,1,2)
sites$basin <- factor( sites$basin, levels=1:4, labels=c("Pacific","Atlantic", "Mediterranean", "Indian") )
sites$coast <- c( 1,1,1,1,1,1,
                  1,1,1,1,1,2,
                  1,2,2,2,2,2,
                  2,1,2,1,1,2,
                  1,2,2,2,2,1,
                  1,1,1,1,2,1,
                  1,1,2,2,2)
sites$coast <- factor( sites$coast, levels=1:2, labels=c("West","East") )

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
with(sites, plot( parmean ~ abs(meanLat)) )
with(sites, plot( sstmean ~ abs(meanLat)))
with(sites, plot( sstmean ~ parmean ) )
dev.off()


m1 <- lm(sstrange~sstmean,sites)
m2 <- lm(sstrange~sstmax,sites)
m3 <- lm(sstrange~sstmin,sites)
AICctab( m1,m2,m3,nobs=nrow(na.omit(sites)) )
# temperature range is much better predicted by minimum temperature that max or mean


windows(3,3)
par(mar=c(5,4,2,2)+0.1)
plot( sstmin~sstmean, data=sites, type='n', ylab="Min or max SST", xlab="Mean annual SST",
      ylim=c(0,35))
  points( sstmax~sstmean, data=sites, col='darkorange' )
  points( sstmin~sstmean, data=sites, col='blue' )
  abline( lm(sstmax~sstmean, data=sites), col='darkorange' )
  abline( lm(sstmin~sstmean, data=sites), col='blue' )
  
  



# read temperature and salinity data collected during predation assays
directenv <- read.csv( '../Data/Environmental Data/BiteMap_EnvironmentalData_Compilation_20180313.csv' )
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
popenv <- left_join( pops.pred.hab, directenvsite )
allenv <- left_join( popenv, sites )



# ### ----------------------------------------------------------------------------------
# ### Figures of environmental data (latitude, temperature, etc.)
# 
# ggplot( allenv, aes(x=meanLat,y=sstmax) ) + geom_point(size=2.5) + 
#   # geom_smooth( method='lm', formula = y~x^2, se=T ) +
#   geom_smooth( aes(group=1), col='black', se=T ) +
#   xlab("Latitude") + 
#   ylab(expression(paste("Maximum annual sea surface temperature (",degree,"C)",sep="")))
# 
# # latitude on y axis for nice map comparison
# windows(1.75,3)
# ggplot( allenv, aes(x=meanLat,y=sstmean) ) +  
#   geom_smooth( aes(group=1),col='black', se=F, lwd=0.5 ) +
#   geom_point(bg="slateblue", size=2, pch=21) +
#   xlab("Degrees latitude") + 
#   ylab("Mean annual\nSST (°C)") +
#   theme_classic() +
#   scale_x_continuous(breaks=c(-30,-15,0,15,30,45,60)) +
#   # Flip axes to let match up with a map
#   coord_flip() 
# 
# 
# # range
# ggplot( allenv, aes(x=abs(meanLat),y=sstrange,fill=hemi) ) + geom_point(size=2,pch=21) + 
#   # geom_smooth( method='lm', formula = y~x^2, se=T ) +
#   xlab("Degrees from equator") + 
#   ylab(expression(paste("Range annual sea surface temperature (",degree,"C)",sep=""))) +
#   scale_fill_manual(values=c("black","white"))
# 
# ggplot( allenv, aes(x=sstmean,y=sstrange,col=basin) ) + geom_point(size=3) + 
#   geom_line(alpha=0.5,size=2)
# 
# 
# # instantaneous temperature as a function of SST
# windows(1.75,3)
# ggplot( allenv, aes(x=sstmean,y=temp) ) + geom_smooth(se=F,col='black',lwd=0.5) + 
#   geom_point( size=2,bg='slateblue', pch=21 ) + 
#   ylab("in situ water temperature (°C)") + xlab("Mean annual\nSST (°C)") + 
#   theme_classic() +
#   scale_y_continuous(position = "right") + 
#   theme(axis.title.y.right  = element_text(angle = 90, vjust=0.5, hjust=0))
# 
# # site labels
# windows(7,5)
# allenv.mean <- ddply( allenv, .(Country,Site, sstmean), summarise, temp=mean(temp) )
# ggplot( allenv.mean, aes(x=sstmean,y=temp) ) + geom_smooth(se=F,col='black',lwd=0.5) + 
#   geom_point( size=2,bg='slateblue', pch=21 ) + 
#   geom_text_repel( aes(label=Site) ) +
#   ylab("in situ water temperature (°C)") + xlab("Mean annual SST (°C)") #+ 
#   # theme_classic() +
#   # scale_y_continuous(position = "right") + 
#   # theme(axis.title.y.right  = element_text(angle = 90, vjust=0, hjust=0.5))
# 
# 
# ###
# 
# 
# 
# 
# 
# 
# 
# ####-------------------------------------------------------------------------
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

# windows(12,7)
# ggplot( popall0, aes(x=hour,y=noteaten,col=habitat, group=ID ) )+ facet_wrap(~Country,ncol=8) + geom_point() +
#   geom_smooth( method='glm', method.args=list(family=poisson), se=F, lwd=0.5 ) +
#   xlab("Hour of predation assay") + ylab("Number of prey remaining")
# dev.off()

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
# log transform fish abundance
rate.env$logfish <- log10(rate.env$TotalFish.std*25+1)
# log transform biomass
rate.env$logbio <- log10(rate.env$biomass.std*1000+1)

# absolute value for latitude
rate.env$abLat <- abs(rate.env$Lat)

# make habitat labels pretty
rate.env <- rate.env %>% 
  mutate( habitat = factor( habitat, c("Unveg", "Seagrass"), c("Unvegetated","Seagrass") ) )

# write these data to disk
# write.csv( rate.env, "Output Data/Bitemap_DATA_analysis_20180825.csv", row.names=FALSE )

##

# 
# # Add in OLD fish data
# fish <- read.csv("Bitemap_predation+fish+env_summary.csv", stringsAsFactors = FALSE)
# fish <- fish %>%
#   dplyr::select( Country, habitat, ENSPIE,ENSPIEdiff, Total, biomass )
# rate.fish <- left_join(rate.env,fish)

## For modeling,
# restrict dataset to rows with both in situ and satellite measurements
rateclean <- rate.env[ !is.na(rate.env$temp) & !is.na(rate.env$sstmean), ]
# which sites predator data, sst, temperature, and both habitats?
sort(unique( rate.env$Site ))                              # 41
sort(unique( rate.env$Site[ !is.na(rate.env$sstmean) ] ))  # 39
sort(unique( rate.env$Site[ !is.na(rate.env$temp) ] ))     # 36
sort(unique( rate.env$Site[ !is.na(rate.env$Total) ] ))    # 30
sort(unique( rate.env$Site[ !is.na(rate.env$MDS1) ] ))     # 35 Yucatan and Florida and India and Italy

# which sites do not have in situ or satellite measurements
unique( rate.env$Country[ is.na(rate.env$temp) ] )
unique( rate.env$Country[ is.na(rate.env$sstmean) ] )
rate.env[ rate.env$Country=="Italy", ]

# write.csv( rate.env[,c("Country","logbio")], "../../../Desktop/biomass.csv", row.names=FALSE )

## Compare rates inside and outside of seagrass across all sites
# boxplot of rates by country and habitat type
# windows(12,7)
# ggplot( rate.env, aes(x=habitat,y=rate, col=habitat)) + facet_wrap(~Country, ncol=8) + geom_boxplot() + geom_point()
# dev.off()

# Sites where rates noticeably higher in seagrass
c( "Australia (NSW2)", "Australia (QLD3)", "Mexico (BN)", "USA (CA)", "USA (CA3)", "Wales" )

# Sites where rates outside of seagrass higher
c( "Australia (QLD4)", "Belize", "USA (FL)", "Korea" )

# Sites with high predation across the board
c( "Australia (NSW)", "Italy", "USA (NC)", "USA (NC2)" )

# Sites with low predation across the board
c( "Australia (QLD)", "Australia (VIC)", "Brazil", "Canada (BC)", "Canada (QC)", "Chile", "Croatia", 
   "France", "India", "Ireland", "Mexico (ICML)", "Mexico (ICML2)", "Norway", "Panama", "USA (AK)", 
   "USA (CA2)", "USA (MA)", "USA (OR)", "USA (VA)", "USA (VA2)", "USA (WA)", "USA (WA2)" )







#### ---------------------------------------------------------------------------------------
#### look at predation data without averaging at site level

# point colors
seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
unveg.color    <- "#d8b365"

# predation rate data
# windows(3,3)
ggplot( rate.env[rate.env$rate>0.005,], aes(x=meanLat, y=rate, col=habitat, fill=habitat ))  +
  # geom_smooth(se=F)  +
  geom_smooth(aes(group=1),method='glm',method.args=list(family=quasibinomial), 
              formula=y~poly(x,2), se=F, col='black', lwd=0.5, lty=2 ) +
  geom_smooth(aes(group=hemi),method='glm',method.args=list(family=quasibinomial), 
              formula=y~poly(x,2), se=F, col='black', lwd=0.5 ) +
  geom_point(pch=21,size=1.5,alpha=0.5) +
  scale_fill_manual(values=c(unveg.color,seagrass.color)) +
  scale_color_manual(values=c(unveg.color,'black')) +
  theme_classic() + 
  scale_x_continuous(breaks=c(-30,-15,0,15,30,45,60)) +
  ylab( expression(paste('Predation rate (',hr^-1,')')) ) +
  xlab( "Latitude" ) +
  theme(legend.justification=c(1,0), legend.position=c(1,0.8), 
        legend.title = element_blank(), legend.background = element_blank() ) +
  # guides(fill = guide_legend(override.aes = list(linetype = 0)),
  #        color = guide_legend(override.aes = list(linetype = 0))) +
  coord_flip() 

##

ggplot(rate.env, aes(x=abLat, y=rate, group=hemi)) + facet_grid(hemi+coast~basin) +  
  geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2) )  +
  geom_point(alpha=0.3) +
  xlab('Latitude') + ylab(expression(paste('Predation rate (',hr^-1,')')))

rate.site <- ddply( rate.env, .(Country,Site,hemi,basin,coast), summarise, rate=mean(rate), abLat=mean(abLat), sstmean=mean(sstmean) )
# windows(9,5)
ggplot(rate.site, aes(x=abLat, y=rate, group=hemi )) + facet_grid(hemi+coast~basin) +  
  geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2), se=F, lwd=0.5 )  +
  geom_text_repel( aes(label=Site), size=3 ) + geom_point(alpha=0.3) +
  xlab('Latitude') + ylab(expression(paste('Predation rate (',hr^-1,')')))
# dev.off()
rate.site$Country[ rate.site$abLat < 22 & rate.site$rate < 0.25 ]
rate.site$Country[ rate.site$rate < 0.25 ]

ggplot(rate.env, aes(x=sstmean, y=rate, group=hemi )) + facet_grid(hemi~basin) + geom_point() + 
  geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2) )  


# run separate models for Southwest Pacific and Northwest Atlantic. Test for quadratic term
# which data?
d <- rate.env[ !is.na(rate.env$sstmean),]
dna <- d[ d$hemi == "North" & d$basin == "Atlantic" & d$coast == "West", ]
  range(dna$abLat)
  ggplot(dna,aes(x=sstmean,y=rate)) + geom_point() + geom_smooth()
dsp <- d[ d$hemi == "South" & d$basin == "Pacific" & d$coast == "West", ]
  range(dsp$abLat)
ggplot(dsp,aes(x=sstmean,y=rate)) + geom_point() + geom_smooth()

m1na <- glmer( rate ~ poly(abLat,1) + (1|Country), dna, family="binomial" ) 
m2na <- glmer( rate ~ poly(abLat,2) + (1|Country), dna, family="binomial" )     
m1na <- glmer( rate ~ poly(sstmean,1) + (1|Country), dna, family="binomial" ) 
m2na <- glmer( rate ~ poly(sstmean,2) + (1|Country), dna, family="binomial" ) 
anova(m1na,m2na)
AICctab(m1na,m2na)

m1sp <- glmer( rate ~ poly(abLat,1) + (1|Country), dsp, family="binomial" ) 
m2sp <- glmer( rate ~ poly(abLat,2) + (1|Country), dsp, family="binomial" ) 
m1sp <- glmer( rate ~ poly(sstmean,1) + (1|Country), dsp, family="binomial" ) 
m2sp <- glmer( rate ~ poly(sstmean,2) + (1|Country), dsp, family="binomial" ) 
anova(m1sp,m2sp)
AICctab(m1sp,m2sp)
#

# evidence for peak at warm sites, but not the warmest sites
# peak is more extreme in southern hemisphere
# ggplot(rate.env, aes(x=sstmean, y=rate, group=hemi )) + geom_point() + 
#   geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2))  
# ggplot(rate.env, aes(x=sstmean, y=rate )) + geom_point() + facet_wrap(~hemi) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2))
# 
# ggplot(rate.env, aes(x=sstrange, y=rate )) + geom_point() +geom_smooth()
# 
# ggplot(rate.env, aes(x=temp, y=rate )) + geom_point() + facet_wrap(~hemi) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial))
# 
# ggplot(rate.env[rate.env$temp > 15, ], aes(x=temp, y=rate )) + geom_point() +  facet_wrap(~hemi) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial))
# 
# # no difference between habitat types at this scale of analysis
# ggplot(rate.env, aes(x=temp, y=rate, group=habitat )) + geom_point() + 
#   geom_smooth(method='glm',method.args=list(family=quasibinomial))
# ggplot(rate.env, aes(x=temp, y=rate, group=habitat )) + geom_point() + facet_wrap(~hemi) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial))


# temperature range
windows(4,3)
ggplot( rate.env, aes(x=abs(meanLat),y=sstrange,size=rate) ) + geom_point(pch=21) + 
   # geom_smooth( data=rate.env[rate.env$rate<0.1 & abs(rate.env$Lat) < 40,], aes(group=1), method='lm', se=F ) +
   # geom_smooth( data=rate.env[rate.env$rate>0.5 & abs(rate.env$Lat) < 40,], aes(group=1), method='lm', se=F ) +
  xlab( "Degrees from equator" ) + 
  ylab( "Range annual sea surface\ntemperature (°C)" ) +
  theme_classic()


# use all fish data available from seines to look at effect on predation
# ggplot(rate.env, aes(x=ENSPIE, y=rate, group=habitat )) + geom_point() + #facet_wrap(~habitat) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic() +
#   xlab("Consumer diversity (ENSPIE)") + ylab("Predation rate")
# ggplot(rate.env, aes(x=ENSPIEfish, y=rate, group=1 )) + geom_point() + #facet_wrap(~habitat) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic() +
#   xlab("Fish diversity (ENSPIE)") + ylab("Predation rate")
# ggplot(rate.env, aes(x=log10(Total), y=rate, group=habitat )) + geom_point() + #facet_wrap(~habitat) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial))
# ggplot(rate.env, aes(x=log10(TotalFish), y=rate, group=1 )) + geom_point() + #facet_wrap(~habitat) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial))
# 
# ggplot(rate.env, aes(x=logbio, y=rate )) + geom_point(alpha=0.1) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()
# # ggplot(rate.env, aes(x=log10(biomass.std), y=rate )) + geom_point(pch=1) +
# #   # geom_smooth(method='glm',method.args=list(family=quasibinomial)) +
# #   geom_smooth( data = rate.env[log10(rate.env$biomass) > c(-4),],
# #                method='glm',method.args=list(family=quasibinomial)) + theme_classic()
# ggplot(rate.env, aes(x=S.warm, y=rate )) + geom_point(alpha=0.1) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()
# ggplot(rate.env, aes(x=S.warm, y=logbio )) + geom_point(alpha=0.1) +
#   geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()

m1 <-  glmer( rate~log10(TotalFish.std*1000+1) + (1|Country), rate.env, family="binomial" ) 
m2 <-   glmer( rate~ENSPIEfish + (1|Country), rate.env, family="binomial" )
m3 <- glmer( rate~logbio + (1|Country), rate.env, family="binomial" ) 
AICctab(m1,m2,m3, nobs=nrow(rate.env))


rate.env[ log10(rate.env$biomass) < -4 & !is.na(rate.env$biomass), ]







###--------------------------------------------------------------------------------------------------
### MODELS



# Final model set as of 13 September 2018
# restrict dataset to rows with in situ and satellite measurements of temperature AND predator surveys
d <- rate.env[ !is.na(rate.env$sstmean) & !is.na(rate.env$temp) & !is.na(rate.env$ENSPIEfish), ]
d <- rate.env %>% 
  filter( !is.na(rate.env$sstmean) & !is.na(rate.env$temp) & !is.na(rate.env$ENSPIEfish) ) %>%
  select( Site, habitat, rate, sstmean,  temp, ENSPIEfish, S.warm, logbio, logfish, MDS1, MDS2, MDS3 ) %>% 
  distinct()
m.sst                    <- glmer( rate~scale(sstmean) + (1|Site), d, family='binomial' )
m.temp                   <- glmer( rate~scale(temp) + (1|Site), d, family='binomial' )
m.temp.habitat           <- glmer( rate~scale(temp) + habitat + (1|Site), d, family='binomial' )
m.temp.logbio            <- glmer( rate~scale(temp) + scale(logbio) + (1|Site), d, family='binomial' )
m.S.warm                 <- glmer( rate~scale(S.warm) + (1|Site), d, family='binomial' )
m.temp.S.warm            <- glmer( rate~scale(temp) + scale(S.warm) + (1|Site), d, family='binomial' )
m.temp.logbio.S.warm     <- glmer( rate~scale(temp) + scale(logbio)+scale(S.warm) + (1|Site), d, family='binomial' )
m.ENSPIE                 <- glmer( rate~scale(ENSPIEfish) + (1|Site), d, family='binomial' )
m.temp.ENSPIE            <- glmer( rate~scale(temp) + scale(ENSPIEfish) + (1|Site), d, family='binomial' )
m.temp.logbio.ENSPIE     <- glmer( rate~scale(temp) + scale(logbio)+scale(ENSPIEfish) + (1|Site), d, family='binomial' )
m.logbio                 <- glmer( rate~scale(logbio) + (1|Site), d, family='binomial' )
m.temp.logabund          <- glmer( rate~scale(temp) + scale(logfish) + (1|Site), d, family='binomial' )
m.composition            <- glmer( rate~ MDS1 +MDS2 +MDS3 + (1|Site), d, family='binomial' )
m0                       <- glmer( rate ~ 1 + (1|Site), d, family="binomial" )

aictable <- AICctab( m.sst, m.composition, m.temp.logabund,m.logbio,m.temp,m.temp.habitat,
                     m.temp.logbio,m.temp.S.warm,m.temp.logbio.S.warm, m.temp.ENSPIE,
                     m.temp.logbio.ENSPIE, m.S.warm, m.ENSPIE, m0, 
                     nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )
aic.df <- as.data.frame(print(aictable))
aic.df$model <- rownames(aic.df)
# write the AIC table to disk
write.csv( aic.df, "Output Data/Bitemap_AICtable_final_model_set.csv", row.names = FALSE )

# diversity and composition
rate.predator <- rate.env %>%
  filter( !is.na(rate.env$ENSPIEfish) & !is.na(rate.env$MDS1) ) %>%
  distinct()
mfull.ENSPIE                 <- glmer( rate~scale(ENSPIEfish) + (1|Site), rate.predator, family='binomial' )
mfull.composition            <- glmer( rate~ MDS1 +MDS2 +MDS3 + (1|Site), rate.predator, family='binomial' )
AICctab( mfull.ENSPIE,mfull.composition,   
         nobs=nrow(rate.env), weights=TRUE, delta = TRUE, base = TRUE )

r.squaredGLMM(mfull.ENSPIE)
r.squaredGLMM(mfull.composition)

rate.env2 <- rate.env %>% distinct()
ggplot( data=rate.env, aes(x=ENSPIEfish, y=rate)) + geom_point()
ggplot( data=rate.env, aes(x=MDS1, y=rate)) + geom_point()

ggplot( rate.env, aes(x=MDS1,y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family=binomial))
ggplot( rate.env, aes(x=MDS2,y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family=binomial))



### MEDIATION MODEL
d <- rate.env %>%
  filter( !is.na(sstmean) & !is.na(MDS1) ) %>%
  select( Site, rate, sstmean, MDS1 ) %>%
  distinct()
dsite <- ddply( d, .(Site), summarise, rate=mean(rate), MDS1=mean(MDS1), sstmean=mean(sstmean))

# step 1
M1 <- glmer( rate~ sstmean + (1|Site), d, family="binomial" )
summary(M1)
# step 2
# M1 <- lm( rate~ sstmean , dsite )
M2 <- lm( MDS1~ sstmean , dsite )
# M3 <- lm( rate~sstmean+MDS1, dsite)
summary(M2)
# step 3
M3 <- glmer( rate~  sstmean + MDS1 + (1|Site), d, family="binomial" )
summary(M3)
# effect of sstmean goes away. This is statistical evidence for full mediation
M3p <- glmer( rate~ poly(sstmean,2) + MDS1 + (1|Site), d, family="binomial" )
summary(M3p)
# it doesn't go away when the polynomial effects are included


# create a composite for the polynomial effect of sst on predation
Mp <- glmer( rate~ poly(sstmean,2) + (1|Site), d, family="binomial" )
summary(Mp)
sstpoly <- data.frame( poly(d$sstmean,2) )

d$SST2 <- with( sstpoly, 112.808*X1 - 114.172*X2 )
dsite <- ddply( d, .(Site), summarise, rate=mean(rate), MDS1=mean(MDS1), 
                sstmean=mean(sstmean), SST2 = mean(SST2))
plot( SST2 ~ sstmean, d )
plot( rate ~ SST2, d )
plot( MDS1 ~ SST2, d )
M2p <- lm( MDS1~ SST2 , dsite )
M3p <- glmer( rate~  SST2 + MDS1 + (1|Site), d, family="binomial" )
with(d, cor.test( MDS1, sstmean) )

## mediation analysis using "Causal Mediation Analysis"
med <- mediation::mediate( M2,M3, sims=1000, treat="sstmean", mediator="MDS1" )
summary(med)
med2 <- mediation::mediate( M2p,M3p, sims=1000, treat="SST2", mediator="MDS1" )
summary(med2)
# significant mediation either way

# just use site means to show patterns of predation against temperature and composition
rate.mean <- ddply( rate.env, .(Country), summarise,
                    sstmean=mean(sstmean), temp=mean(temp), 
                    rate=mean(rate), MDS1=mean(MDS1) )
a <- ggplot( rate.mean, aes(x=sstmean,y=rate)) + 
  geom_point( col='slateblue', alpha=0.8 ) + 
  geom_smooth(method='glm', method.args=list(family=binomial), se=F, col="black", lwd=0.5 ) + 
  geom_smooth(  se=FALSE, col="black", lty=2, lwd=0.5 ) +
  ylim( c(0,1) ) +
  theme_classic() + ylab("Predation rate") + xlab("Mean annual SST") 
b <- ggplot( rate.mean, aes(x=temp,y=rate)) + 
  geom_point( col='slateblue', alpha=0.8 ) + 
  geom_smooth(method='glm', method.args=list(family=binomial), se=F, col="black", lwd=0.5) + 
  geom_smooth(se=FALSE, col="black", lty=2, lwd=0.5 ) +
  ylim( c(0,1) ) +
  theme_classic() + ylab("") + xlab("in situ Temperature")
c <- ggplot( rate.mean, aes(x=MDS1,y=rate)) + 
  geom_point( col='slateblue', alpha=0.8 ) + 
  geom_smooth(method='glm', method.args=list(family=binomial), se=F, col="black", lwd=0.5 ) + 
  geom_smooth( se=FALSE, col="black", lty=2, lwd=0.5 ) +
  ylim( c(0,1) ) +
  theme_classic() + ylab("")
d <- ggplot( rate.mean, aes(x=sstmean,y=MDS1)) + 
  geom_point( col='slateblue', alpha=0.8 ) + 
  geom_smooth(method='lm', se=F, col="black", lwd=0.5 ) + 
  geom_smooth( se=F, col="black", lty=2, lwd=0.5 ) +
  theme_classic() + ylab("MDS1") + xlab("Mean annual SST")
library(cowplot)
windows(5,5)
plot_grid( d,a,b,c, ncol=4 )




#### eight sites have predation rates >0.25. In Fig. 3 there is separation between these points and others
## what is it about these sites?
highpred <- rate.env %>%
  filter( Site %in% c("FL","NC2","NSW1","NSW2","NC","Italy","QLD3","QLD4") ) 
highcomm <- highpred[,78:138]

highagg  <- aggregate( highcomm, by=list(highpred$Site), sum, na.rm=TRUE )
