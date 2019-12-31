#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Preliminary data analysis and exploration
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2019.04.18
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################

## UPDATES
# 2018.10.15: Matt adds Bitemap project to GitHub. Futher changes tracked there, and old changes tracked in archived script files

###   REMOVE NSW1.1 AND 1.2 AND USA(DE) BECAUSE THESE WERE NOT CONDUCTED IN SEAGRASS?
#     NO, THESE WERE KEPT BECAUSE THEY HELP US ESTIMATE HOW RATES CHANGE OVER SPACE

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
library(vegan) # distances, ordination, mantel, etc.

# plotting
library(lattice)
library(viridis)

# point colors
seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
unveg.color    <- "#d8b365"

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
pops   <- read.csv( '../Data/OceanBitemap_Squidpop_Data_20190322.csv', stringsAsFactors = FALSE )
# rename pops columns
names(pops) <- c("timeStamp","name","email","Lat","Long","dateRetrieved","timeRetrieved",
                 "dateDeployed","timeDeployed","habitat","habitatDescription","numberDeployed","missing1",
                 "numberRetrieved","missing24","Country","Institution","prop1","prop24","notes")
# convert vegetated and unvegetated sites to common categories
pops$habitat[ pops$habitat %in% c("Seagrass","Seagrass Meadow","seagrass")] <- "Seagrass"
pops$habitat[ pops$habitat %in% c("Muddy Bottom","Sandy Bottom", "unvegetated ","unvegetated", "Unveg")] <- "Unvegetated"
pops$habitat[ pops$habitat %in% c("Artificial Habitat (dock, breakwater, weir, etc.)","Rocky Reef")] <- NA
pops <- pops[ !is.na(pops$habitat), ]
with(pops, table(Country,habitat))
pops$habitat <- factor(pops$habitat)

## reduce the number of columns
pops <- pops %>%
  dplyr::select( Country, Lat, Long, dateDeployed, habitat, N1=numberDeployed, N24=numberRetrieved, prop1, prop24 )


# get rid of unvegetated site in Belize, because it was actually a coral rubble field
pops <- pops[ pops$Country != "Belize" | pops$habitat != "Unvegetated", ]

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
  dplyr::select( Country, siteName, contact=contactName, Lat=latDecimal, Long=longDecimal, density=shootDensity.m2, 
                 cover=percentCover, canopy=canopyHeight, richness,
                 14:29)
head(sea)  
  
  

## PREDATOR DATA
predator <- read.csv( "Output Data/Bitemap_SEINE_summaries_20190322.csv", stringsAsFactors = FALSE )



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



# summarize predator data to site and habitat
pred.site <-  ddply( predator, .(Country, habitat), summarize, 
       cpua=mean(total.cpua,na.rm=T),    cpua.fish=mean(total.cpua.fish,na.rm=T),
       ENSPIE=mean(ENSPIE,na.rm=T), ENSPIEfish=mean(ENSPIEfish,na.rm=T),
       richness=mean(richness,na.rm=T), richness.fish=mean(richness.fish),
       biomass.area=mean(biomass.area,na.rm=T) ) 

# 
# merge predators
pops.pred <- left_join( pops, pred.site, by=c("Country","habitat") )
# # merge seagrass bed characteristics
# pops.pred.hab <- left_join( pops.pred, sea, by=c("Country","Lat","Long" ) )


##


# Map of Sites
# Winkel-triple projection
library(maps)
library(sf)
library(rnaturalearth) # for countries...see https://bhaskarvk.github.io/user2017.geodataviz/notebooks/02-Static-Maps.nb.html
library(sp)

sites <- siteGPS
coordinates(sites) <- ~meanLong+meanLat
proj4string(sites) <- '+init=epsg:4326'

world <- rnaturalearth::countries110
world <- world[world$name != 'Antarctica',]

newproj = "+proj=tissot"

grid.lines.mj <- sp::gridlines(world,easts = seq(-180,180,by=30), norths = seq(-90,90,by=30))
grid.lines.mi <- sp::gridlines(world,easts = seq(-165,195,by=15), norths = seq(-90,90,by=15))
sites <- spTransform(sites, CRS(newproj))
world <- spTransform(world, CRS(newproj))
grid.lines.mj <- spTransform(grid.lines.mj,CRS(newproj))
grid.lines.mi <- spTransform(grid.lines.mi,CRS(newproj))

# make the figure
windows(12,7)
par(mar = c(1, 1,0,0) + 0.1, bg=NA )
plot( methods::as(world, 'Spatial'), expandBB=c(-0.1,-0.1,-0.1,-0.05) )

# plot gridlines
plot(grid.lines.mi, col=grey(0.95), add=T)
plot(grid.lines.mj, col=grey(0.9), add=T)
# text(labels(grid.lines.mj, side=c(1,4), labelCRS = CRS("+init=epsg:4326")),
     # col = grey(.6), offset=0.3)

# plot land masses and points
plot(world, add=TRUE, border="white", col=grey(0.8))   # previous border grey(0.2)
points(sites, pch=21, bg=scales::alpha("slateblue",0.75), cex=2)

# add labels
sites$pos <- c( 4,4,2,2,4,2,4,4,2,
                4,4,2,2,2,4,2,2,2,
                1,4,1,4,2,2,2,2,2,
                4,2,3,4,4,2,4,2,4,
                2,4,2,4,1,4 )
text(sites, labels=sites$Site, cex=1.2, pos=sites$pos )

# save the plot
library(gridGraphics)
library(grid)
grid.echo()
B <- grid.grab()
dev.off()

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
sites <- read.csv( "Output Data/Bitemap_BioORACLE_20190322.csv", stringsAsFactors = FALSE )
# sites <- left_join( sites,siteGPS[,1:2], by=c("Site") )

# # add ocean basin information here
# sites$Country
# sites$basin <- c( 1,1,1,1,1,1,1,1,1,2,2,2,2,1,3,3,4,2,3,1,1,2,2,2,2,1,1,1,1,2,2,2,2,2,1,2,2,2,1,1,2)
# sites$basin <- factor( sites$basin, levels=1:4, labels=c("Pacific","Atlantic", "Mediterranean", "Indian") )
# sites$coast <- c( 1,1,1,1,1,1,
#                   1,1,1,1,1,2,
#                   1,2,2,2,2,2,
#                   2,1,2,1,1,2,
#                   1,2,2,2,2,1,
#                   1,1,1,1,2,1,
#                   1,1,2,2,2)
# sites$coast <- factor( sites$coast, levels=1:2, labels=c("West","East") )

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
directenv <- read.csv( '../Data/Environmental Data/BiteMap_EnvironmentalData_Compilation_20190228.csv', stringsAsFactors = FALSE )
# combine in situ temperature and aqua
directenv$Temp[ is.na(directenv$Temp) ] <- directenv$aqua[ is.na(directenv$Temp) ]
names(directenv)[names(directenv)=="Seagrass.Unveg"] <- "habitat"
# directenv$habitat[ directenv$habitat=="seagrass"] <- "Seagrass"
directenv$habitat[ directenv$habitat=="Unveg"] <- "Unvegetated"
directenv <- droplevels(directenv)

# get average temp and salinity for each site on each day
directenvmean <- ddply( directenv, .(Site.Name,Lat,Long,habitat,Date), summarize,
                        temp=mean(Temp), sal=mean(Sal) )

# get average temp and salinity for each Country
directenvsite <- ddply( directenv, .(Country,habitat), summarize,
                        temp=mean(Temp), sal=mean(Sal) )

# combine satellite and in situ measurements
popenv <- left_join( pops.pred, directenvsite )
allenv <- left_join( popenv, sites )



# add fishing pressure data from the Sea Around Us Project
sau <- read.csv( "../Fishing Pressure/output/catchmin.csv", stringsAsFactors = FALSE )
allenv <- left_join( allenv, sau, by=c("Site"="site") )

# add human population density data from NASA
hpopdens <- read.csv( "../Human Population Density/NASA_human_pop.csv", stringsAsFactors = FALSE )
hpopdens <- hpopdens %>%
  select( Site, pop, pop30 )
allenv <- left_join( allenv, hpopdens )

# ### ----------------------------------------------------------------------------------
# ### Figures of environmental data (latitude, temperature, etc.)
# 
# ggplot( allenv, aes(x=meanLat,y=sstmax) ) + geom_point(size=2.5) + 
#   # geom_smooth( method='lm', formula = y~x^2, se=T ) +
#   geom_smooth( aes(group=1), col='black', se=T ) +
#   xlab("Latitude") + 
#   ylab(expression(paste("Maximum annual sea surface temperature (",degree,"C)",sep="")))
# 
# latitude on y axis for nice map comparison
windows(1.75,3)
( C <- ggplot( allenv, aes(x=meanLat,y=sstmean) ) +
  geom_smooth( aes(group=1),col='black', se=F, lwd=0.5 ) +
  geom_point(bg="slateblue", size=2, pch=21) +
  xlab("Degrees latitude") +
  ylab("Mean annual\nSST (°C)") +
  theme_classic() +
  scale_x_continuous(breaks=c(-30,-15,0,15,30,45,60),position="top") +
    theme( panel.background = element_rect(fill = "transparent"), # bg of the panel
           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
           panel.grid.major = element_blank(), # get rid of major grid
           panel.grid.minor = element_blank(), # get rid of minor grid
           legend.justification=c(1,0), legend.position=c(1,0.8), 
           legend.title = element_blank(),
           legend.background = element_blank() ) +
  # Flip axes to let match up with a map
  coord_flip() )
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
# site labels
allenv.site <- allenv %>%
  dplyr::group_by(Country,Site, sstmean) %>%
  dplyr::summarise( temp=mean(temp, na.rm=T) )

windows(7,5)
ggplot( allenv.site, aes(x=sstmean,y=temp) ) + 
  geom_smooth(se=F,col='black',lwd=0.5) +
  geom_abline( intercept=0, col='red') + 
  geom_point( size=2,bg='slateblue', pch=21 ) +
  geom_text_repel( aes(label=Site), size=3 ) +
  ylab("in situ water temperature (°C)") + xlab("Mean annual SST (°C)") #+
  # theme_classic() +
  # scale_y_continuous(position = "right") +
  # theme(axis.title.y.right  = element_text(angle = 90, vjust=0, hjust=0.5))

# 
# ###
# 
# 
# 
# 
# 
# 
# 

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
pop01<- rbind( pop, zeros )

# make factor for Country, date, and habitat
pop0$ID <- with( pop0, paste( Country,dateDeployed,round(Lat,3),habitat, sep="." ))
pop01$ID <- with( pop01, paste( Country,dateDeployed,round(Lat,3),habitat, sep="." ))

# windows(12,7)
# ggplot( popall0, aes(x=hour,y=noteaten,col=habitat, group=ID ) )+ facet_wrap(~Country,ncol=8) + geom_point() +
#   geom_smooth( method='glm', method.args=list(family=poisson), se=F, lwd=0.5 ) +
#   xlab("Hour of predation assay") + ylab("Number of prey remaining")
# dev.off()

### Get estimates for individual Poisson models
mods <- dlply( pop0, .(ID), glm, formula=noteaten~hour, family=poisson )
coefs <- ldply( mods, stats::coef )


## merge these back with environment data
popmerge <- pop0[,c(1,4,5,7,11)]
rates <- left_join( popmerge, coefs, by="ID" )
# for some reason I get a lot of duplicates this way
rates <- rates[ !duplicated( rates ), ] 
# merge environment
allenv$ID <- with( allenv, paste( Country,dateDeployed,round(Lat,3),habitat, sep="." ))
rate.env  <- left_join( rates, allenv )
rate.env <- rate.env[ !duplicated( rate.env ), ] 

# convert rate to 1-hour so that rate makes more sense (bigger equals faster)
rate.env$rate <- 1-exp(rate.env$hour)
# log transform fish abundance
rate.env$logabund <- log10(rate.env$cpua+0.01)
rate.env$logfish <- log10(rate.env$cpua.fish+0.01)
# log transform biomass
rate.env$logbio <- log10(rate.env$biomass.area+0.01)


# absolute value for latitude
rate.env$abLat <- abs(rate.env$Lat)

# missing site label for Bodega Bay
rate.env$Site[ rate.env$Country == "USA (CA4)"] <- "CA4"

# # make habitat labels pretty
# rate.env <- rate.env %>% 
#   mutate( habitat = factor( habitat, c("Unveg", "Seagrass"), c("Unvegetated","Seagrass") ) )

# write these data to disk
write.csv( rate.env, "Output Data/Bitemap_rate.env.20190423.csv", row.names=FALSE )

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
windows(12,7)
ggplot( rate.env, aes(x=habitat,y=rate, col=habitat)) + facet_wrap(~Site, ncol=7, scales="free_y") + 
  stat_summary(fun.y = mean, geom = "point", pch=1, size=4, col='black') + geom_point(size=3, alpha=0.5) +   
  scale_color_manual( values=c(seagrass.color,unveg.color) ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Habitat") + ylab("Predation rate")
dev.off()

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


# sites with potential problems
# Victoria -- multiple sites per Macreadie. During left_join, extra rows created
# FIXED
rates[rates$Country=="Australia (VIC)",]
allenv[allenv$Country=="Australia (VIC)",]
rate.env[rate.env$Country=="Australia (VIC)",]
# Croatia  -- same rates for each habitat, multiple times?     YES
rates[rates$Country=="Croatia",]
rate.env[rate.env$Country=="Croatia",]
rate.env[rate.env$Country=="Croatia",]
# VA2      -- same rates for each habitat, multiple times?     FIXED
rates[rates$Country=="USA (VA2)",]
rate.env[rate.env$Country=="USA (VA2)",]
rate.env[rate.env$Country=="USA (VA2)",]





#### ---------------------------------------------------------------------------------------
#### look at predation data without averaging at site level

# point colors
seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
unveg.color    <- "#d8b365"

# predation rate data
windows(2.5,3)
( A <- ggplot( rate.env, aes(x=meanLat, y=rate, col=habitat, fill=habitat, group=hemi ))  +
  geom_smooth( se=F, col='black', lwd=0.5, lty=2 )  +
  # geom_smooth(aes(group=1),method='glm',method.args=list(family=quasibinomial),
  #             formula=y~poly(x,2), se=F, col='black', lwd=0.5, lty=2 ) +
  geom_smooth(aes(group=hemi),method='glm',method.args=list(family=quasibinomial),
              formula=y~poly(x,2), se=F, col='black', lwd=0.5 ) +
  geom_point(pch=21,size=1.5,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',unveg.color)) +
  theme_classic() + 
  scale_x_continuous(breaks=c(-30,-15,0,15,30,45,60)) +
  ylab( expression(paste('Predation rate (',hr^-1,')')) ) +
  xlab( "Latitude" ) +
  # theme(legend.justification=c(1,0), legend.position=c(1,0.8), 
        # legend.title = element_blank(), legend.background = element_blank() ) +
  theme( panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.justification=c(1,0), legend.position=c(1,0.8), 
      legend.title = element_blank(),
      legend.background = element_blank() ) +# get rid of legend panel bg
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0))) +
  coord_flip() )

# plot all together
library( cowplot )
windows(12,4)
plot_grid( A,B,C, ncol=3, labels = "AUTO", rel_widths = c(1,2,0.75),
           align = "v", axis="b", scale=c(1,1.1,1) )


## show latitudinal pattern, separated by major regions
hemicol = c("slateblue","darkorange")
# remove Indian and Mediterranean sites here, and remove Chile + Brazil
rate.plot <- rate.env %>%
  filter( !(basin %in% c("Indian","Mediterranean")) ) %>%
  filter( Site != "Chile", Site != "Brazil" ) 

windows(6,3)
ggplot(rate.plot, aes(x=abLat, y=rate, group=hemi, col=hemi)) + 
  facet_wrap( ~basin, ncol=4 ) +  
  geom_smooth( method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2), se=F ) +
  geom_smooth( method='glm',method.args=list(family=quasibinomial), se=F, lty=2 ) +
  geom_point( alpha=0.3 ) +
  xlab('Latitude') + ylab(expression(paste('Predation rate (',hr^-1,')'))) +
  scale_color_manual( values=hemicol ) +
  guides( color=guide_legend(title="Hemisphere") )

ggplot(rate.plot, aes(x=sstmean, y=rate, group=hemi, col=hemi)) + 
  facet_wrap( ~basin, ncol=4 ) +  
  geom_smooth( method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2), se=F ) +
  geom_smooth( method='glm',method.args=list(family=quasibinomial), se=F, lty=2 ) +
  geom_point( alpha=0.3 ) +
  xlab('Latitude') + ylab(expression(paste('Predation rate (',hr^-1,')'))) +
  scale_color_manual( values=hemicol ) +
  guides( color=guide_legend(title="Hemisphere") )



# site-level rates
rate.site <- ddply( rate.plot, .(Country,Site,hemi,basin,coast), summarise, rate=mean(rate), abLat=mean(abLat), sstmean=mean(sstmean) )
# windows(9,5)
windows(6,3)
ggplot(rate.site, aes(x=abLat, y=rate, group=hemi, col=hemi)) + 
  facet_wrap( ~basin, ncol=4 ) +  
  geom_smooth( method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2), se=F, size=0.5 ) +
  geom_smooth( method='glm',method.args=list(family=quasibinomial), se=F, lty=2, size=0.5 ) +
  # geom_text_repel( aes(label=Site), label.padding = 0.25, box.padding=0.29, size=3, show_guide=FALSE ) +
  geom_point( alpha=0.5 ) +
  xlab('Degrees from equator') + ylab(expression(paste('Predation rate (',hr^-1,')'))) +
  scale_color_manual( values=hemicol ) +
  guides( color=guide_legend(title="Hemisphere") )
ggplot(rate.site, aes(x=sstmean, y=rate, group=hemi, col=hemi)) + 
  facet_wrap( ~basin, ncol=4 ) +  
  geom_smooth( method='glm',method.args=list(family=quasibinomial), formula=y~poly(x,2), se=F, size=0.5 ) +
  geom_smooth( method='glm',method.args=list(family=quasibinomial), se=F, lty=2, size=0.5 ) +
  # geom_text_repel( aes(label=Site), size=3, show_guide=FALSE ) +
  geom_point( alpha=0.5 ) +
  xlab('Mean Annual SST (°C)') + ylab(expression(paste('Predation rate (',hr^-1,')'))) +
  scale_color_manual( values=hemicol ) +
  guides( color=guide_legend(title="Hemisphere") )

# dev.off()
rate.site$Country[ rate.site$abLat < 22 & rate.site$rate < 0.25 ]
rate.site$Country[ rate.site$rate < 0.25 ]


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
ggplot(rate.env, aes(x=ENSPIE, y=rate, group=habitat )) + geom_point() + #facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic() +
  xlab("Consumer diversity (ENSPIE)") + ylab("Predation rate")
ggplot(rate.env, aes(x=ENSPIEfish, y=rate, group=1 )) + geom_point() + #facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic() +
  xlab("Fish diversity (ENSPIE)") + ylab("Predation rate")

ggplot(rate.env, aes(x=log10(cpua), y=rate )) + geom_point() + #facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + geom_smooth(se=F)
ggplot(rate.env, aes(x=log10(cpua+0.011), y=rate )) + geom_point() + #facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + geom_smooth(se=F)



windows(4,4)
ggplot(rate.env, aes(x=biomass.area, y=rate )) + geom_point(alpha=0.1) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()
ggplot(rate.env, aes(x=logbio, y=rate )) + geom_point(alpha=0.1) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()

ggplot(rate.env, aes(x=log10(biomass.area+1), y=rate )) + geom_point(alpha=0.1) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()

ggplot(rate.env, aes(x=S.warm, y=rate )) + geom_point(alpha=0.1) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()
ggplot(rate.env, aes(x=S.warm, y=logbio )) + geom_point(alpha=0.1) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + theme_classic()

m1 <-  glmer( rate~cpua + (1|Country), rate.env, family="binomial" ) 
m2 <-   glmer( rate~ENSPIEfish + (1|Country), rate.env, family="binomial" )
m3 <- glmer( rate~logbio + (1|Country), rate.env, family="binomial" ) 
AICctab(m1,m2,m3, nobs=nrow(rate.env))

# Biomass and Density both 





## independent models of habitat effects within sites

tib <- as_tibble(rate.env)
# which sites have both habitat types
with(tib, table(Site,habitat))
tb <- tib %>%
  group_by(Site) %>%
  dplyr::summarize( seagrass=length(habitat[habitat=="Seagrass"]),
                    unveg=length(habitat[habitat=="Unvegetated"]) ) %>%
  filter( seagrass==0 | unveg==0 )
  
tib <- tib %>%
  filter( !(Site %in% tb$Site) )

tib$habitat <- factor( tib$habitat, levels=c("Unvegetated","Seagrass"))

mods2 <- tib %>%
  nest(-Site) %>% 
  mutate(
    # fit = purrr::map(data, ~ glm(rate ~ habitat, data = .x, family='quasibinomial')),
    fit = purrr::map(data, ~ lm(rate ~ habitat, data = .x)),
    tidied = purrr::map(fit, broom::tidy),
    augmented = purrr::map(fit, broom::augment)
  ) 

mods2 %>% 
  unnest(tidied) %>%
  filter(term=="habitatSeagrass") %>%
  arrange( p.value )
(m2 <- mods2 %>% 
  unnest(tidied) %>%
  filter(term=="habitatSeagrass") %>%
  arrange( -estimate ))
windows(2,2)
ggplot(m2, aes(x=estimate) ) + geom_histogram( col="black", fill="white", bins=25 ) +
  theme_classic() +
  xlim(c(-0.6,0.6))


mods2 %>% 
  unnest(augmented)


###--------------------------------------------------------------------------------------------------
### MODELS



# Final model set as of 20 February 2018
# restrict dataset to rows with in situ and satellite measurements of temperature AND predator surveys
# note there are still NA values for catch, so we shoudn't use this variable
d <- rate.env %>% 
  filter( !is.na(temp) &  !is.na(logfish)  ) %>% #& !is.na(sstmean)  & !is.na(pop)   )  %>%  # only allow data with SST  and seines #& !is.na(temp) -- include this if we want all temperature data
  select( Site, abLat, habitat, rate, sstmean,  temp, chlomax, 
          ENSPIE, ENSPIEfish, richness, richness.fish, 
          S.warm, logbio, logfish, MDS1, MDS2, MDS3,
          catch_min, catch_max, catch_mean, pop, pop30,
          hemi, coast, abLat ) %>% 
  mutate( logcatch=log10(catch_mean), logpop=log10(pop), 
          logpop30=log10(pop30), logchl=log10(chlomax) ) %>%
  distinct()

# splot
windows()
psych::pairs.panels( d[,c('abLat','sstmean','temp','logchl','logcatch','logpop',
                          'logbio', 'logfish', 'ENSPIEfish','richness.fish','MDS1')],
                      scale=F )

# # remove Korea
# d <- d %>%
#   filter( Site!="Korea" ) # to see how catch_mean behaves

# set the family for all glms
fam <- "binomial"

# intercept-only model
m0                       <- glmer( rate ~ 1 + (1|Site), d, family=fam )


# individual predictors
m.sst                    <- glmer( rate~ scale(sstmean) + (1|Site), d, family=fam )
m.sst2                   <- glmer( rate~ poly(sstmean,2) + (1|Site), d, family=fam )
m.temp                   <- glmer( rate~ scale(temp) + (1|Site), d, family=fam )
m.prod                   <- glmer( rate~ scale(logchl) + (1|Site), d, family=fam )
m.catch                  <- glmer( rate~ scale(logcatch) + (1|Site), d, family=fam )
m.hab                    <- glmer( rate~ as.factor(habitat) + (1|Site), d, family=fam )
m.comp                   <- glmer( rate~ scale(MDS1) +scale(MDS2) +scale(MDS3) + (1|Site), d, family=fam )
m.comp1                   <- glmer( rate~ scale(MDS1)  + (1|Site), d, family=fam )
m.bio                    <- glmer( rate~ scale(logbio) + (1|Site), d, family=fam )
# m.abund                  <- glmer( rate~ scale(logabund) + (1|Site), d, family=fam )
m.fish                   <- glmer( rate~ scale(logfish) + (1|Site), d, family=fam )
m.ENSPIE                 <- glmer( rate~ scale(ENSPIE) + (1|Site), d, family=fam )
m.ENSPIEfish             <- glmer( rate~ scale(ENSPIEfish) + (1|Site), d, family=fam )
m.richness.fish          <- glmer( rate~ scale(richness.fish) + (1|Site), d, family=fam )
m.richness               <- glmer( rate~ scale(richness) + (1|Site), d, family=fam )
m.pop                    <- glmer( rate~ scale(logpop) + (1|Site), d, family=fam )
# m.pop30                  <- glmer( rate~ scale(logpop30) + (1|Site), d, family=fam )
# m.S.warm                 <- glmer( rate~ scale(S.warm) + (1|Site), d, family=fam )

AICctab( m.sst2,m.sst, m.temp, m.prod, m.catch, m.hab, m.comp, m.bio, m.fish, 
         m.ENSPIE, m.ENSPIEfish, m.pop, m.pop30, m0, 
         nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )

AICctab( m.sst, m.temp, m.prod, m.catch, m.hab, m.comp, m.bio, m.fish, 
         m.ENSPIEfish, m.richness.fish, m.pop, m.pop30, 
         m0, 
         nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )

# among single predictors, temperature and composition and diversity are strongest. 
# unless quadrtic temperature model included

# nothing can compete with the model with quadratic effects of temperature
# 
# # when just comparing models with temperature data and composition data, composition even better (this includes composition data from videos and transects)
# d <- rate.env %>% 
#   filter( !is.na(sstmean) & !is.na(temp)  & !is.na(MDS1)  )  %>%  #
#   select( Site, habitat, rate, sstmean,  temp, chlomax, ENSPIEfish, S.warm, logbio, logfish, MDS1, MDS2, MDS3, catch_min, catch_max, catch_mean ) %>% 
#   distinct()
# m.temp                   <- glmer( rate~ scale(temp) + (1|Site), d, family=fam )
# m.nmds                   <- glmer( rate~ MDS1 +MDS2 +MDS3 + (1|Site), d, family=fam )
# AICctab( m.temp, m.nmds, m0, 
#          nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )


# combinations with temperature or sst
m.temp.prod             <- glmer( rate~ scale(temp) + scale(logchl) + (1|Site), d, family=fam )
m.temp.habitat          <- glmer( rate~ scale(temp) + habitat + (1|Site), d, family=fam )
m.temp.bio              <- glmer( rate~ scale(temp) + scale(logbio) + (1|Site), d, family=fam )
# m.temp.S.warm           <- glmer( rate~ scale(temp) + scale(S.warm) + (1|Site), d, family=fam )
m.temp.ENSPIE           <- glmer( rate~ scale(temp) + scale(ENSPIEfish) + (1|Site), d, family=fam )
m.temp.richness         <- glmer( rate~ scale(temp) + scale(richness.fish) + (1|Site), d, family=fam )
m.temp.abund            <- glmer( rate~ scale(temp) + scale(logfish) + (1|Site), d, family=fam )
m.temp.catch            <- glmer( rate~ scale(temp) + scale(logcatch) + (1|Site), d, family=fam )
m.temp.comp             <- glmer( rate~ scale(temp) + scale(MDS1) +scale(MDS2) +scale(MDS3) + (1|Site), d, family=fam )
m.temp.pop              <- glmer( rate~ scale(temp) + scale(logpop) + (1|Site), d, family=fam )
AICctab( m.temp.prod, m.temp.habitat, m.temp.bio, 
         m.temp.ENSPIE, m.temp.richness, 
         m.temp.abund,m.temp.catch,m.temp.comp, m.temp.pop,
         m0, m.comp,
         nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )

m.sstmean.prod             <- glmer( rate~ scale(sstmean) + scale(logchl) + (1|Site), d, family=fam )
m.sstmean.habitat          <- glmer( rate~ scale(sstmean) + habitat + (1|Site), d, family=fam )
m.sstmean.bio              <- glmer( rate~ scale(sstmean) + scale(logbio) + (1|Site), d, family=fam )
# m.sstmean.S.warm           <- glmer( rate~ scale(sstmean) + scale(S.warm) + (1|Site), d, family=fam )
m.sstmean.ENSPIE           <- glmer( rate~ scale(sstmean) + scale(ENSPIEfish) + (1|Site), d, family=fam )
m.sstmean.richness         <- glmer( rate~ scale(sstmean) + scale(richness.fish) + (1|Site), d, family=fam )
m.sstmean.abund            <- glmer( rate~ scale(sstmean) + scale(logfish) + (1|Site), d, family=fam )
m.sstmean.catch            <- glmer( rate~ scale(sstmean) + scale(logcatch) + (1|Site), d, family=fam )
m.sstmean.comp             <- glmer( rate~ scale(sstmean) + scale(MDS1) +scale(MDS2) +scale(MDS3) + (1|Site), d, family=fam )
m.sstmean.pop              <- glmer( rate~ scale(sstmean) + scale(logpop) + (1|Site), d, family=fam )
AICctab( m.sstmean.prod, m.sstmean.habitat, m.sstmean.bio, 
         m.sstmean.ENSPIE, m.sstmean.richness, 
         m.sstmean.abund,m.sstmean.catch,m.sstmean.comp, m.sstmean.pop,
         m0, m.comp,
         nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )

# generally models with temp have higher posterior probability. more linear predictor
# but relative probability of models shifts with choice of thermal variable
ggplot(d, aes(x=therm,y=rate)) + geom_point() + geom_smooth()
ggplot(d, aes(x=MDS1,y=rate)) + geom_point() + geom_smooth()
ggplot(d, aes(x=logchl,y=rate)) + geom_point() + geom_smooth()
ggplot(d, aes(x=logpop,y=rate)) + geom_point() + geom_smooth()
ggplot(d, aes(x=logbio,y=rate)) + geom_point() + geom_smooth()
ggplot(d, aes(x=logcatch,y=rate)) + geom_point() + geom_smooth()
ggplot(d, aes(x=therm,y=MDS1)) + geom_point() + geom_smooth()


# composition and catch are in the top when combined with temperature

m.full                   <- glmer( rate~ scale(MDS1) +scale(MDS2) +scale(MDS3) + scale(ENSPIEfish) +scale(temp) + (1|Site), d, family=fam )
m.full2                  <- glmer( rate~ MDS1 +MDS2 +MDS3 + scale(logbio) + scale(logchl) +scale(therm) + factor(habitat) + (1|Site), d, family=fam )
m.full3                  <- glmer( rate~ scale(MDS1) + scale(logcatch) +scale(therm)  + (1|Site), d, family=fam )

AICctab( m.full, m.full2, m.full3,m.sst,m.therm.comp,m.comp,
         nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )


m.prune1              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)  +scale(logbio) +  
                                 scale(temp) + scale(logchl) + 
                                 (1|Site), d, family=fam )
m.prune2              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)+ scale(logbio) + 
                                 scale(therm)  + (1|Site), d, family=fam )
m.prune3              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)+ scale(logbio) + 
                                  scale(logfish) + (1|Site), d, family=fam )
m.prune4              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)+ scale(logbio) +
                                  scale(logchl) + (1|Site), d, family=fam )
m.prune5              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)+ scale(logbio) +
                                  scale(logpop) + (1|Site), d, family=fam )
m.biology             <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)  +scale(logbio) +  
                                 scale(logchl) +  scale(logcatch) +
                                 (1|Site), d, family=fam )
m.nocomp              <- glmer( rate~ scale(logbio) + scale(logchl) +
                                 scale(therm) + scale(logcatch) +
                                 (1|Site), d, family=fam )
m.prune6                 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)  +
                                    scale(therm) + scale(logbio) + scale(logcatch) +
                                  (1|Site), d, family=fam )
m.prune7                 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)  +
                                     scale(therm) + scale(logfish) + scale(logcatch) +
                                     (1|Site), d, family=fam )
m.prune8              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)  +scale(logfish) +  
                                  scale(therm) + scale(logchl) + 
                                  (1|Site), d, family=fam )
m.prune9              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)  +scale(logfish) +  
                                  scale(therm) + 
                                  (1|Site), d, family=fam )
m.prune10              <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3)  +scale(richness.fish) +  
                                  scale(therm) + 
                                  (1|Site), d, family=fam )


aictable <- AICctab( m.sst, m.temp, m.prod, m.hab, 
                     m.comp, m.bio, m.fish, m.catch, m.pop,
                     m.ENSPIEfish, m.richness.fish, 
                     m.therm.prod, m.therm.habitat, m.therm.bio, 
                     m.therm.ENSPIE, m.therm.richness, m.therm.pop,
                     m.therm.abund,m.therm.catch,m.therm.comp,
                     m.prune1, m.prune2,
                     m.prune3, m.prune4, m.prune5, #m.prune6, m.prune7,
                     m.prune8, m.prune9,m.prune10,
                     m.biology, m.nocomp,
                     m.full,m.full2,
                     m0,
                     nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )
aic.df <- as.data.frame(print(aictable))
aic.df$model <- rownames(aic.df)
# write the AIC table to disk
write.csv( aic.df, "Output Data/Bitemap_AICtable_final_model_set.csv", row.names = FALSE )


# # interaction test
# mbest <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
#                 scale(temp) +  scale(logbio) +
#                 (1|Site), d, family=fam ) 
# minter <- glmer( rate~ scale(MDS1)*scale(ENSPIEfish) +  scale(MDS2) + scale(MDS3) +
#                 scale(temp) +  scale(logbio) +
#                 (1|Site), d, family=fam )


## Brute force all model combinations of interest 
options( warn = 1 )
m00 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp) +   
                (1|Site), d, family=fam )
m01 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
               scale(temp) +  scale(logbio) + 
               (1|Site), d, family=fam ) # failed to converge
m02 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
               scale(temp) +  scale(logbio) + scale(logfish) +
               (1|Site), d, family=fam ) # failed to converge
m03 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
               scale(temp) +  scale(logbio) + scale(ENSPIEfish) +
               (1|Site), d, family=fam )
m04 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
               scale(temp) +  scale(logbio) + scale(richness.fish) +
               (1|Site), d, family=fam ) # failed to converge
m05 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
               scale(temp) +  scale(logbio) + scale(logchl) +
               (1|Site), d, family=fam )
# m06 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
#                scale(temp) +  scale(logbio) + scale(logcatch) +
#                (1|Site), d, family=fam )
m07 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
               scale(temp) +  scale(logbio) + scale(logpop) +
               (1|Site), d, family=fam ) # failed to converge
m08 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
               scale(temp) +  scale(logbio) + factor(habitat) +
               (1|Site), d, family=fam ) # failed to converge
m09 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp) +
                (1|Site), d, family=fam )
m10 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp)  + scale(logfish) +
                (1|Site), d, family=fam )
m11 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp)  + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m12 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp)  + scale(richness.fish) +
                (1|Site), d, family=fam )
m13 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp)  + scale(logchl) +
                (1|Site), d, family=fam )
# m14 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
#                 scale(temp)  + scale(logcatch) +
#                 (1|Site), d, family=fam )
m15 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp)  + scale(logpop) +
                (1|Site), d, family=fam )
m16 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(temp)  + factor(habitat) +
                (1|Site), d, family=fam )
m17 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  scale(logbio) + 
                (1|Site), d, family=fam )
m18 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  scale(logbio) + scale(logfish) +
                (1|Site), d, family=fam ) # failed to converge
m19 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  scale(logbio) + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m20 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  scale(logbio) + scale(richness.fish) +
                (1|Site), d, family=fam ) # failed to converge
m21 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  scale(logbio) + scale(logchl) +
                (1|Site), d, family=fam ) 
# m22 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
#                 scale(sstmean) +  scale(logbio) + scale(logcatch) +
#                 (1|Site), d, family=fam )
m23 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  scale(logbio) + scale(logpop) +
                (1|Site), d, family=fam )# failed to converge
m24 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  scale(logbio) + factor(habitat) +
                (1|Site), d, family=fam )
m25 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) +  
                (1|Site), d, family=fam )
m26 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean)  + scale(logfish) +
                (1|Site), d, family=fam )
m27 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean)  + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m28 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean)  + scale(richness.fish) +
                (1|Site), d, family=fam )
m29 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean) + scale(logchl) +
                (1|Site), d, family=fam )
# m30 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
#                 scale(sstmean)  + scale(logcatch) +
#                 (1|Site), d, family=fam )
m31 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean)  + scale(logpop) +
                (1|Site), d, family=fam )
m32 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(sstmean)  + factor(habitat) +
                (1|Site), d, family=fam )
m33 <- glmer( rate~ 
                scale(temp) +  scale(logbio) + 
                (1|Site), d, family=fam )
m34 <- glmer( rate~ 
                scale(temp) +  scale(logbio) + scale(logfish) +
                (1|Site), d, family=fam )
m35 <- glmer( rate~ 
                scale(temp) +  scale(logbio) + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m36 <- glmer( rate~ 
                scale(temp) +  scale(logbio) + scale(richness.fish) +
                (1|Site), d, family=fam )
m37 <- glmer( rate~ 
                scale(temp) +  scale(logbio) + scale(logchl) +
                (1|Site), d, family=fam )
# m38 <- glmer( rate~ 
#                 scale(temp) +  scale(logbio) + scale(logcatch) +
#                 (1|Site), d, family=fam )
m39 <- glmer( rate~ 
                scale(temp) +  scale(logbio) + scale(logpop) +
                (1|Site), d, family=fam )
m40 <- glmer( rate~ 
                scale(temp) +  scale(logbio) + factor(habitat) +
                (1|Site), d, family=fam )
m41 <- glmer( rate~ 
                scale(temp) +  
                (1|Site), d, family=fam )
m42 <- glmer( rate~ 
                scale(temp)  + scale(logfish) +
                (1|Site), d, family=fam )
m43 <- glmer( rate~ 
                scale(temp)  + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m44 <- glmer( rate~ 
                scale(temp)  + scale(richness.fish) +
                (1|Site), d, family=fam )
m45 <- glmer( rate~ 
                scale(temp)  + scale(logchl) +
                (1|Site), d, family=fam )
m46 <- glmer( rate~ 
                scale(temp)  + scale(logcatch) +
                (1|Site), d, family=fam )
m47 <- glmer( rate~ 
                scale(temp)  + scale(logpop) +
                (1|Site), d, family=fam )
m48 <- glmer( rate~ 
                scale(temp)  + factor(habitat) +
                (1|Site), d, family=fam )
m49 <- glmer( rate~ 
                scale(sstmean) +  scale(logbio) + 
                (1|Site), d, family=fam )
m50 <- glmer( rate~ 
                scale(sstmean) +  scale(logbio) + scale(logfish) +
                (1|Site), d, family=fam )
m51 <- glmer( rate~ 
                scale(sstmean) +  scale(logbio) + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m52 <- glmer( rate~ 
                scale(sstmean) +  scale(logbio) + scale(richness.fish) +
                (1|Site), d, family=fam )
m53 <- glmer( rate~ 
                scale(sstmean) +  scale(logbio) + scale(logchl) +
                (1|Site), d, family=fam )
# m54 <- glmer( rate~ 
#                 scale(sstmean) +  scale(logbio) + scale(logcatch) +
#                 (1|Site), d, family=fam )
m55 <- glmer( rate~ 
                scale(sstmean) +  scale(logbio) + scale(logpop) +
                (1|Site), d, family=fam )
m56 <- glmer( rate~ 
                scale(sstmean) +  scale(logbio) + factor(habitat) +
                (1|Site), d, family=fam )
m57 <- glmer( rate~ 
                scale(sstmean) +  
                (1|Site), d, family=fam )
m58 <- glmer( rate~ 
                scale(sstmean)  + scale(logfish) +
                (1|Site), d, family=fam )
m59 <- glmer( rate~ 
                scale(sstmean)  + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m60 <- glmer( rate~ 
                scale(sstmean)  + scale(richness.fish) +
                (1|Site), d, family=fam )
m61 <- glmer( rate~ 
                scale(sstmean)  + scale(logchl) +
                (1|Site), d, family=fam )
# m62 <- glmer( rate~ 
#                 scale(sstmean)  + scale(logcatch) +
#                 (1|Site), d, family=fam )
m63 <- glmer( rate~ 
                scale(sstmean)  + scale(logpop) +
                (1|Site), d, family=fam )
m64 <- glmer( rate~ 
                scale(sstmean)  + factor(habitat) +
                (1|Site), d, family=fam )
m65 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + 
                (1|Site), d, family=fam )
m66 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + scale(logfish) +
                (1|Site), d, family=fam ) # failed to converge
m67 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m68 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + scale(richness.fish) +
                (1|Site), d, family=fam )
m69 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + scale(logchl) +
                (1|Site), d, family=fam )
m70 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + scale(logcatch) +
                (1|Site), d, family=fam )
m71 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + scale(logpop) +
                (1|Site), d, family=fam )
m72 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logbio) + factor(habitat) +
                (1|Site), d, family=fam )
m73 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 
                (1|Site), d, family=fam )
m74 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(logfish) +
                (1|Site), d, family=fam )
m75 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m76 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(richness.fish) +
                (1|Site), d, family=fam )
m77 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                 scale(logchl) +
                (1|Site), d, family=fam )
# m78 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
#                 scale(logcatch) +
#                 (1|Site), d, family=fam )
m79 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                scale(logpop) +
                (1|Site), d, family=fam )
m80 <- glmer( rate~ scale(MDS1) + scale(MDS2) + scale(MDS3) +
                factor(habitat) +
                (1|Site), d, family=fam )
m81 <- glmer( rate~ 
                scale(logbio) + 
                (1|Site), d, family=fam )
m82 <- glmer( rate~ 
                scale(logbio) + scale(logfish) +
                (1|Site), d, family=fam )
m83 <- glmer( rate~ 
                scale(logbio) + scale(ENSPIEfish) +
                (1|Site), d, family=fam )
m84 <- glmer( rate~ 
                scale(logbio) + scale(richness.fish) +
                (1|Site), d, family=fam )
m85 <- glmer( rate~ 
                scale(logbio) + scale(logchl) +
                (1|Site), d, family=fam )
# m86 <- glmer( rate~ 
#                 scale(logbio) + scale(logcatch) +
#                 (1|Site), d, family=fam )
m87 <- glmer( rate~ 
                scale(logbio) + scale(logpop) +
                (1|Site), d, family=fam )
m88 <- glmer( rate~ 
                scale(logbio) + factor(habitat) +
                (1|Site), d, family=fam )

# AIC comparison
# remove all model with catch from Sea Around Us
aictable <- AICctab( m01,m02,m03,m04,m05,    m07,        m10, # several models did not converge
                     m11,m12,m13,    m15,m16,m17,m18,m19,m20,
                     m21,    m23,m24,m25,m26,m27,m28,m29,
                     m31,m32,m33,m34,m35,m36,m37,    m39,m40,
                     m41,m42,m43,m44,m45,    m47,m48,m49,m50,
                     m51,m52,m53,    m55,m56,m57,m58,m59,m60,
                     m61,    m63,m64,m65,m66,m67,m68,m69,
                     m71,m72,m73,m74,m75,m76,m77,    m79,m80,
                     m81,m82,m83,m84,m85,    m87,m88,
                     m0,
                     nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )
# remove models with 7 fixed effects parameters
aictable <- AICctab( m01,                            m09,m10, # several models did not converge
                     m11,m12,m13,    m15,m16,m17,
                                     m25,m26,m27,m28,m29,
                     m31,m32,m33,m34,m35,m36,m37,    m39,m40,
                     m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,
                     m51,m52,m53,    m55,m56,m57,m58,m59,m60,
                     m61,    m63,m64,m65,m66,m67,m68,m69,
                     m71,m72,m73,m74,m75,m76,m77,    m79,m80,
                     m81,m82,m83,m84,m85,    m87,m88,
                     m0,m00,
                     nobs=nrow(d), weights=TRUE, delta = TRUE, base = TRUE )
aic.df <- as.data.frame(print(aictable))
aic.df$model <- rownames(aic.df)
# write to disk
write.csv( aic.df, "Output Data/Bitemap_AICtable_all_models.csv", row.names = FALSE )
aic.df$weight
# best models
m03
m01
m04
m02
m05
m07
m67
m65
m66
m19
m69
m16
m71
m17
m68
m72
m15


# build a smaller table to show off particluar models of interest



# combinations of composition, biomass, abundance, catch, and temperature mediate the quadratic effect
dotplot( ranef(m01,condVar=F) )
param.ci <- confint( m01 , devtol=0.1 )
param.ci2 <- confint( m65 )#, devtol=0.1 )

mod <- m65
par.use <- param.ci2

# make a figure
windows(3.5,4)
dp <- data.frame( mean=fixef(mod), par.use[-c(1),] )
lab <- c("Intercept","MDS1","MDS2","MDS3","in situ Temp", "log10(Biomass)")
# lab <- rownames(dp)
names(dp) <- c("mean", "lcl95","ucl95")
dl <- exp(dp)
dplot <- dp
par( mar=c(5,7,2,2)+0.1, las=1 )
yloc <- nrow(dplot):1
plot( x=rep(0,nrow(dplot) ), y=1:nrow(dplot), type="n", axes=F, 
      xlab="", ylab="", xlim=range(dplot) )
abline(v=0)
axis(1, at=c(-8,-6,-4,-2,0,2,4,6))
mtext( "Estimate", side=1, line=3 )
axis(2, at=yloc, labels = lab, tick = "F", line = -1 )
segments( x0=dplot$lcl95, x1=dplot$ucl95, y0=yloc, y1=yloc)
points( x=dplot$mean, y=yloc, pch=21, bg="white", cex=1.2 )

MuMIn::r.squaredGLMM(m01)
MuMIn::r.squaredGLMM(m65)
rsq::rsq(m01)




# diversity and composition
rate.predator <- rate.env %>%
  filter( !is.na(rate.env$ENSPIEfish) & !is.na(rate.env$MDS1) ) %>%
  distinct()
mfull.ENSPIE                 <- glmer( rate~scale(ENSPIEfish) + (1|Site), rate.predator, family='binomial' )
mfull.composition            <- glmer( rate~ MDS1 +MDS2 +MDS3 + (1|Site), rate.predator, family='binomial' )
AICctab( mfull.ENSPIE,mfull.composition,   
         nobs=nrow(rate.env), weights=TRUE, delta = TRUE, base = TRUE )

MuMIn::r.squaredGLMM(mfull.ENSPIE)
MuMIn::r.squaredGLMM(mfull.composition)

rate.env2 <- rate.env %>% distinct()
ggplot( data=rate.env, aes(x=ENSPIEfish, y=rate)) + geom_point()
ggplot( data=rate.env, aes(x=MDS1, y=rate)) + geom_point()

# site averages
rate.site <- rate.env %>%
  dplyr::group_by( Site, MDS1, MDS2, MDS3, temp, logbio, ENSPIEfish ) %>%
  dplyr::summarise( rate=mean(rate) )

ggplot( rate.site, aes(x=MDS1,y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family=binomial))
ggplot( rate.site, aes(x=MDS2,y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family=binomial))
ggplot( rate.site, aes(x=logbio,y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family=binomial))
ggplot( rate.site, aes(x=temp,y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family=binomial))
ggplot( rate.site, aes(x=ENSPIEfish,y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family=binomial))











### MEDIATION MODEL
d <- rate.env %>%
  filter( !is.na(sstmean) & !is.na(MDS1) ) %>%
  select( Site, rate, sstmean, temp, MDS1, MDS2, MDS3, logbio ) 
  # distinct()
dsite <-  d %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise( rate=mean(rate), MDS1=mean(MDS1),MDS2=mean(MDS2),MDS2=mean(MDS2), sstmean=mean(sstmean), logbio=mean(logbio, na.rm=T), temp=mean(temp) )

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

d$SST2 <- with( sstpoly, 204.30*X1 - 179.64*X2 )
dsite <- ddply( d, .(Site), summarise, rate=mean(rate), MDS1=mean(MDS1), 
                sstmean=mean(sstmean), SST2 = mean(SST2))
plot( SST2 ~ sstmean, d )
plot( rate ~ SST2, d )
plot( MDS1 ~ SST2, d )
M2p <- lm( MDS1~ SST2 , dsite )
M3p <- glmer( rate~  SST2 + MDS1 + (1|Site), d, family="binomial" )
with(d, cor.test( MDS1, sstmean) )

# create a composite for effects of MDS1, 2, and 3
Mm <- glmer( rate~ MDS1 + MDS2 + MDS3 + (1|Site), d, family="binomial" )
summary(Mm)
d$multiv <- with( d, 5.8799*MDS1 - 2.4766*MDS2 - 0.1891*MDS3 )
dsite <- ddply( d, .(Site), summarise, rate=mean(rate), MDS1=mean(MDS1), MDS2=mean(MDS2), MDS3=mean(MDS3), 
                sstmean=mean(sstmean), SST2 = mean(SST2), multiv = mean(multiv) )
M2m <- lm( multiv ~ SST2 , dsite )
M3m <- glmer( rate~  SST2 + multiv + (1|Site), d, family="binomial" )

## mediation analysis using "Causal Mediation Analysis"
  med <- mediation::mediate( M2,M3, sims=100, treat="sstmean", mediator="MDS1" )
summary(med)
med2 <- mediation::mediate( M2p,M3p, sims=500, treat="SST2", mediator="MDS1" )
summary(med2)
med3 <- mediation::mediate( M2m,M3m, sims=1000, treat="SST2", mediator="multiv" )
summary(med3)
# significant mediation either way










# just use site means to show patterns of predation against temperature and composition
rate.mean <- ddply( rate.env, .(Country), summarise,
                    sstmean=mean(sstmean), temp=mean(temp), 
                    rate=mean(rate), MDS1=mean(MDS1) )
# get R2 values for binomial models
r1 <- glm( rate~sstmean, data=rate.mean, family=quasibinomial )
r2 <- glm( rate~temp,    data=rate.mean, family=quasibinomial )
r3 <- glm( rate~MDS1,    data=rate.mean, family=quasibinomial )
library(rsq)
rsq(r1,adj=TRUE);rsq(r2,adj=TRUE);rsq(r3,adj=TRUE)

a <- ggplot( rate.mean, aes(x=sstmean,y=rate)) + 
  geom_smooth(method='glm', method.args=list(family=quasibinomial), se=T, col="black", lwd=0.5 ) + 
  geom_smooth(  se=FALSE, col="black", lty=2, lwd=0.5 ) +
  geom_point( col='slateblue', alpha=0.8 ) + 
  geom_text( y=0.7,x=min(rate.mean$sstmean,na.rm=T) + diff(range(rate.mean$sstmean,na.rm=T))*0.25, 
             label= expression(paste( R^2," = 0.07"  )), size=3, fontface="plain" ) +
  ylim( c(0,1) ) +
  theme_classic() + ylab("Predation rate") + xlab("Mean annual SST") +
  theme(plot.margin = unit(c(0, 0.1, 0, 0)+0.1, "cm"))
b <- ggplot( rate.mean, aes(x=temp,y=rate)) + 
  geom_smooth(method='glm', method.args=list(family=quasibinomial), se=T, col="black", lwd=0.5) + 
  geom_smooth(se=FALSE, col="black", lty=2, lwd=0.5 ) +
  geom_point( col='slateblue', alpha=0.8 ) + 
  geom_text( y=0.7,x=min(rate.mean$temp,na.rm=T) + diff(range(rate.mean$temp,na.rm=T))*0.25, 
             label= expression(paste( R^2," = 0.09"  )), size=3 ) +
  ylim( c(0,1) ) +
  theme_classic() + ylab("") + xlab("in situ Temperature") +
  theme(plot.margin = unit(c(0, 0.1, 0, 0)+0.1, "cm"))
c <- ggplot( rate.mean, aes(x=MDS1,y=rate)) + 
  geom_smooth(method='glm', method.args=list(family=quasibinomial), se=T, col="black", lwd=0.5 ) + 
  geom_smooth( se=FALSE, col="black", lty=2, lwd=0.5 ) +
  geom_point( col='slateblue', alpha=0.8 ) + 
  geom_text( y=0.7,x=min(rate.mean$MDS1,na.rm=T) + diff(range(rate.mean$MDS1,na.rm=T))*0.25, 
             label= expression(paste( R^2," = 0.44"  )), size=3) +
  ylim( c(0,1) ) +
  theme_classic() + ylab("") +
  theme(plot.margin = unit(c(0, 0.1, 0, 0)+0.1, "cm"))
d <- ggplot( rate.mean, aes(x=sstmean,y=MDS1)) + 
  geom_smooth(method='lm', se=F, col="black", lwd=0.5 ) + 
  geom_smooth( se=F, col="black", lty=2, lwd=0.5 ) +
  geom_point( col='firebrick', alpha=0.8 ) + 
  theme_classic() + ylab("MDS1") + xlab("Mean annual SST") +
  theme(plot.margin = unit(c(0, 0.75, 0, 0)+0.1, "cm"))


library(cowplot)
windows(7.5,1.9)
plot_grid( d,a,b,c, ncol=4, labels="AUTO", label_size=10, rel_widths = c(1.12,1,1,1)  )




#### eight sites have predation rates >0.25. In Fig. 3 there is separation between these points and others
## what is it about these sites?
highpred <- rate.env %>%
  filter( Site %in% c("FL","NC2","NSW1","NSW2","NC","Italy","QLD3","QLD4") ) 
highcomm <- highpred[,78:138]

highagg  <- aggregate( highcomm, by=list(highpred$Site), sum, na.rm=TRUE )


# what is the proportion of site with zero or near zero predation rate?


## great circle distance for sites with predation rate data
d <- rate.env %>%
  filter( !is.na(sstmean) & !is.na(MDS1) ) %>%
  select( Site, rate, hemi, basin, coast,sstmean, MDS1 ) %>%
  distinct()
dsite <- ddply( d, .(Site,hemi,basin,coast), summarise, 
                rate=mean(rate), MDS1=mean(MDS1), sstmean=mean(sstmean))

dgcd <- left_join( dsite, siteGPS )
dgcd$group <- with(dgcd,  paste(hemi,coast,basin,sep="_") )
dist.list <- split( dgcd, dgcd$group )
dlist <- dist.list[ unlist(lapply( dist.list, function(z) nrow(z)>1 )) ]

# latitudinal signal
ggplot( dgcd, aes(x=abs(meanLat),y=rate,col=group )) + geom_point() + geom_smooth(se=F)


## matrix of geographic distances (gret circle distance) between sites
source( "VincentryInverseFunction.R" )
# matrix to hold all distances
gcd <- matrix( nrow=nrow(dgcd), ncol=nrow(dgcd) )
# loop over all sites to get distances 
for( i in 1:nrow(dgcd)) {
  for( j in 1:nrow(dgcd)) {
    gcd[i,j] <- gcd.vif( dgcd$meanLong[i], dgcd$meanLat[i], dgcd$meanLong[j], dgcd$meanLat[j] )
  }
}

gcd.df <- data.frame( dgcd, round(gcd/1000,2) )
names( gcd.df )[-c(1:11)] <- dgcd$Site

# get distances for  predation
rate.dist <- vegdist( dgcd$rate, "euclid" )

vegan::mantel(rate.dist, as.dist(gcd), permutations = 1000 )
plot( as.vector(rate.dist) ~ as.vector(gcd[lower.tri(gcd)]) )
distdf <- data.frame(rate.diff=as.vector(rate.dist),geodesic=as.vector(gcd[lower.tri(gcd)]))

windows(3,3)
ggplot( distdf, aes(y=rate.diff, x=geodesic)) + geom_point(alpha=0.2, col='slateblue') + 
  geom_smooth(method='glm', method.args=list(family=poisson), col='black') +
  geom_smooth(se=F, col='black', lty=2) +
  theme_classic() +
  ylab("Difference in predation rate") + xlab("Great circle distance (km)")


# run function across list of grouped data
