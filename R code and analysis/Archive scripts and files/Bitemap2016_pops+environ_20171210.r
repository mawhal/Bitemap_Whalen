#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Preliminary data analysis and exploration
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2017.10.23
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

###   REMOVE NSW1.1 AND 1.2 BECAUSE THESE WERE NOT CONDUCTED IN SEAGRASS

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
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
# geospatial data
library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(rfishbase)
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
pops   <- read.csv( '../Data/OceanBitemap_Squidpop_Data_20180122.csv' )
# rename pops columns
names(pops) <- c("timeStamp","name","email","Lat","Long","dateRetrieved","timeRetrieved",
                 "dateDeployed","timeDeployed","habitat","habitatDescription","numberDeployed","missing1",
                 "numberRetrieved","missing24","Country","Institution","prop1","prop24","notes")
# convert vegetated and unvegetated sites to common categories
pops$habitat[ pops$habitat %in% c("Seagrass","Seagrass Meadow","seagrass")] <- "Seagrass"
pops$habitat[ pops$habitat %in% c("Muddy Bottom","Sandy Bottom", "unvegetated ")] <- "unvegetated"
pops$habitat[ pops$habitat %in% c("Artificial Habitat (dock, breakwater, weir, etc.)","Rocky Reef")] <- NA
pops <- pops[ !is.na(pops$habitat), ]
pops <- droplevels(pops)
levels(pops$habitat) <- c("Seagrass","Unveg")
# coerce proportion after 24 hours to be a numbers
pops$prop24[ pops$prop24=="#DIV/0!" ] <- NA
pops$prop24 <- as.numeric(as.character(pops$prop24))
# factor for hemisphere
pops$hemi <- factor(ifelse( pops$Lat>0,"North","South"))

# FISH SEINING DATA
seines <- read.csv( '../Data/Bitemap_Seine_ALL-DATA_20180122.csv', strip.white = TRUE)
names(seines)[4] <- "habitat"
# convert all vegetated and unvegetated sites to common categories
seines$habitat[ seines$habitat %in% c("Seagrass " )] <- "Seagrass"
seines$habitat[ seines$habitat %in% c("unveg","Unveg","Unvegetated" )] <- "Unveg"
seines <- droplevels(seines)
# for all organisms from France, UNC2, Wales, QLD2, QLD3, multiple lengths by 10 to convert from cm to mm
seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] <- 
  seines$Length[ seines$Country %in% c('France', 'USA (NC2)', 'Wales', 'Australia (QLD2)', 'Australia (QLD3)')] * 10


# merging data
# this is difficult because the replicate seines do not pair up with squidpop reps
# so, in order to merge the two we need to use summaries at the site level
# try to match things based on Lat/Long info
popGPS   <- unique( pops[,c('Lat','Long','Institution','Country')])
seineGPS <- unique( seines[,c('Lat','Long','Country','Site.Name','habitat')])
# compare sites used in seining and predation assays
GPSjoin <- full_join( popGPS,seineGPS )
GPSjoin <- GPSjoin[ with(GPSjoin, order(Lat,Long)), ]
# write.csv( GPSjoin, "Output Data/pops_sienes_GPSmatch.csv", row.names=F )
# to complete merge, we need to average predation and fish metrics at site level

# combine all squidpop and seine replicates from a given location (e.g. NSW2)
siteGPS <- ddply( GPSjoin, .(Country), 
                  summarise, meanLat=mean(Lat), meanLong=mean(Long) )
# for each site within a country
siteGPS.site <- ddply( GPSjoin, .(Country,Site.Name), 
                  summarise, meanLat=mean(Lat), meanLong=mean(Long) )
# write.csv( siteGPS, "Bitemap_sites.csv", row.names=FALSE )






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
files <- list.files( "C:/Users/mawha/Dropbox/Global Databases/BioOracle Data/",pattern = ".asc")
# # read in the raster data files (will assume lat,long for projection)
# r <- lapply( paste0("C:/Users/mawha/Dropbox/Global Databases/BioOracle Data/",files), raster )
# # crop all rasters to same extent, keeping southern hemisphere
# e <- extent(-180,180,-70,70)
# r2 <- lapply( r, function(rast) crop(rast,e) )

########### WorldClim Precipitation ###############
# p <- raster( "C:/Users/mawha/Dropbox/Global Databases/WorldClim/WorldClim_precip_2-5.tif" )

################## input data ######################
# extract lat and long from input
# input <- siteGPS[,c("meanLong","meanLat")]


############ extract environmental data #############
# # average values in all raster cells within a given radius from the GPS pointS
# # Use a "buffer" or radius over which to look for raster cells surrounding each GPS point
# # that has data (note there is no data on land or freshwater for Bio-ORACLE, and no ocean data for WorldClim
# # This will take a while for many variables
# buffer <- 10000  # this is in meters if the map is projected correctly
# # if raster coordinate reference system (CRS) is undefined, it will assume lat/long, which is correct in this case


# use the buff function (defined above) to all rasters
# oracle <- lapply( r2, buff ) # this takes a very long time, average Lat and Long by site first
#   write.csv( oracle, "Output Data/Bitemap_BioOracle.csv")
oracle <- read.csv( "Output Data/Bitemap_BioOracle.csv" )[,-1] 
# precip <- buff( p )
#   write.csv( precip, "Output Data/Bitemap_WorlClimPrecip.csv")
precip <- read.csv( "Output Data/Bitemap_WorlClimPrecip.csv" )[,-1]

# MISSING VALUES FROM Bio-ORACLE
# 7, 110, 130
# GPSjoin[ c(7,110,130), ]
# no data for 
# New South Wales, NSW1.1
# University of Delaware
# Croatia

# combine all of these into a data.frame and give them names
Environmentals <- data.frame( do.call( cbind, oracle ) )
names(Environmentals) <- unlist( strsplit( files, ".asc") )
Environmentals$precip <- precip

sites <- cbind(siteGPS,Environmentals)
# factor for hemisphere
sites$hemi <- factor(ifelse( sites$meanLat>0,"North","South"))


with(sites, plot( sstmean~abs(meanLat), col=hemi) )
with(sites, plot( sstmean~(meanLat)) )

with(sites, plot( sstmin~abs(meanLat)) )
with(sites, plot( sstmax~abs(meanLat)) )
with(sites, plot( sstrange~abs(meanLat)) )
with(sites, plot( sstrange~sstmean,
     xlab="Mean annual SST", ylab="Range annual SST") )
with(sites, plot( sstrange~sstmin) )
plot( sstmin~sstmean, data=sites, type='n', ylab="Min or max SST", xlab="Mean annual SST",
      ylim=c(0,35))
  points( sstmax~sstmean, data=sites, col='darkorange' )
  points( sstmin~sstmean, data=sites, col='blue' )
  abline( lm(sstmax~sstmean, data=sites), col='darkorange' )
  abline( lm(sstmin~sstmean, data=sites), col='blue' )
  
ggplot( data=sites, aes(x=abs(meanLat),y=sstmean)) + facet_wrap(~hemi) + geom_point()
ggplot( data=sites, aes(x=sstmean,y=sstrange)) + facet_wrap(~hemi) + geom_point()


# read temperature and salinity data collected during predation assays
directenv <- read.csv( '../Data/Environmental Data/BiteMap_EnvironmentalData_Compilation_20180109.csv' )
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

# names(pops)[8] <- "Date"
popenv <- left_join(pops,directenvsite)
popall <- left_join( popenv, sites )

# restrict dataset to rows with both in situ and satellite measurements
popclean <- popall[ !is.na(popall$temp) & !is.na(popall$sstmean), ]
popclean$abLat <- abs(popclean$Lat)
# which sites do not have in situ or satellite measurements
unique( popall$Country[ is.na(popall$temp) ] )
unique( popall$Country[ is.na(popall$sstmean) ] )



##### Look at raw predation rates bewtween habitats at 1 and at 24 hours across all sites
ggplot( popall, aes(x=as.numeric(habitat),y=prop1)) + geom_point() + facet_wrap(~Country) +
  geom_smooth( method='glm', formula=y~as.numeric(x), method.args=list(family=binomial), col='red' ) +
  geom_smooth( method='lm', formula=y~as.numeric(x), se=F, alpha=0.5 )

ggplot( popall, aes(x=as.numeric(habitat),y=prop24)) + geom_point() + facet_wrap(~Country) +
  geom_smooth( method='glm', formula=y~as.numeric(x), method.args=list(family=binomial), col='red' ) +
  geom_smooth( method='lm', formula=y~as.numeric(x), se=F, alpha=0.5 )

##### Average data by Country
popsite <- ddply( popall, .(Country), summarise, Lat=mean(Lat), 
                  sstmean=mean(sstmean), temp=mean(temp), salinity=mean(salinity), precip=mean(precip),
                  prop24=mean(prop24) )

# apply a site label
country_split <- strsplit( popsite$Country, split = '[()]' )
popsite$site_label <- unlist(lapply( country_split, function(z){
  if(length(z)==2) z[2] else z[1]
}
))

# which sites had crabs present
crabs <- read.csv( "Output Data/crabsum.csv", stringsAsFactors = FALSE )
names(crabs) <- c("Country","Total")

# merge crab data with predation data
popcrab <- left_join( popsite, crabs )
popcrab$Total[ is.na(popcrab$Total) ] <- 0
popcrab$crabs <- ifelse( popcrab$Total>0,1,0)
popcrab$crabs <- factor(popcrab$crabs,labels=c("no","yes"))

# Crab presence-absence
ggplot( popcrab, aes(x=sstmean,y=prop24, label=site_label, col=crabs) ) + 
  geom_smooth(aes(group=1), method='glm', formula=y~poly(x,2),method.args=list(family=quasibinomial),se=FALSE) +
  geom_text(size=2.5) +  
  scale_color_manual(values=c("black","red")) +
  ylab("Proportion squid missing after 24hr") + xlab("Annual mean SST") +
  guides(color = guide_legend(override.aes = list(linetype = 0)))

ggplot( popcrab, aes(x=temp,y=prop24, label=site_label, col=crabs) ) + 
  geom_smooth(aes(group=1), method='glm',method.args=list(family=quasibinomial),se=TRUE) +
  geom_text(size=2.5) +  
  scale_color_manual(values=c("black","red")) +
  ylab("Proportion squid missing after 24hr") + xlab("in situ temperature") +
  guides(color = guide_legend(override.aes = list(linetype = 0)))


# Total crab abundance
ggplot( popcrab, aes(x=sstmean,y=prop24, label=site_label, size=Total) ) + 
  geom_smooth(aes(group=1), method='glm', formula=y~poly(x,2),method.args=list(family=quasibinomial),se=FALSE) +
  geom_point() +  
  scale_color_manual(values=c("black","red")) +
  ylab("Proportion squid missing after 24hr") + xlab("Annual mean SST") +
  guides(color = guide_legend(override.aes = list(linetype = 0)))

ggplot( popcrab, aes(x=temp,y=prop24,  col=crabs) ) + 
  geom_smooth( method='glm',method.args=list(family=quasibinomial),se=TRUE) +
  geom_point() +  
  scale_color_manual(values=c("slateblue","black")) #+
  ylab("Proportion squid missing after 24hr") + xlab("in situ temperature") +
  # guides(color = guide_legend(override.aes = list(linetype = 0)))





##### look at predation data without averaging at site level
ggplot( popclean, aes(x=temp,y=prop24) ) + geom_point() +  
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + facet_wrap(~habitat)

ggplot( popclean, aes(x=temp,y=prop24) ) + geom_point() +  
  geom_smooth(method='glm',method.args=list(family=quasibinomial)) + facet_grid(hemi~habitat)

ggplot(pops, aes(x=Lat, y=prop24 )) + geom_point() + facet_grid(hemi~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial))  + 
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab('Latitude') 

ggplot(pops, aes(x=Lat, y=prop24,group=hemi )) + geom_point() + facet_wrap(~habitat) +
  geom_smooth(method='glm',method.args=list(family=quasibinomial))  + 
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab('Latitude') 

## model temp replicates using mixed effects model
tempm1 <- glmer( prop24~temp + (1|Country), popclean, family='binomial' )
plot( y=resid(mtemp), x=abs(popclean$Lat[!is.na(popclean$temp)]) )

tempm2 <- glmer( prop24~temp+habitat + (1|Country), popclean, family='binomial' )
summary(tempm2)
vif.mer(tempm2)

tempm3 <- glmer( prop24~temp*habitat+I(temp^2)*habitat + (1|Country), popclean, family='binomial' )
summary(tempm3)
vif.mer(tempm3)

tempm4 <- glmer( prop24~temp+I(temp^2) + (1|Country), popclean, family='binomial' )
summary(tempm4)
vif.mer(tempm4)

tempm5 <- glmer( prop24~temp*habitat + (1|Country), popclean, family='binomial' )
summary(tempm5)
vif.mer(tempm5)

# model comparison
AICctab( tempm1,tempm2,tempm5, nobs=nrow(popclean[!is.na(popclean$temp),]) )



## responses to annual temperature


ggplot(popclean, aes(x=sstmean, y=prop24 )) + geom_point() + facet_wrap(~habitat,ncol=2) +
  geom_smooth(method='glm',  method.args=list(family=quasibinomial) ) + 
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab('in situ temperature') 


## model sst replicates using mixed effects model
sstm2 <- glmer( prop24~sstmean+habitat + (1|Country), popclean, family='binomial' )
summary(sstm2)
summary( glm( prop24~temp+habitat, popclean, family='binomial' ) )

sstm3 <- glmer( prop24~sstmean*habitat+I(sstmean^2)*habitat + (1|Country), popclean, family='binomial' )
summary(sstm3)
vif.mer(sstm3)


sstm4 <- glmer( prop24~sstmean + (1|Country), popclean, family='binomial' )
summary(sstm4)
vif.mer(sstm4)

sstm5 <- glmer( prop24~sstmean+sstrange + (1|Country), popclean, family='binomial' )
summary(sstm5)
vif.mer(sstm5)

sstm6 <- glmer( prop24~sstmin*sstrange + (1|Country), popclean, family='binomial' )
summary(sstm6)
vif.mer(sstm6)

sstm7 <- glmer( prop24~sstmin:sstrange + (1|Country), popclean, family='binomial' )
summary(sstm7)
vif.mer(sstm7)

sstm8 <- glmer( prop24~sstmean+I(sstmean^2) + (1|Country), popclean, family='binomial' )
summary(sstm8)
vif.mer(sstm8)

m0 <- glmer( prop24~1 + (1|Country), popclean, family='binomial' )
summary(m0)

AICctab( m0,sstm2,sstm3,sstm4,sstm5,sstm6,sstm7,sstm8, nobs=nrow(popclean[!is.na(popclean$sstmean),]))
 

# effects of latitude and temp
latm1 <- glmer( prop24~abLat + (1|Country), popclean, family='binomial' )
summary(latm1)

latm2 <- glmer( prop24~abLat+abLat:temp + (1|Country), popclean, family='binomial' )
summary(latm2)
vif.mer(latm2)

latm3 <- glmer( prop24~abLat*sstmean+ (1|Country), popclean, family='binomial' )
summary(latm3)
vif.mer(latm3)

latm4 <- glmer( prop24~abLat+abLat:sstmin + (1|Country), popclean, family='binomial' )
summary(latm4)
vif.mer(latm4)

latm5 <- glmer( prop24~abLat+abLat:sstrange + (1|Country), popclean, family='binomial' )
summary(latm5)
vif.mer(latm5)

latm6 <- glmer( prop24~abLat+I(abLat^2)+ (1|Country), popclean, family='binomial' )
summary(latm6)
vif.mer(latm6)

latm7 <- glmer( prop24~abLat+abLat:sstmax + (1|Country), popclean, family='binomial' )
summary(latm7)
vif.mer(latm7)

latm8 <- glmer( prop24~abLat+hemi + (1|Country), popclean, family='binomial' )
summary(latm8)
vif.mer(latm8)

# compare latitude models
AICctab( latm1,latm2,latm3,latm4,latm5,latm6,latm7,latm8, nobs=nrow(popclean[!is.na(popclean$temp),]))


# compare all models
AICctab( m0,latm1,latm2,latm3,sstm6,sstm3,sstm8,sstm2,tempm2,tempm1, nobs=nrow(popclean) )





### model predictions

# Annual mean sea surface temperature with quadratic
pred.frame <- with(popclean, expand.grid( sstmean=seq(min(sstmean,na.rm=T),max(sstmean,na.rm=T),len=20) ) )
X <- model.matrix(~sstmean+I(sstmean^2),data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
pred <- data.frame(pred.frame,mean=(X%*%fixef(sstm8)))  # these are the point predictions
# Calculate variance between observations within each site (for each combination of fixed effects)
V <- vcov(sstm8)
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T

# Attach to dataframe, calculate standard errors and confidence intervals (using 1.96*sigma may be anti-conservative...)
predictions <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictions <- with(predictions, data.frame(predictions,
                                            p1.lo = mean-1.96*pred.se1,
                                            p1.hi = mean+1.96*pred.se1))

ggplot(popclean, aes(x=sstmean, y=prop24 )) + geom_point(alpha=0.2,col='slateblue') + facet_wrap(~habitat,ncol=2) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(mean) )) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.lo)),lty=2) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.hi)),lty=2) +
  # geom_smooth( method='glm', method.args=list(family=binomial), se=F ) +
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab(expression(paste('Annual mean temperature (',degree,'C)'))) + 
  theme_bw() + theme( panel.grid = element_blank()  )

ggplot(popclean, aes(x=sstmean, y=prop24 )) + geom_point(alpha=0.2,col='slateblue') + 
  geom_line( data=predictions, aes(x=sstmean,y=logistic(mean) )) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.lo)),lty=2) +
  geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.hi)),lty=2) +
  # geom_smooth( method='glm', method.args=list(family=binomial), se=F ) +
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab(expression(paste('Annual mean temperature (',degree,'C)'))) + 
  theme_bw() + theme( panel.grid = element_blank()  )




# Latitude with sstmean
pred.frame <- with(popclean, expand.grid( abLat=seq(min(abLat),max(abLat),len=20), 
                                          # sstmean=seq(min(sstmean),max(sstmean),len=5) ) )
                                          sstmean=mean(sstmean) ) )
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

ggplot(popclean, aes(x=abLat, y=prop24 )) + geom_point(alpha=0.2,col='slateblue') +  
  geom_line( data=predictions, aes(x=abLat,y=logistic(mean), group=sstmean )) +
  geom_line( data=predictions, aes(x=abLat,y=logistic(p1.lo),group=sstmean),lty=2) +
  geom_line( data=predictions, aes(x=abLat,y=logistic(p1.hi),group=sstmean),lty=2) +
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab("Degrees from equator") + 
  theme_bw() + theme( panel.grid = element_blank()  )


# Latitude alone
pred.frame <- with(popclean, expand.grid( abLat=seq(min(abLat,na.rm=T),max(abLat,na.rm=T),len=20) ))
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

ggplot(popclean, aes(x=abLat, y=prop24 )) + geom_point(alpha=0.2,col='slateblue') + 
  geom_line( data=predictions, aes(x=abLat,y=logistic(mean))) +
  # geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.lo)),lty=2) +
  # geom_line( data=predictions, aes(x=sstmean,y=logistic(p1.hi)),lty=2) +
  # geom_smooth( method='glm', method.args=list(family=binomial), se=F ) +
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab("Degrees from equator") + 
  theme_bw() + theme( panel.grid = element_blank()  )



# in situ temperature and habitat
pred.frame <- with(popclean, expand.grid( temp=seq(min(temp,na.rm=T),max(temp,na.rm=T),len=20), 
                                        habitat=unique(habitat) ) )
X <- model.matrix(~temp+habitat,data=pred.frame) # calculate model matrix (formula needs to be same as the model fitted)
pred <- data.frame(pred.frame,mean=(X%*%fixef(tempm2)))  # these are the point predictions
# Calculate variance between observations within each site (for each combination of fixed effects)
V <- vcov(tempm2)
# # model predictions
# pred.frame <- with(popclean, expand.grid( temp=seq(min(temp,na.rm=T),max(temp,na.rm=T),len=20), habitat=unique(habitat) ) )
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

ggplot(popclean, aes(x=temp, y=prop24 )) + geom_point(alpha=0.2,col='slateblue') + facet_wrap(~habitat,ncol=2) +
  geom_line( data=predictions, aes(x=temp,y=logistic(mean) )) +
  geom_line( data=predictions, aes(x=temp,y=logistic(p1.lo)),lty=2) +
  geom_line( data=predictions, aes(x=temp,y=logistic(p1.hi)),lty=2) +
  # geom_smooth( method='glm', method.args=list(family=binomial), se=F ) +
  ylab('Predation intensity\n(prop. squid missing after 24hr)') +
  xlab(expression(paste('in situ temperature (',degree,'C)'))) + 
  theme_bw() + theme( panel.grid = element_blank()  )



