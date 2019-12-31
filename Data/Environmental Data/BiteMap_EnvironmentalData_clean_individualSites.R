#############################################################
### MarineGEO Ocean Bitemap Project
### Script to clean environmental data provided by partners
### Focus is on Temperature and Salinity during predation assays
#############################################################


## Notes about environmental data
# most sites provided point estimates of Temperature and Salinity
# from their sites at the times when squid pops were deployed or collected
# so, most of these are day time Temperatures and do not consider variation at diel scales
# so, when summarizing data, consider averaging over days and also hour surrounding predation assays

# Script updates and accomplishments
# 20171208 - Clean Panama Data provided by Janina Seeman


# load libraries
library(tidyverse)
library( ncdf4 )

# read predation assay data - this will help us find dates and times 
p <- read.csv( "../OceanBitemap_Squidpop_Data_20180109_MW.csv" )
# select relevant columns
p <- p %>%
  select( Country, Institution, Name, Contact.Email=Email, Lat, Long, 
          Date=Date.Squidpops.Deployed..yyyymmdd., Time=Time.Squidpops.Deployed..hhmm..24.hour.clock.,
          #dateRetreive=Date.Squidpops.Retrieved..yyyymmdd., timeRetrieve=Time.Squidpops.Retrieved..hhmm..24.hour.clock.,
          Seagrass.Unveg=Type.Of.Habitat )

# convert vegetated and unvegetated sites to common categories
p$Seagrass.Unveg[ p$Seagrass.Unveg %in% c("Seagrass","Seagrass Meadow")] <- "Seagrass"
p$Seagrass.Unveg[ p$Seagrass.Unveg %in% c("Muddy Bottom","Sandy Bottom", "unvegetated ")] <- "unvegetated"
p$Seagrass.Unveg[ p$Seagrass.Unveg %in% c("Artificial Habitat (dock, breakwater, weir, etc.)","Rocky Reef")] <- NA
p <- p[ !is.na(p$Seagrass.Unveg), ]
p <- droplevels(p)
levels(p$Seagrass.Unveg) <- c("Seagrass","Unveg")

# define dates
p$Date <- ymd(p$Date)
# define times
p$Time
# hours
substrStart <- function(x, n){
  substr(x, 1, nchar(x)-n)
}
hr <- substrStart(p$Time,2)
# minutes
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
mn <- substrRight(p$Time,2)
p$Time <- hm( paste(hr,mn,sep=":") )

# Because data from different sites are heterogeneous, handle each one individually



##################################
###### ITALY
# this one I got from AquaMODIS for the days when the experiment was conducted 
# web link: https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/Daily/4km/sst/2016/
files <- list.files( "Derived Data/Italy/", ".nc" )

italy.temp <- c()
for( i in 1:3 ){
  filename <- paste0("Derived Data/Italy/",files[i])
  nc <- nc_open( filename )
  sst <- ncvar_get( nc, "sst" )
  lon <- ncvar_get( nc, "lon" )
  lat <- ncvar_get( nc, "lat" )
  
  # define bounds for temperature estimates
  latu <- 43.3
  latl <- 42.9
  lonu <- 10
  lonl <- 9.8
  lats <- which( lat>latl & lat<latu )
  lons <- which( lon>lonl & lon<lonu )
  
  italy.temp[i] <- mean(sst[ lons, lats ],na.rm=T)
  nc_close( nc )
}


###############################
###### PANAMA
dpan <- read.csv( "Original Files from Partners/Bitemap_PANAMA_EnvData.csv" )
# select relevant columns
dpan <- dpan %>%
  select( date=Date..MM.DD.YYYY., time=Time..HH.MM.SS., Temp=Temp..C, Sal=Sal.psu ) %>%
  filter( Temp > 28.3, Sal>30 )
# define dates
dpan$dt <- ymd_hms( paste(dpan$date, dpan$time) )
dpan$date <- ymd( dpan$date )


## match up times and dates in environmental data
# isolate Site of interest
ppan <- p[ p$Country=="Panama", ]
dayspan <- unique(ppan$Date)


## mean Temp+Sal by day
dpanday <- dpan %>%
  filter( date %in% dayspan ) %>%
  group_by(date) %>%
  summarize( Temp=mean(Temp,na.rm=T), Sal=mean(Sal,na.rm=T) )
dpanday$dt <- ymd_hms( paste(dpanday$date,"00:00:00") )

## mean Temp+Sal by 24 hour period starting at deployment
ppan$dt <- ymd_hms( paste(ppan$Date,ppan$Time) )
panint <- interval( ppan$dt, ppan$dt + dhours(24) )

dspans <- list()
for( i in 1:length(panint) ){
  dspans[[i]] <- dpan %>%
    filter( dt >= int_start(panint[i])-hours(1), dt <= int_end(panint[i]) ) %>%
    summarize( Temp=mean(Temp,na.rm=T), Sal=mean(Sal,na.rm=T) )
}

paninterval <- data.frame( interval=panint, do.call(rbind,dspans) )

# Get Point estimates closest to deployment times
dpoints <- dpan %>%
  filter( round_date( dt, unit="10minutes" ) %in% round_date( ppan$dt, unit="5minutes" ) )
 

## Plot all estimates together
# temperature
ggplot( paninterval, aes(x=int_start(interval),y=Temp) ) + 
  geom_line( data=dpan, aes(x=dt,y=Temp), col='white') +
  geom_point( data=dpoints, aes(x=dt,y=Temp), col='gray' ) +
  geom_point( data=dpanday, aes(x=dt,y=Temp), col='blue' ) +
  geom_point(col='black') +
  xlim( range(c(dpanday$dt,dpoints$dt)) )

# salinity
ggplot( paninterval, aes(x=int_start(interval),y=Sal) ) + 
  geom_line( data=dpan, aes(x=dt,y=Sal), col='white') +
  geom_point( data=dpoints, aes(x=dt,y=Sal), col='gray' ) +
  geom_point( data=dpanday, aes(x=dt,y=Sal), col='blue' ) +
  geom_point(col='black') +
  xlim( range(c(dpanday$dt,dpoints$dt)) ) + ylim(33.75,34)
        
# Using the interval method
paninterval
env.pan <- data.frame( ppan, do.call(rbind,dspans) )
env.pan 

# --- END OF PANAMA ---
###############


###############################
###### USA (VA) Lefcheck site
dva <- read.csv( "Original Files from Partners/Bitemap_Virginia_Lefcheck_CHE019.38_8122016.csv" )

# select relevant columns
dva <- dva %>%
  select( datetime=SAMPLE_DATETIME, Temp=WTEMP, Sal=SALINITY, Depth=TOTAL_DEPTH ) 
# define dates
dva$dt <- mdy_hms( dva$datetime )

## match up times and dates in environmental data
# isolate Site of interest
pva <- p[ p$Country=="USA (VA)", ]
daysva <- unique(pva$Date)


## mean Temp+Sal by day
dvaday <- dva %>%
  filter( date(dt) %in% daysva ) %>%
  # group_by(date(dt)) %>% # only one date
  summarize( Temp=mean(Temp,na.rm=T), Sal=mean(Sal,na.rm=T) )
dvaday$dt <- ymd_hms( paste(daysva,"00:00:00") )

## mean Temp+Sal by 24 hour period starting at deployment
pva$dt <- ymd_hms( paste(pva$Date,pva$Time) )
vaint <- interval( pva$dt, pva$dt + dhours(24) )

dsvas <- list()
for( i in 1:length(vaint) ){
  dsvas[[i]] <- dva %>%
    filter( dt >= int_start(vaint[i])-hours(1), dt <= int_end(vaint[i]) ) %>%
    summarize( Temp=mean(Temp,na.rm=T), Sal=mean(Sal,na.rm=T) )
}

vainterval <- data.frame( interval=vaint, do.call(rbind,dsvas) )

# Get Point estimates closest to deployment times
dpoints <- dva %>%
  filter( round_date( dt, unit="10minutes" ) %in% round_date( pva$dt, unit="5minutes" ) )


## Plot all estimates together
# temperature
ggplot( vainterval, aes(x=int_start(interval),y=Temp) ) + 
  geom_line( data=dva, aes(x=dt,y=Temp), col='white') +
  geom_point( data=dpoints, aes(x=dt,y=Temp), col='gray' ) +
  geom_point( data=dvaday, aes(x=dt,y=Temp), col='blue' ) +
  geom_point(col='black')  

# salinity
ggplot( vainterval, aes(x=int_start(interval),y=Sal) ) + 
  geom_line( data=dva, aes(x=dt,y=Sal), col='white') +
  geom_point( data=dpoints, aes(x=dt,y=Sal), col='gray' ) +
  geom_point( data=dvaday, aes(x=dt,y=Sal), col='blue' ) +
  geom_point(col='black') 

# Using the interval method
vainterval
env.va <- data.frame( pva, do.call(rbind,dsvas) )
env.va 

# --- END OF USA (VA) ---
###############





###############################
###### USA (CA) San Diego
dsd <- read.csv( "Original Files from Partners/Bitemap_USA (CA)_Mission Bay temperature data Bitemap.csv" )

# select relesdnt columns
dsd <- dsd %>%
  select( year=X.YY, month=MM, day=DD, hour=hh, min=mm, Temp=WTMP ) 
# define dates
dsd$dt <- ymd_hm( with(dsd,paste(year,month,day,hour,min)) )

## match up times and dates in environmental data
# isolate Site of interest
psd <- p[ p$Country=="USA (CA)", ]
dayssd <- unique(psd$Date)


## mean Temp+Sal by day
dsdday <- dsd %>%
  filter( date(dt) %in% dayssd ) %>%
  group_by( date(dt) ) %>% 
  summarize( Temp=mean(Temp,na.rm=T) )
dsdday$dt <- ymd_hms( paste(dayssd,"00:00:00") )

## mean Temp+Sal by 24 hour period starting at deployment
psd$dt <- ymd_hms( paste(psd$Date,psd$Time) )
sdint <- interval( psd$dt, psd$dt + dhours(24) )

dssds <- list()
for( i in 1:length(sdint) ){
  dssds[[i]] <- dsd %>%
    filter( dt >= int_start(sdint[i])-hours(1), dt <= int_end(sdint[i]) ) %>%
    summarize( Temp=mean(Temp,na.rm=T) )
}

sdinterval <- data.frame( intersdl=sdint, do.call(rbind,dssds) )

# Get Point estimates closest to deployment times
dpoints <- dsd %>%
  filter( round_date( dt, unit="10minutes" ) %in% round_date( psd$dt, unit="30minutes" ) )


## Plot all estimates together
# temperature
ggplot( sdinterval, aes(x=int_start(intersdl),y=Temp) ) + 
  geom_line( data=dsd, aes(x=dt,y=Temp), col='white') +
  geom_point( data=dpoints, aes(x=dt,y=Temp), col='gray' ) +
  geom_point( data=dsdday, aes(x=dt,y=Temp), col='blue' ) +
  geom_point(col='black')   +
  xlim( range(c(dsdday$dt,dpoints$dt)) ) + ylim(c(20,22))


# Using the intersdl method
sdinterval
env.sd <- data.frame( psd, do.call(rbind,dssds) )
env.sd$Sal <- 35.5 

# --- END OF USA (CA) San Diego ---
###############



###############################
###### INDIA
# NOTE: environmental data do not overlap experimental data
di1 <- read.csv( "Original Files from Partners/Bitemap_India_Agatti_Lagoon-1.csv" )
di1$Site <- "Agatti"
di2 <- read.csv( "Original Files from Partners/Bitemap_India_Kavaratti_Lagoon.csv" )
di2$Site <- "Kavaratti"
dibind  <- rbind( di1, di2 ) 

# select relevant columns
di <- dibind %>%
  select( Site, Date, Time=Time..GMT.05.30, Temp=Temp..Â.C ) %>%
  filter( Temp>27 ) # some temperature data early in the timeseries are suspect
# define dates
di$dt <- dmy_hms( with(di,paste(Date,Time)) )

## match up times and dates in environmental data
# isolate Site of interest
pi <- p[ p$Country=="India", ]
daysi <- unique(pi$Date)


## mean Temp+Sal over time 
ggplot( di, aes(x=dt,y=Temp,col=Site)) + geom_line()


# Mean temp in January from both datasets
di %>%
  filter( month(dt)==1) %>%
  group_by( Site ) %>%
  summarize( Temp=mean(Temp) )

di %>%
  filter( month(dt)==1) %>%
  summarize( Temp=mean(Temp) )

# Best guess Based on the previous year is 29.5 degrees

# --- END OF India ---
###############


######



##############################
# put all of the data together

# read environmental data
directenv <- read.csv( 'BiteMap_EnvironmentalData_Compilation.csv' )
directenv$Date <- ymd( directenv$Date )
directenv$Time <- hms( directenv$Time )
names( directenv )

directenv2 <- full_join( directenv, env.pan )
directenv2 <- full_join( directenv2, env.va )
directenv2 <- full_join( directenv2, env.sd )

# write the data to file
write.csv( directenv2, "BiteMap_EnvironmentalData_Compilation_20171210.csv" )
