#############################################################
### MarineGEO Ocean Bitemap Project
### Script to extract sst data from AquaMODIS for sll sites 
### not all sites measured temperature
#############################################################


## Notes about environmental data
# most sites provided point estimates of Temperature and Salinity
# from their sites at the times when squid pops were deployed or collected
# so, most of these are day time Temperatures and do not consider variation at diel scales
# so, when summarizing data, consider averaging over days and also hour surrounding predation assays

# Script updates and accomplishments
# 20171208 - Clean Panama Data provided by Janina Seeman
# 

# load libraries
library( tidyverse )
library( lubridate )
library( ncdf4 )

# read predation assay data - this will help us find dates and times 
p <- read.csv( "data/raw/squidpop_assays/Bitemap_Squidpop_Data_20190322.csv" )
# select relevant columns
p <- p %>%
  dplyr::select( Country, Institution, Lat, Long, 
          Date=Date.Squidpops.Deployed..yyyymmdd., Time=Time.Squidpops.Deployed..hhmm..24.hour.clock.,
          Seagrass.Unveg=Type.Of.Habitat )

# convert vegetated and unvegetated sites to common categories
p$Seagrass.Unveg[ p$Seagrass.Unveg %in% c("Seagrass","Seagrass Meadow")] <- "Seagrass"
p$Seagrass.Unveg[ p$Seagrass.Unveg %in% c("Muddy Bottom","Sandy Bottom", "unvegetated ","Unveg")] <- "unvegetated"
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


# Get a breakdown of dates for each site - use these to select AquaMODIS data
p <- p %>%
  mutate( yearmonth = paste(year(Date),month(Date)),
          year = year(Date), julian = format(Date,"%j") )

# get unique dates
dayuni <- sort(unique( p$Date ))
data.frame( dayuni, year=year(dayuni), julian=format(dayuni,"%j") )


p %>%
  group_by( Country ) %>%
  summarise( list(unique(Date)) )

with(p, table( Country, yearmonth ) )

p[p$Country=="Italy",]


##################################################################
### for each site (Country), pull the dates, lats, and longs (ranges?)
## add a column that contains year and julian day like the AquaMODIS file names
p$aqua <- with(p, paste0(year,julian))

## 

psum <- p %>%
  group_by( Country, Seagrass.Unveg ) %>%
  summarise( days=list(unique(aqua)), latran=list(range(Lat)), lonran=list(range(Long)) )

psum$days


### get all netcdf files
files <- list.files( "data/Environmental Data/Derived Data/", ".nc" )
# pull out dates from file names
fileday <- substr( files, 2, 8 )


# 
# # test with first site
# test <- psum[1,]
# # get files
# testfile <- files[ fileday %in% unlist(test$days) ]
# # get lats and longs and add a buffer
# buffer <- 0.2
# lonr <- unlist(test$lonran) + c(-buffer,buffer)
# latr <- unlist(test$latran) + c(-buffer,buffer)
# 
# 
# nc <- nc_open( paste0("Derived Data/",testfile[1]) )
# sst <- ncvar_get( nc, "sst" )
# lon <- ncvar_get( nc, "lon" )
# lat <- ncvar_get( nc, "lat" )
# 
# lats <- which( lat>latr[1] & lat<latr[2] )
# lons <- which( lon>lonr[1] & lon<lonr[2] )
# 
# sstchoose <- mean( sst[ lons, lats ], na.rm=T )
# nc_close(nc)
# 
# # apply set of functions over all files
# lapply( testfile, function(z) {
#   # open netCDF file
#   nc <- nc_open( paste0("Derived Data/",testfile[1]) )
#   sst <- ncvar_get( nc, "sst" )
#   lon <- ncvar_get( nc, "lon" )
#   lat <- ncvar_get( nc, "lat" )
#   
#   lats <- which( lat>latr[1] & lat<latr[2] )
#   lons <- which( lon>lonr[1] & lon<lonr[2] )
#   # extract the SST values and take an average
#   sstchoose <- mean( sst[ lons, lats ], na.rm=T )
#   # close netCDF file
#   nc_close(nc)
#   return(sstchoose) # if NA likely means that it was too cloudy that day for good readings
# })





# wrap lapply for netCDF into another lapply function that repeats the extractions for all sites
aqua.all <- apply( psum, 1, function(z, buffer = 0.2) {
  # get file names associated with a site's assays
  fileuse <- files[ fileday %in% unlist(z['days']) ]

  # get lats and longs and add a buffer
  lonr <- unlist(z['lonran']) + c(-buffer,buffer)
  latr <- unlist(z['latran']) + c(-buffer,buffer)
  
  # apply set of functions over all files
  lapply( fileuse, function(z) {
    # open netCDF file
    nc <- nc_open( paste0("Derived Data/",z ) )
    sst <- ncvar_get( nc, "sst" )
    lon <- ncvar_get( nc, "lon" )
    lat <- ncvar_get( nc, "lat" )
    
    lats <- which( lat>latr[1] & lat<latr[2] )
    lons <- which( lon>lonr[1] & lon<lonr[2] )
    # extract the SST values and take an average
    sstchoose <- mean( sst[ lons, lats ], na.rm=T )
    # close netCDF file
    nc_close(nc)
    return(sstchoose) # if NA likely means that it was too cloudy that day for good readings
  })
  
})


# clean up aqua.all to get a nice looking data.frame
aqua.list <- lapply( aqua.all, unlist )




psum$aqua <- unlist( lapply( aqua.list, mean, na.rm=T ) )
psum$aqua[ psum$aqua=="-Inf"] <- NA


##############################
# put all of the data together

# read environmental data
directenv <- read.csv( 'data/Environmental Data/BiteMap_EnvironmentalData_Compilation_20180313.csv' )
# directenv$Date <- ymd( directenv$Date )
# directenv$Time <- hms( directenv$Time )
names( directenv )

directenv2 <- full_join( directenv, psum )

directenv3 <- directenv2 %>%
  group_by( Country, Seagrass.Unveg ) %>%
  summarise( Temp=mean(Temp,na.rm=T), aqua=mean(aqua,na.rm=T) )

ggplot( directenv3, aes(x=aqua,y=Temp)) + geom_point() + geom_smooth(method='lm')

lm1 <- lm ( Temp ~ aqua, data=directenv3 )
summary(lm1)

directenvsave <- directenv2 %>%
  select( Site.Name, Country, Lat, Long, Seagrass.Unveg, Date, Time, Temp, Sal, Measurment.Number, Data.Collector, Notes, Participant.List, Contact.Email, aqua )

# write the data to file
write.csv( directenvsave, "BiteMap_EnvironmentalData_Compilation_20190228.csv" )
