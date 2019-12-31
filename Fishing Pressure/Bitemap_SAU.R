# Obtaining Sea Around Us data
# Citation policy: http://www.seaaroundus.org/citation-policy/
# script created 13 December 2019
# by Matt Whalen
# for catch, units should be tonnes (thousands of kilograms)?

library(seaaroundus)

## Bitemap data
# get a list of site locations for all sites
sites <- read.csv( "../Data/Bitemap_sites.csv", stringsAsFactors = F )

# generate small polygons for each set of lats and longs

## function to get cells for each site
bitecells <- function( x, cols=4:3, buffer=0.24 ) {
  buff     <- rbind( x[,cols] + buffer, x[,cols] - buffer )
  poly     <- sp::Polygon(buff) 
  polychar <- paste( poly@coords[,1], poly@coords[,2], 
                 sep=" ", collapse=", " )
  wkt      <- paste0( "POLYGON ((",polychar, "))" )
  return( getcells(wkt) )
}

# test
getcelldata( cells = bitecells( sites[21,] ) )

# apply to every site to get a list of all cells we want to investigate
cells <- lapply( split(sites, seq(nrow(sites)) ), bitecells )

# get all data for each list element
sau <- lapply( cells, function(z) getcelldata(cells=z) )
# rename list elements for easier tracking
names(sau) <- sites$Site
# write this output to disk so we don't have to make multiple API calls
# currently writing this as multiple files to make it easier to track
mapply( write.csv, sau, file=paste0( 'output/', names(sau), '.csv'))


# summarize workflow
library( tidyverse )
testdata %>%
  group_by( cell_id ) %>%
  filter( sector_type_name != "Industrial" ) %>%
  filter( functional_group_name %in% keep ) %>%
  summarise( catch_sum=sum(catch_sum) )
x <- testdata %>%
  group_by( cell_id ) %>%
  filter( functional_group_name %in% keep ) %>%
  summarize( functional_groups=list(unique(functional_group_name)) )
# yikes this makes a big difference cell-to-cell! See India and USA (CA)
sort(unique( testdata$functional_group_name ))
testdata[ testdata$functional_group_name == "demersalmd", ]


# select/omit particular groups (e.g. we don't want large pelagics)
keep <- c( "demersallg", "demersalmd", "demersalsm",
           "benthopelagicmd", "benthopelagicsm",
           "reef-associatedlg", "reef-associatedmd", "reef-associatedsm",
           "pelagicmd", "pelagicsm", "lobsters,crab", "flatfishsm md" )


## for each site, and all cells for that site, find the cell with minimum catch (perhaps this is the most conservative vieW?)
sau %>%
  filter( sector_type_name != "Industrial" )
library( data.table )
saudf <- bind_rows(sau, .id = "site")
# save this to disk
write.csv( saudf, "output/ALL.csv", row.names = FALSE )


cellsum <- saudf %>%
  group_by( site, cell_id ) %>%
  filter( sector_type_name != "Industrial" ) %>%
  filter( functional_group_name %in% keep ) %>%
  summarise( catch_sum=sum(catch_sum) )

cellsummary <- cellsum %>%
  group_by( site ) %>%
  summarise( catch_min = min(catch_sum), catch_mean = mean(catch_sum), catch_max = max(catch_sum) )

# save min catch to disk
write.csv( cellsummary, "output/catchmin.csv", row.names = FALSE )
