#########################
## Bitemap Author List
## merge two documents containing author information to compile a complete list
## identify missing information for authors
##
## Matt Whalen
## 20 March 2019
#############################


# load libraries
library( tidyverse )


# read data
d1 <- read.csv( "bitemap_authors.csv", stringsAsFactors = F )
d2 <- read.csv( "Bitemap Authors_googledoc.csv", stringsAsFactors = F )




names(d1)
names(d2)

# merge the data
d <- full_join( d1,d2, by=c("last","first"))

# write to disk
write.csv( d, "bitemap_authors_merge.csv", row.names=F )
