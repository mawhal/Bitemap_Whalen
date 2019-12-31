#########################
## Bitemap Author List
## manipulate the current co-author list from google drive
##
## Matt Whalen
## 20 March 2019
#############################


# load libraries
library( tidyverse )


# read data
d <- read.csv( "bitemap_authors_email_orcid_20190424.csv", stringsAsFactors = FALSE)



# select the columns we want
d <- d %>%
  select( name, email )


# split up the names to extract first and last names only
nsplit <- strsplit( d$name, " " )

# first name
nfirst <- unlist( lapply( nsplit, function(z) z[1] ) )

# last name
nwhich <- unlist( lapply( nsplit, length ) )
nlast <- mapply( function(x,y) x[y], nsplit, nwhich )

# put it together
science <- data.frame( nfirst, nlast, d$email )

# write to disk
write.csv( science, "../Manuscript/Science submission/author_list_prep.csv", row.names = FALSE )
