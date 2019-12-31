##########################################################
### Reef Life Survey Traits
### Based on work from Rick Stuart-Smith and Graham Edgar
### assess trait types from their large list
### list provided to Matt Whalen by Ross Whippo
###

### read the data
d <- read.csv( "Traits_all-species_edit.csv" )


# what are the unique options for each traits?
names(d)

x <- apply( d[,9:12], 2, unique )
x
data.frame(x$Trophic.group)
