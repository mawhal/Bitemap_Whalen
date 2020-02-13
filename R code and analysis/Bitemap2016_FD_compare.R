## combine all FD measures from different analyses

# libararies
library( tidyverse )


# read data
fd_pa <- read_csv( "Output Data/FunctionalDiversity_indices_PA.csv")
fd_pa$method <- "presence-absence"
fd_pa_filter <- read_csv( "Output Data/FunctionalDiversity_indices_PA_filtered.csv")
fd_pa_filter$method <- "presence-absence filtered"
fd_abun <- read_csv( "Output Data/FunctionalDiversity_indices.csv" )
fd_abun$method <- "abundance"
fd_abun_length <- read_csv( "Output Data/FunctionalDiversity_indices_length.csv" )
fd_abun_length$method <- "abundance with length"
fd_abun_filter <- read_csv( "Output Data/FunctionalDiversity_indices_filtered.csv" )
fd_abun_filter$method <- "abundance filtered"


fds <- rbind( fd_pa, fd_pa_filter, fd_abun, fd_abun_length, fd_abun_filter )


## read site data
sites <- read_csv( "Output Data/Bitemap_BioORACLE_20190107.csv" )

## select columns
sites <- sites %>%
  select( Country, Site, basin, coast, meanLat, meanLong, chlomax, chlomean, chlomin, chlorange, cloudmax, cloudmean, cloudmin,
          parmax, parmean, salinity, sstmax, sstmean, sstmin, sstrange )

## consumption rate data
rate.env <- read_csv( "Output Data/Bitemap_rate.env.20190423.csv" ) %>% 
  select( Country, Site, habitat, rate )
rate.hab <- rate.env %>% 
  group_by( Country, Site, habitat ) %>% 
  summarize( rate=mean(rate,na.rm=T) )


## merge
ss <- left_join( rate.hab, sites )
fdss <- left_join( fds, ss )




# some plots
ggplot( fdss, aes(x=FRic, y=rate, col=method) ) + geom_point() + geom_smooth(se=F,span=0.66)



# correlations with rate
cor.na <- function(x,y) cor(x,y, use = "complete.obs")
fd.cor <- fdss %>% 
  group_by(method) %>% 
  summarize( FRic=cor.na(FRic,rate), FEve=cor.na(FEve,rate),
             FDis=cor.na(FDis,rate),RaoQ=cor.na(RaoQ,rate),
             FGR=cor.na(FGR,rate))
write_csv( fd.cor, "Output Data/FunctionalDiversity_COMPARE.csv")

fd.cors <- fd.cor %>% 
  gather(key = "index", value="cor",-method)

ggplot(fd.cors, aes(x=method,y=cor)) +facet_wrap(~index) +geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
