#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Cleaning fish community data, calculating summaries
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 17 August 2019
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 



########################################################################################
########################################################################################



# library(profvis)  - can look at places where the code is slow, bottlenecks, etc.
# profvis({

# load libraries
# data.frame manipulation and joining
library(tidyverse)



## read data

# seine abundance + biomass
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )


# environmental data from remote sensing and oceanographic expeditions
oracle <- read_csv( "Output Data/Bitemap_BioORACLE_20190107.csv")[,1:36]

# squidpop consumption rates
pop <- read_csv( "Output Data/Bitemap_rate.env.20190423.csv" )
pop <- pop %>% 
  dplyr::group_by( Country,Site,habitat ) %>% 
  dplyr::summarise( sstmean=mean(sstmean), temp=mean(temp), rate=mean(rate) )

# trait information that includes video data
video <- read_csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv" )
video <- video %>% 
  select( Country, habitat, Genus ) %>% 
  mutate( habitat=tolower(habitat) ) %>% 
  distinct()
video$habitat[ video$habitat == "seagrass" ] <- "Seagrass"
video$habitat[ video$habitat == "unveg" ] <- "Unvegetated"
video$Genus[ video$Genus=="Athrina" ] <- "Atherina"

genfam <- seine %>% select(family,Genus) %>% distinct()
video1 <- left_join(video, genfam)

# lookup taxa with taxize
taxa <- video1 %>% 
  filter( is.na(family) ) %>% 
  select( Genus ) %>% 
  distinct() %>%
  arrange( Genus ) %>% 
  unlist()
taxa <- na.omit(taxa)
taxa <- taxa[ !(taxa %in% c("NO ID","Not identified")) ]
# do the lookup
# library(taxize)
# options(ENTREZ_KEY = "a06da45de96c4b7d680015d4c5b46694d908")
# taxclass <- taxize::classification( taxa, db="ncbi", 
#                                     ENTREZ_KEY = "a06da45de96c4b7d680015d4c5b46694d908" )
# taxfam <- rbind(taxclass) %>% 
#   filter(rank=='family') %>% 
#   select( family=name, Genus=query )
# # write taxfam to disk so we don't have to look it up every time
# write_csv( taxfam, "Output Data/families_taxize_RDA.csv" )
taxfam <-  read_csv( "Output Data/families_taxize_RDA.csv" )

famfam <- bind_rows( genfam, taxfam )
video <- left_join(video, famfam)
video <- filter( video, !is.na(family) )


## summarize seine data
# sum of catch per seine
seine.totals <- seine %>% 
  dplyr::group_by(Country,Site.Name,Lat,Long,habitat,Date,Time,SPECIES_NAME) %>% 
  dplyr::summarize( cpua=sum(cpua,na.rm=T) )
# because we want taxon specific catch we need to include zeros when the taxon was not caught
seine.all <- seine.totals %>% 
  dplyr::group_by(Country,Site.Name,Lat,Long,habitat,Date,Time) %>% 
  tidyr::spread( SPECIES_NAME, cpua, fill=0 ) %>% 
  tidyr::gather( SPECIES_NAME, cpua, -Country,-Site.Name,-Lat,-Long,-habitat,-Date,-Time )
# avearage by site and habitat
mean.catch <- seine.all %>% 
  dplyr::group_by(Country,Site.Name,habitat,SPECIES_NAME) %>% 
  dplyr::summarize( cpua=mean(cpua) )
  
# spread out seine data
mean.catch %>% 
  dplyr::group_by(Country,Site.Name,habitat) %>% 
  tidyr::spread( SPECIES_NAME, cpua, fill=0 )

# merge conusmer data with site information and taxonomy
tax <- seine %>% 
  select( phylum, family, SPECIES_NAME ) %>% 
  distinct()
mean.catch.tax <- left_join( mean.catch, tax )

site <- select( oracle, Site, Country, meanLat, meanLong, chla=chlomean )

mean.catch.site <- right_join( site, mean.catch.tax )

## summarize data by family
catch.fam <- mean.catch.site %>% 
  dplyr::group_by( Site, habitat, family ) %>% 
  dplyr::summarize( P=sum(cpua), meanLat=mean(meanLat), meanLong=mean(meanLong) ) %>% 
  dplyr::mutate( P=ifelse(P>0,1,0) )


# merge with video data
video <- right_join(site,video)
video <- select(video, Site, habitat, meanLat, meanLong, family)
video <- video %>% filter( !is.na(family) )
video$P <- 1

catch.video <- bind_rows( catch.fam, video )
catch.video <- catch.video %>% 
  distinct()

# because we are merging seine and video data in some cases, we get discrepancies
# for example: Sparids were caught on video but not in the seine in Croatia
catch.video %>% 
  dplyr::group_by(Site,habitat,family,meanLat,meanLong) %>% 
  dplyr::summarise( count=length(P) ) %>% 
  filter( count>1 )

# since there are only ones and zeros and no repeats, we can add up the presence values
catch.video <- catch.video %>% 
  dplyr::group_by(Site,habitat,meanLat,meanLong,family) %>% 
  dplyr::summarise( P=sum(P) ) 

# spread out again
fam.catch.wide <- catch.video %>% 
  group_by( Site, habitat,meanLat,meanLong ) %>% 
  spread( family, P, fill=0 )


# can average across habitat types, too

fam.data <- fam.catch.wide[,-c(1:4)]

fam.meta <- fam.catch.wide[,1:4]
fam.meta <- left_join( fam.meta, pop )

fam.meta <- fam.meta %>%  tidyr::unite( SH, Site, habitat, remove=FALSE )

# rownames(fam.data) <- fam.meta$SH

# write to disk
fam.all <- bind_cols(fam.meta,fam.data)
write_csv( fam.all, "Output Data/consumer_presence_wide.csv")

## add ordination results, traits, functional diversity, unfiltered abundance
rdaunc <- read_csv( "Output Data/multivar_unconstr_sites.csv" )
rdaunc <- rdaunc %>% separate( SH, c("Site", "habitat") )
brda <- read_csv( "Output Data/biomass_RDAselected.csv" ) # 30 sites
brda <- brda %>%
  select( Site, habitat, biomass, cpua )
active <- read_csv( "Output Data/consumer_active_ratio_PA.csv" )
active <- active %>% 
  select( Site, habitat, act.ratio=act.ratio.ind )
fds  <- read_csv( "Output Data/FunctionalDiversity_indices_PA.csv" ) # fewer sites
fds <- as.data.frame(fds)
fds[is.na(fds)] <- 0
cwm <- read_csv( "Output Data/FunctionalDiversity_CWM_PA.csv" )
cwm_length <- read_csv( "Output Data/FunctionalDiversity_CWM_length.csv")
cwm_length <- cwm_length %>% select( Site, habitat, length )
# abundance
meancpua <- mean.catch %>% 
  dplyr::group_by( Country, habitat, Site.Name ) %>% 
  dplyr::summarize( tot.cpua = sum(cpua,na.rm=T)) %>% 
  dplyr::group_by( Country, habitat ) %>% 
  dplyr::summarize( mean.cpua = mean(tot.cpua,na.rm=T)) %>% 
  dplyr::ungroup()
meancpua <- left_join( meancpua, select(site, Site,Country) )
meancpua <- meancpua %>% 
  select( Site, habitat, mean.cpua ) %>% 
  filter( complete.cases(meancpua) )


addon <- full_join( full_join(full_join(full_join(full_join(full_join(rdaunc,brda),active),fds),cwm), cwm_length), meancpua)
addon <- addon 



rate.mean <- left_join(fam.meta, addon)
rate.mean <- left_join(rate.mean, select(site,Country,Site,chla) )

rate.mean.pairs <- rate.mean %>% 
  mutate( log.cpua = log10(cpua+0.001),log.mean.cpua=log10(mean.cpua+0.001),
          logit.rate = car::logit(rate),abLat=abs(meanLat) ) %>% 
  select( Site, habitat,`degrees\nfrom\nequator`=abLat,meanLong, `mean\nannual\nSST`=sstmean, 
          PCoA1=MDS1, PCoA2=MDS2,chla,
          `proportion\nactive\nforagers`=act.ratio, `selected\nabundance`=log.cpua, 
          `consumption\nrate`=logit.rate, FRic, act, length, `total\nabundance`=log.mean.cpua )

# bivariate correlations
chart.Correlation <- source("chart.correlation.r")$value

dim1 = 8.66142
windows(dim1/2,dim1/2)
svg( "Fig3b.svg", width = dim1/2, height = dim1/2 )

chart.Correlation( rate.mean.pairs[,c("mean\nannual\nSST",
                                      "proportion\nactive\nforagers",
                                      "FRic",
                                      "PCoA1",
                                      "selected\nabundance",
                                      "consumption\nrate") ],
                   histogram = F, method="spearman")
dev.off()


## bigger one for the supplement
# inlcude Latitude, total abundance, FRic based on video data, length, CWMs from FD
# also include total abundance from seine data
windows(8,8)
chart.Correlation( rate.mean.pairs[,c("degrees\nfrom\nequator",
                                      "mean\nannual\nSST",
                                      "chla",
                                      "total\nabundance",
                                      "length",
                                      "proportion\nactive\nforagers",
                                      "FRic",
                                      "PCoA1","PCoA2",
                                      "selected\nabundance",
                                      "consumption\nrate") ],
                   histogram = F, method="spearman" )

rate.mean.pairs$color <- ifelse( rate.mean.pairs$`consumption\nrate`>1,"black","gray" )
pairs( rate.mean.pairs[,c("degrees\nfrom\nequator",
                                      "mean\nannual\nSST",
                                      "PCoA1","PCoA2",
                                      "consumption\nrate") ],
                   histogram = F, method="spearman", col=rate.mean.pairs$color, xpd=T )


# look at realtionships with CWMs
rate.cwm <- rate.mean
theme_set( theme_classic())
# reorder factors
rate.cwm$troph <- factor(rate.cwm$troph, levels=c("omnivore","benthic invertivore;planktivore","benthic invertivore","benthic invertivore;higher carnivore","higher carnivore"))
rate.cwm$act <- ifelse( rate.cwm$act==0,1,0)
# rate.cwm <- rate.cwm %>% filter( !is.na(meanLat) )
ggplot( rate.cwm, aes(y=rate,x=act) ) + geom_point() +geom_smooth(method='lm',se=T)
# ggplot( rate.cwm, aes(y=rate),x=length) ) + geom_point() +geom_smooth(se=T)
a <- ggplot( rate.cwm, aes(y=rate,x=reorder(feed, rate, FUN = mean)) ) + geom_boxplot() + geom_point() +
  xlab("Feeding mode") + ylab("Consumption rate") + coord_flip()
  # theme( axis.text.x = element_text(angle=-45,hjust=0))
b <- ggplot( rate.cwm, aes(y=rate,x=reorder(troph, rate, FUN = median)) ) + geom_boxplot() + geom_point() +
  xlab("Trophic group") + ylab("Consumption rate") + coord_flip()
  # theme( axis.text.x = element_text(angle=-45,hjust=0))
c <- ggplot( rate.cwm, aes(y=rate,x=reorder(watercol, rate, FUN = mean)) ) + geom_boxplot() + geom_point() +
  xlab("Water column use") + ylab("Consumption rate") + coord_flip()
  # theme( axis.text.x = element_text(angle=-45,hjust=0))
d <- ggplot( rate.cwm, aes(y=rate,x=reorder(body, rate, FUN = mean)) ) + geom_boxplot() + geom_point() +
  xlab("Lateral body shape") + ylab("Consumption rate") + coord_flip()
  # theme( axis.text.x = element_text(angle=-45,hjust=0))
e <- ggplot( rate.cwm, aes(y=rate,x=length) )  + geom_point() + geom_smooth(method='glm',method.args=list(family=quasibinomial),se=F) +
  xlab("Body size (mm)") + ylab("Consumption rate") 
f <- ggplot( rate.cwm, aes(y=rate,x=act) )  + geom_point() + geom_smooth(method='glm',method.args=list(family=quasibinomial),se=F) +
  xlab("Dominant activity") + ylab("Consumption rate") 

t2 <- cowplot::plot_grid( a,b,c,d, ncol=2, labels=c("C","D","E","F"), align='hv')
t1 <- cowplot::plot_grid( f,e, ncol=2, labels="AUTO", align='hv')

windows(9,7)
cowplot::plot_grid( t1,t2, ncol=1, rel_widths = c(1,2),  rel_heights = c(1,2), align='hv' )
cowplot::plot_grid( f,e,a,b,c,d, ncol=2, align='hv' )



a <- ggplot( rate.cwm, aes(y=FRic,x=reorder(feed, rate, FUN = mean)) ) + geom_boxplot() + geom_point() +
  xlab("Feeding mode") + ylab("Functional Richness") + coord_flip()
# theme( axis.text.x = element_text(angle=-45,hjust=0))
b <- ggplot( rate.cwm, aes(y=FRic,x=reorder(troph, rate, FUN = median)) ) + geom_boxplot() + geom_point() +
  xlab("Trophic group") + ylab("Functional Richness") + coord_flip()
# theme( axis.text.x = element_text(angle=-45,hjust=0))
c <- ggplot( rate.cwm, aes(y=FRic,x=reorder(watercol, rate, FUN = mean)) ) + geom_boxplot() + geom_point() +
  xlab("Water column use") + ylab("Functional Richness") + coord_flip()
# theme( axis.text.x = element_text(angle=-45,hjust=0))
d <- ggplot( rate.cwm, aes(y=FRic,x=reorder(body, rate, FUN = mean)) ) + geom_boxplot() + geom_point() +
  xlab("Lateral body shape") + ylab("Functional Richness") + coord_flip()
# theme( axis.text.x = element_text(angle=-45,hjust=0))
e <- ggplot( rate.cwm, aes(y=FRic,x=length) )  + geom_point() + 
  geom_smooth(method='lm',se=F) +
  xlab("Body size (mm)") + ylab("Functional Richness") 
f <- ggplot( rate.cwm, aes(y=FRic,x=act) )  + geom_point() + 
  geom_smooth(method='lm',se=F) +
  xlab("Dominant activity") + ylab("Functional Richness") 

t2 <- cowplot::plot_grid( a,b,c,d, ncol=2, labels=c("C","D","E","F"), align='hv')
t1 <- cowplot::plot_grid( f,e, ncol=2, labels="AUTO", align='hv')

windows(9,7)
cowplot::plot_grid( f,e,a,b,c,d, ncol=2, labels="AUTO", align='hv' )




#### why does consumption rate dip at the equator?
## patterns of presence-absence and abundance
## patterns of family presence-absence across latitude
# 
# get capscale values to select and order families of interest
taxa.corr <- read_csv( "Output Data/multivar_constr_taxa.csv" ) %>% arrange(-CAP1)
# get the twenty families with highest correlation with consumption rate
# add more families based on video data?
fam.choose <- taxa.corr$family[taxa.corr$family != c("Caesionidae")]
fam.choose <- c(taxa.corr$family[1:18],"Labridae","Ostraciidae")
fam.choose <- c(taxa.corr$family[1:13],"Labridae",taxa.corr$family[99:104])

pa.choose <- fam.data %>% select( fam.choose )
pa <- bind_cols(fam.meta,pa.choose)
# gather
pa.gather <- pa %>% 
  group_by(SH,Site,habitat,meanLat,meanLong,Country,sstmean,temp,rate) %>% 
  gather(family,PA,fam.choose)
pa.gather <- left_join(pa.gather,taxa.corr)
pa.gather$family <- factor(pa.gather$family)
pa.gather$family <- fct_reorder(.f=pa.gather$family, .x=pa.gather$CAP1, .desc=TRUE)

# visualize
windows(7,6)
ggplot( pa.gather, aes(x=abs(meanLat),y=PA)) + facet_wrap(~family) +
  geom_point(alpha=0.5) + geom_smooth(se=F)
ggplot( pa.gather, aes(x=abs(meanLat),y=PA)) + facet_wrap(~family) +
  geom_point(col="slateblue",alpha=0.5,size=2) + 
  geom_smooth(method='glm',formula=y~poly(x,2),method.args=list(family="binomial"), se=F, col="black", lwd=0.5) +
  geom_smooth(se=F,lty=3,col="black",lwd=0.5) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  ylab("Family presence") + xlab("Degrees from equator")





#### how do pairwise differences in consumption correlate with differences in composition across site pairs (seagrass vs unvegetated)
##
# remove sites with both habitat types
nhab <- rate.mean %>% group_by(Site) %>% summarize( n=length(rate) ) %>% filter(n>1)
nhabseine <- rate.mean %>% filter(!is.na(cpua)) %>% group_by(Site) %>% summarize( n=length(cpua) ) %>% filter(n>1)
diffs <- rate.mean %>% 
  ungroup() %>% 
  filter( Site %in% nhab$Site ) %>% 
  arrange(Site, habitat) %>% 
  group_by(Site) %>% 
  summarize( MDS1=diff(MDS1), MDS2=diff(MDS2), abund=diff(cpua), rate=diff(rate) )


psych::pairs.panels(diffs[,-1])
