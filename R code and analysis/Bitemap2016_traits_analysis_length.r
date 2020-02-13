#############################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Use traits of consumers in analyses of community stucture, consumption, etc.
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 30 August 2019
#############################################################################


## load libraries
library(tidyverse) 
library(FD)



## SEINE DATA FOR ABUNDANCE DATA AND SIZE DATA
## read seine data
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )
seine <- seine %>% 
  select( Country, habitat, Lat, Long, Date, Time, sciName=SPECIES_NAME, length=Length, biomass, abun=Abundance )
# change some taxon names
seine$sciName[ seine$sciName=="Pelates sexlineatus"] <- "Helotes sexlineatus"
unique(seine$sciName[ !(seine$sciName %in% trait.final$sciName) ]) # still a few missing taxa
# summarize median per capita length and biomass per seine
seine_per <- seine %>% 
  group_by( Country, habitat, Lat, Long, Date, Time, sciName ) %>% 
  summarize( length=median(length,na.rm=T), biomass=median(biomass,na.rm=T), abun=sum(abun,na.rm=T) )
# summarize mean across multiple seines within habitats
seine_hab <- seine_per %>% 
  group_by( Country, habitat, sciName ) %>% 
  summarize( present=1 )
  # summarize( present=1 )
# summarize mean across species
seine_species <- seine_per %>% 
  group_by( sciName ) %>% 
  summarize( length=mean(length,na.rm=T), biomass=mean(biomass,na.rm=T), abun=mean(abun,na.rm=T) )



## TRAITS
## read trait data - see traits_imput.R and traits_clean.r for details
trait.final <- read_csv( "Output Data/traits_clean_final.csv" )

# merge length data from seines
trait.final <- left_join( trait.final, seine_species )



### switch which traits to use
trait.use <- trait.final %>% 
  select( sciName, name, act, feed, troph, watercol, body, length  )



## mice package 
md <- trait.use %>%
  mice::md.pattern()






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
s1 <- left_join( seine_hab, ss )
s1 <- left_join( s1, trait.final )


# filter out taxa not in trait.final and vice versa
s1 <- s1 %>% 
  filter( sciName %in% trait.use$sciName )
trait.use <- trait.use %>% 
  filter( sciName %in% s1$sciName )


# duplicate entries?
s1 %>%
  arrange( Site, habitat, sciName ) %>%
  group_by( Site, habitat, sciName) %>%
  summarise( entries=length(sciName) ) %>%
  filter( entries>1 )


### summarize trait data by site
# filter out only fish?
suse <- s1 #filter( s1, phylum=="Chordata" )

# indvidual and taxon based ratios below are currently the same because only one individual per site*habitat combination
# ratio of active to ambush at each site by individuals 
rat.ind <- suse %>%
  dplyr::group_by( Site, habitat ) %>%
  dplyr::summarize(act.ratio.ind = length(act[act %in% c("active","ambush.active")]) / length(sciName)  ) 
# and by taxon
rat.tax <- suse %>%
  dplyr::select( Site, habitat, sciName, act ) %>%
  dplyr::distinct() %>%
  dplyr::group_by( Site, habitat ) %>%
  dplyr::summarize( act.ratio.tax = length(act[act %in% c("active","ambush.active")]) / length(sciName),
                    richness=length(sciName)) 
# merge
rat <- left_join( rat.ind, rat.tax )
rat.site <- left_join( rat, ss )

# write to disk
# write_csv( rat, "Output Data/consumer_active_ratio2.csv" )


# plot
ggplot( rat.site, aes(x=abs(meanLat), y=act.ratio.ind)) + geom_point() + geom_smooth()
ggplot( rat.site, aes(x=act.ratio.ind,  y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))

# nice figure
# windows(3,3)
ggplot( rat.site, aes(x=act.ratio.tax,  y=rate)) + 
  geom_point(pch=16,alpha=0.5) +
  xlab("Ratio of active to restful taxa") + ylab("Consumption rate") +
  theme_classic()
summary(glm(rate~act.ratio.tax,data=rat.site,family="quasibinomial"))


# repeat with length
len.ind <- suse %>%
  dplyr::group_by( Site, habitat ) %>%
  dplyr::summarize(len.mean = mean(length)) 
len.site <- left_join( len.ind, ss )
# plot
ggplot( len.site, aes(x=abs(meanLat), y=len.mean)) + geom_point() + geom_smooth()
ggplot( len.site, aes(x=len.mean,  y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial' ))















## community weighted mean traits and functional diversity indices
# FD package

## trait matrix
xuse <- data.frame( trait.use[ trait.use$sciName %in% s1$sciName,-c(1,2)] )
xuse[xuse==""] <- NA
xuse$act[ xuse$act == "active;ambush" ] <- "ambush"
# to numeric
# xuse[,c(2)] <- apply(xuse[,c(2)],2,as.numeric)
xuse$act <- as.numeric( as.factor(xuse$act))-1
summary(xuse)
row.names(xuse) <- trait.use$name
xuse <- xuse[ order(rownames(xuse)), ]
# xuse$squidpop <- as.numeric(as.character(xuse$squidpop))
mice::md.pattern( xuse )
# seine.keep <- s1[ s1$sciName %in% trait.keep$sciName, ]
  



## abundance matrix
comm <- s1 %>% 
  ungroup() %>% 
  unite( "SH", Site, habitat  ) %>% 
  select( SH, name, abun ) %>% 
  spread( name,abun,fill=0  )
ause <- as.matrix( comm[,-1] )
rownames(ause) <- comm$SH
  
# put in another matrix
xmat <- list( trait=xuse, abun=ause )



gd <- gowdis( xmat$trait, asym.bin = 1 )
m1 <- dbFD( gd, xmat$abun, w.abun = T, 
            corr="lingoes", print.pco=T, calc.FGR = T, clust.type="kmeans",
            km.inf.gr=3,km.sup.gr=10,
            calc.CWM = FALSE, m="min", calc.FDiv = FALSE, stand.FRic = TRUE )

# Community weight means
cwm <- functcomp( xmat$trait, xmat$abun )
cwm$SH <- comm$SH
cwm <- cwm %>% 
  separate( SH, c("Site","habitat"))
cwm.ss <- left_join( cwm, ss )

write_csv( cwm, "Output Data/FunctionalDiversity_CWM_length.csv")

# quantitative traits
ggplot( cwm.ss, aes( x=abs(meanLat), y=as.numeric(as.character(act)))) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))
ggplot( cwm.ss, aes( x=as.numeric(as.character(act)), y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))
# qualitative traits
ggplot( cwm.ss, aes( x=abs(meanLat), y=feed)) + geom_point() 
ggplot( cwm.ss, aes( x=abs(meanLat), y=troph)) + geom_point() 
ggplot( cwm.ss, aes( x=abs(meanLat), y=watercol)) + geom_point() 

ggplot( cwm.ss, aes( x=feed, y=rate)) + geom_point() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
ggplot( cwm.ss, aes( x=troph, y=rate)) + geom_point() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
ggplot( cwm.ss, aes( x=watercol, y=rate)) + geom_point() +
  theme(axis.text.x = element_text(angle=45,hjust=1))



## Other dimensions of functional diversity
names(m1)
# number of functional groups
fd <- with(m1, data.frame(FRic,qual.FRic,FEve,FDis, RaoQ,FGR) )
fd$SH <- rownames(fd)
fd <- fd %>% 
  separate( SH, c("Site", "habitat"))

fds <- left_join(fd,ss)
psych::pairs.panels( fds[,c("FRic","FEve","FDis","RaoQ","FGR","rate")])
#RaoQ, FGR, Fric
# from the help page: Rao's quadratic entropy (Q) is computed from the uncorrected species-species distance matrix via divc. See Botta-Dukát (2005) for details. Rao's Q is conceptually similar to FDis, and simulations (via simul.dbFD) have shown high positive correlations between the two indices (Laliberté and Legendre 2010). Still, one potential advantage of FDis over Rao's Q is that in the unweighted case (i.e. with presence-absence data), it opens possibilities for formal statistical tests for differences in FD between two or more communities through a distance-based test for homogeneity of multivariate dispersions (Anderson 2006); see betadisper for more details.
# psych::pairs.panels( fds[,c("FRic","RaoQ","FGR","rate")])
write_csv( fd, "Output Data/FunctionalDiversity_indices_length.csv" )


