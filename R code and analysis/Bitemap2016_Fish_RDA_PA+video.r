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



## Seine Data
# calculate different summaries of fish community that can be related to squidpop consumption intensity
# summaries: total abundance, functional group abundance, richness, diversity metrics
#            fish biomass (need to calculate based on size and taxonomy)

########################################################################################
###  Seines were pulled for different distances. We take this into account,
###  but we still lack information on how wide seines were!
########################################################################################


## Video Data
# only use video data to identify which fish families were present


########################################################################################
########################################################################################



# library(profvis)  - can look at places where the code is slow, bottlenecks, etc.
# profvis({

# load libraries
# data.frame manipulation and joining
library(tidyverse)
library(plyr)
library(ggrepel) # for plotting with text
# library(reshape2)
library(cowplot) # for arranging multiple plots
# geospatial data
# library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(vegan)
library(viridis)
library(here)



## read data

# seine abundance + biomass
seine <- read_csv( "Output Data/Bitemap_seine_abundance_biomass.csv" )

# environmental data from remote sensing and oceanographic expeditions
oracle <- read_csv( "Output Data/Bitemap_BioORACLE_20190107.csv")[,1:36]

# squidpop consumption rates
pop <- read_csv( "Output Data/Bitemap_rate.env.20190423.csv" )
pop <- pop %>% 
  dplyr::group_by( Site,habitat ) %>% 
  dplyr::summarise( sstmean=mean(sstmean), temp=mean(temp), rate=mean(rate) )

# trait information that includes video data
video <- read_csv( "../Data/Video Data/Bitemap_Video_Data_ALL.csv" )
video <- video %>% 
  select( Country, habitat, Genus ) %>% 
  mutate( habitat=tolower(habitat) ) %>% 
  distinct()
video$habitat[ video$habitat == "seagrass" ] <- "Seagrass"
video$habitat[ video$habitat == "unveg" ] <- "Unvegetated"
# fix some names
video$Genus <- gsub( "Athrina", "Atherina", video$Genus )
video$Genus <- gsub( "Eucinostramus", "Eucinostomus", video$Genus )
video$Genus <- gsub( "Lactophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Lactrophris", "Lactophrys", video$Genus )
video$Genus <- gsub( "Alters", "Aluterus", video$Genus )
video$Genus <- gsub( "Sphyraema", "Sphyraena", video$Genus )

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

# vet taxa without family info
video %>% filter( Country %in% c("Mexico (ICML)","Mexico (ICML2)")) %>% 
  arrange( Country, family, Genus )
filter( video, is.na(family) )
video$family[ video$Genus %in% "Aluterus"] <- "Monacanthidae"

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
tax$family[ tax$SPECIES_NAME == "Menticirrhus ophicephalus" ] <- "Sciaenidae"
mean.catch.tax <- left_join( mean.catch, tax )

site <- select( oracle, Site, Country )

mean.catch.site <- right_join( site, mean.catch.tax )

## summarize data by family
catch.fam <- mean.catch.site %>% 
  dplyr::group_by( Site, habitat, family ) %>% 
  dplyr::summarize( P=sum(cpua) ) %>% 
  dplyr::mutate( P=ifelse(P>0,1,0) )


# merge with video data
video <- left_join(site,video)
video <- select(video, Site, habitat, family)
video <- video %>% filter( !is.na(family) )
video$P <- 1

catch.video <- bind_rows( catch.fam, video )
catch.video <- catch.video %>% 
  distinct()



# because we are merging seine and video data in some cases, we get discrepancies
# for example: Sparids were caught on video but not in the seine in Croatia
catch.video %>% 
  dplyr::group_by(Site,habitat,family) %>% 
  dplyr::summarise( count=length(P) ) %>% 
  filter( count>1 )

# since there are only ones and zeros, we can add up the presence values
catch.video <- catch.video %>% 
  dplyr::group_by(Site,habitat,family) %>% 
  dplyr::summarise( P=sum(P) ) 

# spread out again
fam.catch.wide <- catch.video %>% 
  group_by( Site, habitat ) %>% 
  spread( family, P, fill=0 )


# can average across habitat types, too

fam.data <- fam.catch.wide[,-c(1:2)]

fam.meta <- fam.catch.wide[,1:2]
fam.meta <- left_join( fam.meta, pop )

fam.meta <- fam.meta %>%  tidyr::unite( SH, Site, habitat, remove=FALSE )

rownames(fam.data) <- fam.meta$SH

# write to disk
fam.all <- bind_cols(fam.meta,fam.data)
write_csv( fam.all, "Output Data/consumer_presence_wide.csv")



##### MULTIVARIATE SECTION


# RDA, CAP
j1 <- capscale( fam.data~1,fam.meta, distance="raup", scale=T )
j2 <- capscale( fam.data~scale(rate),fam.meta, distance="raup", scale=T )
j <- dbrda( fam.data~scale(rate),fam.meta, distance="raup", scale=T )
j2.sst <- capscale( fam.data~scale(sstmean),fam.meta, distance="raup", scale=T, na.action=na.exclude )
j.sst <- dbrda( fam.data~scale(sstmean),fam.meta, distance="raup", scale=T, na.action=na.exclude )
j3 <- capscale( fam.data~scale(rate)+scale(temp),fam.meta, distance="raup", scale=T, na.action=na.exclude )
# r  <- rda( y, scale=T )
# summaries
custum <- function(model){
  list( coef(model),
        RsquareAdj(model),
        summary(model)$cont$importance[,1:6] )
}
custum(j2)
custum(j)
custum(j2.sst)
custum(j.sst)
custum(j3)

# plots
par(mar=c(5,4,2,2))
plot(j1,choices=c(1,2), scaling=2)
plot(j2,choices=c(1,2), scaling=2)
plot(j2,choices=c(1,2), scaling=3)
plot(j2,choices=c(1,2), scaling='symmetric', correlation=TRUE )
plot(j,choices=c(1,2), scaling=2)
plot(j3,choices=c(1,2), scaling=0 )

## get species vectors (can be thought of correlations with CAP axes?)
j2.v <- data.frame( family=row.names(j2$CCA$v), j2$CCA$v ) %>% arrange(CAP1)
j3.v <- data.frame( family=row.names(j3$CCA$v), j3$CCA$v ) %>% arrange(CAP2)
# families associated with extreme rates and temperatures
j2.v[ abs(j2.v$CAP1) > sd(j2.v$CAP1)*1.75, ]  
j3.v[ abs(j3.v$CAP1) > sd(j3.v$CAP1)*1.5 , ]  

# extract axes
nax <- 1:2
juse <- j2
s4 <- scores(juse,choices = nax, scaling = 0)$sites
# colnames(s2)[1] <- paste( colnames(s2)[1],"2",sep="_" )
# extract taxa
t4 <- scores(juse,choices = nax, scaling=0)$species




# also write to disk the vector for the constrained axes
vec.dir <- scores(juse, display="bp")[1]



# customize an RDA figure
# points make sure rate points to the right
capraup <- s4
capraup[,1] <- capraup[,1] * vec.dir
capspec <- t4
capspec[,1] <- capspec[,1] * vec.dir
# rates
capraup <- data.frame(capraup,fam.meta)
coef(juse) # Function coef will give the regression coefficients from centred environmental variables (constraints and conditions) to linear combination scores. The coefficients are for unstandardized environmental variables. The coefficients will be NA for aliased effects.
# proportion explained
RsquareAdj(juse)  # explains 20-25% of variation in consumption rate?
R2 <- eigenvals(juse)/sum(eigenvals(juse))
R2
summary(juse)
cummr2 <- R2
for(i in 2:length(eigenvals(juse))){
  cummr2[i] <- cummr2[i]+cummr2[i-1]  
}

# join with rate.mean
sr <- data.frame( fam.meta, s4 )


ggplot( data=sr, aes(x=CAP1,y=rate) ) + geom_point(size=5) +
  geom_smooth( aes(group=1), method='glm', method.args=list(family='quasibinomial')) 



asite <- ggplot( capraup, aes(x=CAP1,y=MDS1)) + 
  # geom_hline(yintercept = 0, col='gray', lty=2 ) +
  # geom_vline(xintercept = 0, col='gray', lty=2 ) +
  geom_point(aes(size=rate, col=temp),  pch=16, alpha=0.5) +
  xlab(paste0(names(capraup)[1],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(capraup)[2],' (',round(R2[2],3)*100, '%)')) +
  scale_color_viridis() 
  # scale_fill_manual(values=c("green","gray25"))#+ 

# add species
capspec <- data.frame( capspec, family=rownames(capspec) )
capspec <- capspec %>% 
  mutate( dir = atan2(MDS1,-CAP1)/pi*180 ) %>% 
  mutate( angle= (dir %% 360 + 90) %% 360 ) %>% 
  arrange( -angle )
capspec2 <- capspec#[ abs(capspec$CAP1) > sd(capspec$CAP1)*1, ]  
capspec2 %>% 
  mutate( image_length = sqrt(abs(CAP1/max(CAP1))) ) %>% 
  arrange( CAP1 )

####### write constrained ordination sites and taxa to disk
# write_csv( sr, "Output Data/multivar_constr_sites.csv" )
# write_csv( capspec, "Output Data/multivar_constr_taxa.csv" )



# custom
# capspec2$yadj <- c(0.1, 0.1, 0.05, 0, -0.025, -0.025, 0, 0 )
# capspec2 <- capspec2 %>% 
  # mutate( ynew = MDS1+yadj )

aspec <- asite + 
  # geom_text_repel( data=capspec2, aes(x=CAP1,y=MDS1,label=family), col='slateblue',
  #                  box.padding = 0.1 ) +
  # geom_text( data=capspec2, aes(x=CAP1,y=ynew,label=family), col='slateblue',
             # hjust = c(0,0,0,1,0,0,0,0), nudge_x = c(0.1,0.1,0.1,-0.1,0.1,0.1,0.1,0.1 ) ) +
  # geom_segment( data=capspec2, aes(x=0,y=0,xend=CAP1*1,yend=CAP2*1), col='slateblue',
  #               arrow = arrow(length = unit(0.2, "cm")) ) +
  geom_point( data=capspec2, aes(x=CAP1*1,y=MDS1*1), 
              fill='white', col='black', alpha=0.75,
              pch=23, size=3 ) 
  # geom_segment( data=data.frame(CAP1=1.8,MDS1=0), aes(x=0,y=0,xend=CAP1,yend=MDS1),
  #               arrow=arrow(length = unit(0.2,"cm")))
  
# ## all taxa correlations
# aspec2 <- ggplot( capraup, aes(x=CAP1,y=MDS1)) + 
#   geom_segment( data=capspec, aes(x=0,y=0,xend=CAP1,yend=MDS1, col=angle),
#                 arrow = arrow(length = unit(0.2, "cm")) ) +
#   geom_point(aes(size=rate),  pch=16, alpha=0.25) +
#   xlab(paste0('CAP1 (',round(R2[1],3)*100, '%)')) +
#   ylab(paste0('MDS1 (',round(R2[2],3)*100, '%)')) +
#   scale_fill_manual(values=c("green","gray25")) + 
#   scale_color_viridis() +
#   theme_minimal() +
#   theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) 
# aspec2
b <- aspec + theme_minimal() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.border = element_rect( fill=NA ),
         legend.background = element_rect(linetype = 1) ) 
b



###### UNCONSTRAINED
# CAP
# presence-absence
h <- capscale( fam.data~1,fam.meta, 
                distance="raup", scale=T, na.action = na.omit )
# 

# get scores
model <- h
ts <- scores(model,choices = nax, scaling=2)$sites
# colnames(s2)[1] <- paste( colnames(s2)[1],"2",sep="_" )
# extract taxa
tt <- scores(model,choices = nax, scaling=2)$species
# tv <- scores(model, display="bp")[1]
# rsware
R2 <- eigenvals(model)/sum(eigenvals(model))
R2

# write unconstrained ordination to disk
tsdf <- data.frame( SH=fam.meta$SH, ts )
write_csv( tsdf, "Output Data/multivar_unconstr_sites.csv" )
  

## get colors for sites
# based on sst
fam.enviro <- left_join( catch.video, pop )
fam.enviro <- fam.enviro %>% 
  dplyr::group_by( Site, habitat, family ) %>% 
  dplyr::summarize( P=max(P), sstmean=mean(sstmean,na.rm=T), temp=mean(temp,na.rm=T),
             rate=mean(rate,na.rm=T) ) %>% 
  dplyr::group_by( family,P ) %>% 
  dplyr::summarize( sstmean=mean(sstmean,na.rm=T), temp=mean(temp,na.rm=T),
             rate=mean(rate,na.rm=T) ) %>% 
  filter( P==1 )
  

# make a figure
# ts[,1] <- ts[,1] * tv
# tt[,1] <- tt[,1] * tv

# rates
capunc <- data.frame(ts,fam.meta)
capunc <- capunc %>% 
  select( Site, MDS1, MDS2, SST=sstmean, rate )
xlimits <- c(-1.3,1.3)
ylimits <- c(-1.9,1.1)
a <- ggplot( capunc, aes(x=MDS1,y=MDS2,fill=SST)) + 
  geom_point(aes(size=rate), pch=21, alpha=0.75) +
  # geom_text_repel(aes(label=Site), point.padding = 0.5) +
  xlab(paste0(names(capunc)[2],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(capunc)[3],' (',round(R2[2],3)*100, '%)')) +
  scale_fill_viridis(option = "C", limits=c(5,30)) +
  # guides(size=FALSE) +
  scale_x_continuous(limits = xlimits, breaks = c(-1,0,1)) +
  scale_y_continuous(limits = ylimits, breaks = c(-1,0,1)) +
  theme_minimal() +
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_rect( fill=NA ),
         # legend.background = element_rect(linetype = 1),
         # legend.title = element_blank(),
         legend.key.size = unit(0.35, "cm"),
         legend.key.width = unit(0.5,"cm") ) 

tt <- as.data.frame(tt)
tt$family <- row.names(tt)
spec.score.all <- left_join( tt %>% arrange(family), capspec2 %>% arrange(family), by='family' )
spec.score.all %>% filter(family=="Monacanthidae")
spec.score.all %>% filter(family=="Lethrinidae")
specmax <- spec.score.all[abs(spec.score.all$CAP1) > 0.12 & 
                 (abs(spec.score.all$MDS1.x) > sd(spec.score.all$MDS1.x)*1.2 |
                    abs(spec.score.all$MDS2) > sd(spec.score.all$MDS2)*1.5), ]
specmax %>% arrange(CAP1,MDS2)
tt2 <- tt[ abs(tt$MDS1) > sd(tt$MDS1)*1.5 |
                       abs(tt$MDS2) > sd(tt$MDS2)*1.5   , ] 
selected.families <- c( "Cancridae", "Pleuronectidae", "Cottidae", "Embiotocidae",
                        "Atherinidae", "Tetraodontidae", "Labridae", "Sparidae", 
                        "Mugilidae", "Gerreidae", "Portunidae", "Atherinopsidae", 
                        "Hemiramphidae", "Haemulidae" )
tt2 <- tt[ tt$family %in% selected.families ,]
tt2 <- left_join( tt2, fam.enviro )
aspec <- ggplot( capunc, aes(x=MDS1,y=MDS2)) + 
  # geom_hline(yintercept = 0, col='gray', lty=2 ) +
  # geom_vline(xintercept = 0, col='gray', lty=2 ) +
  # geom_point(aes(size=rate, col=temp),  pch=16, alpha=0.5) +
  xlab(paste0(names(capunc)[2],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(capunc)[3],' (',round(R2[2],3)*100, '%)')) +
  geom_point( data=tt2, aes(x=MDS1*1,y=MDS2*1, fill=temp), 
              col='black', alpha=0.75,
              pch=23, size=3 ) +
  scale_fill_viridis(option = "C", limits=c(5,30)) +
  scale_x_continuous(limits = xlimits, breaks = c(-1,0,1)) +
  scale_y_continuous(limits = ylimits, breaks = c(-1,0,1)) 
b <- aspec + theme_minimal() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.border = element_rect( fill=NA ),
         legend.background = element_rect(linetype = 1) ) 
b


windows(7,2.75)
plot_grid( a, b, ncol=2, labels = "AUTO", align = 'hv', axis="tblr",
          rel_widths = c(1,1) )
ggsave( "Figs/Fig2_capscale.pdf", width=7, height=2.5, dpi=600 )
ggsave( "Figs/Fig2_capscale.jpg", width=7, height=2.5, dpi=600 )


captemp <- data.frame(ts,fam.meta[ !is.na(fam.meta$temp), ])
ggplot( captemp, aes(x=CAP1,y=MDS1,col=habitat)) + 
  geom_hline(yintercept = 0, col='gray', lty=2 ) +
  geom_vline(xintercept = 0, col='gray', lty=2 ) +
  geom_point(aes(size=rate),  alpha=0.5) +
  xlab(paste0('CAP1 (',round(R2[1],3)*100, '%)')) +
  ylab(paste0('MDS1 (',round(R2[2],3)*100, '%)')) +
  scale_color_manual(values=c("green","gray25")) +
  theme_minimal() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) 
  


ggplot( data=capunc, aes(x=temp,y=MDS1) ) + geom_point(size=5) +
  theme_classic()
ggplot( data=capunc, aes(x=temp,y=MDS2) ) + geom_point(size=5) +
  theme_classic()
ggplot( data=capunc, aes(x=temp,y=rate) ) + geom_point(size=5) +
  theme_classic()
ggplot( data=capunc, aes(x=sstmean,y=rate) ) + geom_point(size=5) +
  theme_classic()
ggplot( data=capunc, aes(x=MDS1,y=rate) ) + geom_point(size=5) +
  theme_classic()
ggplot( data=capunc, aes(x=MDS2,y=rate) ) + geom_point(size=5) +
  theme_classic()


# pull in 
psych::pairs.panels( capunc[,c("temp","sstmean","MDS1","MDS2","rate")],
                     ellipses = F                   )

