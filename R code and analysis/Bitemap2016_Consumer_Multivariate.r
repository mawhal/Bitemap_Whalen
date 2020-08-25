#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Cleaning fish community data, calculating summaries
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2019.02.19
#####################################################################

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

#####################################################################

## UPDATES FOR THIS SCRIPT
# 2018.10.15: Matt adds Bitemap project to GitHub. Futher changes tracked there, and old changes tracked in archived script files


## Seine Data
# calculate different summaries of fish community that can be related to predation intensity
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





# load libraries
# data.frame manipulation and joining
library(tidyverse)
library(plyr)
library(ggrepel) # for plotting with text
library(reshape2)
library(cowplot) # for arranging multiple plots
# geospatial data
# library(raster) # note that select() also occurs in dplyr
# library(velox) # for faster extract
# accessing data from FishBase
library(rfishbase)
library(taxize)



## read data
group <- read_csv( "R code and analysis/Output Data/consumer_family_community_data.csv" )
fam.meta <- group[,1:8]
fam.data <- group[,-c(1:8)]

phyla <- read_csv( "R code and analysis/Output Data/predator_families+phyla.csv" )

# add site information
sites <- read_csv( "R code and analysis/Output Data/Bitemap_BioORACLE_20190107.csv" )
sites <- sites[ ,c(1:31)]
# sites <- left_join( sites,siteGPS[,1:2], by=c("Country") )

# # add ocean basin information here
# sites$Country
# sites$basin <- c( 1,1,1,1,1,1,
#                   1,1,1,2,2,1,
#                   2,1,3,3,4,2,
#                   3,1,1,2,2,2,
#                   2,1,1,1,1,2,
#                   2,2,2,2,1,2,
#                   2,2,1,1,2)
# sites$basin <- factor( sites$basin, levels=1:4, labels=c("Pacific","Atlantic", "Mediterranean", "Indian") )
# sites$coast <- c( 1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,2,2,2,2,1,2,1,1,2,1,2,2,2,2,1,1,1,1,1,2,1,1,1,2,2,2)
# sites$coast <- factor( sites$coast, levels=1:2, labels=c("West","East") )





##### MULTIVARIATE SECTION

# NMDS on families
library(vegan)
library(viridis)

# # test a number of dimensions
# mds2 <- metaMDS( fam.data, "raup", k = 2, binary=TRUE, trymax = 100 )
# mds3 <- metaMDS( fam.data, "raup", k = 3, binary=TRUE, trymax = 100 )
# mds4 <- metaMDS( fam.data, "raup", k = 4, binary=TRUE, trymax = 100 )
# mds5 <- metaMDS( fam.data, "raup", k = 5, binary=TRUE, trymax = 100 )
# mds6 <- metaMDS( fam.data, "raup", k = 6, binary=TRUE, trymax = 100 )
# mds7 <- metaMDS( fam.data, "raup", k = 7, binary=TRUE, trymax = 100 )
# mds8 <- metaMDS( fam.data, "raup", k = 8, binary=TRUE, trymax = 100 )
# mds9 <- metaMDS( fam.data, "raup", k = 9, binary=TRUE, trymax = 100 )
# mds10 <- metaMDS( fam.data, "raup", k = 10, binary=TRUE, trymax = 100 )
# mds11 <- metaMDS( fam.data, "raup", k = 11, binary=TRUE, trymax = 100 )
# mds12 <- metaMDS( fam.data, "raup", k = 12, binary=TRUE, trymax = 100 )
# mds13 <- metaMDS( fam.data, "raup", k = 13, binary=TRUE, trymax = 100 )
# mds.comp <- list( mds2, mds3, mds4, mds5, mds6, mds7, mds8, mds9, mds10, mds11, mds12, mds13 )
# stress <- unlist( lapply( mds.comp, function(z) z$stress ) )
# plot( stress ) # no break point

# windows(7,4)
mds <- metaMDS( fam.data, "raup", k = 3, binary=TRUE, trymax = 100 )
mds

stressplot( mds )


# capscale
row.names(fam.data) <- fam.meta$Site
# reduce number of families to most abundant
y <- fam.data[-10, which( colSums(fam.data) > 3 )] # remove Delaware. No SST data
x <- group[-10,]

fam.cap <- capscale( y ~ meanLat + poly(sstmean,2,raw = F) , x, dist="raup")
plot(fam.cap)
points( fam.cap, pch=21, cex=x$rate*10)

fam.rda <- dbrda( y ~ poly(sstmean,2), x, dist="raup")
fam.bioenv <- bioenv( y ~ sstmean, x,index="raup")
# pull out constrained axes and look at relationship with SST
cap1 <- fam.cap$CCA$u
cap2 <- fam.cap$Ybar
plot(x$rate ~ cap1[,1] )
plot(x$rate ~ cap2[,1] )
rda1

## rotate the NMDS to mean annual SST
mds.rot <- MDSrotate( mds, group$sstmean, na.rm=TRUE )

## rotate the original NMDS to mean predation rate
mds.rot.rate <- MDSrotate( mds, group$rate, na.rm=TRUE )

mds.use <- mds

plot( mds.rot.rate, "sites", cex=(fam.meta$rate*3)+1, xlim=c(-1,1) )
vec.fam <- envfit( mds.rot.rate$points, fam.data, perm=1000 )
# add to plot
plot( vec.fam, p.max=0.05, col="blue" )

vec.fam.df <- as.data.frame(vec.fam$vectors$arrows*sqrt(vec.fam$vectors$r))
vec.fam.df$family <- rownames(vec.fam.df)
vec.fam.df$p      <- vec.fam$vectors$pvals

# summarize species richness of families associated with colder or warmer temps
vec.pos <- vec.fam.df %>% filter( MDS1>0 )
vec.neg <- vec.fam.df %>% filter( MDS1<0 )

# significant 
vec.pos[ vec.pos$p < 0.05, ]
# strongly in direction of increasing predation rate
vec.pos.strong <- vec.pos[ vec.pos$MDS1 > 0.25, ]
# write to disk
# write.csv( vec.pos.strong, "Output Data/NMDS_rateRotate_positive.csv", row.names = FALSE )

# plot them all together
mdsdf <- cbind( group, mds.rot.rate$points )
windows(3,3)
ggplot( mdsdf, aes(x=MDS1,y=MDS2) ) + 
  geom_point(  aes(size=rate),shape=1, col='slateblue') +
  geom_segment(data=vec.pos.strong,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.2, "cm")),colour="black", alpha=0.5) + 
  geom_text_repel(data=vec.pos.strong,aes(x=MDS1,y=MDS2,label=family),
                  point.padding=0.9, box.padding = 0.1,
                  size=4, col="black", segment.color="slategray", segment.alpha = 0.5) +
  # xlim(xlimits) + ylim(ylimits) +
  theme_classic() + theme( axis.line = element_blank() )



plot( mds.use, "sites", cex=(fam.meta$rate*3)+1, xlim=c(-1,1) )
vec.fam <- envfit( mds.rot.rate$points, fam.data, perm=1000 )
# add to plot
plot( vec.fam, p.max=0.05, col="blue" )

vec.fam.df <- as.data.frame(vec.fam$vectors$arrows*sqrt(vec.fam$vectors$r))
vec.fam.df$family <- rownames(vec.fam.df)
vec.fam.df$p      <- vec.fam$vectors$pvals
#write to disk
# write.csv( vec.fam.df, "Output Data/family_vectors.csv", row.names = FALSE )

fam.sel <- rownames(vec.fam$vectors[[1]])[ vec.fam$vectors$pvals < 0.05 ]

# # plot species richness within each families
# by( predators, predators$Country, function(z) sort(unique( z$SPECIES_NAME )) )
# # only select relevant columns of predator data.frame
# pred.sel <- predators %>%
#   select( Country, Site.Name, Lat, Long, habitat, SPECIES_NAME, family ) %>%
#   distinct( )
# # aggregate to count the number of rows in each family at each site
# pred.richness <- ddply( pred.sel, .(Country, Lat, Long, habitat, family), 
#                         summarize, S = length(SPECIES_NAME) )
# summary(pred.richness)
# # tend to be very few species within families
# # look at species list for sites with high species richness within particular families
# pred.richness[ pred.richness$S >2, ]
# predators[ predators$Country=="Canada (BC)" & predators$habitat=="Seagrass" & predators$family=="Cottidae",]
# predators[ predators$Country=="USA (CA)" & predators$habitat=="Seagrass" & predators$family=="Serranidae",]
# predators[ predators$Country=="USA (TX)" & predators$habitat=="Seagrass" & predators$family=="Sciaenidae",]
# predators[ predators$family=="Diodontidae",]
# unique( predators$Country[ predators$family=="Sillaginidae" ] )
# # merge with the environmental data
# pred.rich.env <- full_join( pred.richness, sites )
# # select influential families
# pre.sel <- pred.rich.env[ pred.rich.env$family %in% fam.sel, ]
# 
# 
# # pad with zeros
# S <- 0
# zeros <- expand.grid( Site = fam.meta$Site, family=unique(phyla$family), S=S )
# # pre.pad <- full_join( pre.sel,zeros, by=c("Site","family","S") )



# gather the families
group.sel <- group  # cbind( fam.meta, group[,fam.sel] ) to pick families with "signifcant" loadings on NMDS
group.sel.melt <- melt( group.sel, id.vars=1:8, value.name = "S", variable.name="family"  ) 
# add phylum back in
group.sel.melt <- left_join( group.sel.melt, phyla, by="family" )

# 
# # windows( 6,3 )
# ggplot( group.sel.melt, aes(x=sstmean, y=S ) ) + facet_wrap(phylum~family, ncol=5) +
#   geom_smooth(se=F) + 
#   geom_point( pch=1 ) +
#   ylab( "Abundance of each family" ) + xlab("Mean Annual SST (?C)") 
# 
# ggplot( pre.pad, aes(x=sstmean, y=S ) ) + facet_wrap(~family, ncol=10) +
#   geom_smooth(se=F) + 
#   geom_point( pch=1 ) +
#   ylab( "Species richness" ) + xlab("Mean Annual SST (?C)")

# dev.off()


# summarize species richness of families associated with colder or warmer temps
vec.pos <- vec.fam.df %>% filter( MDS1>0 )
vec.neg <- vec.fam.df %>% filter( MDS1<0 )

# aggregate and combine
pos.richness <- plyr::ddply( group.sel.melt[ group.sel.melt$family %in% vec.pos$family, ], .(sstmean), 
                             summarize, S=sum(S) )
pos.richness$thermal <- "S.warm"
neg.richness <- plyr::ddply( group.sel.melt[ group.sel.melt$family %in% vec.neg$family, ], .(sstmean), 
                             summarize, S=sum(S) )
neg.richness$thermal <- "S.cold"
total.richness <- plyr::ddply( group.sel.melt, .(sstmean), summarize, S=sum(S) )
total.richness$thermal <- "S.all"
richness.groups <- rbind(pos.richness,neg.richness,total.richness)

# merge richness and sites information
# first melt and cast richness.groups
richness.melt <- melt( richness.groups, id=c(1,3) )
richness.cast <- dcast( richness.melt, sstmean~thermal )
sites.rich <- full_join( sites, richness.cast )



# plot them
windows( 3.5,3 )
ggplot( richness.groups, aes(x=sstmean,y=S) ) + facet_wrap( ~thermal ) +
  # geom_smooth( se=F, col='black' ) + 
  geom_point( pch=1, col='slateblue', alpha=0.7 ) +
#  stat_smooth(method="glm", method.args=list(family="quasipoisson"), 
 #               formula = y ~ splines::ns(x, 2), col="black") +
  ylab( "Species richness" ) + xlab("Mean Annual SST (?C)")



#
adonis2( fam.data ~ coast:(hemi+basin), data=group )
adonis2( fam.data ~ rate, data=group )
adonis2( fam.data ~ coast:(hemi+basin)+rate, data=group )

# create interaction dummies
group$hemicoast  <- with( group, paste(hemi,coast) )
group$coastbasin <- with( group, paste(coast,basin) )

# combine data and NMDS results
mdsfd <- cbind( group, mds.use$points )
# get rid of some columns
mdsfd <- mdsfd %>%
  select( -sstmean )

# now merge in the rotated nmds results
sites.nmds <- full_join( sites.rich, mdsfd, by=c("Site","meanLat","hemi","basin","coast" ) )
sites.nmds <- sites.nmds[ !is.na(sites.nmds$Country), ]
# strong correlation with mean annual SST
plot( MDS1 ~ sstmean, sites.nmds )



### Ellipses for groups
plot( mds, "sites" )
ord<-ordiellipse( mds, group$coastbasin, display = "sites", 
                  kind = "se", conf = 0.95, label = T )
## new data.frame to show ellipses
# https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
# define function for ellipse...DELETED

## unwrap ordisurf object
# https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/
#ordisurf:
ordi<-ordisurf( mds.use, group$sstmean, display = "sites", 
               kind = "se", conf = 0.95, label = T )
ordi.grid <- ordi$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.pred <- expand.grid(x = c(ordi.grid$x), y = c(ordi.grid$y)) #get x and ys
ordi.pred$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.na <- data.frame(na.omit(ordi.pred)) #gets rid of the nas
# ordi.na #looks ready for plotting!
# # repeat with omnivory
# ordi2<-ordisurf( mds.use, fam.meta2$herbprop, display = "sites", 
#                 kind = "se", conf = 0.95, label = T )
# ordi.grid2 <- ordi2$grid #extracts the ordisurf object
# str(ordi.grid2) #it's a list though - cannot be plotted as is
# ordi.pred2 <- expand.grid(x = c(ordi.grid2$x), y = c(ordi.grid2$y)) #get x and ys
# ordi.pred2$z <- as.vector(ordi.grid2$z) #unravel the matrix for the z scores
# ordi.na2 <- data.frame(na.omit(ordi.pred2)) #gets rid of the nas



# the plots
# windows(6,4)
xlimits <- c(-1.5,1.5)
ylimits <- c(-1.0,1.0)
labelsize <- 12
( sstsurf <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
  stat_contour( data = ordi.na, aes(x = x, y = y, z = z, colour = ..level..),
               binwidth = 1.85, size=0.8 ) + #can change the binwidth depending on how many contours you want
  geom_text_repel( aes(label=Site), point.padding = 0.1 ) +
  geom_point( aes(size=group$rate) ) +
  scale_color_viridis( name = "mean\nannual SST" ) +
  scale_size_continuous( name = "predation\nrate" ) +
  xlim(xlimits) + ylim(ylimits) +
  theme_classic() + theme( axis.line = element_blank(), 
                           axis.text=element_text(size=labelsize),
                           axis.title=element_text(size=labelsize+2) ) ) 

# # omnivory
# windows(6,4)
# xlimits <- c(-1.5,1.5)
# ylimits <- c(-1.0,1.0)
# labelsize <- 12
# ( sstsurf <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
#     stat_contour( data = ordi.na2, aes(x = x, y = y, z = z*100, colour = ..level..),
#                   binwidth = 1, size=0.8 ) + #can change the binwidth depending on how many contours you want
#     geom_text_repel( aes(label=Site), point.padding = 0.1 ) +
#     geom_point( aes(size=rate, col=fam.meta2$herbprop*100) ) +
#     scale_color_viridis( name = "Number omnivorous taxa" ) +
#     scale_size_continuous( name = "predation\nrate" ) +
#     xlim(xlimits) + ylim(ylimits) +
#     theme_classic() + theme( axis.line = element_blank(), 
#                              axis.text=element_text(size=labelsize),
#                              axis.title=element_text(size=labelsize+2) ) ) 


with( mdsfd, cor.test( MDS1, sstmean ) )
cor.test( mds$points[,1], mdsfd$sstmean )
cor.test( mds.rot$points[,1], mdsfd$sstmean )

plot( x=fam.meta$propomn, y=group$rate )

# 
vec.fam.sel <- vec.fam.df[vec.fam.df$p < 0.05,]

(vectors <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
  geom_point(  col='slateblue', size=2 ) +
  geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.2, "cm")),colour="slateblue", alpha=0.5) + 
  geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family),
                  point.padding=0.01, box.padding = 1,
                  size=4, col="black", segment.color="slategray", segment.alpha = 0.5) +
  xlim(xlimits) + ylim(ylimits) +
  theme_classic() + theme( axis.line = element_blank() ))

(centroids <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
    geom_point(  col='slateblue', size=2, alpha=0.3 ) +
    geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2),colour="black", alpha=0.3) +
    geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family),
                    point.padding=0.05, box.padding = 0,
                    size=3.5, col="black") +
    xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank() ))

(centroids2 <- ggplot( mdsdf, aes(x=MDS1,y=MDS2) ) + 
    # geom_point(  col='slateblue', size=2, alpha=0.3 ) +
    geom_point( aes(size=rate), col='slateblue', alpha=0.3 ) +
    geom_segment(data=vec.pos.strong,aes(x=0,xend=MDS1,y=0,yend=MDS2),
                 colour="black", alpha=0.8) +
    geom_text_repel(data=vec.pos.strong,aes(x=MDS1,y=MDS2,label=family),
                    point.padding=0.5, box.padding = 0.5,
                    size=4.5, col="slateblue4") +
    scale_size_continuous( name = "predation\nrate" ) +
    xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank(), 
                             axis.text=element_text(size=labelsize),
                             axis.title=element_text(size=labelsize+2),
                             legend.position = "none") ) 


vec.fam.sel <- vec.fam.df[ vec.fam.df$p < 0.05, ]
percentage <- 0.15
percentage <- 0.2  # for UNC IMS interview
mds.move <- diff(range(mds.rot.rate$points[,1])) * percentage
vec.fam.sel <- vec.fam.df[ vec.fam.df$MDS1 < -mds.move | vec.fam.df$MDS1 > mds.move, ]
diff(range(mds.rot.rate$points[,1]))
# add omnivory from Fishbase script
troph2 <- read_csv( "R code and analysis/Output Data/trophic_family_strong_omni.csv" )
vec.fam.sel <- left_join(vec.fam.sel,troph2)

# make a color for centrouds
vec.fam.sel$color <- factor( ifelse( vec.fam.sel$MDS1 >0, "black", "firebrick" ) )
windows(4,4)
(centroids3 <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
    # geom_point(  col='slateblue', size=2, alpha=0.3 ) +
    geom_point( aes(size=rate), col='slateblue', alpha=0.3 ) +
    # geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2,col=color),
    #             alpha=0.75, arrow = arrow(length = unit(0.1,"cm")), linejoin='mitre' ) +
    # geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family, col=color),
    #                 point.padding=0.1, box.padding = 0.3,
    #                 size=4, #col="sla5teblue4",
    #                 segment.color="slateblue1", segment.size=0.3, alpha=0.8) +
    scale_size_continuous( name = "predation\nrate" ) +
    scale_color_manual( values=c("black","firebrick")) +
    xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank(), 
                             axis.text=element_text(size=labelsize),
                             axis.title=element_text(size=labelsize+2),
                             legend.position = "none") ) 

   windows(10,4)
# plot_grid(centroids,sstsurf, labels=c("A","B"), ncol = 2, nrow =1,label_size=18, rel_widths = c(1,1.25) )
plot_grid(sstsurf,centroids3, labels=c("A","B"), ncol = 2, nrow =1,label_size=18, rel_widths = c(1.25,1) )


# supplementary figure to show all predator families
vec.fam.sel <- vec.fam.df[vec.fam.df$p < 1,]  # to show all families
windows( 6,6 )
(vectors <- ggplot( mdsfd, aes(x=MDS1,y=MDS2) ) + 
    geom_point(  col='slateblue', size=1.5 ) +
    geom_segment(data=vec.fam.sel,aes(x=0,xend=MDS1,y=0,yend=MDS2),
                 arrow = arrow(length = unit(0.2, "cm")),colour="slateblue", alpha=0.5) + 
    geom_text_repel(data=vec.fam.sel,aes(x=MDS1,y=MDS2,label=family),
                    point.padding = 1,
                    size=3, col="black", segment.color="black", segment.alpha = 0.5,
                    force=6) +
     xlim(xlimits) + ylim(ylimits) +
    theme_classic() + theme( axis.line = element_blank() ))
# dev.off()


# what's driving communities?
# SIMPER? Bubble Plots

# try PCA, best fits??
# RDA? Correspondence analysis


## Redundant Descriminant analysis
rda1 <- rda(  fam.data ~ sstmean + meanLat, group, na.action=na.exclude )
plot(rda1)

## Canonical Correspondence analysis
cca1 <- cca(  fam.data ~ sstmean, group, na.action=na.exclude )
cca1 <- cca(  fam.data ~ sstmean, group, na.action=na.exclude )
plot(cca1)

# get rid of NA
fam.na <- na.omit( fam.pa )
# separate community data from site data
fam.meta <- fam.na[,1:6]
fam.data <- fam.na[,-c(1:6)]
## environmental data matrix
vare.cca <- cca( fam.data ~ sstmean + basin, fam.meta )
vare.cca
plot(vare.cca)

## RDA
rda1 <- rda( fam.data ~ sstmean + hemi*basin*coast, fam.meta )
anova(rda1)
summary(rda1)
plot(rda1) 

rda1$Ybar


################################
## END OF MULTIVARIATE SECTION
##
################################

