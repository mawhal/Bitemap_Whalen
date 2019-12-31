###################################################################################
#                                                                                ##
# Ocean Bitemap Seine data: Analyses of global fish communities                  ##
# Data are current as of 20170221                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20170222                                                          ##
#                                                                                ##
################################################################################### 

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the seine to be paired with squidpop data and designed to give rough global metrics of 
# top-down pressure in marine 
# environments. The datasets were delivered via email to the MarineGEO database (SI M Drive)

# RECENT CHANGES

# 20170222 version 4
# make nMDS method 'altGower' used 'Site.Name' as factor
# 20160930 version 3, 
# Version 2 - Expanded Fish Characteristic analyses 
# Version 1 - Everything's new!

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# INSTALL/LOAD PACKAGES                                                           #
# PREPARE DATASET: Read in data                                                   #
# OBSERVED FISH RICHNESS                                                          #
# FISH ABUNDANCE                                                                  #
# FISH SIZE DIST                                                                  #
# SUMMARY FUNCTION                                                                #               
###################################################################################

#

###################################################################################
# INSTALL/LOAD PACKAGES                                                           #
###################################################################################

# install required packages
# install.packages('dplyr')
# install.packages('tidyr')
# install.packages('plyr')
# install.packages('ggplot2')

# load required packages
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
detach(package:plyr)    
library(dplyr)
library(car)
library(vegan)
library(reshape)

# set working directory (if neccessary)
setwd("~/Desktop/Smithsonian/TMON/Dropbox (Smithsonian)/Bitemap/Bitemap Scripts/Seine")

###################################################################################
# PREPARE DATASET: Read in data                                                   #
###################################################################################

# Data files are in folder: smb://SI-NAS1.SMB.US.SINET.SI.EDU/MarineGEO-TMON/MarineGEO Database/Network Activities/Bitemap/2016/Bitemap Fish Survey/

# NOTE: In order to run summarySE(), you must first run SUMMARY FUNCTION at the end of this script
# NOTE: Be sure to update the .csv file name in the read.csv() command

#  read in survey data (last updated 20170221)
bitemap.seine.full <- read.csv('Bitemap_Seine_ALL-DATA_20170221.csv')

# view data set
glimpse(bitemap.seine.full)

# load plyr for function 
library(plyr)
# make factors 'vegetated' 'unvegetated'
bitemap.seine.full$Seagrass.Unveg <- revalue(bitemap.seine.full$Seagrass.Unveg, 
                                             c("Seagrass" = "vegetated", "Unveg" = "unvegetated"))
bitemap.seine.full$Seagrass.Unveg <- as.factor(bitemap.seine.full$Seagrass.Unveg)

# detach plyr
detach(package:plyr)

# create genus-species column
bitemap.seine.full <- unite(bitemap.seine.full, "genus.species", c(Genus,Species), sep = "_", remove = FALSE)

# create unique event ID column
bitemap.seine.full <- unite(bitemap.seine.full, "ID", c(Site.Name,Seagrass.Unveg,Date), sep ="_", remove = FALSE)
bitemap.seine.full$ID <- as.factor(bitemap.seine.full$ID)

# create site/habitat type ID
bitemap.seine.full <- unite(bitemap.seine.full, "Site.Type", c(Site.Name, Seagrass.Unveg), sep = "_", remove= FALSE)
bitemap.seine.full$Site.Type <- as.factor(bitemap.seine.full$Site.Type)

# make Length numeric
bitemap.seine.full$Length <- as.character(bitemap.seine.full$Length)
bitemap.seine.full$Length <- as.numeric(bitemap.seine.full$Length)

# relevel factors for color
bitemap.seine.full$Seagrass.Unveg <- ordered(bitemap.seine.full$Seagrass.Unveg, levels = c("unvegetated", "vegetated"))

# replace NAs in abundance with 1
bitemap.seine.full$Abundance[is.na(bitemap.seine.full$Abundance)] <- 1

# view data set
glimpse(bitemap.seine.full)

###################################################################################
# OBSERVED FISH RICHNESS                                                          #
###################################################################################

#filter out inverts
bitemap.seine.fish <- filter(bitemap.seine.full, Phylum != "Arthropoda"| is.na(Phylum))
bitemap.seine.fish <- filter(bitemap.seine.fish, Phylum != "Mollusca"| is.na(Phylum))
bitemap.seine.fish$Phylum <- factor(bitemap.seine.fish$Phylum)

# Table of raw richness on all seines (plyr must be detached)
Richness <- bitemap.seine.fish %>%
  group_by(ID, Site.Type, round(Lat, 2), Seagrass.Unveg) %>%
  summarise(richness = length(unique(genus.species))) %>%
  arrange(Site.Type)
names(Richness)[names(Richness)=="round(Lat, 2)"] <- "Lat"

# Remove sites that didn't complete entire fish seine protocol
Plot_richness <- Richness[c(1:12, 15:66, 69:98),]

# Plot observed fish richness by latitude/habitat type (need change some "1" values to "0")
ggplot(Plot_richness, aes(factor(Lat), richness)) + geom_boxplot(aes(fill = factor(Seagrass.Unveg))) +  labs(title = "Bitemap Fish Richness by Latitude", y = "Number of Species", x = "Latitude") 

# ANOVA of observed richness
rich_aov <- aov(richness ~ Seagrass.Unveg * Lat, data = Plot_richness)
summary(rich_aov)
# 20170222
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Seagrass.Unveg      1  112.9  112.86  10.816 0.001437 ** 
#  Lat                 1  146.3  146.27  14.017 0.000319 ***
#  Seagrass.Unveg:Lat  1   18.2   18.18   1.742 0.190182    
#Residuals          90  939.2   10.44                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Mean fish richness by habitat type box
ggplot(Richness, aes(x=factor(Seagrass.Unveg), y=richness, fill = Seagrass.Unveg) ) +geom_boxplot()  +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Observed Fish Richness",
       x="Type of Habitat", y="Number of Species")

# Raw richness by latitude (no S Hemi)
ggplot(Plot_richness, aes(Lat, richness, colour = Seagrass.Unveg)) + geom_point() + geom_smooth(method = lm) + xlim(27,68) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title = element_text(hjust = 0.5)) +
  labs(title="Observed Fish Richness",
       x="Latitude", y="Number of Species")


###################################################################################
# FISH ABUNDANCE                                                                  #
###################################################################################

# total fish abundance per seine
Tot.abun.fish <- bitemap.seine.fish %>%
  group_by(ID, Site.Name, Date, round(Lat,2), Seagrass.Unveg) %>%
  summarise(abundance = sum(Abundance)) 
names(Tot.abun.fish)[names(Tot.abun.fish)=="round(Lat, 2)"] <- "Lat"

# total number of fish counted
sum(Tot.abun.fish$abundance)
#20170222
#[1] 16441

# mean observed fish abundance by habitat type
fish.abun.hab.box <- ggplot(Tot.abun.fish, aes(x=Seagrass.Unveg, y=abundance, fill = Seagrass.Unveg))
fish.abun.hab.box + geom_boxplot() +
  geom_point(aes(x=factor(Seagrass.Unveg), y=abundance)) +
  labs(title = "Observed Fish Abundance by Habitat Type", y = "Abundance", x = "Habitat Type") +
  ylim(0,1000) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Observed Fish Abundance",
       x="Type of Habitat", y="Number of Individuals")

# Remove sites that didn't complete entire fish seine protocol
Plot_abun <- Tot.abun.fish[c(1:12, 15:66, 69:98),]

# Plot observed fish richness by latitude/habitat type (need change some "1" values to "0")
ggplot(Plot_abun, aes(factor(Lat), abundance)) + geom_boxplot(aes(fill = factor(Seagrass.Unveg)))

# Raw abundance by latitude
#ggplot(Tot.abun.fish, aes(x=Lat, y=abundance, color = Seagrass.Unveg)) + 
#  geom_point() + stat_smooth(method = "lm") +
#  ylim (0,400)

# mean observed fish abundance per site and habitat
fish.site.box <- ggplot(Tot.abun.fish, aes(x=Site.Type, y=abundance, fill = Seagrass.Unveg))
fish.site.box + geom_boxplot() +
  geom_point(aes(x=Site.Type, y=abundance)) +
  labs(title = "Observed Mean Fish Abundance by Site", y = "Abundance", x = "Site") +
  ylim(0,200)

###################################################################################
# FISH SIZE DIST                                                                  #
###################################################################################

# log transform length data
bitemap.seine$loglength <- log(bitemap.seine$Length)

# length by latitude scatter
scatterplot(Length ~ Lat, data = bitemap.seine)

# heat map of length by latitude
length.lat <- ggplot(bitemap.seine, aes(x=Lat, y=Length))
length.lat + geom_bin2d(bins = c(8,20)) +
  labs(title = "Fish Length Density by Latitude", y = "Length", x = "Latitude") +
  ylim(0,350)

# mean length of fish by habitat type
length.hab <- ggplot(bitemap.seine, aes(x=factor(Seagrass.Unveg), y = Length, fill = Seagrass.Unveg))
length.hab + geom_boxplot() +
  geom_point(aes(x=Seagrass.Unveg, y = Length)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Observed Fish Length",
       x="Habitat Type", y="Fish Length") +
  ylim(0,310)

ggplot(bitemap.seine, aes(x=Lat, y=Length, color = Seagrass.Unveg)) + 
  geom_point() + stat_smooth(method = "lm") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Fish Length by Latitude",
       x="Latitude", y="Length") +
  ylim(0,310)

###################################################################################
# FISH DIVERSITY                                                                  #
###################################################################################

#filter out inverts
bitemap.seine.fish <- filter(bitemap.seine.full, Phylum != "Arthropoda"| is.na(Phylum))
bitemap.seine.fish <- filter(bitemap.seine.fish, Phylum != "Mollusca"| is.na(Phylum))
bitemap.seine.fish$Phylum <- factor(bitemap.seine.fish$Phylum)

#subset data
seine.div <- select(bitemap.seine.fish, Site.Name, round(Lat,2), Seagrass.Unveg, genus.species, Abundance)

# detach plyr
detach(package:plyr)
#sum species totals
seine.div.tot <- seine.div %>%
  group_by(Site.Name, Seagrass.Unveg, genus.species, Lat) %>%
  summarise(sum(Abundance))
 
names(seine.div.tot)[names(seine.div.tot)=="sum(Abundance)"] <- "Abundance"
  
#matrix
seine.mat <- seine.div.tot %>%
  spread(genus.species, Abundance)
seine.mat[is.na(seine.mat)] <- 0

seine.mat <- seine.mat[,c(1:3,5:194)]
seine.mat <- seine.mat[rowSums(seine.mat[,4:193])!=0,]
simpson <- diversity(seine.mat[,4:193], "simpson")
seine.mat$simpson <- simpson
simp_div <- seine.mat[,c(1:3, 194)]

# simpson diversity by habitat type
div.plot <- ggplot(simp_div, aes(x=Seagrass.Unveg, y=simpson, fill = Seagrass.Unveg))
div.plot + geom_boxplot() + geom_point(aes(x=Seagrass.Unveg, y = simpson)) +
  ylim(0.01,1) +
  labs(title = "Diversity by Habitat Type", y = "Diversity", x = "Habitat Type") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) 

# Simpson diversity by latitude
ggplot(simp_div, aes(x=Lat, y=simpson, color = Seagrass.Unveg)) + 
  geom_point() + stat_smooth(method = "lm") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  xlim(27,68) +
  labs(title="Simpson Diversity by Latitude",
       x="Latitude", y="Diversity")

specnumber(seine.mat[,4:193], groups = "Site.Name")
# make matrix readable by metaMDS
seine.mat.runme <- seine.mat[,4:193]

# MDS of community
bitemap.fish.MDS <- metaMDS(seine.mat.runme, dist = "altGower")
plot.data<- cbind(seine.mat, bitemap.fish.MDS$points)
ggplot(plot.data, aes(x=MDS1, y=MDS2, pch = Seagrass.Unveg)) +
  geom_point(aes(x=MDS1, y=MDS2, color = Site.Name), size = 4)

###################################################################################
# SUMMARY FUNCTION                                                                #
###################################################################################

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

one_hour_anova <- anova(lm(Amount.Of.Bait.Missing.After.One.Hour ~ Number.Of.Squidpops.Deployed * Type.Of.Habitat, bitemap.seine))
twentyfour_hour_anova <- anova(lm(Amount.Of.Bait.Missing.After.24.Hours ~ Number.Of.Squidpops.Retrieved.After.24.Hours * Type.Of.Habitat, bitemap.seine))


