###################################################################################
#                                                                                ##
# Ocean Bitemap Seine data: Analyses of global fish communities                  ##
# Data are current as of 20160913                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20160913                                                          ##
#                                                                                ##
###################################################################################

# TO DO
#     

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the seine to be paired with squidpop data and designed to give rough global metrics of 
# top-down pressure in marine 
# environments. The datasets were delivered via email to the MarineGEO database (SI M Drive)

# WHAT'S NEW IN THIS VERSION

# Version 1 - Everything's new!

#

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# INSTALL/LOAD PACKAGES                                                           #
# PREPARE DATASET: Read in data                                                   #
# OBSERVED FISH RICHNESS                                                          #
# FISH ABUNDANCE                                                                  #
# FISH SIZE DIST
# SUMMARY FUNCTION                                                                #                #
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

# set working directory (if neccessary)
setwd('~/Desktop/Smithsonian/TMON/Data/Bitemap/Bitemap Scripts/Seine/CSVs')

###################################################################################
# PREPARE DATASET: Read in data                                                   #
###################################################################################

# Data files are in folder: Amphitrite/Desktop/Smithsonian/TMON/Data/Bitemap/Bitemap Scripts/
#     Seine/CSVs

# NOTE: In order to run summarySE(), you must first run SUMMARY FUNCTION at the end of this script
# NOTE: Be sure to update the .csv file name in the read.csv() command

#  read in survey data (last updated 20160913)
bitemap.seine.full <- read.csv('Bitemap_Seine_ALL-DATA_20160913.csv')

# view data set
glimpse(bitemap.seine.full)

# create genus-species column
bitemap.seine.full <- unite(bitemap.seine.full, "genus.species", c(Genus,Species), sep = "_", remove = FALSE)

# create unique event ID column
bitemap.seine.full <- unite(bitemap.seine.full, "ID", c(Site.Name,Seagrass.Unveg,Date, Time), sep ="_", remove = FALSE)
bitemap.seine.full$ID <- as.factor(bitemap.seine.full$ID)

# create site/habitat type ID
bitemap.seine.full <- unite(bitemap.seine.full, "Site.Type", c(Site.Name, Seagrass.Unveg), sep = "_", remove= FALSE)
bitemap.seine.full$Site.Type <- as.factor(bitemap.seine.full$Site.Type)

#filter out inverts
bitemap.seine <- filter(bitemap.seine.full, bitemap.seine.full$Fish.Invert.Other != "Invert")

# view data set
glimpse(bitemap.seine)

###################################################################################
# OBSERVED FISH RICHNESS                                                          #
###################################################################################

# Table of raw richness on all seines
Richness <- bitemap.seine %>%
  group_by(Site.Name, Seagrass.Unveg, Lat) %>%
  summarise(richness = length(unique(genus.species))) %>%
  arrange(richness)

# Mean fish richness by habitat type bar
Hab_rich <- Richness %>%
  group_by(Seagrass.Unveg) %>%
  summarise(mean(richness), sd(richness)) 
Hab_rich$sterror <- Hab_rich$sd/sqrt(length(unique(Richness$Site.Name)))

names(Hab_rich)[names(Hab_rich) == 'mean(richness)'] <- 'mean.spp'
names(Hab_rich)[names(Hab_rich) == 'sd(richness)'] <- 'sd.spp'

rich_limits <- aes(ymax = mean.spp + sterror, ymin = mean.spp - sterror)

ggplot(Hab_rich, aes(x=factor(Seagrass.Unveg), y=mean.spp)) + 
  labs(title = "Mean Observed Fish Richness by Habitat Type", y = "Number of Species", 
       x = "Habitat Type") +
  stat_summary(fun.y="mean", geom="bar") +
  geom_errorbar(rich_limits, width = 0.2) 

# Mean fish richness by habitat type box
#fish.rich.box <- ggplot(Hab_rich, aes(x=factor(Seagrass.Unveg), y=mean.spp))
#fish.rich.box +geom_boxplot() +
#  geom_point(aes(x=factor(Seagrass.Unveg), y=mean.spp) )

# Raw richness by latitude
ggplot(Richness, aes(x=Lat, y=richness, color = Seagrass.Unveg)) + 
  geom_point() + stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

###################################################################################
# FISH ABUNDANCE                                                                  #
###################################################################################

# replace NAs in abundance with 1
bitemap.seine$Abundance[is.na(bitemap.seine$Abundance)] <- 1

# total fish abundance per seine
Tot.abun.fish <- bitemap.seine %>%
  group_by(ID, Site.Name, Date, Time, Seagrass.Unveg, Site.Type) %>%
  summarise(abundance = sum(Abundance)) 

# mean observed fish abundance by habitat type
fish.abun.hab.box <- ggplot(Tot.abun.fish, aes(x=Seagrass.Unveg, y=abundance))
fish.abun.hab.box + geom_boxplot() +
  geom_point(aes(x=factor(Seagrass.Unveg), y=abundance)) +
  labs(title = "Observed Fish Abundance by Habitat Type", y = "Abundance", x = "Habitat Type") +
  ylim(0,1000)

# mean observed fish abundance per site and habitat
fish.site.box <- ggplot(Tot.abun.fish, aes(x=Site.Type, y=abundance))
fish.site.box + geom_boxplot() +
  geom_point(aes(x=Site.Type, y=abundance, color = Seagrass.Unveg)) +
  labs(title = "Observed Mean Fish Abundance by Site", y = "Abundance", x = "Site") +
  ylim(0,200)

###################################################################################
# FISH SIZE DIST                                                                  #
###################################################################################

# length by latitude scatter
scatterplot(Length ~ Lat, data = bitemap.seine)

length.lat <- ggplot(bitemap.seine, aes(x=Lat, y=Length))

# heat map of length by latitude
length.lat + geom_bin2d(bins = c(6,2))

w#################
#check structure of data
glimpse(bitemap.seine)

#arrange presence of species by habitat type
num_fish <- bitemap.seine %>%
  group_by(Site.Name, Species) %>%
  summarise(total_abd = unique(Species)) 

#calculate observed richness by location
num_fish_table <- num_fish %>%
  group_by(Site.Name) %>%
  summarise(richness = n())


#create rank abundance plot of abundances by habitat type
num_fish %>%
  group_by(Site.Name) %>%
  arrange(desc(log(total_abd))) %>%
  mutate(rank = row_number(desc(log(total_abd)))) %>%
  ggplot(aes(x = rank, y = log(total_abd))) + geom_line() + facet_wrap(~location)

#arrange survey sites as rows, species abundnace sums as columns
site_abd <- bz.fish %>%
  filter(Data.Type == "Q") %>%
  unite(unique_ID, ID, locale, sep = "_") %>%
  unite(common.binomial, Common.Name, Binomial, sep = "_") %>%
  group_by(unique_ID, common.binomial) %>%
  summarize(abund = sum(Number)) %>%
  spread(common.binomial, abund) 



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



iris %>%
  group_by(Species) %>%
  summarise(avg = mean(Sepal.Width)) %>%
  arrange(avg)
irish
iris
