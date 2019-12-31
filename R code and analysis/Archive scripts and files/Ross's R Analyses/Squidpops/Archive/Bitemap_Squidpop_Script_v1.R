###################################################################################
#                                                                                ##
# Ocean Bitemap Squidpop data: Analyses of global consumption rates              ##
# Data are current as of 20160921                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20160921                                                          ##
#                                                                                ##
###################################################################################

# TO DO
#     

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the Squidpop assay designed to give rough global metrics of top-down pressure in marine 
# environments. The datasets are extracted from the Ocean Bitemap online portal
# http://bitemap.wordpress.com via Google Form entries. 

# WHAT'S NEW IN THIS VERSION

# Version 1 - Everything's new!

#

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# INSTALL/LOAD PACKAGES                                                           #
# PREPARE DATASET: Read in data                                                   #
# PLOT DEPLOYMENTS BY HABITAT TYPE                                                #
# PLOT MEAN LOSS OF BAIT BY HABITAT TYPE: One hour, twenty-four hours             #
# SUMMARY FUNCTION                                                                                #
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

# set working directory (if neccessary)
# setwd("~/Desktop/Smithsonian/TMON/Dropbox (Smithsonian)/Bitemap/Bitemap Scripts/Squidpops")

###################################################################################
# PREPARE DATASET: Read in data                                                   #
###################################################################################

# Data files are in folder: setwd("~/Desktop/Smithsonian/TMON/Dropbox (Smithsonian)/Bitemap/Bitemap Scripts/Squidpops")

# NOTE: In order to run summarySE(), you must first run SUMMARY FUNCTION at the end of this script
# NOTE: Be sure to update the .csv file name in the read.csv() command

#  read in survey data (last updated 20170112)
bitemap.squidpop <- read.csv('Ocean Bitemap Squidpop Data - Squidpops 20170112.csv')

# view data set
glimpse(bitemap.squidpop)

# change numbers to numeric
bitemap.squidpop$Proportion.Missing..24.hours. <- as.numeric(as.character(bitemap.squidpop$Proportion.Missing..24.hours.))

# replace long factor names
#levels(bitemap.squidpop$Type.Of.Habitat) <- c(levels(bitemap.squidpop$Type.Of.Habitat), 
#                                              "Artificial Habitat")
#bitemap.squidpop$Type.Of.Habitat[bitemap.squidpop$Type.Of.Habitat == 
#              'Artificial Habitat (dock, breakwater, weir, etc.)'] <- 'Artificial Habitat'
#bitemap.squidpop$Type.Of.Habitat <- factor(bitemap.squidpop$Type.Of.Habitat)

#subset and rename habitat factors
bitemap.squidpop <- bitemap.squidpop[bitemap.squidpop$Type.Of.Habitat %in% c("Seagrass Meadow", "Sandy Bottom", "Muddy Bottom", "Rocky Reef"), ]
bitemap.squidpop$Type.Of.Habitat <- revalue(bitemap.squidpop$Type.Of.Habitat, c("Seagrass Meadow" = "vegetated", "Sandy Bottom" = "unvegetated", "Muddy Bottom" = "unvegetated", "Rocky Reef" = "unvegetated"))
bitemap.squidpop$Type.Of.Habitat <- as.factor(bitemap.squidpop$Type.Of.Habitat)

# summarize 1 hour missing bait
summary.1hr <- summarySE(bitemap.squidpop, measurevar="Proportion.Missing..1.hour.", groupvars="Type.Of.Habitat")

# summarize 24 hour missing bait
summary.24hr <- summarySE(bitemap.squidpop, measurevar="Proportion.Missing..24.hours.", groupvars="Type.Of.Habitat", na.rm = TRUE)

###################################################################################
# PLOT DEPLOYMENTS BY HABITAT TYPE                                                #
###################################################################################

# histogram of deployments by habitat type
qplot(factor(Type.Of.Habitat), data=bitemap.squidpop, geom="bar", fill=factor(Type.Of.Habitat)) +
       labs(title="Number Of Deployments By Habitat", x="Type of Habitat", y="Number of Deployments") + 
  theme(axis.text.x = element_text(colour="grey20",size=20, angle=20, vjust=0.5),
        axis.text.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20", size=20, vjust=1),
        axis.title.y = element_text(colour="grey20", size=20, vjust=0.75),
        title = element_text(colour="grey20", size=20), legend.position='none')

###################################################################################
# PLOT MEAN LOSS OF BAIT BY HABITAT TYPE: One hour, twenty-four hours             #
###################################################################################

# barchart of 1hour mean missing by habitat type w/ 95% conf int
ggplot(summary.1hr, aes(x=Type.Of.Habitat, y=Proportion.Missing..1.hour., fill=Type.Of.Habitat)) + 
  geom_bar(position=position_dodge(), stat="identity") + labs(title="Squidpop Bait Loss, 1 hr",
                x="Type of Habitat", y="Proportion Missing") + 
  theme(axis.text.x = element_text(colour="grey20",size=20, angle=20, vjust=0.5),
        axis.text.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20", size=20, vjust=1),
        axis.title.y = element_text(colour="grey20", size=20, vjust=0.75),
        title = element_text(colour="grey20", size=20), legend.position='none') + ylim(0,1) + guides(colour=FALSE) +
  geom_errorbar(aes(ymin=Proportion.Missing..1.hour.-ci, ymax=Proportion.Missing..1.hour.+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

# barchart of 24hour mean missing by habitat type w/ 95% conf int
ggplot(summary.24hr, aes(x=Type.Of.Habitat, y=Proportion.Missing..24.hours., fill=Type.Of.Habitat)) + 
  geom_bar(position=position_dodge(), stat="identity") + labs(title="Squidpop Bait Loss, 24 hrs",
                                                              x="Type of Habitat", y="Proportion Missing") + 
  theme(axis.text.x = element_text(colour="grey20",size=20, angle=20, vjust=0.5),
        axis.text.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20", size=20, vjust=1),
        axis.title.y = element_text(colour="grey20", size=20, vjust=0.75),
        title = element_text(colour="grey20", size=20), legend.position='none') + ylim(0,1) + guides(colour=FALSE) +
  geom_errorbar(aes(ymin=Proportion.Missing..24.hours.-ci, ymax=Proportion.Missing..24.hours.+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

# boxplot of 1hr missing by habitat type
one.hour.box <- ggplot(bitemap.squidpop, aes(x=Type.Of.Habitat, y=Proportion.Missing..1.hour., fill=Type.Of.Habitat))
one.hour.box +geom_boxplot() +
  geom_point(aes(x=Type.Of.Habitat, y=Proportion.Missing..1.hour.) ) +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
  labs(title="Squidpop Bait Loss, 1 hr",
       x="Type of Habitat", y="Proportion Missing")

# boxplot of 24hr missing by habitat type
tf.hour.box <- ggplot(bitemap.squidpop, aes(x=Type.Of.Habitat, y=Proportion.Missing..24.hours., fill=Type.Of.Habitat))
tf.hour.box +geom_boxplot() +
  geom_point(aes(x=Type.Of.Habitat, y=Proportion.Missing..24.hours.) ) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Squidpop Bait Loss, 24 hrs",
       x="Type of Habitat", y="Proportion Missing")

###################################################################################
# LATITUDINAL PATTERNS                                                            #
###################################################################################

ggplot(bitemap.squidpop, aes(x = Latitude..GPS..decimal., y = Proportion.Missing..1.hour., color = Type.Of.Habitat)) +
  geom_point() + geom_smooth(method = "lm") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Squidpop Bait Loss, 1 hr",
       x="Latitude", y="Proportion Missing")

ggplot(bitemap.squidpop, aes(x = Latitude..GPS..decimal., y = Proportion.Missing..24.hours., color = Type.Of.Habitat)) +
  geom_point() + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Squidpop Bait Loss, 24 hrs",
       x="Latitude", y="Proportion Missing")


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

one_hour_anova <- anova(lm(Amount.Of.Bait.Missing.After.One.Hour ~ Number.Of.Squidpops.Deployed * Type.Of.Habitat, bitemap.squidpop))
twentyfour_hour_anova <- anova(lm(Amount.Of.Bait.Missing.After.24.Hours ~ Number.Of.Squidpops.Retrieved.After.24.Hours * Type.Of.Habitat, bitemap.squidpop))
twentyfour_hour_anova 
