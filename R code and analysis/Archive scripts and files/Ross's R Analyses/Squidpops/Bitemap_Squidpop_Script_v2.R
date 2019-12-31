###################################################################################
#                                                                                ##
# Ocean Bitemap Squidpop data: Analyses of global consumption rates              ##
# Data are current as of 20170406                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20170406                                                          ##
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

# 20170406

#

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# INSTALL/LOAD PACKAGES                                                           #
# PREPARE DATASET: Read in data                                                   #
# PLOT DEPLOYMENTS BY HABITAT TYPE                                                #
# PLOT MEAN LOSS OF BAIT BY HABITAT TYPE: One hour, twenty-four hours             #
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

# set working directory 
setwd("~/Dropbox (Personal)/Bitemap Manuscript 2017/R Analyses/Squidpops")

###################################################################################
# PREPARE DATASET: Read in data                                                   #
###################################################################################

# Data files are in folder: setwd("~/Desktop/Smithsonian/TMON/Dropbox (Smithsonian)/Bitemap/Bitemap Scripts/Squidpops")

# NOTE: In order to run summarySE(), you must first run SUMMARY FUNCTION at the end of this script
# NOTE: Be sure to update the .csv file name in the read.csv() command

#  read in survey data (last updated 20170406)
bitemap.squidpop <- read.csv('20170405_Bitemap_Squidpop.csv')

# view data set
glimpse(bitemap.squidpop)

# change numbers to numeric
bitemap.squidpop$Proportion.Missing..24.hours. <- as.numeric(as.character(bitemap.squidpop$Proportion.Missing..24.hours.))
bitemap.squidpop$Proportion.Missing..1.hour. <- as.numeric(as.character(bitemap.squidpop$Proportion.Missing..1.hour.))

#subset and rename habitat factors
bitemap.squidpop <- bitemap.squidpop[bitemap.squidpop$Type.Of.Habitat %in% c("Seagrass Meadow", "Sandy Bottom", "Muddy Bottom", "Rocky Reef"), ]
bitemap.squidpop$Type.Of.Habitat <- revalue(bitemap.squidpop$Type.Of.Habitat, c("Seagrass Meadow" = "vegetated", "Sandy Bottom" = "unvegetated", "Muddy Bottom" = "unvegetated", "Rocky Reef" = "unvegetated"))
bitemap.squidpop$Type.Of.Habitat <- as.factor(bitemap.squidpop$Type.Of.Habitat)


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

# slope of consumption by habitat type

####FIRST HOUR ONLY FOR SLOPES#########
#######################################

bitemap.squidpop$initial <- rep(0)
bitemap.squidpop$hour <- rep(1)
ggplot(bitemap.squidpop, aes(hour, Proportion.Missing..1.hour.)) +
  geom_point() +
  xlim(0,1) +
  geom_abline(aes(yintercept = 0), onehr$'c(mean(Proportion.Missing..1.hour., n...)' +
  facet_grid(.~Type.Of.Habitat) 

all_1 <- sp_binom %>%
  filter(Time == "1") %>%
  group_by(Site,Type.Of.Habitat, Year, Time) %>%
  summarise(sum(Detached.Not))
names(all_1)[names(all_1) == "sum(Detached.Not)"]<- "detached"
Site <- all_1$Site
Year <- all_1$Year
Type.Of.Habitat <- all_1$Type.Of.Habitat
Time <- rep(0, 29)
detached <- rep(0,29)
all_0 <- data.frame(Site, Year, Type.Of.Habitat, Time, detached)
all_all <- bind_rows(all_1, all_0)


spslopes <- all_all %>% 
  group_by(Type.Of.Habitat) %>% 
  do({
    mod = lm(detached ~ Time, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

spslopes

ggplot(all_all, aes(x=Time, y=detached, colour = Type.Of.Habitat)) +
  geom_point(position=position_jitter(width = 0.1)) +
  facet_grid(.~Type.Of.Habitat) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_text(data=data.frame(x=0.5, y=22, label=c("22.83", "12.80","19.67","0.83","5.17"), Type.Of.Habitat = c("fore reef", "mangrove", "patch reef", "sand", "seagrass")), aes(x,y,label=label), inherit.aes=FALSE)










# mean 1hr detachment values
onehr <- bitemap.squidpop %>%
  group_by(Type.Of.Habitat) %>%
  summarise(c(mean(Proportion.Missing..1.hour., na.rm = TRUE)))

# mean 24hr detachment values
tfhr <- bitemap.squidpop %>%
  group_by(Type.Of.Habitat) %>%
  summarise(c(mean(Proportion.Missing..24.hours., na.rm = TRUE)))


# Dotplot of 1hour mean missing by habitat type w/ 95% conf int
ggplot(summary.1hr, aes(x=Type.Of.Habitat, y=Proportion.Missing..1.hour., fill=Type.Of.Habitat)) + 
  geom_errorbar(aes(ymin=Proportion.Missing..1.hour.-ci, ymax=Proportion.Missing..1.hour.+ci), width=.2,                    
                # Width of the error bars
                position=position_dodge(.9)) +
  geom_point(aes(colour = Type.Of.Habitat), cex = 5) + labs(title = "Squidpop Bait Loss, 1 Hour", x="Type of Habitat", y="Proportion Missing") + 
  theme(axis.text.x = element_text(colour="grey20",size=20, angle=20, vjust=0.5),
        axis.text.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20", size=20, vjust=1),
        axis.title.y = element_text(colour="grey20", size=20, vjust=0.75),
        title = element_text(colour="grey20", size=20), legend.position='none') + ylim(0,1) + guides(colour=FALSE) 


# Dotplot of 24hour mean missing by habitat type w/ 95% conf int
ggplot(summary.24hr, aes(x=Type.Of.Habitat, y=Proportion.Missing..24.hours., fill=Type.Of.Habitat)) + 
  geom_errorbar(aes(ymin=Proportion.Missing..24.hours.-ci, ymax=Proportion.Missing..24.hours.+ci), width=.2,                    
# Width of the error bars
position=position_dodge(.9)) +
  geom_point(aes(colour = Type.Of.Habitat), cex = 5) + labs(title = "Squidpop Bait Loss, 24 Hours", x="Type of Habitat", y="Proportion Missing") + 
  theme(axis.text.x = element_text(colour="grey20",size=20, angle=20, vjust=0.5),
        axis.text.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20", size=20, vjust=1),
        axis.title.y = element_text(colour="grey20", size=20, vjust=0.75),
        title = element_text(colour="grey20", size=20), legend.position='none') + ylim(0,1) + guides(colour=FALSE) 
  

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
# T-TESTS                                                                         #
###################################################################################

onehr_ttest <- t.test(Proportion.Missing..1.hour. ~ Type.Of.Habitat, data = bitemap.squidpop)
onehr_ttest

tfourhr_ttest <- t.test(Proportion.Missing..24.hours. ~ Type.Of.Habitat, data = bitemap.squidpop)
tfourhr_ttest

###################################################################################
# LATITUDINAL PATTERNS GRAPHS                                                     #
###################################################################################

ggplot(bitemap.squidpop, aes(x = abs(Latitude..GPS..decimal.), y = Proportion.Missing..1.hour., color = Type.Of.Habitat)) +
  geom_point() + geom_smooth(method = "lm") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Squidpop Bait Loss, 1 hr",
       x="Latitude", y="Proportion Missing")

ggplot(bitemap.squidpop, aes(x = abs(Latitude..GPS..decimal.), y = Proportion.Missing..24.hours., color = Type.Of.Habitat)) +
  geom_point() + geom_smooth(method = "lm") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Squidpop Bait Loss, 24 hrs",
       x="Latitude", y="Proportion Missing")


###################################################################################
# LINEAR MODEL OF RESULTS BY LAT                                                  #
###################################################################################

onehr_lm <- lm(Proportion.Missing..1.hour. ~ Latitude..GPS..decimal. + Type.Of.Habitat, data = bitemap.squidpop)
summary(onehr_lm)

tfourhr_lm <- lm(Proportion.Missing..24.hours. ~ Latitude..GPS..decimal. + Type.Of.Habitat, data = bitemap.squidpop)
summary(tfourhr_lm)

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
