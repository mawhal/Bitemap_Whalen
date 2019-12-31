###################################################################################
#                                                                                ##
# Ocean Bitemap Seine data: Analyses of global fish communities                  ##
# Data are current as of 20170221                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20170515                                                          ##
#                                                                                ##
################################################################################### 

# METADATA:

# This script analyzes data from the Ocean Bitemap program run by the MarineGEO program, 
# Tennenbaum Marine Observatories Network - Smithsonian, focusing on 
# the seine to be paired with squidpop data and designed to give rough global metrics of 
# top-down pressure in marine 
# environments. The datasets were delivered via email to the MarineGEO database (SI M Drive)

# RECENT CHANGES

# 20170712 updated data file to new version, added functional diversity analyses
# 20170515 Analyses from "final" seine data added
# 20170406 All folders and files moved from SI Dropbox
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
# FUNCTIONAL DIVERSITY                                                            #
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
setwd("~/Dropbox (Personal)/Bitemap Manuscript 2017/R Analyses/Seine")

###################################################################################
# PREPARE DATASET: Read in data                                                   #
###################################################################################

# Data files are in folder: smb://SI-NAS1.SMB.US.SINET.SI.EDU/MarineGEO-TMON/MarineGEO Database/Network Activities/Bitemap/2016/Bitemap Fish Survey/

# NOTE: In order to run summarySE(), you must first run SUMMARY FUNCTION at the end of this script
# NOTE: Be sure to update the .csv file name in the read.csv() command

#  read in survey data (last updated 20170221)
bitemap.seine.full <- read.csv('Bitemap_Seine_ALL-DATA_20170712_RW.csv')

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
bitemap.seine.full <- unite(bitemap.seine.full, "genus.species", c(Genus,Species), sep = " ", remove = FALSE)

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

# view species names for errors
unique(bitemap.seine.full$genus.species)

# Correct species spelling mistakes (from Duffy)
# Rename misspelled or non-standard variables
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Albula  vulpes"] <- "Albula vulpes"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherina  hepsetus"] <- "Atherina hepsetus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherinops  affinis"] <- "Atherinops affinis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherinops affins"] <- "Atherinops affinis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherinosoma Atherinosoma microstoma"] <- "Atherinosoma microstoma"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Cheilodactylus  variegatus"] <- "Cheilodactylus variegatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Carcinus Carcinus maenas"] <- "Carcinus maenas"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Crangon  septemspinosa"] <- "Crangon septemspinosa"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ctenogobius boteosoma"] <- "Ctenogobius boleosoma"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Cymatogaster  aggregata"] <- "Cymatogaster aggregata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Cynoscion  nebulosus"] <- "Cynoscion nebulosus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ditema  temminckii"] <- "Ditema temminckii"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Elops  saurus"] <- "Elops saurus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Embiotoca  jacksoni"] <- "Embiotoca jacksoni"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Embiotocidae sp "] <- "Embiotocidae spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Engraulis  ringens"] <- "Engraulis ringens"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Epinephelus sp. 1"] <- "Epinephelus spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Etropus  longimanus"] <- "Etropus longimanus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Farfantepenaus spp."] <- "Farfantepenaeus spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Farfantepenaeus spp."] <- "Farfantepenaeus spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Fundulus  similis"] <- "Fundulus similis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gasterosteus aculeatus aculeatus"] <- "Gasterosteus aculeatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Girella Girella zebra"] <- "Girella zebra"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobies species complex"] <- "Goby spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae sp. 1"] <- "Goby sp1"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae sp. 2"] <- "Goby sp2"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae sp. 3"] <- "Goby_sp3"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae sp. 4"] <- "Goby sp4"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma  bosci"] <- "Gobiosoma bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma sp."] <- "Gobiosoma sp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma  bosci"] <- "Gobiosoma bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma  bosci"] <- "Gobiosoma bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma bosc"] <- "Gobiosoma bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobius  geniporus"] <- "Gobius geniporus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gymnapistes Gymnapistes marmoratus"] <- "Gymnapistes marmoratus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Heteroclinus Heteroclinus perspicillatus"] <- "Heteroclinus perspicillatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hippocampus  sp."] <- "Hippocampus sp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hipsoblennius  sordidus"] <- "Hypsoblennius sordidus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hypsoblennius  jenkinsi"] <- "Hypsoblennius jenkinsi"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Lagodon  rhomboides"] <- "Lagodon rhomboides"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Leiostomus  xanthurus"] <- "Leiostomus xanthurus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Leptocottus  armatus"] <- "Leptocottus armatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Liza  argentea"] <- "Liza argentea"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Lutjanus  russellii"] <- "Lutjanus russellii"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menidia  menidia"] <- "Menidia menidia"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menidia beryllinia"] <- "Menidia beryllina"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menidia  beryllina"] <- "Menidia beryllina"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menticirrhus   ophicephalus"] <- "Menticirrhus ophicephalus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Metacarcinus  magister"] <- "Metacarcinus magister"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Micrometrus  minimus"] <- "Micrometrus minimus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Nesogobius Nesogobius maccullochi"] <- "Nesogobius maccullochi"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Odontesthes  regia"] <- "Odontesthes regia"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Orthopistis chrysoptera"] <- "Orthopristis chrysoptera"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Palaemon  elegans"] <- "Palaemon elegans"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Palaemonetes  spp."] <- "Palaemonetes spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Palaemonetes spp."] <- "Palaemonetes spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Paralichthys  californicus"] <- "Paralichthys californicus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Paralichthys  microps"] <- "Paralichthys microps"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Petroscirtis  variabilis"] <- "Petroscirtis variabilis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Phanerodon  furcatus"] <- "Phanerodon furcatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pholis  ornata"] <- "Pholis ornata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Platycephalidae sp. 1"] <- "Platycephalidae sp1"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Platycephalus Platycephalus aurimaculatus"] <- "Platycephalus aurimaculatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pomatoschistus sp."] <- "Pomatoschistus sp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Prionotus  punctatus"] <- "Prionotus punctatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Processa  macrophthalma"] <- "Processa macrophthalma"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Psettichthy melanostictus"] <- "Psettichthys melanostictus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Rhinobatos  productus"] <- "Rhinobatos productus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Romaleon  polyodon"] <- "Romaleon polyodon"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Scartichthys  viridis"] <- "Scartichthys viridis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "scorpaenichthys  marmoratus"] <- "Scorpaenichthys marmoratus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "sculpin sp."] <- "sculpin spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sebastes sp"] <- "Sebastes spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Siganus  fuscescens"] <- "Siganus fuscescens"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sillago  ciliata"] <- "Sillago ciliata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sillago maculatus"] <- "Sillago maculata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Solea  solea"] <- "Solea solea"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Stigmatopora Stigmatopora argus"] <- "Stigmatopora argus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sygnathus  floridae"] <- "Syngnathus floridae"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sygnathus floridae"] <- "Syngnathus floridae"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sygnathus griseolineatus"] <- "Syngnathus griseolineatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus  floridae"] <- "Syngnathus floridae"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus  leptorhynchus\xe6"] <- "Syngnathus leptorhynchus_xe6"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus  leptorhynchus\xca"] <- "Syngnathus leptorhynchus_xca"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus sp."] <- "Syngnathus spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Tetractenos Tetractenos glaber"] <- "Tetractenos glaber"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Torquiginer sp.r"] <- "Torquiginer spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Umbrina  roncador"] <- "Umbrina roncador"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Unknown goby sp "] <- "Goby sp6"

# From Paul York's work
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ambassis marinus"] <- "Ambassis marianus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Arothron hipidus"] <- "Arothron hispidus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Centropogon  australis"] <- "Centropogon australis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ditema temminckii"] <- "Ditrema temminckii"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ephinephelus coioides"] <- "Epinephelus coioides"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gasterosteus aculateus"] <- "Gasterosteus aculeateus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gerres subfasciatus "] <- "Gerres subfasciatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Herring sp"] <- "Clupeidae spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hyporhamphus regularis ardelio"] <- "Hyporhamphus regularis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Monacanthis chinensis"] <- "Monacanthus chinensis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Oncorhyncus tshawytscha"] <- "Oncorhynchus tshawytscha"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Paralichthys albiguttata"] <- "Paralichthys albigutta"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pelates sexlineatus"] <- "Helotes sexlineatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pipefish sp"] <- "Syngnathidae spp."
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pleuronectes platessa "] <- "Pleuronectes platessa"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pseudoblennius dottoides"] <- "Pseudoblennius cottoides"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Rhombosolea Rhombosolea tapirina"] <- "Rhombosolea tapirina"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sillago temminck"] <- "Sillago japonica"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Suarida nebulosa"] <- "Saurida nebulosa"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus griseolineatus"] <- "Ambassis marianus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus leptorhynchus xe6"] <- "Syngnathus leptorhynchus"

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
#Plot_richness <- Richness[c(1:12, 15:66, 69:98),]
# OR RUN AS IS
Plot_richness <- Richness

# Plot observed fish richness by latitude/habitat type (need change some "1" values to "0")
ggplot(Plot_richness, aes(factor(abs(Lat)), richness)) + geom_boxplot(aes(fill = factor(Seagrass.Unveg))) +  labs(title = "Bitemap Fish Richness by Latitude", y = "Number of Species", x = "Latitude") 

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
#20170515
#[1] 18759

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
fish.site.box <- ggplot(Tot.abun.fish, aes(x=Seagrass.Unveg, y=abundance, fill = Seagrass.Unveg))
fish.site.box + geom_boxplot() +
  geom_point(aes(x=Seagrass.Unveg, y=abundance)) +
  labs(title = "Observed Mean Fish Abundance by Site", y = "Abundance", x = "Site") +
  ylim(0,200)

###################################################################################
# FISH SIZE DIST                                                                  #
###################################################################################

# log transform length data
bitemap.seine.fish$loglength <- log(bitemap.seine.fish$Length)

# length by latitude scatter
scatterplot(Length ~ abs(Lat), data = bitemap.seine.fish)

# heat map of length by latitude
length.lat <- ggplot(bitemap.seine.fish, aes(x=abs(Lat), y=loglength), na.rm = TRUE)
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
# FUNCTIONAL DIVERSITY                                                            #
###################################################################################

# import character list
characters <- read.csv("Traits_all-species_edit.csv")

# fix level duplication
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="planktivore"] <- "Planktivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="higher carnivore"] <- "Higher carnivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="omnivore"] <- "Omnivore"
levels(characters$Trophic.group)

# import length bins
#Lbins <- read.csv("lengthbins.csv")
# import these bins for soler sizes
# Lbins <- read.csv("lengthbins.csv")

# gather data so every fish is a row
#Seine_func <- bitemap.seine.full %>%
#  gather("size", "total", 24:51)
#filter out inverts
Seine_func <- filter(bitemap.seine.full, Phylum != "Arthropoda"| is.na(Phylum))
Seine_func <- filter(Seine_func, Phylum != "Mollusca"| is.na(Phylum))
Seine_func$Phylum <- factor(Seine_func$Phylum)


# match characters to RLS data
names(characters)[names(characters)=="CURRENT_TAXONOMIC_NAME"] <- "genus.species"
Seine_func <- left_join(Seine_func, characters, by="genus.species")

# identify species not in character list
unique(Seine_func$genus.species[!(Seine_func$genus.species %in% characters$genus.species)])

# create MaxLength bin column
RLS_func <- left_join(RLS_func, Lbins, by="size")

# create functional entity column 
RLS_func <- RLS_func %>%
  unite_("entity", c("Lbin", "Trophic.group", "Water.column"), sep = "_", remove = FALSE)

# make functional entity a factor
RLS_func$entity <- as.factor(RLS_func$entity)

# remove extraneous columns
RLS_func <- RLS_func %>%
  select(Site.Name, Habitat.x, Year, Depth, Method, Block, Species, total, Family, Class, entity)

# remove 0 occurrences
RLS_func <- RLS_func %>%
  filter(total != "0")

# extract data to analyse pelagic and cryptic seperately
RLS_func1 <- RLS_func %>%
  filter(Method == "1")

RLS_func2 <- RLS_func %>%
  filter(Method == "2")

# Continue with M1/M2 analysis

detach(package:plyr)
# sum functional groups for each site M1/M2
Func_anal <- RLS_func %>%
  group_by(Site.Name, Habitat.x, entity) %>%
  summarise(sum(total))

# change column name
names(Func_anal)[names(Func_anal)=="sum(total)"] <- "total"

# spread for diversity analysis
Func_anal <- Func_anal %>%
  spread(entity, total, fill = 0)

# simpson diversity of functional groups
func_simp <- Func_anal[,3:40]
func_cats <- Func_anal[,1:2]
attach(func_cats)


simpson <- diversity(func_simp, index = "simpson")
func_cats$simpson <- simpson


# functional diversity simpson boxplot M1 & M2
div.box <- ggplot(func_cats, aes(x=Habitat.x, y=simpson, fill = Habitat.x) )
div.box + geom_boxplot() +
  geom_point(aes(x=Habitat.x, y=simpson)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Pelagic/Benthic Functional Simpson Diversity (2015-16)",
       x="Type of Habitat", y="Simpson Value")

# extract values
func_simp_sum <- func_cats%>%
  group_by(Habitat.x) %>%
  summarise(median(simpson))

# M1 functional diversity

# sum functional groups for each site M1 
Func_anal1 <- RLS_func1 %>%
  group_by(Site.Name, Habitat.x, entity) %>%
  summarise(sum(total))

# change column name
names(Func_anal1)[names(Func_anal1)=="sum(total)"] <- "total"

# spread for diversity analysis
Func_anal1 <- Func_anal1 %>%
  spread(entity, total, fill = 0)

# simpson diversity of functional groups
func_simp1 <- Func_anal1[,3:54]
func_cats1 <- Func_anal1[,1:2]
attach(func_cats1)


simpson <- diversity(func_simp1, index = "simpson")
func_cats1$simpson <- simpson


# functional diversity simpson boxplot M1 
div.box <- ggplot(func_simp1, aes(x=Habitat.x, y=simpson, fill = Habitat.x) )
div.box + geom_boxplot() +
  geom_point(aes(x=Habitat.x, y=simpson)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Pelagic Functional Simpson Diversity (2015-16)",
       x="Type of Habitat", y="Simpson Value")

# M2 functional diversity

# sum functional groups for each site M2
Func_anal2 <- RLS_func2 %>%
  group_by(Site.Name, Habitat.x, entity) %>%
  summarise(sum(total))

# change column name
names(Func_anal2)[names(Func_anal2)=="sum(total)"] <- "total"

# spread for diversity analysis
Func_anal2 <- Func_anal2 %>%
  spread(entity, total, fill = 0)

# simpson diversity of functional groups
func_simp2 <- Func_anal2[,3:22]
func_cats2 <- Func_anal2[,1:2]
attach(func_cats2)


simpson <- diversity(func_simp2, index = "simpson")
func_cats2$simpson <- simpson


# functional diversity simpson boxplot M1 
div.box <- ggplot(func_simp2, aes(x=Habitat.x, y=simpson, fill = Habitat.x) )
div.box + geom_boxplot() +
  geom_point(aes(x=Habitat.x, y=simpson)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Cryptic Functional Simpson Diversity (2015-16)",
       x="Type of Habitat", y="Simpson Value")

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


