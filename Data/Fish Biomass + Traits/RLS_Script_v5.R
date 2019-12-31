###################################################################################
#                                                                                ##
# Reef Life Survey: Analyses of fish communities                                 ##
# Data are current as of 20170119                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20170413                                                          ##
#                                                                                ##
###################################################################################

# TO DO

# fix "# calculate biomass per observation" section
# add method 0 back into richness
# measures of trophic positions
# aesthetics of all figures

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# RICHNESS                                                                        #
# FUNCTIONAL DIVERSITY                                                            #
# BIODIVERSITY DATA                                                               #
# FISH DENSITY                                                                    #
# FISH BIOMASS                                                                    #
# SPECIES LISTS                                                                   #
# SIZE SPECTRA                                                                    #
# MULTIVARIATE DATA                                                               #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 20170413 Added combined M1/M2 analyses and removed to new file
# 20170405 Script updated to v5, log10 for density, size spectra added
# 20170404 Added biomass analyses and MDS
# 20170403 Added famililal diversity analyses
# 20170317 Added functional diversity analyses 
# 20170227 Separated old code from new code, brought in species list code. 
# 20170119 Script taken from Bioacoustics_RLS_Script_v4.R for Carrie Bow manuscript RLS analyses

# Deveoped from Bioacoustics RLS script v4 from Leray fish analysis script

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(nlme)
library(plyr)
library(car)
library(ggplot2)
library(reshape2) 
library(dplyr)
library(tidyr)
library(rfishbase)
library(vegan)
library(grid)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# Before import, be sure to change first column header "ID" to "Line_ID" in .csv. 
# Delete lat/long columns, and first empty row if present. Add Habitat Column.

# Change working directory
setwd("~/Dropbox (Personal)/CBC Manuscript 2017/R Analyses/RLS")

# Import the survey data (version 2 has correct Block values)
RLS_survey_data <- read.csv("CBCRLS_CORE.csv")

glimpse(RLS_survey_data)

# Create unique ID for each survey
RLS_survey_data <- RLS_survey_data %>%
  unite_("event", c("Site.Name", "Date", "Block"), sep = "_", remove = FALSE)

# make unique ID a factor
RLS_survey_data$event <- as.factor(RLS_survey_data$event)

# replace NAs with 0 
RLS_survey_data[is.na(RLS_survey_data)] <- 0

#remove method 1 with no given abundance
RLS_survey_data <- RLS_survey_data %>%
  subset(Total !="0")

# fix incorrect taxa
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Sphoeroides spenglerii"] <- "Sphoeroides spengleri"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Xyricthys splendens"] <- "Xyrichtys splendens"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Ctenogonius stigmaturus"] <- "Ctenogobius stigmaturus"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Ophistignathus sp."] <- "Opistognathus sp."
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Ctenogonius stigmaturus"] <- "Ctenogobius stigmaturus"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Hyppolytidae"] <- "Hippolytidae"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Paradiplogrammus bairdi"] <- "Callionymus bairdi"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Cryptotoums roseus"] <- "Cryptotomus roseus"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Eucinostoums gula"] <- "Eucinostomus gula"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Atherinidae"] <- "Atherinid spp."
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Clupeoid spp."] <- "Clupeidae"


# remove unknown and freshwater species
RLS_survey_data <- RLS_survey_data %>%
  filter(Species != c("Gambusia sp.")) %>%
  filter(Species != c("Scarine spp."))

# all fish community all years 
RLS_1_2 <- RLS_survey_data %>%
  filter(Inverts < 1)
# remove incidental inverts
RLS_1_2 <- RLS_1_2 %>%
  filter(Species != "Hippolytidae") %>%
  filter(Species != "Mithracid spp.")

#method 1 community 
RLS_1 <- RLS_survey_data %>%
  filter(Method == "1")

#method 1 community 2015
RLS_1_2015 <- RLS_survey_data %>%
  filter(Method == "1") %>%
  filter(Year == 2015)

#method 1 community 2016
RLS_1_2016 <- RLS_survey_data %>%
  filter(Method == "1") %>%
  filter(Year == 2016)

#method 2 fish
RLS_2 <- RLS_survey_data %>%
  filter(Method == "2") %>%
  filter(Inverts == 0)

#method 2 fish 2015
RLS_2_2015 <- RLS_survey_data %>%
  filter(Method == "2") %>%
  filter(Inverts == 0) %>%
  filter(Year == 2015)

#method 2 fish 2016
RLS_2_2016 <- RLS_survey_data %>%
  filter(Method == "2") %>%
  filter(Inverts == 0) %>%
  filter(Year == 2016)

#inverts
RLS_inverts <- RLS_survey_data %>%
  filter(Inverts > 0)

#inverts 2015
RLS_inverts_2015 <- RLS_survey_data %>%
  filter(Inverts > 0) %>%
  filter(Year == 2015)

#inverts 2016
RLS_inverts_2016 <- RLS_survey_data %>%
  filter(Inverts > 0) %>%
  filter(Year == 2016)

# remove method 0 
#RLS_survey_data_12 <- RLS_survey_data %>%
#  subset(Method != "0")
# remove method 1 with no given abundance
#RLS_survey_data_12 <- RLS_survey_data_12 %>%
# subset(Total !="0")

###################################################################################
# RICHNESS                                                                        #
###################################################################################

#set dataset to analyse
#year, method
x <- RLS_1_2016 %>%
  spread(Species, Total, fill = 0)
x <- x[,c(6,9, 50:158)]
x <- x %>%
  group_by(Site.Name, Habitat) %>%
  summarise_each(funs(sum))
y <- x[,3:111]
z <- x[,1:2]
attach(z)

# Estimated richness using Chao2
pool <- y %>%
  specpool(Habitat)

# STOP, CHANGE NAME
# date_CBC_method_measure_year.csv
#Save as CSV
write.csv(pool, "20172224_CBC_invertRichness.csv")

# chao plot
pool$habitat <- as.factor(c("Fore Reef", "Mangrove", "Patch Reef", "Sand", "Seagrass"))
boot.box <- ggplot(pool, aes(x=habitat, y=boot, fill = habitat) )
limits <- aes(ymax = boot + boot.se, ymin=boot - boot.se)
boot.box + geom_point() +
  geom_errorbar(limits, width=0.2) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Chao Richness 2016 (M1)",
       x="Type of Habitat", y="Richness")

plot(chao ~ habitat, data = pool)

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
Lbins <- read.csv("lengthbins.csv")

# gather data so every fish is a row
RLS_func <- RLS_1_2 %>%
 gather("size", "total", 24:51)

# match characters to RLS data
names(characters)[names(characters)=="CURRENT_TAXONOMIC_NAME"] <- "Species"
RLS_func <- left_join(RLS_func, characters, by="Species")

# identify species not in character list
unique(RLS_1_2$Species[!(RLS_1_2$Species %in% characters$Species)])

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
# BIODIVERSITY DATA                                                               #
###################################################################################

#set dataset to analyse
#year, method
x <- RLS_1 %>%
  spread(Species, Total, fill = 0)
x <- x[,c(6,9, 50:193)]
x <- x %>%
  group_by(Site.Name, Habitat) %>%
  summarise_each(funs(sum))
y <- x[,3:146]
z <- x[,1:2]
attach(z)


# simpson diversity
simpson <- diversity(y, index = "simpson")
z$simpson <- simpson

# shannon diversity
shannon <- diversity(y, index = "shannon")
z$shannon <- shannon

# effective number of species
z$ENS <- exp(shannon)

#STOP, RENAME

# date_CBC_method_measure_year.csv
# Write diversity values to csv 
write.csv(z, "20170224_CBC_M1_div_2016.csv")

# simpson box plot
div.box <- ggplot(x, aes(x=Habitat, y=simpson, fill = Habitat) )
div.box +geom_boxplot() +
  geom_point(aes(x=Habitat, y=simpson)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Simpson Diversity 15-16 (M1)",
       x="Type of Habitat", y="Simpson Value")



# Familial diversity

# import character list
characters <- read.csv("Traits_all-species.csv")

# remove inverts
RLS_1_2 <- RLS_survey_data %>%
  filter(Inverts < 1)

# gather data so every fish is a row
RLS_fam <- RLS_1_2 %>%
  gather("size", "total", 24:51)

# match characters to RLS data
names(characters)[names(characters)=="CURRENT_TAXONOMIC_NAME"] <- "Species"
RLS_fam <- left_join(RLS_fam, characters, by="Species")

# identify species not in character list
unique(RLS_1_2$Species[!(RLS_1_2$Species %in% characters$Species)])


# remove extraneous columns
RLS_fam <- RLS_fam %>%
  select(Site.Name, Habitat.x, Year, Depth, Method, Block, Species, total, Family, Class)

# remove 0 occurrences
RLS_fam <- RLS_fam %>%
  filter(total != "0")

# extract data to analyse pelagic and cryptic seperately
RLS_fam1 <- RLS_fam %>%
  filter(Method == "1")

RLS_fam2 <- RLS_fam %>%
  filter(Method == "2")

# sum familial groups for each site M1 and M2
Fam_anal <- RLS_func %>%
  group_by(Site.Name, Habitat.x, Family) %>%
  summarise(sum(total))

# change column name
names(Fam_anal)[names(Fam_anal)=="sum(total)"] <- "total"

# spread for diversity analysis
Fam_anal <- Fam_anal %>%
  spread(Family, total, fill = 0)

# simpson diversity of familial groups
fam_simp <- Fam_anal[,3:47]
fam_cats <- Fam_anal[,1:2]
attach(fam_cats)


simpson <- diversity(fam_simp, index = "simpson")
fam_cats$simpson <- simpson


# familial diversity simpson boxplot M1 & M2
div.box <- ggplot(fam_simp, aes(x=Habitat.x, y=simpson, fill = Habitat.x) )
div.box + geom_boxplot() +
  geom_point(aes(x=Habitat.x, y=simpson)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Pelagic/Benthic Familial Simpson Diversity (2015-16)",
       x="Type of Habitat", y="Simpson Value")

# M1 familial diversity

# sum family groups for each site M1 
Fam_anal1 <- RLS_fam1 %>%
  group_by(Site.Name, Habitat.x, Family) %>%
  summarise(sum(total))

# change column name
names(Fam_anal1)[names(Fam_anal1)=="sum(total)"] <- "total"

# spread for diversity analysis
Fam_anal1 <- Fam_anal1 %>%
  spread(Family, total, fill = 0)

# simpson diversity of familial groups
fam_simp1 <- Fam_anal1[,3:38]
fam_cats1 <- Fam_anal1[,1:2]
attach(fam_cats1)


simpson <- diversity(fam_simp1, index = "simpson")
fam_cats1$simpson <- simpson


# familial diversity simpson boxplot M1 
div.box <- ggplot(fam_simp1, aes(x=Habitat.x, y=simpson, fill = Habitat.x) )
div.box + geom_boxplot() +
  geom_point(aes(x=Habitat.x, y=simpson)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Pelagic Familial Simpson Diversity (2015-16)",
       x="Type of Habitat", y="Simpson Value")

# M2 familial diversity

# sum functional groups for each site M2
Fam_anal2 <- RLS_fam2 %>%
  group_by(Site.Name, Habitat.x, Family) %>%
  summarise(sum(total))

# change column name
names(Fam_anal2)[names(Fam_anal2)=="sum(total)"] <- "total"

# spread for diversity analysis
Fam_anal2 <- Fam_anal2 %>%
  spread(Family, total, fill = 0)

# simpson diversity of familial groups
fam_simp2 <- Fam_anal2[,3:19]
fam_cats2 <- Fam_anal2[,1:2]
attach(fam_cats2)


simpson <- diversity(fam_simp2, index = "simpson")
fam_cats2$simpson <- simpson


# familial diversity simpson boxplot M2 
div.box <- ggplot(fam_simp2, aes(x=Habitat.x, y=simpson, fill = Habitat.x) )
div.box + geom_boxplot() +
  geom_point(aes(x=Habitat.x, y=simpson)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Cryptic Familial Simpson Diversity (2015-16)",
       x="Type of Habitat", y="Simpson Value")

###################################################################################
# FISH DENSITY                                                                    #
###################################################################################

#set dataset to analyse
#year, method
x <- RLS_1 %>%
  spread(Species, Total, fill = 0)

x <- x[,c(6,9, 50:188)] #remove col 11 for single year analyses

x <- x %>%
  group_by(Site.Name,Habitat) %>% # remove "Year" for single year analyses
  summarise_each(funs(sum))

y <- x[,3:141] # start col 3-single year, start col 4-all year
y$sums <- rowSums(y)
z <- x[,1:2] #1:2 for single year, c(1,3) for all years
z$sums <- log10(1 + y$sums)

# Add back in 0 occurrence sites for M2
# Twin Cays Mangrove 15, Curlew Sand 15, Twin Cays Seagrass 16, Bluground Mangrove 16, Tobacco Mangrove 16

  # For 15-16

  Site.Name <- c("Twin Cays Mangrove", "Curlew Sand", "Twin Cays   Seagrass", "Blueground Mangrove", "Tobacco Mangrove")
  Habitat <- c("Mangrove", "Sand", "Seagrass", "Mangrove", "Mangrove")
  sums <- c("0", "0", "0", "0", "0")
  M2zeros <- data.frame(Site.Name, Habitat, sums)
  M2zeros$sums <- as.numeric(as.character(M2zeros$sums))
  z <- bind_rows(z, M2zeros)
  
  # For 15
  
  Site.Name <- c("Twin Cays Mangrove", "Curlew Sand")
  Habitat <- c("Mangrove", "Sand")
  sums <- c("0", "0")
  M2zeros <- data.frame(Site.Name, Habitat, sums)
  M2zeros$sums <- as.numeric(as.character(M2zeros$sums))
  z <- bind_rows(z, M2zeros)
  
  # For 16

  Site.Name <- c("Twin Cays   Seagrass", "Blueground Mangrove", "Tobacco Mangrove")
  Habitat <- c("Seagrass", "Mangrove", "Mangrove")
  sums <- c("0", "0", "0")
  M2zeros <- data.frame(Site.Name, Habitat, sums)
  M2zeros$sums <- as.numeric(as.character(M2zeros$sums))
  z <- bind_rows(z, M2zeros)  
  

#STOP, RENAME

# date_CBC_method_measure_year.csv
# Write diversity values to csv 
write.csv(z, "20170301_densM1_CBC16.csv")

# Plot of fish density
# density box plot
dens.box <- ggplot(z, aes(x=Habitat, y=sums, fill = Habitat) )
dens.box +geom_boxplot() +
  geom_point(aes(x=Habitat, y=sums)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Log Fish Density per 2500m^3 (M1) 15-16",
       x="Type of Habitat", y="Density")

###################################################################################
# FISH BIOMASS                                                                    #
###################################################################################

# import biomass coefficients
coef <- read.csv("20160710_RLS_biomass_coefs.csv")

# remove inverts
RLS_1_2 <- RLS_survey_data %>%
  filter(Inverts < 1)

# gather data so every fish is a row
RLS_biom <- RLS_1_2 %>%
  gather("size", "total", 24:51)

# turn lengths into numbers
RLS_biom$size <-  gsub("X([0-9])", "\\1", RLS_biom$size)
RLS_biom$size <- as.numeric(RLS_biom$size)

# match characters to RLS data
names(coef)[names(coef)=="SPECIES_NAME"] <- "Species"
RLS_biom <- full_join(RLS_biom, coef, by="Species")

# identify species not in character list
unique(RLS_1_2$Species[!(RLS_1_2$Species %in% coef$Species)])

# create biomass column 
RLS_biom$biomass <- RLS_biom$total*(RLS_biom$A*RLS_biom$size^RLS_biom$B)

# remove extraneous columns
RLS_biom <- RLS_biom %>%
  select(Site.Name, Habitat, Year, Depth, Method, Block, Species, size, total, biomass)

# make year a factor
RLS_biom$Year <- as.factor(RLS_biom$Year)

# remove 0 occurrences
RLS_biom <- RLS_biom %>%
  filter(total != "0")

# extract data to analyse pelagic 
RLS_biom1 <- RLS_biom %>%
  filter(Method == "1")

RLS_biom2 <- RLS_biom %>%
  filter(Method == "2")

# sum biomass per site and year pelagic
M1_biom <- RLS_biom1 %>%
  group_by(Site.Name, Habitat, Year) %>%
  summarise(sum(biomass, na.rm = TRUE))

# change column name
names(M1_biom)[names(M1_biom)=="sum(biomass, na.rm = TRUE)"] <- "biomass"

# make biomass into kg and log
M1_biom$biomass <- log10(M1_biom$biomass*(0.001))



# graph biomass by habitat type
bm.box <- ggplot(M1_biom, aes(x=Habitat, y=biomass, fill = Habitat))
bm.box + geom_boxplot() +
  geom_point(aes(x=Habitat, y=biomass)) +
  facet_grid(.~Year) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(title="Fish Biomass By Habitat Type",
       x = "Type of Habitat", y="Biomass") 


###################################################################################
# SPECIES LISTS                                                                   #
###################################################################################

# create complete species list
allspecies <- unique(RLS_survey_data$Species)

# invert-only species list
invertspp <- RLS_survey_data %>%
  subset(Inverts != "0")
invertspp <-  unique(invertspp$Species)

# fish-only species list
fishspp <- RLS_survey_data %>%
  subset(Inverts == "0")
fishspp <- unique(fishspp$Species)

# write fish species list to csv
write.csv(fishspp, "20170227_fishspplist_CBC.csv")

###################################################################################
# SIZE SPECTRA                                                                    #
###################################################################################

# import length bins
Lbins <- read.csv("lengthbins.csv")
# import these bins for soler sizes
#Lbins <- read.csv("lengthbins_soler.csv")

# gather data so every fish is a row
RLS_spec <- RLS_1_2 %>%
  gather("size", "total", 24:51)

# create size bin column
RLS_spec <- left_join(RLS_spec, Lbins, by="size")

# remove extraneous columns
RLS_spec <- RLS_spec %>%
  select(Site.Name, Habitat, Year, Depth, Method, Block, Species, total, Lbin)

# sum occurrences of size classes
RLS_dist <- RLS_spec %>%
  group_by(Site.Name, Habitat, Year, Lbin) %>%
  summarise(sum(total))

# rename summarised column
names(RLS_dist)[names(RLS_dist)=="sum(total)"] <- "total"

#reorder soler factors
#RLS_dist$Lbin <- factor(RLS_dist$Lbin, levels = c("<7.5","7.5-30",">30"))
#reorder 10 bin factors
RLS_dist$Lbin <- factor(RLS_dist$Lbin, levels = c("2.5-4.9","5-9.9","10-19.9","20-39.9","40-74.9","70-124.5","125-187.4","187.5-249.9","250-399.9","400"))

ggplot(RLS_dist, aes(Habitat, log10(total + 1))) +
  geom_bar(aes(fill = Lbin), position = "dodge", stat="identity") +
  facet_grid(.~Year)

###################################################################################
# MULTIVARIATE DATA                                                               #
###################################################################################

# Diversity ordination

# TAXONOMIC

# spread RLS data into community matrix
RLS_full_mat <- RLS_1_2 %>%
  spread(Species, Total, fill = 0)
RLS_mat <- RLS_full_mat[,c(6, 9, 11, 50:213)]
RLS_mat <- RLS_mat %>%
  group_by(Site.Name, Habitat, Year) %>%
  summarise_each(funs(sum))

# run matrix 
RLS_mat_run <- RLS_mat[,4:167]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
plot_data<- data.frame(RLS_mat[,1:3], RLS_mds_points)
ggplot(plot_data, aes(x=MDS1, y=MDS2, fill = Habitat)) +
  geom_point(aes(x=MDS1, y=MDS2, color = Habitat), size = 4)

# FUNCTIONAL 

# import character list
characters <- read.csv("Traits_all-species_edit.csv")

# fix level duplication
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="planktivore"] <- "Planktivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="higher carnivore"] <- "Higher carnivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="omnivore"] <- "Omnivore"
levels(characters$Trophic.group)

# match characters to RLS data
names(characters)[names(characters)=="CURRENT_TAXONOMIC_NAME"] <- "Species"
RLS_func_1_2 <- left_join(RLS_1_2, characters, by="Species")

# replace species names with trophic group
RLS_func_1_2$Species <- RLS_func_1_2$Trophic.group

# identify species not in character list
unique(RLS_1_2$Species[!(RLS_1_2$Species %in% characters$Species)])



# all fishes

# spread RLS data into community matrix
RLS_fullfunc_mat <- RLS_func_1_2 %>%
  spread(Species, Total, fill = 0)
RLS_func_mat <- RLS_fullfunc_mat[,c(6, 9, 11, 61:69)]
RLS_func_mat <- RLS_func_mat %>%
  group_by(Site.Name, Habitat.x, Year) %>%
  summarise_each(funs(sum))

# run matrix 
RLS_mat_run <- RLS_func_mat[,4:12]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
vec.sp<-envfit(RLS_mds$points, RLS_mat_run, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$group<-rownames(vec.sp.df)
plot_data<- data.frame(RLS_func_mat[,1:3], RLS_mds_points)
ggplot(plot_data, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(x=MDS1, y=MDS2, color = plot_data$Habitat, size = 3)) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
      arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit.aes=FALSE) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=group),size=5) +
  coord_fixed()


# Pelagic fish (start from before all fishes)

# filter to method 1
RLS_func_1 <- RLS_func_1_2 %>%
  filter(Method == "1")

# spread RLS data into community matrix
RLS_fullfunc_mat1 <- RLS_func_1 %>%
  spread(Species, Total, fill = 0)
RLS_func_mat1 <- RLS_fullfunc_mat1[,c(6, 9, 11, 61:68)]
RLS_func_mat1 <- RLS_func_mat1 %>%
  group_by(Site.Name, Habitat.x, Year) %>%
  summarise_each(funs(sum))

# run matrix 
RLS_mat_run1 <- RLS_func_mat1[,4:11]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run1)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
vec.sp<-envfit(RLS_mds$points, RLS_mat_run1, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$group<-rownames(vec.sp.df)
plot_data<- data.frame(RLS_func_mat1[,1:3], RLS_mds_points)
ggplot(plot_data, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(x=MDS1, y=MDS2, color = plot_data$Habitat, size = 3)) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit.aes=FALSE) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=group),size=5) +
  coord_fixed()

# Cryptic Fishes


# filter to method 2
RLS_func_2 <- RLS_func_1_2 %>%
  filter(Method == "2")

# spread RLS data into community matrix
RLS_fullfunc_mat2 <- RLS_func_2 %>%
  spread(Species, Total, fill = 0)
RLS_func_mat2 <- RLS_fullfunc_mat2[,c(6, 9, 11, 61:67)]
RLS_func_mat2 <- RLS_func_mat2 %>%
  group_by(Site.Name, Habitat.x, Year) %>%
  summarise_each(funs(sum))

# run matrix 
RLS_mat_run2 <- RLS_func_mat2[,4:10]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run2)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
vec.sp<-envfit(RLS_mds$points, RLS_mat_run2, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$group<-rownames(vec.sp.df)
plot_data<- data.frame(RLS_func_mat2[,1:3], RLS_mds_points)
ggplot(plot_data, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(x=MDS1, y=MDS2, color = plot_data$Habitat, size = 3)) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit.aes=FALSE) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=group),size=5) +
  coord_fixed()




# FUNCTIONAL REEF ONLY

# constrain to reef
RLS_reef <- RLS_1_2 %>%
  filter(Habitat == c("Fore Reef", "Patch Reef"))

# import character list
characters <- read.csv("Traits_all-species_edit.csv")

# fix level duplication
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="planktivore"] <- "Planktivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="higher carnivore"] <- "Higher carnivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="omnivore"] <- "Omnivore"
levels(characters$Trophic.group)

# match characters to RLS data
names(characters)[names(characters)=="CURRENT_TAXONOMIC_NAME"] <- "Species"
RLS_func_1_2 <- left_join(RLS_reef, characters, by="Species")

# replace species names with trophic group
RLS_func_1_2$Species <- RLS_func_1_2$Trophic.group

# identify species not in character list
unique(RLS_1_2$Species[!(RLS_1_2$Species %in% characters$Species)])

# spread RLS data into community matrix M1/M2
RLS_fullfunc_mat <- RLS_func_1_2 %>%
  spread(Species, Total, fill = 0)
RLS_func_mat <- RLS_fullfunc_mat[,c(6, 9, 11, 61:69)]
RLS_func_mat <- RLS_func_mat %>%
  group_by(Site.Name, Habitat.x, Year) %>%
  summarise_each(funs(sum))

# run matrix 
RLS_mat_run <- RLS_func_mat[,4:12]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
vec.sp<-envfit(RLS_mds$points, RLS_mat_run, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$group<-rownames(vec.sp.df)
plot_data<- data.frame(RLS_func_mat[,1:3], RLS_mds_points)
ggplot(plot_data, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(x=MDS1, y=MDS2, color = plot_data$Habitat, size = 3)) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit.aes=FALSE) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=group),size=5) +
  coord_fixed()



# spread RLS data into community matrix M1 only
RLS_fullfunc_mat1 <- RLS_func_1_2 %>%
  filter(Method == "1") %>%
  spread(Species, Total, fill = 0)
RLS_func_mat1 <- RLS_fullfunc_mat1[,c(6, 9, 11, 61:68)]
RLS_func_mat1 <- RLS_func_mat1 %>%
  group_by(Site.Name, Habitat.x, Year) %>%
  summarise_each(funs(sum))

# run matrix 
RLS_mat_run1 <- RLS_func_mat1[,4:11]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run1)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
vec.sp<-envfit(RLS_mds$points, RLS_mat_run, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$group<-rownames(vec.sp.df)
plot_data<- data.frame(RLS_func_mat1[,1:3], RLS_mds_points)
ggplot(plot_data, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(x=MDS1, y=MDS2, color = plot_data$Habitat, size = 3)) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit.aes=FALSE) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=group),size=5) +
  coord_fixed()


# spread RLS data into community matrix M1 only
RLS_fullfunc_mat2 <- RLS_func_1_2 %>%
  filter(Method == "2") %>%
  spread(Species, Total, fill = 0)
RLS_func_mat2 <- RLS_fullfunc_mat2[,c(6, 9, 11, 61:66)]
RLS_func_mat2 <- RLS_func_mat2 %>%
  group_by(Site.Name, Habitat.x, Year) %>%
  summarise_each(funs(sum))

# run matrix 
RLS_mat_run2 <- RLS_func_mat2[,4:9]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run2)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
vec.sp<-envfit(RLS_mds$points, RLS_mat_run, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$group<-rownames(vec.sp.df)
plot_data<- data.frame(RLS_func_mat2[,1:3], RLS_mds_points)
ggplot(plot_data, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(x=MDS1, y=MDS2, color = plot_data$Habitat, size = 3)) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit.aes=FALSE) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=group),size=5) +
  coord_fixed()


###################################################################################
# END OF SCRIPT                                                                   #
###################################################################################
