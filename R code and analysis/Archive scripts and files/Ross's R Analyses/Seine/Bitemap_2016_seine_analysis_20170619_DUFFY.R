###################################################################################
#                                                                                ##
# Ocean Bitemap 2016 Seine data: Analyses                                        ##
# Data are current as of 2017-06-19                                              ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# Based on script originated by Ross Whippo, modified by Emmett Duffy            ##
# Last updated 2017-06-22                                                        ##
#                                                                                ##
################################################################################### 

# TO DO

# draft first model using logistic regression (see former bitemap scripts from Jon)
# get taxonomic data for each fish species
# extract bio-oracle data on temp, salinity

# METADATA:

# This script analyzes data on fish assemblages from the Ocean Bitemap 2016 project led 
# by the Smithsonian's MarineGEO program, Tennenbaum Marine Observatories Network. It focuses 
# fishes sampled by (somewhat) standardized seine hauls paired with squidpop deployments.  
# The datasets were delivered by participating sceintists via email to the MarineGEO 
# database, currently housed on the SI M Drive. 

# RECENT CHANGES

# 20170621 Added Paul York's spreadsheet on fish taxonomy and functional traits
# 20170614 New script started by selecting chunks of file: Bitemap_2016_seine_analysis_20170613.R
# 20170613 First version edited by JED


###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# LOAD PACKAGES                                                                   #
# PREPARE DATASET: READ IN DATA: FISH                                             #
# PREPARE DATASET: READ IN DATA: SQUIDPOPS                                        #
# OBTAIN SUMMARY STATISTICS BY SEINE HAUL                                         #
# ASSEMBLE DATA FRAME: FISH SPECIES AND SUMMARY METRICS PER SEINE                 #
# ASSEMBLE DATA FRAME: COMBINED FISH DATA AND SQUIDPOP DATA                       #
# PLOTS                       #
# EXTRACT BIO-ORACLE DATA                                                         #
# ADD GEOGRAPHIC VARIABLES                                                        #
#                                                                                 #
#                                                                                 #
###################################################################################


###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# load required packages
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
# detach(package:plyr)    
library(dplyr)
library(car)
library(vegan)
library(reshape)
library(reshape2)
library(psych)


###################################################################################
# PREPARE DATASET: READ IN DATA: FISH                                             #
###################################################################################

# Data files are in folder: 
#  smb://SI-NAS1.SMB.US.SINET.SI.EDU/MarineGEO-TMON/MarineGEO Database/Network Activities/Bitemap/2016/Bitemap Fish Survey/

# Read in survey data (last updated 20170221)
bitemap.seine.full <- read.csv('Bitemap_Seine_ALL-DATA_20170616_JED.csv')

# create genus-species column
# NOTE: unite command requires library tidyr
bitemap.seine.full <- unite(bitemap.seine.full, "genus.species", c(Genus,Species), sep = "_", remove = FALSE)

# create unique event ID column
bitemap.seine.full <- unite(bitemap.seine.full, "ID", c(Site.Name,Seagrass.Unveg,Date), sep ="_", remove = FALSE)
bitemap.seine.full$ID <- as.factor(bitemap.seine.full$ID)

# create site/habitat type ID
bitemap.seine.full <- unite(bitemap.seine.full, "Site.Type", c(Site.Name, Seagrass.Unveg), sep = "_", remove= FALSE)
bitemap.seine.full$Site.Type <- as.factor(bitemap.seine.full$Site.Type)

# Create overarching site name to group paired sites 
# NOTE: revalue requires detaching 
bitemap.seine.full$site.paired <- bitemap.seine.full$Site.Name
bitemap.seine.full$site.paired <- revalue(bitemap.seine.full$site.paired, 
  c("Cymodocea_Moreton_Bay" = "Moreton_Bay", "Unvegetated_Moreton_Bay" = "Moreton_Bay",
  "Choked Pocket" = "Choked", "Choked Sandspit" = "Choked",
  "Redondo N-3" = "Redondo", "Redondo N-4" = "Redondo", "Redondo N1" = "Redondo", 
  "Redondo N2" = "Redondo", "Redondo S-5" = "Redondo", "Redondo S-6" = "Redondo",
  "Urunga Lagoon 1" = "Urunga Lagoon", "Urunga Lagoon 2" = "Urunga Lagoon", "Urunga Lagoon 3" = "Urunga Lagoon", 
  "Urunga Lagoon 4" = "Urunga Lagoon", "Urunga Lagoon 5" = "Urunga Lagoon", "Urunga Lagoon 6" = "Urunga Lagoon",     
  "Urunga Lagoon 7" = "Urunga Lagoon", "Urunga Lagoon 8" = "Urunga Lagoon", "Urunga Lagoon 9" = "Urunga Lagoon" 
))

# make Length numeric
bitemap.seine.full$Length <- as.character(bitemap.seine.full$Length)
bitemap.seine.full$Length <- as.numeric(bitemap.seine.full$Length)

# Rename habitat factors
# NOTE: library plyr is needed for 'revalue'
bitemap.seine.full$Seagrass.Unveg <- revalue(bitemap.seine.full$Seagrass.Unveg, c("Seagrass" = "vegetated", "Seagrass " = "vegetated", "unveg" = "unvegetated", "Unveg" = "unvegetated", "Unvegetated" = "unvegetated", "Sand" = "unvegetated"))
bitemap.seine.full$Seagrass.Unveg <- as.factor(bitemap.seine.full$Seagrass.Unveg)

# relevel factors for Seagrass/Unveg
bitemap.seine.full$Seagrass.Unveg <- ordered(bitemap.seine.full$Seagrass.Unveg, levels = c("vegetated", "unvegetated"))

# replace NAs in abundance with 1
bitemap.seine.full$Abundance[is.na(bitemap.seine.full$Abundance)] <- 1

# Review species names for syntactical glitches
str(bitemap.seine.full)
# bitemap.seine.full$genus.species <- as.factor(bitemap.seine.full$genus.species)
# levels(bitemap.seine.full$genus.species)

# Rename misspelled or non-standard variables
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Acanthaluteres_Acanthaluteres spilomelanurus"] <- "Acanthaluteres_spilomelanurus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Albula _vulpes"] <- "Albula_vulpes"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherina _hepsetus"] <- "Atherina_hepsetus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherinops _affinis"] <- "Atherinops_affinis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherinops_affins"] <- "Atherinops_affinis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Atherinosoma_Atherinosoma microstoma"] <- "Atherinosoma_microstoma"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Cheilodactylus _variegatus"] <- "Cheilodactylus_variegatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Carcinus_Carcinus maenas"] <- "Carcinus_maenas"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Crangon _septemspinosa"] <- "Crangon_septemspinosa"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ctenogobius_boteosoma"] <- "Ctenogobius_boleosoma"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Cymatogaster _aggregata"] <- "Cymatogaster_aggregata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Cynoscion _nebulosus"] <- "Cynoscion_nebulosus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ditema _temminckii"] <- "Ditema_temminckii"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Elops _saurus"] <- "Elops_saurus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Embiotoca_ jacksoni"] <- "Embiotoca_jacksoni"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Embiotocidae_sp "] <- "Embiotocidae_spp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Engraulis _ringens"] <- "Engraulis_ringens"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Epinephelus_sp. 1"] <- "Epinephelus_sp1"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Etropus _longimanus"] <- "Etropus_longimanus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Farfantepenaus_spp."] <- "Farfantepenaeus_spp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Farfantepenaeus_spp."] <- "Farfantepenaeus_spp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Fundulus _similis"] <- "Fundulus_similis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gasterosteus_aculeatus aculeatus"] <- "Gasterosteus_aculeatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Girella_Girella zebra"] <- "Girella_zebra"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobies_species complex"] <- "Gobies_species_complex"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae_sp. 1"] <- "Gobiidae_sp1"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae_sp. 2"] <- "Gobiidae_sp2"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae_sp. 3"] <- "Gobiidae_sp3"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiidae_sp. 4"] <- "Gobiidae_sp4"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma _bosci"] <- "Gobiosoma_bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma_sp."] <- "Gobiosoma_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma _bosci"] <- "Gobiosoma_bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma _bosci"] <- "Gobiosoma_bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobiosoma_bosc"] <- "Gobiosoma_bosci"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gobius _geniporus"] <- "Gobius_geniporus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gymnapistes_Gymnapistes marmoratus"] <- "Gymnapistes_marmoratus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Heteroclinus_Heteroclinus perspicillatus"] <- "Heteroclinus_perspicillatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hippocampus _sp."] <- "Hippocampus_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hipsoblennius _sordidus"] <- "Hypsoblennius_sordidus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hypsoblennius_ jenkinsi"] <- "Hypsoblennius_jenkinsi"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Lagodon _rhomboides"] <- "Lagodon_rhomboides"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Leiostomus _xanthurus"] <- "Leiostomus_xanthurus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Leptocottus _armatus"] <- "Leptocottus_armatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Liza _argentea"] <- "Liza_argentea"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Lutjanus _russellii"] <- "Lutjanus_russellii"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menidia _menidia"] <- "Menidia_menidia"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menidia_beryllinia"] <- "Menidia_beryllina"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menidia _beryllina"] <- "Menidia_beryllina"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Menticirrhus  _ophicephalus"] <- "Menticirrhus_ophicephalus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Metacarcinus _magister"] <- "Metacarcinus_magister"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Micrometrus_ minimus"] <- "Micrometrus_minimus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Nesogobius_Nesogobius maccullochi"] <- "Nesogobius_maccullochi"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Odontesthes_ regia"] <- "Odontesthes_regia"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Orthopistis_chrysoptera"] <- "Orthopristis_chrysoptera"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Palaemon _elegans"] <- "Palaemon_elegans"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Palaemonetes _spp."] <- "Palaemonetes_spp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Palaemonetes_spp."] <- "Palaemonetes_spp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Paralichthys _californicus"] <- "Paralichthys_californicus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Paralichthys_ microps"] <- "Paralichthys_microps"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Petroscirtis _variabilis"] <- "Petroscirtis_variabilis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Phanerodon _furcatus"] <- "Phanerodon_furcatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pholis _ornata"] <- "Pholis_ornata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Platycephalidae_sp. 1"] <- "Platycephalidae_sp1"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Platycephalus_Platycephalus aurimaculatus"] <- "Platycephalus_aurimaculatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pomatoschistus_sp."] <- "Pomatoschistus_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Prionotus _punctatus"] <- "Prionotus_punctatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Processa _macrophthalma"] <- "Processa_macrophthalma"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Psettichthy_melanostictus"] <- "Psettichthys_melanostictus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Rhinobatos _productus"] <- "Rhinobatos_productus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Romaleon _polyodon"] <- "Romaleon_polyodon"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Scartichthys _viridis"] <- "Scartichthys_viridis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "scorpaenichthys _marmoratus"] <- "Scorpaenichthys_marmoratus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "sculpin_sp."] <- "sculpin_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sebastes_sp"] <- "Sebastes_spp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Siganus _fuscescens"] <- "Siganus_fuscescens"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sillago_ ciliata"] <- "Sillago_ciliata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sillago_maculatus"] <- "Sillago_maculata"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Solea _solea"] <- "Solea_solea"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Stigmatopora_Stigmatopora argus"] <- "Stigmatopora_argus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sygnathus _floridae"] <- "Syngnathus_floridae"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sygnathus_floridae"] <- "Syngnathus_floridae"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sygnathus_griseolineatus"] <- "Syngnathus_griseolineatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus _floridae"] <- "Syngnathus_floridae"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus _leptorhynchus\xe6"] <- "Syngnathus_leptorhynchus_xe6"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus _leptorhynchus\xca"] <- "Syngnathus_leptorhynchus_xca"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus_sp."] <- "Syngnathus_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Tetractenos_Tetractenos glaber"] <- "Tetractenos_glaber"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Torquiginer_sp.r"] <- "Torquiginer_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Umbrina _roncador"] <- "Umbrina_roncador"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Unknown goby_sp "] <- "Unknown_goby_sp"

# From Paul York's work
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ambassis_marinus"] <- "Ambassis_marianus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Arothron_hipidus"] <- "Arothron_hispidus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Centropogon _australis"] <- "Centropogon_australis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ditema_temminckii"] <- "Ditrema_temminckii"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Ephinephelus_coioides"] <- "Epinephelus_coioides"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gasterosteus_aculateus"] <- "Gasterosteus_aculeateus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Gerres_subfasciatus "] <- "Gerres_subfasciatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Herring_sp"] <- "Clupeidae_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Hyporhamphus_regularis ardelio"] <- "Hyporhamphus_regularis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Monacanthis_chinensis"] <- "Monacanthus_chinensis"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Oncorhyncus_tshawytscha"] <- "Oncorhynchus_tshawytscha"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Paralichthys_albiguttata"] <- "Paralichthys_albigutta"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pelates_sexlineatus"] <- "Helotes_sexlineatus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pipefish_sp"] <- "Syngnathidae_sp"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pleuronectes_platessa "] <- "Pleuronectes_platessa"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Pseudoblennius_dottoides"] <- "Pseudoblennius_cottoides"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Rhombosolea_Rhombosolea tapirina"] <- "Rhombosolea_tapirina"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Sillago_temminck"] <- "Sillago_japonica"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Suarida_nebulosa"] <- "Saurida_nebulosa"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus_griseolineatus"] <- "Ambassis_marianus"
bitemap.seine.full$genus.species[bitemap.seine.full$genus.species == "Syngnathus_leptorhynchus_xe6"] <- "Syngnathus_leptorhynchus"

# Read in fish families and trophic groups
fish_trophic_groups <- read.csv('bitemap_fish_families_trophic_groups_20170621.csv')
unique(fish_trophic_groups$trophic.group)

# streamline trophic group designations
fish_trophic_groups$trophic.group.simple <- fish_trophic_groups$trophic.group

levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="benthic.invertivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="benthic.invertivore.detritivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="benthic.invertivore.herbivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="benthic.invertivore.herbivore.detritivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="benthic.invertivore.piscivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="benthic.invertivore.planktivore"] <- "carnivore"

levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="detritivore.herbivore"] <- "herbivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="detritivore.planktivore.benthic.invertivore"] <- "planktivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="detritivore.benthic.invertivore"] <- "carnivore"

levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="herbivore.benthic.invertivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="Herbivore.benthic.invertivore.piscivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="herbivore.benthic.invertivore.planktivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="piscivore.benthic.invertivore"] <- "carnivore"

levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="Planktivore"] <- "planktivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="planktivore.benthic.invertivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="planktivore.benthic.invertivore.piscivore"] <- "carnivore"
levels(fish_trophic_groups$trophic.group.simple)[levels(fish_trophic_groups$trophic.group.simple)=="planktivore.detritivore"] <- "planktivore"

levels(fish_trophic_groups$trophic.group.simple)

# match fish families and trophic groups to species
bitemap.seine.full <- left_join(bitemap.seine.full, fish_trophic_groups, by = "genus.species")


###################################################################################
# PREPARE DATASET: READ IN DATA: SQUIDPOPS                                        #
###################################################################################

#  read in survey data (last updated 20170406)
bitemap.squidpop <- read.csv('Ocean Bitemap Squidpop Data (Raw) downloaded 20170616.csv')

str(bitemap.squidpop)

# change numbers to numeric
bitemap.squidpop$Proportion.Missing..24.hours. <- as.numeric(as.character(bitemap.squidpop$Proportion.Missing..24.hours.))
bitemap.squidpop$Proportion.Missing..1.hour. <- as.numeric(as.character(bitemap.squidpop$Proportion.Missing..1.hour.))


# subset and rename habitat factors
bitemap.squidpop <- bitemap.squidpop[bitemap.squidpop$Type.Of.Habitat %in% c("Seagrass", "Seagrass Meadow", "Sandy Bottom", "Muddy Bottom", "Rocky Reef"), ]
# NOTE: library plyr is needed for 'revalue'
bitemap.squidpop$Type.Of.Habitat <- revalue(bitemap.squidpop$Type.Of.Habitat, c("Seagrass Meadow" = "vegetated", "Sandy Bottom" = "unvegetated", "Muddy Bottom" = "unvegetated", "Rocky Reef" = "unvegetated"))
bitemap.squidpop$Type.Of.Habitat <- as.factor(bitemap.squidpop$Type.Of.Habitat)

names(bitemap.squidpop)

# Rename wonky variables
names(bitemap.squidpop)[names(bitemap.squidpop) == "Latitude..GPS..decimal."] <- "latitude"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Longitude"] <- "longitude"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Date.Squidpops.Deployed..yyyymmdd."] <- "date.deployed"
names(bitemap.squidpop)[names(bitemap.squidpop) == "If.You.Chose..Other...Please.Describe.The.Habitat"] <- "other.habitat.describe"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Number.Of.Squidpops.Deployed"] <- "n.squidpops.deployed"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Amount.Of.Bait.Missing.After.One.Hour"] <- "squid.missing.1hr"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Number.Of.Squidpops.Retrieved.After.24.Hours"] <- "n.pops.retrieved.24hr"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Amount.Of.Bait.Missing.After.24.Hours"] <- "squid.missing.24hr"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Time.Squidpops.Deployed..hhmm..24.hour.clock."] <- "time.deployed"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Proportion.Missing..1.hour."] <- "proportion.missing.1hr"
names(bitemap.squidpop)[names(bitemap.squidpop) == "Proportion.Missing..24.hours."] <- "proportion.missing.24hr"

# Export the data frames
write.csv(bitemap.squidpop, "bitemap.squidpop.20170622.csv", row.names = F)


###################################################################################
# OBTAIN SUMMARY STATISTICS BY SEINE HAUL                                         #
###################################################################################

# Filter out inverts
# NOTE: Requires library: dplyr
bitemap.seine.fish <- filter(bitemap.seine.full, Phylum != "Arthropoda"| is.na(Phylum))
bitemap.seine.fish <- filter(bitemap.seine.fish, Phylum != "Mollusca"| is.na(Phylum))
bitemap.seine.fish$Phylum <- factor(bitemap.seine.fish$Phylum)

# Calculate fish richness per seine haul
# detach(package:plyr)
richness_seine <- bitemap.seine.fish %>%
  group_by(ID, site.paired, Site.Type, Site.Name, Seagrass.Unveg) %>%
  summarise(richness.seine = length(unique(genus.species))) %>%
  arrange(Site.Type)

hist(richness_seine$richness.seine) # right-skewed
richness_seine$log10.richness.seine <- log10(richness_seine$richness.seine)
hist(richness_seine$log10.richness.seine) # better
nrow(richness_seine) # 183


# Calculate total fish abundance per seine haul
fish_abundance_seine <- bitemap.seine.fish %>%
  group_by(ID, site.paired, Site.Type, Site.Name, Date, Seagrass.Unveg) %>%
  summarise(fish.abundance.seine = sum(Abundance)) 
# names(fish_abundance_seine)[names(fish_abundance_seine)=="round(Lat, 2)"] <- "Lat"
nrow(fish_abundance_seine) # 183

hist(fish_abundance_seine$fish.abundance.seine) # right-skewed
fish_abundance_seine$log10.fish.abundance.seine <- log10(fish_abundance_seine$fish.abundance.seine + 1)
hist(fish_abundance_seine$log10.fish.abundance.seine) # nice


# Calculate abundances by trophic group per seine haul
# Requires library(plyr)
trophic_group_abundance_seine <- ddply(bitemap.seine.fish, c("ID","trophic.group.simple"), summarise, abundance.trophic.group= sum(Abundance))

# Create a seine x trophic group matrix
# NOTE: dcast requires library: reshape2
trophic_group_per_seine <- dcast(trophic_group_abundance_seine, ID ~ trophic.group.simple,  sum,  value.var = "abundance.trophic.group")
trophic_group_per_seine[is.na(trophic_group_per_seine)] <- 0
nrow(trophic_group_per_seine) # 183
names(trophic_group_per_seine)

# total number of fish counted
sum(fish_abundance_seine$fish.abundance.seine) # 29373


###################################################################################
# ASSEMBLE DATA FRAME: FISH SPECIES AND SUMMARY METRICS PER SEINE                 #
###################################################################################

# Create a dataframe with desired metadata 
seine_meta  <- bitemap.seine.fish[ ,c("ID", "site.paired", "Site.Name", "Date", "Lat", "Seagrass.Unveg", "genus.species", "Abundance") ] 

# detach plyr
# detach(package:plyr)
# Sum abundances to obtain total for each species per site
seine_fish_sp_abunds <- seine_meta %>%
  group_by(ID, Site.Name, Date, Seagrass.Unveg, genus.species) %>%
  summarise(Abundance = sum(Abundance))
nrow(seine_fish_sp_abunds) # 1001


# Create a seine x species matrix
# NOTE: dcast requires library: reshape2
fish_per_seine <- dcast(seine_fish_sp_abunds, ID ~ genus.species,  sum,  value.var = "Abundance")
fish_per_seine[is.na(fish_per_seine)] <- 0
nrow(fish_per_seine) # 183
names(fish_per_seine)

# Add metadata and summary data to the seine x species matrix
fish_per_seine$Site.Name <- bitemap.seine.fish$Site.Name[match(fish_per_seine$ID, bitemap.seine.fish$ID)]
fish_per_seine$site.paired <- bitemap.seine.fish$site.paired[match(fish_per_seine$ID, bitemap.seine.fish$ID)]
fish_per_seine$Site.Type <- bitemap.seine.fish$Site.Type[match(fish_per_seine$ID, bitemap.seine.fish$ID)]
fish_per_seine$Lat <- bitemap.seine.fish$Lat[match(fish_per_seine$ID, bitemap.seine.fish$ID)]
fish_per_seine$Long <- bitemap.seine.fish$Long[match(fish_per_seine$ID, bitemap.seine.fish$ID)]
fish_per_seine$Seagrass.Unveg <- bitemap.seine.fish$Seagrass.Unveg[match(fish_per_seine$ID, bitemap.seine.fish$ID)]

fish_per_seine$fish.abundance.seine <- fish_abundance_seine$fish.abundance.seine[match(fish_per_seine$ID, fish_abundance_seine$ID)]
fish_per_seine$log10.fish.abundance.seine <- fish_abundance_seine$log10.fish.abundance.seine[match(fish_per_seine$ID, fish_abundance_seine$ID)]
fish_per_seine$fish.richness.seine <- richness_seine$richness.seine[match(fish_per_seine$ID, richness_seine$ID)]
fish_per_seine$log10.fish.richness.seine <- richness_seine$log10.richness.seine[match(fish_per_seine$ID, richness_seine$ID)]

fish_per_seine$carnivores.seine <- trophic_group_per_seine$carnivore[match(fish_per_seine$ID, trophic_group_per_seine$ID)]
fish_per_seine$herbivores.seine <- trophic_group_per_seine$herbivore[match(fish_per_seine$ID, trophic_group_per_seine$ID)]
fish_per_seine$planktivores.seine <- trophic_group_per_seine$planktivore[match(fish_per_seine$ID, trophic_group_per_seine$ID)]
fish_per_seine$piscivores.seine <- trophic_group_per_seine$piscivore[match(fish_per_seine$ID, trophic_group_per_seine$ID)]

fish_per_seine$log10.carnivores.seine <- log10(fish_per_seine$carnivores.seine +1)
fish_per_seine$log10.herbivores.seine <- log10(fish_per_seine$herbivores.seine +1)
fish_per_seine$log10.planktivores.seine <- log10(fish_per_seine$planktivores.seine +1)
fish_per_seine$log10.piscivores.seine <- log10(fish_per_seine$piscivores.seine +1)

names(fish_per_seine)

# Reorder columns so metadata and summary data are up front 
fish_per_seine <- fish_per_seine[,c(1, 281:298, 2:280)]

# Export the data frames
write.csv(fish_per_seine, "fish_per_seine_20170622.csv", row.names = F)


###################################################################################
# ASSEMBLE DATA FRAME: COMBINED FISH DATA AND SQUIDPOP DATA                       #
###################################################################################

nrow(fish_per_seine) # 183
nrow(bitemap.squidpop) # 294

# Because seine samples and squidpop deployments did not retain unique identifiers that 
# connected them, it has been a process of slow, painstaking manual detective work to come 
# up with a way to pair them for analysis. I have taken the unique "ID" from the fish_per_seine
# data frame and manually pasted that into a corresponding column in the bitemap.squidpop 
# data frame. Now we have to reimport the doctored dataframe so that we can compile the 
# two together and run analyses. 

# However, note that the two data frames need only be connected by "site", each of which
# should have up to 3 replicate seines and 3 replocate squidpop deployments in each of 
# vegetated and unvegetated habitats. 

# Read in modified squidpop dataframe
bitemap_data_modified <- read.csv('bitemap.squidpop.20170615.modified.csv')
names(bitemap_data_modified)
nrow(bitemap_data_modified) # 167

names(fish_per_seine)

# Create new combined data frame with squidpop and fish seine data
bitemap2016_data <- fish_per_seine

# Add squidpop data to combined dataframe:
bitemap2016_data$Name <- bitemap_data_modified$Name[match(bitemap2016_data$ID, bitemap_data_modified$ID)]
bitemap2016_data$n.squidpops.deployed <- bitemap_data_modified$n.squidpops.deployed[match(bitemap2016_data$ID, bitemap_data_modified$ID)]
bitemap2016_data$squid.missing.1hr <- bitemap_data_modified$squid.missing.1hr[match(bitemap2016_data$ID, bitemap_data_modified$ID)]
bitemap2016_data$n.pops.retrieved.24hr <- bitemap_data_modified$n.pops.retrieved.24hr[match(bitemap2016_data$ID, bitemap_data_modified$ID)]
bitemap2016_data$squid.missing.24hr <- bitemap_data_modified$squid.missing.24hr[match(bitemap2016_data$ID, bitemap_data_modified$ID)]
bitemap2016_data$seine.data <- bitemap_data_modified$seine.data[match(bitemap2016_data$ID, bitemap_data_modified$ID)]
bitemap2016_data$latitude.absolute <- abs(bitemap2016_data$Lat)

# Calculate variables for squidpop loss
bitemap2016_data$proportion.missing.1hr <- bitemap2016_data$squid.missing.1hr / bitemap2016_data$n.squidpops.deployed
bitemap2016_data$proportion.missing.24hr <- bitemap2016_data$squid.missing.24hr / bitemap2016_data$n.pops.retrieved.24hr

# # Create overarching site name to group paired sites 
# # NOTE: revalue requires library(plyr)
# bitemap2016_data$site.paired <- bitemap2016_data$Site.Name
# bitemap2016_data$site.paired <- revalue(bitemap2016_data$site.paired, 
#   c("Cymodocea_Moreton_Bay" = "Moreton_Bay", "Unvegetated_Moreton_Bay" = "Moreton_Bay",
#   "Choked Pocket" = "Choked", "Choked Sandspit" = "Choked",
#   "Redondo N-3" = "Redondo", "Redondo N-4" = "Redondo", "Redondo N1" = "Redondo", 
#   "Redondo N2" = "Redondo", "Redondo S-5" = "Redondo", "Redondo S-6" = "Redondo",
#   "Urunga Lagoon 1" = "Urunga Lagoon", "Urunga Lagoon 2" = "Urunga Lagoon", "Urunga Lagoon 3" = "Urunga Lagoon", 
#   "Urunga Lagoon 4" = "Urunga Lagoon", "Urunga Lagoon 5" = "Urunga Lagoon", "Urunga Lagoon 6" = "Urunga Lagoon",     
#   "Urunga Lagoon 7" = "Urunga Lagoon", "Urunga Lagoon 8" = "Urunga Lagoon", "Urunga Lagoon 9" = "Urunga Lagoon" 
# ))

names(bitemap2016_data)
nrow(bitemap2016_data) # 183

# Reorder columns so relevant data are up front 
bitemap2016_data <- bitemap2016_data[,c(1:19, 299:307, 20:298)]

# Create separate data set with only vegetated
bitemap2016_data.vegetated <- droplevels(subset(bitemap2016_data, Seagrass.Unveg == "vegetated"))
nrow(bitemap2016_data.vegetated) # 89
bitemap2016_data.unvegetated <- droplevels(subset(bitemap2016_data, Seagrass.Unveg == "unvegetated"))
nrow(bitemap2016_data.unvegetated) # 94


# Obtain mean values per site
bitemap2016_site_means <- ddply(bitemap2016_data, c("site.paired", "Seagrass.Unveg"), summarize, 
  fish.abundance.mean = mean(fish.abundance.seine, na.rm = T), 
  log10.fish.abundance.mean = mean(log10.fish.abundance.seine, na.rm = T), 
  fish.richness.seine.mean = mean(fish.richness.seine, na.rm = T),                      
  log10.fish.richness.seine.mean = mean(log10.fish.richness.seine, na.rm = T), 
  carnivore.abundance.mean = mean(carnivores.seine, na.rm = T), 
  herbivore.abundance.mean = mean(herbivores.seine, na.rm = T), 
  planktivore.abundance.mean = mean(planktivores.seine, na.rm = T), 
  piscivore.abundance.mean = mean(piscivores.seine, na.rm = T), 
  log10.carnivore.abundance.mean = mean(log10.carnivores.seine, na.rm = T), 
  log10.herbivore.abundance.mean = mean(log10.herbivores.seine, na.rm = T), 
  log10.planktivore.abundance.mean = mean(log10.planktivores.seine, na.rm = T), 
  log10.piscivore.abundance.mean = mean(log10.piscivores.seine, na.rm = T), 
  proportion.missing.1hr.mean = mean(proportion.missing.1hr, na.rm = T),
  proportion.missing.24hr.mean = mean(proportion.missing.24hr, na.rm = T)
)

names(bitemap2016_site_means)
nrow(bitemap2016_site_means) # 63

# Add metadata back in 
bitemap2016_site_means$latitude <- bitemap.seine.fish$Lat[match(bitemap2016_site_means$site.paired, bitemap.seine.fish$site.paired)]
bitemap2016_site_means$longitude <- bitemap.seine.fish$Long[match(bitemap2016_site_means$site.paired, bitemap.seine.fish$site.paired)]
bitemap2016_site_means$Data.Collector <- bitemap.seine.fish$Data.Collector[match(bitemap2016_site_means$site.paired, bitemap.seine.fish$site.paired)]
bitemap2016_site_means$latitude.absolute <- abs(bitemap2016_site_means$latitude)

write.csv(bitemap2016_site_means, "bitemap2016_site_means_20170622.csv", row.names = F)

# Create separate data set with only vegetated
bitemap2016_site_means_vegetated <- droplevels(subset(bitemap2016_site_means, Seagrass.Unveg == "vegetated"))
nrow(bitemap2016_site_means_vegetated) # 32
bitemap2016_site_means_unvegetated <- droplevels(subset(bitemap2016_site_means, Seagrass.Unveg == "unvegetated"))
nrow(bitemap2016_site_means_unvegetated) # 31


###################################################################################
# PLOTS                       #
###################################################################################

pairs.panels(bitemap2016_site_means[,c("latitude", "latitude.absolute", "proportion.missing.1hr.mean", "proportion.missing.24hr.mean", 
  "log10.fish.richness.seine.mean", "log10.fish.abundance.mean","log10.carnivore.abundance.mean", 
  "log10.herbivore.abundance.mean", "log10.planktivore.abundance.mean", "log10.piscivore.abundance.mean")], 
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 4)


plot(log10.fish.abundance.site ~ Seagrass.Unveg, data = bitemap2016_site_means)
plot(log10.fish.richness.site ~ Seagrass.Unveg, data = bitemap2016_site_means)
plot(proportion.missing.1hr.site ~ Seagrass.Unveg, data = bitemap2016_site_means)
plot(proportion.missing.24hr.site ~ Seagrass.Unveg, data = bitemap2016_site_means)

plot(log10.fish.richness.site ~ latitude.absolute, data = bitemap2016.site.means)
plot(log10.fish.richness.site ~ latitude.absolute, data = bitemap2016.site.means.vegetated)
plot(log10.fish.richness.site ~ latitude.absolute, data = bitemap2016.site.means.unvegetated)

plot(proportion.missing.1hr.site ~ latitude.absolute, data = bitemap2016.site.means.vegetated)
plot(proportion.missing.1hr.site ~ latitude.absolute, data = bitemap2016.site.means.unvegetated)

plot(proportion.missing.24hr.site ~ latitude.absolute, data = bitemap2016.site.means.vegetated)
plot(proportion.missing.24hr.site ~ latitude.absolute, data = bitemap2016.site.means.unvegetated)

plot(proportion.missing.1hr.site ~ log10.fish.abundance.site, data = bitemap2016.site.means.vegetated)
plot(proportion.missing.1hr.site ~ log10.fish.abundance.site, data = bitemap2016.site.means.unvegetated)

plot(proportion.missing.24hr.site ~ log10.fish.abundance.site, data = bitemap2016.site.means.vegetated)
plot(proportion.missing.24hr.site ~ log10.fish.abundance.site, data = bitemap2016.site.means.unvegetated)

plot(proportion.missing.1hr.site ~ log10.fish.richness.site, data = bitemap2016.site.means.vegetated)
plot(proportion.missing.1hr.site ~ log10.fish.richness.site, data = bitemap2016.site.means.unvegetated)

plot(proportion.missing.24hr.site ~ log10.fish.richness.site, data = bitemap2016.site.means.vegetated)
plot(proportion.missing.24hr.site ~ log10.fish.richness.site, data = bitemap2016.site.means.unvegetated)


###################################################################################
# EXTRACT BIO-ORACLE DATA                                                         #
###################################################################################

# This script adapted from Matt Whalen to extract Bio-ORACLE data 
library(raster)
# read in the raster data
sst.min   <- raster("sstmin.asc")
sst.max   <- raster("sstmax.asc")
sst.mean  <- raster("sstmean.asc")
salinity.BO <- raster("salinity.asc")
nitrate.BO <- raster("nitrate.asc")
# phos.BO <- raster("phosphate.asc")
chlomean.BO <- raster("chlomean.asc")

# stack the rasters into a single object (sort of like a list or array)
oracle.data <- stack(sst.min, sst.max, sst.mean, salinity.BO, 
  nitrate.BO, chlomean.BO)
# For some reason phosphate.BO doesn;t work. Message: "Error in compareRaster(x) : different extent" 

oracle.data  # note that you can see the names of all the constituent rasters in the "names" field

# Cteate a dataframe of bitemap latitude and longitude
bitemap.geo <- fish.per.seine[ ,c("Long", "Lat") ] 
write.csv(bitemap.geo, "bitemap.geo.csv", row.names = F)


# Rename variables
names(bitemap.geo)[names(bitemap.geo)=="Long"] <- "longitude"
names(bitemap.geo)[names(bitemap.geo)=="Lat"] <- "latitude"

## Extract data for each site and add to the ZEN data we already have
# zen2014.data$sst.min <- extract( sst.min, zen.geo, method='bilinear' )
# No Bio-ORACLE data available for Croatia, France subsite B
# fish.per.seine$sst.min   <- unlist(lapply(extract( sst.min, bitemap.geo, method='bilinear', buffer=10000 ),mean,na.rm=T))
# fish.per.seine$sst.max   <- unlist(lapply(extract( sst.max, bitemap.geo, method='bilinear', buffer=10000 ),mean,na.rm=T))
fish.per.seine$sst.mean  <- unlist(lapply(extract( sst.mean, bitemap.geo, method='bilinear', buffer=10000 ),mean,na.rm=T))
fish.per.seine$salinity.BO  <- unlist(lapply(extract( salinity.BO, bitemap.geo, method='bilinear', buffer=10000 ),mean,na.rm=T))
fish.per.seine$nitrate.BO   <- unlist(lapply(extract( nitrate.BO, bitemap.geo, method='bilinear', buffer=10000 ),mean,na.rm=T))
# fish.per.seine$chlomean.BO  <- unlist(lapply(extract( chlomean.BO, bitemap.geo, method='bilinear', buffer=10000 ),mean,na.rm=T))

# write the extracted BioOracle environmental data to new data frame:
write.csv( fish.per.seine, "fish.per.seine_with_BioOracle_20170614.csv", row.names = F )
# Few data for Texas sites. wtf?


###################################################################################
# ADD GEOGRAPHIC VARIABLES                                                        #
###################################################################################

library(sp)
library(rworldmap)
library(maptools)
# For a short introduction type : 	 vignette('rworldmap')

# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees

coords2country = function(points)
{  
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}  

points <- bitemap.geo
coords2country(points) # wtf? Texas = India ???


###################################################################################
# ???                 #
###################################################################################

plot(fish.per.seine$fish.richness ~ fish.per.seine$sst.mean)
summary(lm(fish.per.seine$fish.richness ~ fish.per.seine$sst.mean))

plot(fish.per.seine$log10.fish.richness ~ fish.per.seine$sst.mean)
summary(lm(fish.per.seine$log10.fish.richness ~ fish.per.seine$sst.mean ))

fish.per.seine$sst.mean.sqrd <- (fish.per.seine$sst.mean)^2

plot(fish.per.seine$log10.fish.abundance.total ~ fish.per.seine$sst.mean)
summary(lm(fish.per.seine$log10.fish.abundance.total ~ fish.per.seine$sst.mean + fish.per.seine$sst.mean.sqrd))

