#####################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Preliminary data analysis and exploration
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# updated 2020.01.10
#####################################################################

# load libraries
library(tidyverse)


# set colors
seagrass.color <- "#5ab4ac"    # see http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3 ALSO http://www.color-hex.com/color-palette/31668
unveg.color    <- "#d8b365"


### read  data
predator <- read.csv( "Output Data/Bitemap_SEINE_summaries_20190322.csv", stringsAsFactors = FALSE )
predator <- predator %>% select( -sstmean )
sites <- read.csv( "Output Data/Bitemap_BioORACLE_20190107.csv", stringsAsFactors = FALSE )
sites <- sites %>% select(Country:rate)


# get site means
pred.graph <- left_join( predator, sites )
pred.means <- pred.graph %>%
  group_by( Site, Site.Name, habitat, ENSPIE, richness, meanL, biomass.area,sstmean ) %>% 
  summarize( tot.cpua = sum(cpua,na.rm=T)) %>% 
  group_by( Site, habitat ) %>%
  summarize( sstmean=mean(sstmean,na.rm=T), mean.cpua=mean(tot.cpua,na.rm=T),
                    ENSPIE=mean(ENSPIE,na.rm=T), 
                    richness=mean(richness,na.rm=T), 
                    biomass.area=mean(biomass.area,na.rm=T),meanL=mean(meanL,na.rm=T) )




# axis label size needs to be smaller
labsize = 12

# graph fish abundance by habitat type
a <- ggplot(pred.means, aes(x=sstmean, y=(mean.cpua+0.01), col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',"gold4")) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  scale_y_log10(breaks=c(0.001,0.01,0.1,1)) +
  # scale_x_discrete( labels=c("Seagrass", "Unvegetated")) +
  # scale_fill_manual(values=c(seagrass.color,unveg.color)) + 
  labs( x = "Mean annual SST (°C)", y=expression(paste("Catch per ", m^2)) ) 

# graph predator diversity by habitat type
b <- ggplot(pred.means, aes(x=sstmean, y=ENSPIE, col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',"gold4")) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Mean annual SST (°C)", y=expression(paste("Effective number species")) ) 

# graph fish size (80th percentile) by habitat type
c <- ggplot(pred.means, aes(x=sstmean, y=meanL, col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',"gold4")) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  labs( x = "Mean annual SST (°C)", y=expression(paste("Fish total length (cm)")) ) 

# graph fish biomass by habitat type
d <- ggplot(pred.means, aes(x=sstmean, y=(biomass.area+0.01), col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',"gold4")) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  scale_y_log10(breaks=c(0.1,1,10)) +
  labs( x = "Mean annual SST (°C)", y=expression(paste("Fish biomass (g per ",m^2,")")) ) 

# graph fish richness by habitat type
e <- ggplot(pred.means, aes(x=sstmean, y=richness, col=habitat, fill = habitat)) + 
  geom_smooth(aes(group=habitat), se=T,  lwd=0.5 ) +
  geom_point(pch=21,size=2,alpha=0.5) +
  scale_fill_manual(values=c(seagrass.color,unveg.color)) +
  scale_color_manual(values=c('black',"gold4")) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=labsize),
        legend.position = "none") +
  scale_y_continuous(breaks=c(0,2,4,6,8)) +
  labs( x = "Mean annual SST (°C)", y=expression(paste("Species richness")) ) 

# plot all of it together
windows(12,3)
cowplot::plot_grid( e,b,a,d, labels="AUTO", ncol=4,
                    align="v",
                    label_x = c(0.05,0.10,0.05,0.10) )

