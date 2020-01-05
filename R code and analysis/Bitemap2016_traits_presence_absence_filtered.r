#############################################################################
# MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian
# Ocean Bitemap 2016
# Use traits of consumers in analyses of community stucture, consumption, etc.
# Code by Matt Whalen, Ross Whippo, and Emmett Duffy
# created 30 August 2019
#############################################################################


## load libraries
library(tidyverse) 
library(FD)

# read trait data and community data
s1 <- read_csv( "Output Data/Bitemap_presence+absence_trait.csv" )
trait.use <- read_csv( "Output Data/Bitemap_presence+absence_trait_use.csv" )

# filter taxa based on multivariate analysis
multivar <- read_csv( "Output Data/multivar_constr_taxa.csv" )
fam.pos <- multivar %>% filter( CAP1>0 )
suse <- s1 %>% filter( family %in% fam.pos$family)

# indvidual and taxon based ratios below are currently the same because only one individual per site*habitat combination
# ratio of active to ambush at each site by individuals 
rat.ind <- suse %>%
  dplyr::group_by( Site, habitat ) %>%
  dplyr::summarize(act.ratio.ind = length(act[act %in% c("active","ambush.active")]) / length(sciName)  ) 
# and by taxon
rat.tax <- suse %>%
  dplyr::select( Site, habitat, sciName, act ) %>%
  dplyr::distinct() %>%
  dplyr::group_by( Site, habitat ) %>%
  dplyr::summarize( act.ratio.tax = length(act[act %in% c("active","ambush.active")]) / length(sciName),
                    richness=length(sciName)) 
# merge
rat <- left_join( rat.ind, rat.tax )
rat.site <- left_join( rat, suse )

# write to disk
write_csv( rat, "Output Data/consumer_active_ratio_PA_filtered.csv" )


# plot
ggplot( rat.site, aes(x=abs(meanLat), y=act.ratio.ind)) + geom_point() + geom_smooth()
ggplot( rat.site, aes(x=abs(meanLat), y=act.ratio.tax)) + geom_point() + geom_smooth()
ggplot( rat.site, aes(x=act.ratio.tax,  y=act.ratio.ind)) + geom_point() + geom_smooth()
ggplot( rat.site, aes(x=act.ratio.ind,  y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))
ggplot( rat.site, aes(x=act.ratio.tax,  y=rate)) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial' ))
# there is some relationship here, but not very strong -- just one trait that is not very well identified across all taxa
hilo <- rat.site %>% 
  filter( act.ratio.tax==1, rate<0.5 ) %>% 
  arrange(rate) %>% 
  select( Site ) %>% 
  unlist()
suse %>% 
  filter( Site %in% hilo ) %>% 
  select( Site,  family, sciName ) %>% 
  distinct()

# nice figure
# windows(3,3)
ggplot( rat.site, aes(x=act.ratio.tax,  y=rate)) + 
  geom_point(pch=16,alpha=0.5) +
  xlab("Ratio of active to restful taxa") + ylab("Consumption rate") +
  theme_classic()
summary(glm(rate~act.ratio.tax,data=rat.site,family="quasibinomial"))

## community weighted mean traits
# FD package

# matrix of traits
# do not allow mutliple values per species (take the mode)
# mode function (https://www.r-bloggers.com/computing-the-mode-in-r/)
Mode = function(x){ 
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = names(ta)[ta == tam]
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(paste(mod,collapse = ";"))
}
Mode( c('active','ambush','active') )
Mode( c('active','ambush') )
Mode( c('active','ambush', NA))
Mode( c(0,1, NA))

# trait.mode <- trait %>%
#   dplyr::mutate( phylum=ifelse(family %in% c("Portunidae","Penaeidae","Palaemonidae", 'Cancridae'),
#                         "Arthropoda","Chordata") ) %>%
#   dplyr::group_by( sciName ) %>%
#   dplyr::summarize( phylum = Mode(phylum), activty = Mode(active.ambush), feeding = Mode(feedingImpute), 
#              trophic = Mode(trophicImpute), position = Mode(waterColumnImpute),
#              squidpop = round(mean(eat.squidImpute),0), BodyShapeI=Mode(BodyShapeI),
#              BodyShapeII= Mode(BodyShapeII), Forehead=Mode(Forehead), TypeofMouth=Mode(TypeofMouth),
#              PosofMouth=Mode(PosofMouth), CShape=Mode(CShape), Attributes=Mode(Attributes),
#              TL=mean(TL,na.rm=T), HL=mean(HL,na.rm=T),BD=mean(BD,na.rm=T), SnoutTipX=mean(SnoutTipX,na.rm=T),
#              SnoutTipY=mean(SnoutTipY,na.rm=T), AspectRatio=mean(AspectRatio,na.rm=T), EncIndex=mean(EncIndex,na.rm=T),
#              DietTrop=mean(DietTroph,na.rm=T),DietSeTroph=mean(DietSeTroph,na.rm=T), 
#              FoodTroph=mean(FoodTroph,na.rm=T),FoodSeTroph=mean(FoodSeTroph,na.rm=T),
#              # SoftBottom=Mode(SoftBottom), Sand=Mode(Sand), Mud=Mode(Mud), HardBottom=Mode(HardBottom),
#              # SeaGrassBeds=Mode(SeaGrassBeds),CoralReefs=Mode(CoralReefs),Estuaries=Mode(Estuaries),
#              # Mangroves=Mode(Mangroves),Intertidal=Mode(Intertidal),Saltmarshes=Mode(Saltmarshes),
#              Herbivory2=Mode(Herbivory2),Climate=Mode(Climate), FoodI=Mode(FoodI),FoodII=Mode(FoodII) )
# trait.mode[,26:35] <- apply(trait.mode[,25:34],2,as.numeric)



# # add size as a trait
# size <- pred_biom %>%
#   dplyr::group_by( sciName ) %>%
#   dplyr::summarize( length = mean(Length, na.rm=T) )
# 
# trait.size <- left_join( trait.mode.clean, size )

# summarize seine data to get abundances by site
# only retain site and taxa that are in trait - more in trait because based on videos + seine
# trait.keep <- trait.size[ trait.size$sciName %in% s1$sciName, ]
# trait.keep <- trait.keep[ trait.keep$squidpop == 1, ]
trait.use <- data.frame( trait.use[ trait.use$sciName %in% suse$sciName, ] )
xuse <- trait.use[,-c(1,2)]
xuse[xuse==""] <- NA
xuse$act[ xuse$act == "active;ambush" ] <- "ambush"
# to numeric
# xuse[,c(2)] <- apply(xuse[,c(2)],2,as.numeric)
xuse$act <- as.numeric( as.factor(xuse$act))-1
summary(xuse)
row.names(xuse) <- trait.use$name
xuse <- xuse[ order(rownames(xuse)), ]
# xuse$squidpop <- as.numeric(as.character(xuse$squidpop))
mice::md.pattern( xuse )
# seine.keep <- s1[ s1$sciName %in% trait.keep$sciName, ]
  


### Ignore abundance patterns for Presence-Absence analysis,
# but can fill a matrix with ones and zeros
# abundance patterns
scomm <- suse %>% mutate( abun=1 )
comm <- scomm %>%
  ungroup() %>%
  unite( "SH", Site, habitat  ) %>%
  select( SH, name, abun ) %>%
  spread( name,abun,fill=0  )
ause <- as.matrix( comm[,-1] )
rownames(ause) <- comm$SH


# sort(unique(seine.keep$sciName))
# sort(unique(trait.keep$sciName))
# a <- seine.keep %>%
#   dplyr::mutate( area=Distance*width.transect ) %>%
#   dplyr::mutate( cpua = Abundance/area ) %>%
#   dplyr::group_by( Site, sciName ) %>%
#   dplyr::summarise( n = sum(cpua,na.rm=T) ) %>%
#   # group_by( Site) %>%
#   spread( sciName, n, fill=0 )
# a <- data.frame(a)
# rownames(a) <- a$Site
# a <- a[,-1]
# ause <- a
# 
# b <- seine.keep %>%
#   dplyr::mutate( area=Distance*width.transect ) %>%
#   dplyr::mutate( bpua = biomass/area ) %>%
#   dplyr::group_by( Site, sciName ) %>%
#   dplyr::summarise( n = sum(bpua,na.rm=T) ) %>%
#   # group_by( Site) %>%
#   spread( sciName, n, fill=0 )
# b <- data.frame(b)
# rownames(b) <- b$Site
# b <- b[,-1]
# buse <- b

# use dbFD
# m1 <- dbFD( xuse, ause, w.abun = T, corr="cailliez", calc.FGR=T, clust.type='ward.D' )
# cwm <- functcomp( as.matrix(xuse), as.matrix(ause) )
# xmat <- list( trait=xuse, abun=ause )
# debugonce(FD::dbFD)
# x <- as.list( body(dbFD) )
# as.list( x[[54]] )
# as.list( x[[c(54,3,2,3,2)]] )[[3]]
# trace(dbFD, browser, 
#       at = list( c(54,3,2,3,2) )
#       )
# 
# dbFD( xmat$trait, xmat$abun, w.abun = T, 
#       corr="lingoes", print.pco=T, calc.FGR = T, clust.type="kmeans",
#       km.inf.gr=3,km.sup.gr=20 )
# untrace(dbFD)



gd <- gowdis( xuse, asym.bin = 1 )
m1 <- dbFD( gd, ause, w.abun = T, 
            corr="lingoes", print.pco=T, calc.FGR = T, clust.type="kmeans",
            km.inf.gr=3,km.sup.gr=20,
            calc.CWM = FALSE, m="min", calc.FDiv = FALSE, stand.FRic = FALSE )

# Community weight means
cwm <- functcomp( xuse, ause )
cwm$SH <- comm$SH
cwm <- cwm %>% 
  separate( SH, c("Site","habitat"))
# combine with environmental data
ss <- suse %>% 
  select( Site, habitat, meanLat, meanLong, sstmean, rate ) %>% 
  distinct()
cwm.ss <- left_join( cwm, ss  )

write_csv( cwm, "Output Data/FunctionalDiversity_CWM_PA_filtered.csv")

# quantitative traits
ggplot( cwm.ss, aes( x=abs(meanLat), y=as.numeric(as.character(act)))) + geom_point() + geom_smooth(method='glm', method.args=list(family='quasibinomial'))
# ggplot( cwm.ss, aes( x=abs(meanLat), y=food)) + geom_point() + geom_smooth()
# ggplot( cwm.ss, aes( x=abs(meanLat), y=asp)) + geom_point() + geom_smooth()
# ggplot( cwm.ss, aes( x=abs(meanLat), y=tl)) + geom_point() + geom_smooth()
# ggplot( cwm.ss, aes( x=abs(meanLat), y=snoutx)) + geom_point() + geom_smooth()
# ggplot( cwm.ss, aes( x=abs(meanLat), y=snouty)) + geom_point() + geom_smooth()
# qualitative traits
ggplot( cwm.ss, aes( x=abs(meanLat), y=feed)) + geom_point() 
ggplot( cwm.ss, aes( x=abs(meanLat), y=troph)) + geom_point() 
ggplot( cwm.ss, aes( x=abs(meanLat), y=watercol)) + geom_point() 



# 
# fdm <- do.call( cbind, m1 )
# # functional groups 
# fdm$Site <- rownames(cwm)
# fdm.site <- left_join( fdm, sites )
# group abundance
# gr.abun <- m1$gr.abun
# gr.abun$SH <- comm$SH 
# gr.abun <- gr.abun %>% 
#   separate( SH, c("Site", "habitat") )
# gr.abun <- left_join( gr.abun, ss )
# # can we define groups based on trait bundles?
# # taxa in each group
# groups <- data.frame( name=names(m1$spfgr), functgroup=m1$spfgr )
# seine.group <- left_join( s1, groups )
# ggplot( seine.group, aes(x=functgroup,y=rate) ) + geom_point()

# summarize groups by site (distribution based on presence-absence and abundance)
# functional group richness
# without abundance
# group.site <- seine.group %>%
#   dplyr::group_by( Site, habitat, sciName, functgroup ) %>%
#   dplyr::summarize( n=sum(abun,na.rm=T))
# group.long <- group.site %>%
#   dplyr::ungroup() %>%
#   dplyr::select( Site, habitat, functgroup ) %>%
#   dplyr::group_by( Site, habitat, functgroup ) %>%
#   dplyr::summarize( count=length(functgroup) ) #%>%
# # spread( functgroup, count, fill=0 )
# # names(group.wide)[-1] <- paste0("group",names(group.wide)[-1])
# group.rate <- right_join( ss, group.long )
# group.rate <- group.rate %>% arrange(rate)
# # size of point is number in each group, location of points based on group number and site (rate?)
# site.ord <- unique(group.rate$Site)
# reps <- table(group.rate$Site)[site.ord]
# group.rate$index <- rep( 1:length(site.ord), as.vector(reps) )
# # functional group richness
# ggr <- ggplot( group.rate, aes( y=factor(functgroup), x=index, size=count, col=rate)) + 
#   geom_point() + scale_color_continuous(guide=FALSE)
# 
# # based on abundance
# group.n <- group.site %>%
#   dplyr::ungroup() %>%
#   dplyr::select( Site, habitat, functgroup, abun ) %>%
#   dplyr::group_by( Site, habitat, functgroup ) %>%
#   dplyr::summarize( count=length(functgroup), cpua=sum(cpua) )
# group.rate.n <- right_join( sites, group.n )
# group.rate.n <- group.rate.n %>% arrange(rate)
# 
# # size of point is number in each group, location of points based on group number and site (rate?)
# group.rate.n$index <- rep( 1:length(site.ord), as.vector(reps) )
# # richness within functional groups 
# ggn <- ggplot( group.rate.n, aes( y=factor(functgroup), x=index, size=cpua, col=rate)) + 
#   geom_point()
# 
# windows(5,5)
# plot_grid( ggr,ggn,ncol=1, align = 'hv')
#   # group4, group6, group8
# 
# par(mfrow=c(3,3))
# for(i in 1:9) {
#   plot( formula(paste0('rate~group',i)),gr.abun )
# }
# dev.off()
# plot(m1$x.axes[,1:2])

# 
# 
# # remaining questions:
# # do all families associated with squidpop consumption fall into the same functional group - NO
# # what is it about the taxa strongly associated with consumption?
# # which families are represented by each functional group?
# spec.group <- data.frame( sciName=names(m1$spfgr), group=m1$spfgr )
# fam.group <- left_join( spec.group, select(trait, sciName, family ))
# fam.group %>%
#   dplyr::group_by( group ) %>%
#   dplyr::summarize( n=length(unique(family)), fams = paste(sort(unique(family)),collapse=";") )
# # CWM for functional groups
# b2 <- seine.group %>%
#   dplyr:: mutate( area=Distance*width.transect ) %>%
#   dplyr::mutate( bpua = biomass/area ) %>%
#   dplyr::group_by( functgroup, sciName ) %>%
#   dplyr::summarise( n = sum(bpua,na.rm=T) )
# b2s <- b2 %>%
#   # group_by( Site) %>%
#   spread( sciName, n, fill=0 )
# b2s <- data.frame(b2s)
# rownames(b2s) <- b2s$functgroup
# b2s <- b2s[,-1]
# buse <- b2s
# trait.keep2 <- trait.keep[ trait.keep$sciName %in% b2$sciName, ]
# # trait.keep <- trait.keep[ trait.keep$squidpop == 1, ]
# xuse2 <- trait.keep2[,-c(1,6)]
# row.names(xuse2) <- trait.keep2$sciName
# cwm2 <- functcomp( as.matrix(xuse2), as.matrix(buse) )
# cwm2$functgroup <- (rownames(cwm2))
# plot( activty ~ functgroup, cwm2 )
# choose <- 12:18
# par( mfcol=c( 3,3  ))
# for( i in choose ){
#   plot(  cwm2$functgroup, as.numeric(cwm2[,i] ),
#          xlab="Functional group", ylab=names(cwm2)[i] )
# }
# 
# plot( as.numeric(HL) ~ (functgroup), cwm2 )
# plot( as.numeric(BD) ~ (functgroup), cwm2 )


## Other dimensions of functional diversity
names(m1)
# number of functional groups
fd <- with(m1, data.frame(FRic,qual.FRic,FEve,FDis, RaoQ,FGR) )
fd$SH <- rownames(fd)
fd <- fd %>% 
  separate( SH, c("Site", "habitat"))

fds <- left_join(fd,ss)
psych::pairs.panels( fds[,c("FRic","FEve","FDis","RaoQ","FGR","rate")])
#RaoQ, FGR, Fric
# from the help page: Rao's quadratic entropy (Q) is computed from the uncorrected species-species distance matrix via divc. See Botta-Dukát (2005) for details. Rao's Q is conceptually similar to FDis, and simulations (via simul.dbFD) have shown high positive correlations between the two indices (Laliberté and Legendre 2010). Still, one potential advantage of FDis over Rao's Q is that in the unweighted case (i.e. with presence-absence data), it opens possibilities for formal statistical tests for differences in FD between two or more communities through a distance-based test for homogeneity of multivariate dispersions (Anderson 2006); see betadisper for more details.
psych::pairs.panels( fds[,c("FRic","RaoQ","FGR","rate")])
write_csv( fd, "Output Data/FunctionalDiversity_indices_PA_filtered.csv" )

with(fds,plot(FRic~cwm$act))
