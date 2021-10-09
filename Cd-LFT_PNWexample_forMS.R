############## Code for analyses in Anderegg et al. ##################
# "Representing plant functional diversity in land models: 
# An evolutionary approach to making ‘Functional Types’ functional"
######################################################################3

# last updated 10.8.21 by ldla
# previous files: Cd-LFT_testVardecompd.R; CutCode_forLFTanalysis.R

##### Load Packages
require(tidyr)
require(lme4)
require(lmerTest)
require(stringr)
require(RColorBrewer)
require(dplyr)
require(MuMIn)
require(ggplot2)
source("/Users/leeanderegg/Desktop/R functions, general/ggplot_helpers.R")
##### Variance Decomp of LES dataset for LFT analysis ######

# getting baad data as of 10.25.19
# devtools::install_github("ropenscilabs/datastorr")
# install.packages("remotes")
# remotes::install_github("traitecoevo/baad.data")



mypal <- brewer.pal(n=9, "Set1")
palette(mypal)
# set the colors that wil denote within-species, within-genus, within family and across CWMs
colchoices <- c(1,2,4,3,6)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### load BAAD data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

baad.data::baad_data_version_current()
  # currently 1.0.1 on 10/02/17
  # still 1.0.1 on 098/30/21
baad.all <- baad.data::baad_data()
baad.all$dictionary[which(baad.all$dictionary$variable=="growingCondition"),c("variable", "label", "units")]

baad <- baad.all$data
baad$genus <- sapply(X=baad$speciesMatched, FUN = function(x){strsplit(x, split = " ")[[1]][1]})
baad$studyName <- factor(baad$studyName)
baad$genus <- factor(baad$genus)
baad$family <- factor(baad$family)
baad$Al.As <- baad$a.lf/baad$a.ssbh
baad$Ab.Bl <- baad$m.so/baad$m.rt
baad$RSR <- baad$m.rt/baad$m.so # making root:shoot ratio for matching with Ledo 2017 below
baad$lmf <- baad$m.lf / (baad$m.lf + baad$m.st)
baad$taper <- baad$d.ba/baad$h.t
baad$d.bh[which(is.na(baad$d.bh))] <- 4/5 * baad$d.ba[which(is.na(baad$d.bh))]



#colnames(baad)[which(colnames(baad) %in% c("family"))] <- c("Family", "a.lf")

# a.lf = leaf area
# a.ssbh = sapwood area at breast height
# h.t = height
# m.lf = leaf mass
# m.st = stem mass
# m.so = aboveground mass
# m.rf/c = mass of fine/coarse roots
# m.rt = root mass
# r.st = wood density
# r.ss = sapwood density
# r.sh = heartwood desnity


baad$log.Al.As <- log(baad$Al.As, base=10)
baad$log.r.st <- log(baad$r.st, base=10)
baad$log.Ab.Bl <- log(baad$Ab.Bl, base=10)
baad$log.RSR <- log(baad$RSR, base=10)
baad$log.h.t <- log(baad$h.t, base=10)
baad$lm.dbh <- with(baad, m.lf/(3.1415 *(d.bh/2)^2))
baad$log.lm.dbh <- log(baad$lm.dbh, base=10)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
####### Global Wood Density Database ######
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GWD <- read.csv("/Users/leeanderegg/Dropbox/WD project/GlobalWoodDensityDatabase.csv", header=T)
GWD$Genus <- factor(sapply(X=as.character(GWD$Binomial), FUN = function(x){strsplit(x, split = " ")[[1]][1]}))
# bunch of replicated genera, and 100+ families.

colnames(GWD)[4] <- "WD"

# #WDvar2 <- lmer(WD~ 1 + (1|Genus/Binomial), GWD)
# WDvar2 <- lmer(WD~ 1 + (1|Family) + (1|Genus) + (1|Binomial), GWD)
# logWDvar2 <- lmer(log(WD, base=10)~ 1 + (1|Family) + (1|Genus) + (1|Binomial), GWD)
# 
# rawWDvariance2 <- data.frame(VarCorr(WDvar2))
# WDvariance2 <- data.frame(VarCorr(logWDvar2))




##### ####### Al:As & WD variance decomp #############
colnames(GWD)
colnames(baad)
# combine WD data in GWD and baad
wdbaad <- baad %>% filter(r.st>0) %>% select(WD=r.st, Binomial=speciesMatched, Family=family, Genus=genus, Study=studyName) %>% mutate(WD=WD/1000)
wdGWD <- GWD %>% select(WD, Binomial, Family, Genus,Study=Reference.Number)
wdbaad$Study <- as.character(wdbaad$Study)
wdGWD$Study <- as.character(wdGWD$Study)
allwd <- rbind(wdbaad, wdGWD)

 WDvar3 <- lmer(WD~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)
 logWDvar3 <- lmer(log(WD, base=10)~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)
 
 rawWDvariance3 <- data.frame(VarCorr(WDvar3))
 rawWDvariance3$scaledVar <- rawWDvariance3$vcov/ sum(rawWDvariance3$vcov)
 WDvariance3 <- data.frame(VarCorr(logWDvar3))


## Calculate whole plant Al.As for obs without sapwood area but with DBH or basal area
## make a flag for how Al.As was calcualted
baad$Al.As_type <- NA
baad$Al.As_type[which(!is.na(baad$Al.As))] <- "direct.dbh"
baad$Al.As.others <- NA

baad$Al.As.others[which(is.na(baad$a.ssbh) & !is.na(baad$a.ssba))] <- with(baad[which(is.na(baad$a.ssbh) & !is.na(baad$a.ssba)),], a.lf/a.ssba)
baad$Al.As_type[which(is.na(baad$a.ssbh) & !is.na(baad$a.ssba))] <- "direct.basal"

baad$Al.As_type[which(is.na(baad$a.ssbh) & !is.na(baad$a.stbh) & !is.na(baad$a.lf) & is.na(baad$Al.As.others))] <- "nosw.dbh"
baad$Al.As.others[which(is.na(baad$a.ssbh) & !is.na(baad$a.stbh) & !is.na(baad$a.lf) & is.na(baad$Al.As.others))] <- with(baad[which(is.na(baad$a.ssbh) & !is.na(baad$a.stbh) & !is.na(baad$a.lf) & is.na(baad$Al.As.others)),], a.lf/(a.stbh))

baad$Al.As_type[which(is.na(baad$Al.As) &
                        !is.na(baad$a.stba) &
                        !is.na(baad$a.lf) &
                        is.na(baad$Al.As.others))] <- "nosw.basal"
baad$Al.As.others[which(is.na(baad$Al.As) 
                        & !is.na(baad$a.stba) 
                        & !is.na(baad$a.lf) 
                        & baad$Al.As_type=="nosw.basal")] <- with(baad[which(is.na(baad$Al.As) & !is.na(baad$a.stba) & !is.na(baad$a.lf) & baad$Al.As_type=="nosw.basal"),], a.lf/(a.stba))

baad$Al.As.all <- baad$Al.As
baad$Al.As.all[which(is.na(baad$Al.As))] <- baad$Al.As.others[which(is.na(baad$Al.As))]


# this new Al.As.all has 3798 measurements, with 234 species, 68 studies, 136 genera, 62 families
# updated: new database has 10409 measurements (9585 with heights)
# AlAsvar <- lmer(log(Al.As.all)~ Al.As_type + (1|family/genus/speciesMatched), baad)
# AlAsvariance <- data.frame(VarCorr(AlAsvar))
# AlAsvariance$scaledVar <- AlAsvariance$vcov/ sum(AlAsvariance$vcov)


   
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## R:S database (2017) ###########
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


RtS <- read.csv("/Users/leeanderegg/Desktop/LFT_analysis/nph14863-sup-0002-NotesS1.csv", header=T)
RtS$vegetationTYPE[which(RtS$vegetationTYPE=="TropMF ")] <- "TropMF"
RtS$vegetationTYPE[which(RtS$vegetationTYPE=="PlantTem ")] <- "PlantTem"
RtS$vegetationTYPE <- factor(RtS$vegetationTYPE)
RtS$log.AGB <- log(RtS$AGB_kg, base=10)
RtS$BGB_kg[which(RtS$BGB_kg==0)] <- NA # there are 94 measurements with 0 for BGB
RtS$log.BGB <- log(RtS$BGB_kg, base=10)
RtS$RS_RATIO[which(RtS$RS_RATIO==0)] <- NA # there are 94 measurements with 0 for BGB
RtS$log.RSR <- log(RtS$RS_RATIO, base=10)
RtS$taper <- RtS$DBH_cm/RtS$H_m



######## Xylem Functional Traits Database XFT #################
  # note: exported from Xylem functional traits database Master 13 March 2015.xls but had to change a bunch of column names
xft <- read.csv("/Users/leeanderegg/Desktop/LFT_analysis/XFT_traits_13_March_2015.csv", header=T,na.strings = "")
  # turns out there are 300+ uncleaned species
xft$Cleaned.family <- as.character(xft$Cleaned.family)
xft$Cleaned.family[which(is.na(xft$Cleaned.binomial))] <- xft$Family[which(is.na(xft$Cleaned.binomial))]
xft$Cleaned.genus <- as.character(xft$Cleaned.genus)
xft$Cleaned.genus[which(is.na(xft$Cleaned.binomial))] <- xft$Genus[which(is.na(xft$Cleaned.binomial))]
xft$Species <- as.character(xft$Species)
xft$Species[which(xft$Species==" limon")] <- "limon"
xft$Species.cleaning <- xft$Species
xft$Species.cleaning[grep(" x ", xft$Species)] <- "hybrid"
xft$Species.cleaning <- str_replace(xft$Species.cleaning, " \\(shade\\)", "")
xft$Species.cleaning[grep("�", xft$Species.cleaning)] <- "tremuloides"
xft$Species.cleaning <- str_replace(xft$Species.cleaning, " ","")
xft$Cleaned.binomial <- as.character(xft$Cleaned.binomial)
xft$Cleaned.binomial[which(is.na(xft$Cleaned.binomial))] <- paste(xft$Genus[which(is.na(xft$Cleaned.binomial))], xft$Species.cleaning[which(is.na(xft$Cleaned.binomial))])
# replacing 's' with 'S'
xft$Plant.organ[which(xft$Plant.organ=="s")] <- "S"


xft$Cleaned.family <- factor(xft$Cleaned.family)
xft$Cleaned.genus <- factor(xft$Cleaned.genus)
xft$Cleaned.binomial <- factor(xft$Cleaned.binomial)
xft$P50 <- str_replace(xft$P50, "_","-")
xft$P50 <- as.numeric(xft$P50)

## turns out there are 300+ rows that don't have 'Cleaned' family, genus, spp




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### Combine all the variance decomps together
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# the plan: includ - LMA, Nmass, LL, WD, RSR, P50, Ks

##### . LMA, Nmass & LL ###########
data.all <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample/DerivedData/AllTraitData_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)

## with all spp used for the hierarchical analysis
logLMAvar <- lmer(log.LMA~ 1 + (1|Project) + (1|Family) + (1|Genus) + (1|Species), data.all)
logLLvar <- lmer(log.LL~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all[-which(data.all$Project=="CO"),])
logNmassvar <- lmer(log.Nmass~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all)
logNareavar <- lmer(log.Narea~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all)
# technically, Genus and Species random effects should be nested within Family, but the nested model never converged.
# tests on subsets yeild extremely similar variance estimates for nested and non-nested effects.

# create dataframes with variance parameter estimates
LMAvariance <- data.frame(VarCorr(logLMAvar))
LLvariance <- data.frame(VarCorr(logLLvar))
Nmassvariance <- data.frame(VarCorr(logNmassvar))
Nareavariance <- data.frame(VarCorr(logNareavar))



####### . WD  ######
wdbaad <- baad %>% filter(r.st>0) %>% select(WD=r.st, Binomial=speciesMatched, Family=family, Genus=genus, Study=studyName) %>% mutate(WD=WD/1000)
wdGWD <- GWD %>% select(WD, Binomial, Family, Genus,Study=Reference.Number)
wdbaad$Study <- as.character(wdbaad$Study)
wdGWD$Study <- as.character(wdGWD$Study)
allwd <- rbind(wdbaad, wdGWD)

WDvar3 <- lmer(WD~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)
logWDvar3 <- lmer(log(WD, base=10)~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)

rawWDvariance3 <- data.frame(VarCorr(WDvar3))
#rawWDvariance3$scaledVar <- rawWDvariance3$vcov/ sum(rawWDvariance3$vcov)
#WDvariance3 <- data.frame(VarCorr(logWDvar3))



######## . P50, Ks from XFT ###################

logp50var <- lmer(log(-1*P50)~ Plant.organ + (1|Cleaned.family) + (1|Cleaned.genus) + (1|Cleaned.binomial), xft)
logKsvar <- lmer(log(Ks)~ Plant.organ + (1|Cleaned.family) + (1|Cleaned.genus) + (1|Cleaned.binomial), xft)

# create dataframes with variance parameter estimates
P50variance <- data.frame(VarCorr(logp50var))
Ksvariance <- data.frame(VarCorr(logKsvar))



### . R:S from baad (supplemented with Ledo 2017) ######

RSbaad <- baad[which(baad$RSR>0),] %>% select(latitude, longitude, vegetationTYPE=vegetation, Binomial=speciesMatched, Family=family, Genus=genus, RSR, log.RSR)
baadspp <- str_replace(unique(RSbaad$speciesMatched),pattern = " ",replacement = "_")
RtSunique <- RtS[which(!RtS$SPECIES %in% baadspp),] %>% select(latitude, longitude, vegetationTYPE, Binomial=SPECIES, Family=FAMILY, Genus=GENUS, RSR=RS_RATIO, log.RSR)
allRSR <- rbind(RSbaad, RtSunique)
allRSR$Binomial <- str_replace(allRSR$Binomial, "_"," ")

logRSRvar <- lmer(log.RSR ~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allRSR)

logRSRvariance <- data.frame(VarCorr(logRSRvar))



###### . Al:As from Trugman et al 2019 ######


# This code was originally developed for doing taxonomic variance decomposition and also
# estimating the effect of plant size on Al.As
# I used full code here but we only want the pure taxonomic variance decomposition 
# (where size-related variation is lumped int residual or within-species variation)


# 1) Restrictive analysis with only measurements with direct sapwood area values (at DBH or basal)
AlAsvar1 <- lmer(log(Al.As.all)~ Al.As_type + (1|family/genus/speciesMatched), baad[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0),])
AlAsvar2 <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + (1|family/genus/speciesMatched), baad[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0),])
AlAsvar3 <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + log(h.t):family - family + (1|family/genus/speciesMatched), baad[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0),])
#AlAsvar4 <- lmer(log(Al.As.all)~ log(h.t) * Al.As_type + (log(h.t)|family/genus/speciesMatched), baad[which(baad$Al.As_type %in% c("direct.dbh","direct.basal")),])
AIC(AlAsvar1, AlAsvar2, AlAsvar3)

AlAsvar.pure <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + (1|family/genus/speciesMatched), baad[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0),])
AlAsvariance.pure <- data.frame(VarCorr(AlAsvar.pure))
AlAsvariance.pure$scaledVar <- AlAsvariance.pure$vcov/ sum(AlAsvariance.pure$vcov)

# 2) Analysis with both direct and indirect (using basal area or DBH)
AlAsvar1 <- lmer(log(Al.As.all)~ Al.As_type + (1|family/genus/speciesMatched), baad[which(baad$h.t>0),])
AlAsvar2 <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + (1|family/genus/speciesMatched), baad[which(baad$h.t>0),])
AlAsvar3 <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + log(h.t):family - family + (1|family/genus/speciesMatched), baad[which(baad$h.t>0),])
AlAsvar4 <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + log(h.t):speciesMatched - speciesMatched + (1|family/genus/speciesMatched), baad[which(baad$h.t>0),])
#AlAsvar4 <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + (log(h.t)|family/genus/speciesMatched), baad[which(baad$Al.As_type %in% c("direct.dbh","direct.basal")),])
AIC(AlAsvar1, AlAsvar2, AlAsvar3, AlAsvar4)
#seperate lines per species are actually the best (could do this with random slopes, but then I'm not sure how to interpret the variance decomp)

# 08/19/19 version with slopes per species
AlAsvar.all <- lmer(log(Al.As.all)~ log(h.t) + Al.As_type + log(h.t):speciesMatched - speciesMatched + (1|family/genus/speciesMatched), baad[which(baad$h.t>0),])
AlAsvariance.all <- data.frame(VarCorr(AlAsvar.all))
AlAsvariance.all$scaledVar <- AlAsvariance.all$vcov/ sum(AlAsvariance.all$vcov)


effects <- AlAsvar.all@beta[-c(1:5)] # pull out the betas for species, minus the intercept, ref and effects of type
ref <- AlAsvar.all@beta[2]
effects.real <- effects + ref
AlAsvar.all.effects <- c(ref, effects.real)


## Just the taxonomic variance decomp
#### THIS IS WHAT IS USED FOR ANDEREGG ET AL. 2021
AlAsvar.all.noh <- lmer(log(Al.As.all)~ Al.As_type  + (1|family/genus/speciesMatched), baad[which(baad$h.t>0),])
AlAsvariance.all.noh <- data.frame(VarCorr(AlAsvar.all.noh))
AlAsvariance.all.noh$scaledVar <- AlAsvariance.all.noh$vcov/ sum(AlAsvariance.all.noh$vcov)

# determining how much variation is due to Height, and whether it primarily comes from w/in species component
r.squaredGLMM(AlAsvar.all)
r.squaredGLMM(AlAsvar.all.noh)
# so 36.298 % is explained in .all
# 7.65 % is explained in .all.noh (so 7% is Al.As_type)
# so 36.3 - 7.65 = 28.65% explained by height 
# (maybe some of that 7 is joint height and Al.As_type, but I think we should attribute it to type for methodological concerns)



Al.As_pure_n <- c(length(unique(baad$family[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0)])),length(unique(baad$genus[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0)])),length(unique(baad$speciesMatched[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0)])),length(baad$family[which(baad$Al.As_type %in% c("direct.dbh","direct.basal") & baad$h.t>0)]) )
Al.As_all_n <- c(length(unique(baad$family[which(!is.na(baad$Al.As.all) & baad$h.t>0)])),length(unique(baad$genus[which(!is.na(baad$Al.As.all)& baad$h.t>0)])),length(unique(baad$speciesMatched[which(!is.na(baad$Al.As.all)& baad$h.t>0)])),length(baad$family[which(!is.na(baad$Al.As.all)& baad$h.t>0)]) )

#new.vardecomps.forann <- data.frame(rawWDvariance3[c(3,2,1,4),c("grp","scaledVar")], AlAsvariance$scaledVar[c(3,2,1,4)])
# note: as of 4/29/19 i just copied the new Al_As column to the old spreadsheet
#new.vardecomps.forann <- data.frame(rawWDvariance3[c(3,2,1,4),c("grp","scaledVar")], AlAsvariance$scaledVar[c(3,2,1,4)], AlAsvariance.pure$scaledVar[c(3,2,1,4)], Al.As_pure_n, AlAsvariance.all$scaledVar[c(3,2,1,4)], Al.As_all_n, AlAsvariance.all.noh$scaledVar[c(3,2,1,4)])
# note: as of 5/31/19 I added Al_As.pure and Al_As.all
new.vardecomps.forann <- data.frame(rawWDvariance3[c(3,2,1,4),c("grp","scaledVar")], AlAsvariance$scaledVar[c(3,2,1,4)], AlAsvariance.pure$scaledVar[c(3,2,1,4)], Al.As_pure_n, Al.As_all_n, AlAsvariance.all$scaledVar[c(3,2,1,4)], AlAsvariance.all.noh$scaledVar[c(3,2,1,4)], AlAsvariance.all$vcov[c(3,2,1,4)], AlAsvariance.all.noh$vcov[c(3,2,1,4)])
# note: as of 8/19/19 I added unscaled variances and updated the wheight to have a slope per species

colnames(new.vardecomps.forann) <- c("TaxoLeve","WD","Al.As_pure", "Al.As_pure_wheight","pure_wheight_n", "all_wheight_n","Al.As_all_wheight", "Al.As_all_noheight", "Al.As_all_wheight_raw", "Al.As_all_noheight_raw" )
#write.csv(new.vardecomps.forann, file = "/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Al_As-vardecom-forAnna-20190603.csv" )
#new.vardecoms.forann <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Al_As-vardecom-forAnna-20190603.csv")


# combine all variance estimates, leaving out "Project"
# updated 06.17.21 to bring in Al:As
traitvars <- data.frame(LMAvariance[which(LMAvariance$grp !="Project"),4], LLvariance[which(LLvariance$grp !="Project"),4], Nmassvariance[which(Nmassvariance$grp !="Project"),4], Nareavariance[which(Nareavariance$grp !="Project"),4],
                        rawWDvariance3[,4], P50variance[,4],Ksvariance[,4], logRSRvariance[,4], new.vardecoms.forann[new.vardecoms.forann$X,"Al.As_all_noheight"])
colnames(traitvars) <- c("logLMA", "logLL", "logNmass", "logNarea", "WD","logP50","logKs","logRSR","logAlAs")
rownames(traitvars) <- c("BtwSpecies", "BtwGenera", "BtwFamilies", "WtinSpecies")


traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}
# re-order rows
traitvars_scaled2 <- traitvars_scaled[c(3,2,1,4),]






#_______________________________________________________________________________
############## FIGURE 1: Variance Decomp #####
#quartz(width=3.75, height=4)
quartz(width=6.81, height=3.5)
par(mgp=c(3,.7,0), cex.lab=1.3, cex.axis=1.1, mfrow=c(1,1), mar=c(6,2,3,6), oma=c(0,3,0,0))
cols <- rev(c(mypal[colchoices[c(1,2,3)]],"#CCCCCC"))#brewer.pal(11, "RdBu")[c(11,9,1, 6)]
# With old colors, alpha=CC, with new colors, alpha=99
bp <-barplot(as.matrix(traitvars_scaled2[,c("logNmass","logLL","logLMA","WD","logP50","logKs","logNarea","logAlAs","logRSR")]),beside=F,legend.text = F,xpd = T,las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
             , names.arg = c(expression(paste(log(N[mass]))),"log(LeafLife)","log(LMA)", "WD",expression(paste(log(P[50]))),expression(paste(log(K[s]))), expression(paste(log(N[area]))), expression(paste(log(A[L]:A[S]))),"log(R:S)")#log(Height)","AboveGrnd:\nBelowGrnd\nbiomass")
             , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
             , xlab="", mgp=c(2.5,.8,0))
#legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
#text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.4, font=1)
#mtext("a)", side=3, adj=-.1, line=1.3)
mtext("Proportion of Total Variance", side=2, line=2.8)
#mtext("BAAD Allometry data", side=3, line=0.3)
#quartz(width=5, height=4)
#cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
#legend(xpd=T, x = 11, y=0.7, legend=rev(c("btw Fams","btw Gen","btw Spp")), fill=rev(paste0(cols[1:3],"99")), ncol=1, bty="n",  cex=1.2)
# mtext(text="28 spp w/ replication, 1000+ spp,\n500+ genera, 150+ families",side = 1,line = 3.3)
# legend(xpd=T, x = 0, y=1.3, legend=rownames(rtraitvars_scaled2), fill=paste0(cols,"CC"), ncol=2, bty="n",  cex=1.2)
legend(xpd=NA, x = 11, y=0.7, legend=rev(c("btw Fams","btw Gen","btw Spp","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=1, bty="n",  cex=1)

#_______________________________________________________________________________



# #### Plotting Subset for Presentation ########
# ## making everything the same size as the Ecol Let Vardecmp FIG 1
# quartz(width=4.33, height=4.73)
# par(mfrow=c(2,2), mgp=c(2,.7,0), cex.lab=1.1, cex.axis=1.1, mar=c(4.5,2,1.5,2),oma=c(0,2,3.8,0))
# 
# bp <-barplot(as.matrix(baadtraitvars_scaled3[,c("BAADrawWD","LMF","logHeight","logAGB_BGB"),]),beside=F,legend.text = F,xpd = T,las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
#              , names.arg = c("","","","")#c("Wood Density","Leaf Mass Fraction","log(Height)","AboveGrnd:BelowGrnd\nbiomass") #c("log(WD)","raw WD", "Leaf Mass\nFraction", "log(Height)","AboveGrnd:\nBelowGrnd\nbiomass")
#              , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
#              , xlab="", mgp=c(2.5,.8,0))
# #legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
# #text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.4, font=1)
# #mtext("a)", side=3, adj=-.1, line=1.3)
# mtext("Proportion of Total Variance", side=2, line=2.8)
# text(x = bp, y= par("usr")[3]-.05,labels =  c("Wood Density", "Leaf Mass Fraction",expression(paste(log[10](Height))),expression(paste(M[AbvGr]:M[BlwGr]))), srt=40, adj=1,xpd=NA, cex=1.1, font=1)
# 
# 
# 











#_______________________________________________________________________________
############# Leaf Lifespan Exploration ##############################
#_______________________________________________________________________________

# bring in full GLOPNET+ database from Anderegg et al. 2018 Eco Let

data.all <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample/DerivedData/AllTraitData_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)
# fix some misnamed oak families
data.all$Family[which(data.all$Family=="Facaeae")] <- "Fagaceae"
# calculate spp averages
allspp <- data.all %>% group_by(Species, Genus, Family) %>% summarise( lslog.LL = mean(log.LL, na.rm=T), lslog.LMA = mean(log.LMA, na.rm=T), lslog.Nmass = mean(log.Nmass, na.rm=T), lslog.Narea = mean(log.Narea, na.rm=T)
                                                                       ,slog.LL = log(mean(10^log.LL, na.rm=T),base=10), slog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), slog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), slog.Narea = log(mean(10^log.Narea, na.rm=T),base=10)
                                                                       ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
colnames(allspp) <- gsub("slog", "log", colnames(allspp))



#### Genus means ___________________________________________________________________
allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], lglog.LL = mean(llog.LL, na.rm=T),lglog.LMA = mean(llog.LMA, na.rm=T),lglog.Nmass = mean(llog.Nmass, na.rm=T), lglog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,glog.LL = log(mean(10^log.LL, na.rm=T),base=10), glog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), glog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), glog.Narea = log(mean(10^log.Narea, na.rm=T),base=10), nspp = n() 
                                                    ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
# taking mean of raw values or logged values doesn't matter all that much yet
colnames(allgen) <- gsub("glog", "log", colnames(allgen))
colnames(allgen)[2] <- "Family"
# 939 genera from 211 families
# allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], Needle.Broad = unique(na.omit(NeedleBroad))[1], glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
#                                                     ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() )



#### Family Means _______________________________________________________________________
allfam <- allgen %>% group_by(Family) %>% summarise(lflog.LL = mean(llog.LL, na.rm=T), lflog.LMA = mean(llog.LMA, na.rm=T), lflog.Nmass = mean(llog.Nmass, na.rm=T),lflog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,flog.LL = log(mean(10^log.LL, na.rm=T),base=10), flog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), flog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), flog.Narea = log(mean(10^log.Narea, na.rm=T), base=10)
                                                    , tnspp = sum(nspp), ngen=n() )
colnames(allfam) <- gsub("flog", "log", colnames(allfam))

#### Family Means for fams w/ > 3 species
fam.data <- allfam # new: 211 families (up from 189)
fam.dataclean <- allfam[which(allfam$tnspp>2),] # 101 families, up from 97 families




##### dataset of only genera w/ >5 species _____________________________________________
# for w/in genus analysis
gen.data <- allspp[which(allspp$Genus %in% names(which(xtabs(~Genus, allspp)>=5))),]
# 750 measurements of 73 genera.
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels


##### Spp. in Family level data, for families w/ >5 species in them ________________________
sppinfam.data <- allspp[which(allspp$Family %in% names(which(xtabs(~Family,allspp)>=5))),]
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels


# making it so that logged LL of 12 months = 1 (i.e. log base 12)
data.all$log12.LL <- with(data.all, log(10^log.LL, base=12))
#LES$log12.LL <- with(LES, log(10^log.LL, base=12))
allspp$log12.LL <- with(allspp, log(10^log.LL, base=12))
allgen$log12.LL <- with(allgen, log(10^log.LL, base=12))
allfam$log12.LL <- with(allfam, log(10^log.LL, base=12))

common.gens <- names(which(xtabs(~Genus, allspp[which(allspp$log.LL>0),])>=5))

# figure out range within genera
genLLrange <- data.all %>% group_by (Genus) %>% summarise(minLL = min(log12.LL, na.rm=T), maxLL= max(log12.LL, na.rm=T), minT = min(MAT, na.rm=T), n=n()-length(which(is.na(log12.LL))))

# ID genera that have deciduous and evergreen spp
flipfloppers <- genLLrange[which(genLLrange$minLL<1 & genLLrange$maxLL>1 & genLLrange$n>3),]


quartz(width=4.4, height=3.8)
p2 <- ggplot(allspp, aes(x=MAT,y=log12.LL)) + geom_point(col="grey") +
    #geom_point(data=allspp[which(allspp$Genus %in% c("Quercus", "Salix","Betula","Populus")),], aes(col=Genus)) + geom_smooth(data=allspp[which(allspp$Genus %in% c("Quercus", "Salix","Betula","Populus")),], method="lm",se=F, aes(col=Genus))+
    geom_point(data=allspp[which(allspp$Genus %in% common.gens),], aes(col=Genus)) + geom_smooth(data=allspp[which(allspp$Genus %in% common.gens),], method="lm",se=F, aes(col=Genus))+
    geom_hline(yintercept=1)+
        ggtitle("Within Genera") +
  labs(x="Mean Annual Temp (ºC)", y=expression(paste(log[12],"(Leaf Lifespan in months)"))) + 
  
    theme(legend.text = element_text(size=8), legend.key.height = unit(.4, 'cm'))

p1 <- ggplot(allspp, aes(x=MAT,y=log12.LL)) + geom_point(col="grey") +
    geom_point(data=allspp[which(allspp$Family %in% c("Fabaceae", "Fagaceae","Asteraceae", "Ericaceae","Rosaceae","Myrtaceae", "Pinaceae")),], aes(col=Family)) + 
    geom_smooth(data=allspp[which(allspp$Family %in% c("Fabaceae","Fagaceae", "Asteraceae", "Ericaceae","Rosaceae","Myrtaceae","Pinaceae")),], method="lm", se=F,aes(col=Family)) +
    geom_hline(yintercept=1) +
    labs(x="Mean Annual Temp (ºC)", y=expression(paste(log[12],"(Leaf Lifespan in months)"))) + 
    ggtitle("Within Families")
 # "Acacia","Piper","Eucalyptus",
# (p3 <-ggplot(allspp[which(allspp$MAT>10),], aes(x=VPD, y=log12.LL)) ) + geom_point()+
#   geom_point(data=allspp[which(allspp$Genus %in% common.gens),], aes(col=Genus)) + geom_smooth(data=allspp[which(allspp$Genus %in% common.gens),], method="lm",se=F, aes(col=Genus))+
#   geom_hline(yintercept=1)+
#   ggtitle("Hydraulic Region")
# )
p3 <- ggplot(allspp, aes(x=MAT,y=log12.LL)) + geom_point(col="grey") +
    #geom_point(data=allspp[which(allspp$Genus %in% c("Quercus")),], aes(col=Genus), size=3) +
    geom_point(data=data.all[which(data.all$Species %in% names(which(xtabs(~Species, data.all[which(data.all$log12.LL>1),])>30))),], aes(col=Species), alpha=0.1 ) +
    geom_smooth(data=data.all[which(data.all$Species %in% names(which(xtabs(~Species, data.all[which(data.all$log12.LL>1),])>30))),], aes(col=Species), method="lm", se=F) +
  theme(legend.position="none") + 
  ggtitle("Theory") +
  geom_hline(yintercept=1) +
  labs(x="Mean Annual Temp (ºC)", y=expression(paste(log[12],"(Leaf Lifespan in months)"))) 

  #geom_point(data=traits[which(traits$Family=="Pinaceae"),], aes(col=SP.ID), alpha=0.06) + 
    #geom_smooth(data=traits[which(traits$Family=="Pinaceae"),], aes(col=SP.ID), method="lm", se=F) + theme(legend.position = "none")


# verticle
quartz(width=4, height=6)
multiplot(p1,p2, cols=1)

# single
quartz(width=3.5, height=3.75)
p3

#horizontal
quartz(width=7.5, height=3)
multiplot(p1,p2, cols=2)

















#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###  ***** Phylogeny and LFT Generation   ###########
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

require(stringr)
require(taxize)
require(brranching)
require(ape)
require(dplyr)
require(tidyverse)

### step one: Make a Species List ###############3

biomass <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample/DerivedData/PNW_Biomass_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names=1)
traits <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample/DerivedData/PNW_Trait_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names=1)


spp.codes <- traits %>% group_by(GE.SP, Family) %>% summarise(Code=unique(SP.ID), fullspp = unique(FullSpecies), n_traits=n())
spplist_raw <- unique(traits$GE.SP)
spplist <- str_replace(spplist_raw[-grep("unknown", spplist_raw)], "\\.", " ")

spplist1 <- spplist[order(spplist)]

sppordered <- spplist1[c(10,33,34,7,8,19,30,31,32,13,14,9,12,15,16,29,18,28,20,21,22,27,26,24,25,23,37,36,4,2,1,6,3,5,11,35,17)]
LeafLife <- c("E","E","D","D","D","E","E","D","D","E","D","D","D","D","D","D","D",rep("E", times=20))

Spp.descriptions <- data.frame(sppordered, LeafLife)
Spp.descriptions$SP.ID <- spp.codes$Code[match(sppordered,spp.codes$fullspp)]
Spp.descriptions$n_traits <- spp.codes$n_traits[match(sppordered,spp.codes$fullspp)]
Spp.descriptions$Family <- spp.codes$Family[match(sppordered,spp.codes$fullspp)]

## clean the names with taxize:
#taxize::use_entrez()
#usethis::edit_r_environ()
spplist_phylo <- phylomatic_names(spplist)
spp_phylo <- phylomatic(spplist_phylo)


pinus1 <- read.tree("/Users/leeanderegg/Desktop/LFT_analysis/subtree-node-mrcaott8444ott8454-Pinus--Cathaya.tre")

### summarize species statistics to plot on tips ####

biomass$bio1 <- with(biomass, SPP_O1_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)
biomass$bio2 <- with(biomass, SPP_O2_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)
biomass$bio3 <- with(biomass, SPP_O3_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)
biomass$bio4 <- with(biomass, SPP_O4_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)

sp1 <- data.frame(Species=biomass$SPP_O1_ABBREV, biomass=biomass$bio1, perc_BA = biomass$SPP_O1_BASAL_AREA_FRACTION)
sp2 <- data.frame(Species=biomass$SPP_O2_ABBREV, biomass=biomass$bio2, perc_BA = biomass$SPP_O2_BASAL_AREA_FRACTION)
sp3 <- data.frame(Species=biomass$SPP_O3_ABBREV, biomass=biomass$bio3, perc_BA = biomass$SPP_O3_BASAL_AREA_FRACTION)
sp4 <- data.frame(Species=biomass$SPP_O4_ABBREV, biomass=biomass$bio4, perc_BA = biomass$SPP_O4_BASAL_AREA_FRACTION)

biomasslong <- rbind(sp1[-which(is.na(sp1$biomass)|is.na(sp1$Species)),], sp2[-which(is.na(sp2$biomass) | sp2$Species=="none"),],sp3[-which(is.na(sp3$biomass) | sp2$Species=="none"),],sp4[-which(is.na(sp4$biomass) | sp2$Species=="none"),])
tot_bio <- biomasslong %>% group_by(Species) %>% summarise(n_plots=n(), total_biomass=sum(biomass, na.rm=T), tot_dominance = sum(perc_BA, na.rm=T), max_dominance = max(perc_BA, na.rm=T), mean_dominance=mean(perc_BA, na.rm=T), median_dominance=median(perc_BA, na.rm=T))
tot_bio_sorted <- tot_bio[match(Spp.descriptions$SP.ID,tot_bio$Species),]
# note: this drops CERLED and PRUEMA, both of which don't have any trait measurements. but they're only tiny (<5%) fraction of 1&2 plots, respectively

nacp.sums <- left_join(Spp.descriptions, tot_bio, by=c("SP.ID"="Species"))


xtabs(~GE.SP, traits)
write.csv(nacp.sums, "SpeciesSummaries_forLFTanalysis_V1.csv")

#### ** note: I manually added in some columns about dominance and also the first few LFTs
nacp.sums <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/SpeciesSummaries_forLFTanalysis_V2.csv")




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############## . assigning traits their FT #####################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


traits$PFT <- nacp.sums$PFT[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Deep.LFT <- nacp.sums$Deep[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Mid.LFT <- nacp.sums$Mid[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Shallow.LFT <- nacp.sums$Shallow[match(traits$SP.ID, nacp.sums$SP.ID)]

# with Larix broken out
traits$PFT.larix <- nacp.sums$PFT.larix[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Deep.LFT.larix <- nacp.sums$Deep.larix[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Mid.LFT.larix <- nacp.sums$Mid.larix[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Shallow.LFT.larix <- nacp.sums$Shallow.larix[match(traits$SP.ID, nacp.sums$SP.ID)]






####### . cluster analysis, unsupervised #############

# Ward agglomerative clustering
clust1 <- hclust(d=dist(traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")]),method = "ward.D", members = NULL)
clust1.group9 <-  cutree(clust1, k = 9)
  # deprecated using cutree rather than old rect.hclust that produced a list
#clust1.groups <- rep(NA, nrow(traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),]))
# for(i in 1:length(clust1.groups)){
#   for(j in 1:9){
#     if(i %in% clust1.group9[[j]]) clust1.groups[i] <- j
#   }
#   
# }
traits$clust1 <- NA
traits$clust1[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust1.group9


#K-means clustering
clust.kmeans <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 9,nstart = 100)
traits$kmeans <- NA
traits$kmeans[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans$cluster


## summarizing how many species fall into which categories


####### . Analysis of kmeans consistency with different n ########
clust.kmeans.5 <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 5,nstart = 100)
traits$kmeans.5 <- NA
traits$kmeans.5[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans.5$cluster

clust.kmeans.6 <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 6,nstart = 100)
traits$kmeans.6 <- NA
traits$kmeans.6[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans.6$cluster

clust.kmeans.10 <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 10,nstart = 100)
traits$kmeans.10 <- NA
traits$kmeans.10[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans.10$cluster



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######### ** Univariate explained variance ##############
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Nmass.PFTmod <- lm(log.Nmass~PFT, traits)
Nmass.Deep.LFTmod <- lm(log.Nmass~Deep.LFT, traits)
Nmass.Mid.LFTmod <- lm(log.Nmass~Mid.LFT, traits)
Nmass.Shallow.LFTmod <- lm(log.Nmass~Shallow.LFT.larix, traits)

LMA.PFTmod <- lm(log.LMA~PFT, traits)
LMA.Deep.LFTmod <- lm(log.LMA~Deep.LFT, traits)
LMA.Mid.LFTmod <- lm(log.LMA~Mid.LFT, traits)
LMA.Shallow.LFTmod <- lm(log.LMA~Shallow.LFT.larix, traits)
LMA.clust1mod <- lm(log.LMA~clust1, traits)


LL.PFTmod <- lm(log.LL~PFT, traits)
LL.Deep.LFTmod <- lm(log.LL~Deep.LFT, traits)
LL.Mid.LFTmod <- lm(log.LL~Mid.LFT, traits)
LL.Shallow.LFTmod <- lm(log.LL~Shallow.LFT.larix, traits)
LL.clust1mod <- lm(log.LL~clust1, traits)

log.Shallow <- data.frame(trait = c("LMA","N","LL"),
                      r2 = c(summary(LMA.Shallow.LFTmod)$r.squared,
                             summary(Nmass.Shallow.LFTmod)$r.squared,
                             summary(LL.Shallow.LFTmod)$r.squared)
)


### average variance explained with clusters
ward <- data.frame(trait = c("LMA","N","LL"),
                   r2 = c(summary(lm(LMA~factor(clust1), traits))$r.squared,
                          summary(lm(LEAF_NITROGEN~factor(clust1), traits))$r.squared,
                          summary(lm(LEAF_LIFE~factor(clust1), traits))$r.squared)
)

kmeans.r2 <- data.frame(trait = c("LMA","N","LL"),
                        r2 = c(summary(lm(LMA~factor(kmeans), traits))$r.squared,
                               summary(lm(LEAF_NITROGEN~factor(kmeans), traits))$r.squared,
                               summary(lm(LEAF_LIFE~factor(kmeans), traits))$r.squared)
)

## Note: as of 09.12.2021, everything uses log-traits
log.ward <- data.frame(trait = c("LMA","N","LL"),
                       r2 = c(summary(lm(log.LMA~factor(clust1), traits))$r.squared,
                              summary(lm(log.Nmass~factor(clust1), traits))$r.squared,
                              summary(lm(log.LL~factor(clust1), traits))$r.squared)
)

log.kmeans.r2 <- data.frame(trait = c("LMA","N","LL"),
                            r2 = c(summary(lm(log.LMA~factor(kmeans), traits))$r.squared,
                                   summary(lm(log.Nmass~factor(kmeans), traits))$r.squared,
                                   summary(lm(log.LL~factor(kmeans), traits))$r.squared)
)







#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############## FIG: plotting boxplot of biomass to show hyperdominance ##########
quartz(width=2.5, height=5)
par(mar=c(3.5,1,1,1), mgp=c(2.5,1,0))
barplot(height = rev(nacp.sums$total_biomass/1000), names.arg = NA, horiz = T, las=1, xlab="total biomass (kg C)")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########### Variance of CWMs explained by LFTs ############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# for .v1 and .v2 used PFT, Deep, Mid and Shallow.larix

biomass$PFT1 <- nacp.sums$PFT[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT1 <- nacp.sums$Deep[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT1 <- nacp.sums$Mid[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT1 <- nacp.sums$Shallow.larix[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]

biomass$PFT2 <- nacp.sums$PFT[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT2 <- nacp.sums$Deep[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT2 <- nacp.sums$Mid[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT2 <- nacp.sums$Shallow.larix[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]


biomass$PFT3 <- nacp.sums$PFT[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT3 <- nacp.sums$Deep[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT3 <- nacp.sums$Mid[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT3 <- nacp.sums$Shallow.larix[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]


biomass$PFT4 <- nacp.sums$PFT[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT4 <- nacp.sums$Deep[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT4 <- nacp.sums$Mid[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT4 <- nacp.sums$Shallow.larix[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]


PFT.traitmeans <- traits %>% group_by(PFT) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))
Deep.LFT.traitmeans <- traits %>% group_by(Deep.LFT) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))
Mid.LFT.traitmeans <- traits %>% group_by(Mid.LFT) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))
Shallow.LFT.traitmeans <- traits %>% group_by(Shallow.LFT.larix) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))




#######. LMA CWM by disaggregation method #################3
# calcualte community weighted LMA by PFT
PFT.LMA1 <- PFT.traitmeans$LMA[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.LMA2 <- PFT.traitmeans$LMA[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.LMA3 <- PFT.traitmeans$LMA[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.LMA4 <- PFT.traitmeans$LMA[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.LMA <- apply(data.frame(PFT.LMA1,PFT.LMA2,PFT.LMA3,PFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted LMA by Deep LFT
Deep.LFT.LMA1 <- Deep.LFT.traitmeans$LMA[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.LMA2 <- Deep.LFT.traitmeans$LMA[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.LMA3 <- Deep.LFT.traitmeans$LMA[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.LMA4 <- Deep.LFT.traitmeans$LMA[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.LMA <- apply(data.frame(Deep.LFT.LMA1,Deep.LFT.LMA2,Deep.LFT.LMA3,Deep.LFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted LMA by Mid LFT
Mid.LFT.LMA1 <- Mid.LFT.traitmeans$LMA[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.LMA2 <- Mid.LFT.traitmeans$LMA[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.LMA3 <- Mid.LFT.traitmeans$LMA[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.LMA4 <- Mid.LFT.traitmeans$LMA[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.LMA <- apply(data.frame(Mid.LFT.LMA1,Mid.LFT.LMA2,Mid.LFT.LMA3,Mid.LFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted LMA by Shallow LFT
Shallow.LFT.LMA1 <- Shallow.LFT.traitmeans$LMA[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.LMA2 <- Shallow.LFT.traitmeans$LMA[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.LMA3 <- Shallow.LFT.traitmeans$LMA[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.LMA4 <- Shallow.LFT.traitmeans$LMA[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.LMA <- apply(data.frame(Shallow.LFT.LMA1,Shallow.LFT.LMA2,Shallow.LFT.LMA3,Shallow.LFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)



#######. LL CWM by disaggregation method #################3
# calcualte community weighted LL by PFT
PFT.LL1 <- PFT.traitmeans$LL[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.LL2 <- PFT.traitmeans$LL[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.LL3 <- PFT.traitmeans$LL[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.LL4 <- PFT.traitmeans$LL[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.LL <- apply(data.frame(PFT.LL1,PFT.LL2,PFT.LL3,PFT.LL4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted LL by Deep LFT
Deep.LFT.LL1 <- Deep.LFT.traitmeans$LL[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.LL2 <- Deep.LFT.traitmeans$LL[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.LL3 <- Deep.LFT.traitmeans$LL[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.LL4 <- Deep.LFT.traitmeans$LL[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.LL <- apply(data.frame(Deep.LFT.LL1,Deep.LFT.LL2,Deep.LFT.LL3,Deep.LFT.LL4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted LL by Mid LFT
Mid.LFT.LL1 <- Mid.LFT.traitmeans$LL[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.LL2 <- Mid.LFT.traitmeans$LL[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.LL3 <- Mid.LFT.traitmeans$LL[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.LL4 <- Mid.LFT.traitmeans$LL[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.LL <- apply(data.frame(Mid.LFT.LL1,Mid.LFT.LL2,Mid.LFT.LL3,Mid.LFT.LL4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted LL by Shallow LFT
Shallow.LFT.LL1 <- Shallow.LFT.traitmeans$LL[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.LL2 <- Shallow.LFT.traitmeans$LL[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.LL3 <- Shallow.LFT.traitmeans$LL[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.LL4 <- Shallow.LFT.traitmeans$LL[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.LL <- apply(data.frame(Shallow.LFT.LL1,Shallow.LFT.LL2,Shallow.LFT.LL3,Shallow.LFT.LL4),MARGIN=1,FUN=sum, na.rm=T)







#######. Nmass CWM by disaggregation method #################3
# calcualte community weighted Nmass by PFT
PFT.Nmass1 <- PFT.traitmeans$Nmass[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.Nmass2 <- PFT.traitmeans$Nmass[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.Nmass3 <- PFT.traitmeans$Nmass[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.Nmass4 <- PFT.traitmeans$Nmass[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.Nmass <- apply(data.frame(PFT.Nmass1,PFT.Nmass2,PFT.Nmass3,PFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted Nmass by Deep LFT
Deep.LFT.Nmass1 <- Deep.LFT.traitmeans$Nmass[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.Nmass2 <- Deep.LFT.traitmeans$Nmass[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.Nmass3 <- Deep.LFT.traitmeans$Nmass[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.Nmass4 <- Deep.LFT.traitmeans$Nmass[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.Nmass <- apply(data.frame(Deep.LFT.Nmass1,Deep.LFT.Nmass2,Deep.LFT.Nmass3,Deep.LFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted Nmass by Mid LFT
Mid.LFT.Nmass1 <- Mid.LFT.traitmeans$Nmass[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.Nmass2 <- Mid.LFT.traitmeans$Nmass[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.Nmass3 <- Mid.LFT.traitmeans$Nmass[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.Nmass4 <- Mid.LFT.traitmeans$Nmass[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.Nmass <- apply(data.frame(Mid.LFT.Nmass1,Mid.LFT.Nmass2,Mid.LFT.Nmass3,Mid.LFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)



# calcualte community weighted Nmass by Shallow LFT
Shallow.LFT.Nmass1 <- Shallow.LFT.traitmeans$Nmass[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.Nmass2 <- Shallow.LFT.traitmeans$Nmass[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.Nmass3 <- Shallow.LFT.traitmeans$Nmass[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.Nmass4 <- Shallow.LFT.traitmeans$Nmass[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.Nmass <- apply(data.frame(Shallow.LFT.Nmass1,Shallow.LFT.Nmass2,Shallow.LFT.Nmass3,Shallow.LFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)






#######. Narea CWM by disaggregation method #################3
# calcualte community weighted Narea by PFT
PFT.Narea1 <- PFT.traitmeans$Narea[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.Narea2 <- PFT.traitmeans$Narea[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.Narea3 <- PFT.traitmeans$Narea[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.Narea4 <- PFT.traitmeans$Narea[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.Narea <- apply(data.frame(PFT.Narea1,PFT.Narea2,PFT.Narea3,PFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted Narea by Deep LFT
Deep.LFT.Narea1 <- Deep.LFT.traitmeans$Narea[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.Narea2 <- Deep.LFT.traitmeans$Narea[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.Narea3 <- Deep.LFT.traitmeans$Narea[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.Narea4 <- Deep.LFT.traitmeans$Narea[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.Narea <- apply(data.frame(Deep.LFT.Narea1,Deep.LFT.Narea2,Deep.LFT.Narea3,Deep.LFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted Narea by Mid LFT
Mid.LFT.Narea1 <- Mid.LFT.traitmeans$Narea[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.Narea2 <- Mid.LFT.traitmeans$Narea[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.Narea3 <- Mid.LFT.traitmeans$Narea[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.Narea4 <- Mid.LFT.traitmeans$Narea[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.Narea <- apply(data.frame(Mid.LFT.Narea1,Mid.LFT.Narea2,Mid.LFT.Narea3,Mid.LFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)



# calcualte community weighted Narea by Shallow LFT
Shallow.LFT.Narea1 <- Shallow.LFT.traitmeans$Narea[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.Narea2 <- Shallow.LFT.traitmeans$Narea[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.Narea3 <- Shallow.LFT.traitmeans$Narea[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.Narea4 <- Shallow.LFT.traitmeans$Narea[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.Narea <- apply(data.frame(Shallow.LFT.Narea1,Shallow.LFT.Narea2,Shallow.LFT.Narea3,Shallow.LFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)


biomass$log.PFT.LMA <- log(biomass$PFT.LMA, base=10)
biomass$log.Deep.LFT.LMA <- log(biomass$Deep.LFT.LMA, base=10)
biomass$log.Mid.LFT.LMA <- log(biomass$Mid.LFT.LMA, base=10)
biomass$log.Shallow.LFT.LMA <- log(biomass$Shallow.LFT.LMA, base=10)

biomass$log.PFT.LL <- log(biomass$PFT.LL, base=10)
biomass$log.Deep.LFT.LL <- log(biomass$Deep.LFT.LL, base=10)
biomass$log.Mid.LFT.LL <- log(biomass$Mid.LFT.LL, base=10)
biomass$log.Shallow.LFT.LL <- log(biomass$Shallow.LFT.LL, base=10)

biomass$log.PFT.Nmass <- log(biomass$PFT.Nmass, base=10)
biomass$log.Deep.LFT.Nmass <- log(biomass$Deep.LFT.Nmass, base=10)
biomass$log.Mid.LFT.Nmass <- log(biomass$Mid.LFT.Nmass, base=10)
biomass$log.Shallow.LFT.Nmass <- log(biomass$Shallow.LFT.Nmass, base=10)




################ . CW variance explained by each predictor ###########

log.CWexplans <- data.frame(Type = c("PFT","Deep LFT", "Mid LFT", "Shallow LFT"), 
                            log.LMA = c(summary(lm(log.cw_LMAp_if~log.PFT.LMA, biomass))$r.squared,
                                        summary(lm(log.cw_LMAp_if~log.Deep.LFT.LMA, biomass))$r.squared,
                                        summary(lm(log.cw_LMAp_if~log.Mid.LFT.LMA, biomass))$r.squared,
                                        summary(lm(log.cw_LMAp_if~log.Shallow.LFT.LMA, biomass))$r.squared),
                            log.LL = c(summary(lm(log.cw_LLp_if~log.PFT.LL, biomass))$r.squared,
                                       summary(lm(log.cw_LLp_if~log.Deep.LFT.LL, biomass))$r.squared,
                                       summary(lm(log.cw_LLp_if~log.Mid.LFT.LL, biomass))$r.squared,
                                       summary(lm(log.cw_LLp_if~log.Shallow.LFT.LL, biomass))$r.squared),
                            log.Nmass = c(summary(lm(log.cw_Nmassp_if~log.PFT.Nmass, biomass))$r.squared,
                                          summary(lm(log.cw_Nmassp_if~log.Deep.LFT.Nmass, biomass))$r.squared,
                                          summary(lm(log.cw_Nmassp_if~log.Mid.LFT.Nmass, biomass))$r.squared,
                                          summary(lm(log.cw_Nmassp_if~log.Shallow.LFT.Nmass, biomass))$r.squared))



CWexplans <- data.frame(Type = c("PFT","Deep LFT", "Mid LFT", "Shallow LFT"), 
                        LMA = c(summary(lm(cw_LMAp_if~PFT.LMA, biomass))$r.squared,
                                summary(lm(cw_LMAp_if~Deep.LFT.LMA, biomass))$r.squared,
                                summary(lm(cw_LMAp_if~Mid.LFT.LMA, biomass))$r.squared,
                                summary(lm(cw_LMAp_if~Shallow.LFT.LMA, biomass))$r.squared),
                        LL = c(summary(lm(cw_LLp_if~PFT.LL, biomass))$r.squared,
                               summary(lm(cw_LLp_if~Deep.LFT.LL, biomass))$r.squared,
                               summary(lm(cw_LLp_if~Mid.LFT.LL, biomass))$r.squared,
                               summary(lm(cw_LLp_if~Shallow.LFT.LL, biomass))$r.squared),
                        Nmass = c(summary(lm(cw_Nmassp_if~PFT.Nmass, biomass))$r.squared,
                                  summary(lm(cw_Nmassp_if~Deep.LFT.Nmass, biomass))$r.squared,
                                  summary(lm(cw_Nmassp_if~Mid.LFT.Nmass, biomass))$r.squared,
                                  summary(lm(cw_Nmassp_if~Shallow.LFT.Nmass, biomass))$r.squared))

rowMeans(CWexplans[,-1])







#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############# FD within-communities analysis ########################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#using FD package, I can calcualte a whole bunch of FD metrics real easy.
# I just need to set up the data appropriately.

require(FD)
length(unique(c(unique(biomass$SPP_O1_ABBREV), unique(biomass$SPP_O2_ABBREV), unique(biomass$SPP_O3_ABBREV), unique(biomass$SPP_O4_ABBREV))))
# 31 species have meaningful BA (well, probably 29-30 because some 'unknowns' and some 'none')
length(unique(traits$SP.ID))
# 39 species in trait datasets

rel.abun <- data.frame(matrix(data=0, nrow=nrow(biomass), ncol = length(unique(traits$SP.ID))), row.names = biomass$PLOT_ID)
colnames(rel.abun) <- unique(traits$SP.ID)

for (i in unique(traits$SP.ID)){
  rel.abun[which(biomass$SPP_O1_ABBREV==i),i] <-biomass$SPP_O1_BASAL_AREA_FRACTION[which(biomass$SPP_O1_ABBREV==i)]
  rel.abun[which(biomass$SPP_O2_ABBREV==i),i] <-biomass$SPP_O2_BASAL_AREA_FRACTION[which(biomass$SPP_O2_ABBREV==i)]
  rel.abun[which(biomass$SPP_O3_ABBREV==i),i] <-biomass$SPP_O3_BASAL_AREA_FRACTION[which(biomass$SPP_O3_ABBREV==i)]
  rel.abun[which(biomass$SPP_O4_ABBREV==i),i] <-biomass$SPP_O4_BASAL_AREA_FRACTION[which(biomass$SPP_O4_ABBREV==i)]
}

## Note: 5 plots have "0" relative abund (remove), and a couple have slightly less than 100%
tot.rel.abun <- rowSums(rel.abun, na.rm=T)
# remove those sites with <90 rel abun accounted for (8 sites)
rel.abun.cl <- rel.abun[-which(tot.rel.abun < 90),] 

## Note: a number of species have no meaningful BA
# -> two options: 1)  find the sites where they occur and give them 1%, or 2) remove them

#starting with simply removing them
rel.abun.clean <- rel.abun.cl[,which(colSums(rel.abun.cl)>0)]


### aggregate rel.abun columns to PFTs and LFTs
rel.abun.PFT <- data.frame(matrix(data=0, nrow=nrow(rel.abun.clean), ncol = length(unique(nacp.sums$PFT))), row.names = rownames(rel.abun.clean))
colnames(rel.abun.PFT) <- unique(nacp.sums$PFT)
for (i in unique(nacp.sums$PFT)){
  rel.abun.PFT[,i]  <- rowSums(rel.abun.clean[which(colnames(rel.abun.clean) %in% nacp.sums$SP.ID[which(nacp.sums$PFT==i)])])
}

rel.abun.Deep <- data.frame(matrix(data=0, nrow=nrow(rel.abun.clean), ncol = length(unique(nacp.sums$Deep))), row.names = rownames(rel.abun.clean))
colnames(rel.abun.Deep) <- unique(nacp.sums$Deep)
for (i in unique(nacp.sums$Deep)){
  rel.abun.Deep[,i]  <- rowSums(rel.abun.clean[which(colnames(rel.abun.clean) %in% nacp.sums$SP.ID[which(nacp.sums$Deep==i)])])
}

rel.abun.Mid <- data.frame(matrix(data=0, nrow=nrow(rel.abun.clean), ncol = length(unique(nacp.sums$Mid))), row.names = rownames(rel.abun.clean))
colnames(rel.abun.Mid) <- unique(nacp.sums$Mid)
for (i in unique(nacp.sums$Mid)){
  rel.abun.Mid[,i]  <- rowSums(rel.abun.clean[which(colnames(rel.abun.clean) %in% nacp.sums$SP.ID[which(nacp.sums$Mid==i)])])
}

rel.abun.Shallow.larix <- data.frame(matrix(data=0, nrow=nrow(rel.abun.clean), ncol = length(unique(nacp.sums$Shallow.larix))), row.names = rownames(rel.abun.clean))
colnames(rel.abun.Shallow.larix) <- unique(nacp.sums$Shallow.larix)
for (i in unique(nacp.sums$Shallow.larix)){
  rel.abun.Shallow.larix[,i]  <- rowSums(rel.abun.clean[which(colnames(rel.abun.clean) %in% nacp.sums$SP.ID[which(nacp.sums$Shallow.larix==i)])])
}



################# Species Means (across all sites) ##########
# Note: mLMA = LMA_HSA, and mLMA_PSA = LMA_PSA, I've used LMA_HSA in full dataset creation
# also, mlog.Trait = mean of logged traits
# log.Trait = log of mean traits

spp.traits <- traits %>% group_by(SP.ID) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                       , LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID))
                                                       , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNmass = mean(LEAF_NITROGEN, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                       , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                       , climPC1 = mean(climPC1, na.rm=T), climPC2 = mean(climPC2, na.rm=T), climPC3 = mean(climPC3, na.rm=T)
                                                       , soil_N= mean(soil_N, na.rm=T), soil_pH=mean(soil_pH, na.rm=T), ASA = mean(ASA, na.rm=T), LAI_O=mean(LAI_O, na.rm=T), AG_TGROWTH = mean(AG_TGROWTH, na.rm=T))
# log species mean traits for between species analysis
spp.traits$log.LMA <- log(spp.traits$mLMA, base=10)
spp.traits$log.LL <- log(spp.traits$mLLmonths, base=10)
spp.traits$log.Nmass <- log(spp.traits$mNmass, base=10)
spp.traits$log.Narea <- log(spp.traits$mNarea, base=10)


# trim down to only the traits I care about to match with CWM analysis
spp.tr <- data.frame(spp.traits %>% select(log.LMA, log.LL, log.Nmass))
row.names(spp.tr) <- spp.traits$SP.ID
## fill in unknown LL for missing species
# ARBMEN 14.79108 mo in GLOPNET
# CERLED 11.74 from C. betuloides in GLOPNET. set to 12
# JUNOCC 77.625 from J monosperma in Glopnet
# PINFLE 36.31 in GLOPNET
spp.tr$log.LL[which(row.names(spp.tr)=="ARBMEN")] <- log(14.79, base=10)
spp.tr$log.LL[which(row.names(spp.tr)=="CERLED")] <- log(12, base=10)
spp.tr$log.LL[which(row.names(spp.tr)=="JUNOCC")] <- log(77.625, base=10)
spp.tr$log.LL[which(row.names(spp.tr)=="PINFLE")] <- log(36.31, base=10)

spp.ordered <- spp.tr[match(colnames(rel.abun), rownames(spp.tr)),]

# remove spp with 0  rel abun in rel.abun.clean
spp.ordered.clean <- spp.ordered[which(row.names(spp.ordered) %in% colnames(rel.abun.clean)),]


# do the same thing for PFT/LFT trait means
PFT.ordered <- data.frame(PFT.traitmeans[match(colnames(rel.abun.PFT), PFT.traitmeans$PFT),] %>% select(log.LMA, log.LL, log.Nmass)) # get rid of NA row at bottom
rownames(PFT.ordered) <- names(rel.abun.PFT)

Deep.ordered <- data.frame(Deep.LFT.traitmeans[match(colnames(rel.abun.Deep), Deep.LFT.traitmeans$Deep.LFT),] %>% select(log.LMA, log.LL, log.Nmass)) # get rid of NA row at bottom
rownames(Deep.ordered) <- names(rel.abun.Deep)

Mid.ordered <- data.frame(Mid.LFT.traitmeans[match(colnames(rel.abun.Mid), Mid.LFT.traitmeans$Mid.LFT),] %>% select(log.LMA, log.LL, log.Nmass)) # get rid of NA row at bottom
rownames(Mid.ordered) <- names(rel.abun.Mid)

Shallow.ordered <- data.frame(Shallow.LFT.traitmeans[match(colnames(rel.abun.Shallow.larix), Shallow.LFT.traitmeans$Shallow.LFT.larix),] %>% select(log.LMA, log.LL, log.Nmass)) # get rid of NA row at bottom
rownames(Shallow.ordered) <- names(rel.abun.Shallow.larix)





######### calculate Species means for each plot in which they occur
spp.plot.traits <- traits %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                     , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNmass = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                     , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                     , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                                     , climPC1 = unique(climPC1), climPC2 = unique(climPC2), climPC3 = unique(climPC3)
)
# create unique species-plot tag
spp.plot.traits$SP.PLOT <- paste(spp.plot.traits$SP.ID, spp.plot.traits$PLOT_ID, sep="-")



####### Calculating FD with SPP means #############

FDspp <- dbFD(spp.ordered.clean, a=rel.abun.clean)
# FEVe: Could not be calculated for communities with <3 functionally singular species. 
# FDis: Equals 0 in communities with only one functionally singular species. 
# FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species. 
# FRic: Dimensionality reduction was required. The last PCoA axis (out of 3 in total) was removed. 
# FRic: Quality of the reduced-space representation = 0.9342504 
# FDiv: Could not be calculated for communities with <3 functionally singular species. 

## roughly 120 (50%) of plots have only 1 or 2 sppecies, so FEve, FRic, and FDiv are NA for those plots
# but FDis and RaoQ work

FDPFT <- dbFD(PFT.ordered, a=rel.abun.PFT)

FDDeep <- dbFD(Deep.ordered, a=rel.abun.Deep)
FDMid <- dbFD(Mid.ordered, a=rel.abun.Mid)
FDShallow <- dbFD(Shallow.ordered, a= rel.abun.Shallow.larix)


# summarize the explanatory power of the four Functional Types for Rao Quad Ent and Func Dispersioin
FDsummaries <- data.frame(FT = c("PFT","Deep","Mid","Shallow"), RaoQ = rep(NA, 4), FDis= rep(NA, 4))

FDsummaries$RaoQ[1] <- summary(lm(FDspp$RaoQ~FDPFT$FDis))$r.squared
FDsummaries$RaoQ[2] <- summary(lm(FDspp$RaoQ~FDDeep$FDis))$r.squared
FDsummaries$RaoQ[3] <- summary(lm(FDspp$RaoQ~FDMid$FDis))$r.squared
FDsummaries$RaoQ[4] <- summary(lm(FDspp$RaoQ~FDShallow$FDis))$r.squared

FDsummaries$FDis[1] <- summary(lm(FDspp$FDis~FDPFT$FDis))$r.squared
FDsummaries$FDis[2] <- summary(lm(FDspp$FDis~FDDeep$FDis))$r.squared
FDsummaries$FDis[3] <- summary(lm(FDspp$FDis~FDMid$FDis))$r.squared
FDsummaries$FDis[4] <- summary(lm(FDspp$FDis~FDShallow$FDis))$r.squared




########## . Table of FT ~ clim R2s ##############

climR2 <- data.frame(PFT = rep(NA, times=4), DeepLFT = rep(NA, times=4), MidLFT=rep(NA, times=4), ShallowLFT = rep(NA, times=4), Ward = rep(NA, times=4), Kmeans=rep(NA, times=4))

climR2[1,] <- round(c(summary(lm(tmean.gy.c~factor(PFT), traits))$r.squared,
            summary(lm(tmean.gy.c~factor(Deep.LFT), traits))$r.squared,
            summary(lm(tmean.gy.c~factor(Mid.LFT), traits))$r.squared,
            summary(lm(tmean.gy.c~factor(Shallow.LFT.larix), traits))$r.squared,
            summary(lm(tmean.gy.c~factor(clust1), traits))$r.squared,
            summary(lm(tmean.gy.c~factor(kmeans), traits))$r.squared),2)

climR2[2,] <- round(c(summary(lm(ppt.gy.mm~factor(PFT), traits))$r.squared,
                summary(lm(ppt.gy.mm~factor(Deep.LFT), traits))$r.squared,
                summary(lm(ppt.gy.mm~factor(Mid.LFT), traits))$r.squared,
                summary(lm(ppt.gy.mm~factor(Shallow.LFT.larix), traits))$r.squared,
                summary(lm(ppt.gy.mm~factor(clust1), traits))$r.squared,
                summary(lm(ppt.gy.mm~factor(kmeans), traits))$r.squared), 2)

climR2[3,] <- round(c(summary(lm(vpd.gy.max~factor(PFT), traits))$r.squared,
                summary(lm(vpd.gy.max~factor(Deep.LFT), traits))$r.squared,
                summary(lm(vpd.gy.max~factor(Mid.LFT), traits))$r.squared,
                summary(lm(vpd.gy.max~factor(Shallow.LFT.larix), traits))$r.squared,
                summary(lm(vpd.gy.max~factor(clust1), traits))$r.squared,
                summary(lm(vpd.gy.max~factor(kmeans), traits))$r.squared), 2)

climR2[4,] <- round(c(summary(lm(cmi.gy.mm~factor(PFT), traits))$r.squared,
                summary(lm(cmi.gy.mm~factor(Deep.LFT), traits))$r.squared,
                summary(lm(cmi.gy.mm~factor(Mid.LFT), traits))$r.squared,
                summary(lm(cmi.gy.mm~factor(Shallow.LFT.larix), traits))$r.squared,
                summary(lm(cmi.gy.mm~factor(clust1), traits))$r.squared,
                summary(lm(cmi.gy.mm~factor(kmeans), traits))$r.squared), 2)

rownames(climR2) <- c("MAT","MAP","maxVPD","CMI")
# climR2[5,] <- c(summary(lm(climPC1~factor(PFT), traits))$r.squared,
#                 summary(lm(climPC1~factor(Deep.LFT), traits))$r.squared,
#                 summary(lm(climPC1~factor(Mid.LFT), traits))$r.squared,
#                 summary(lm(climPC1~factor(Shallow.LFT.larix), traits))$r.squared,
#                 summary(lm(climPC1~factor(clust1), traits))$r.squared,
#                 summary(lm(climPC1~factor(kmeans), traits))$r.squared)
# 
# 
# climR2[6,] <- c(summary(lm(climPC2~factor(PFT), traits))$r.squared,
#                 summary(lm(climPC2~factor(Deep.LFT), traits))$r.squared,
#                 summary(lm(climPC2~factor(Mid.LFT), traits))$r.squared,
#                 summary(lm(climPC2~factor(Shallow.LFT.larix), traits))$r.squared,
#                 summary(lm(climPC2~factor(clust1), traits))$r.squared,
#                 summary(lm(climPC2~factor(kmeans), traits))$r.squared)
# 
# climR2[7,] <- c(summary(lm(climPC3~factor(PFT), traits))$r.squared,
#                 summary(lm(climPC3~factor(Deep.LFT), traits))$r.squared,
#                 summary(lm(climPC3~factor(Mid.LFT), traits))$r.squared,
#                 summary(lm(climPC3~factor(Shallow.LFT.larix), traits))$r.squared,
#                 summary(lm(climPC3~factor(clust1), traits))$r.squared,
#                 summary(lm(climPC3~factor(kmeans), traits))$r.squared)




#___________________________________________________________________________
###################### TABLE S1: Parameter Tables #########################

##### . PFTs #############
PFT.params <- traits[-which(is.na(traits$PFT)),] %>% group_by(PFT) %>% summarise (LMA = mean(LMA, na.rm=T), LeafLife=mean(LEAF_LIFE, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T))
PFT.params$P50 <- NA
PFT.params$P50[which(PFT.params$PFT=="BD")] <- mean(xft$P50[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form=="T" & xft$Group=="Angiosperm" & xft$Phenology %in% c("W","D") & xft$Plant.organ =="S")], na.rm=T) 
#temperate deciduous broadleaf = biom=TMR/TMS, growth=T, group=Angiosperm, Phenology = W/D
PFT.params$P50[which(PFT.params$PFT=="BE")] <- mean(xft$P50[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form=="T" & xft$Group=="Angiosperm" & xft$Phenology=="E" & xft$Plant.organ =="S")], na.rm=T) 
#temperate evergreen broadleaf = biom=TMR/TMS, growth=T, group=Angiosperm, Phenology = E
PFT.params$P50[which(PFT.params$PFT=="EN")] <- mean(xft$P50[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form=="T" & xft$Group=="Gymnosperm" & xft$Phenology!="W" & xft$Plant.organ =="S")], na.rm=T) 
#ever needle = biom=TMR/TMS, growth=T, group=Gymnosperm, Phenology = E
PFT.params$Ks <- NA
PFT.params$Ks[which(PFT.params$PFT=="BD")] <- mean(xft$Ks[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form=="T" & xft$Group=="Angiosperm" & xft$Phenology %in% c("W","D") & xft$Plant.organ =="S")], na.rm=T) 
#temperate deciduous broadleaf = biom=TMR/TMS, growth=T, group=Angiosperm, Phenology = W/D
PFT.params$Ks[which(PFT.params$PFT=="BE")] <- mean(xft$Ks[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form=="T" & xft$Group=="Angiosperm" & xft$Phenology=="E" & xft$Plant.organ =="S")], na.rm=T) 
#temperate evergreen broadleaf = biom=TMR/TMS, growth=T, group=Angiosperm, Phenology = E
PFT.params$Ks[which(PFT.params$PFT=="EN")] <- mean(xft$Ks[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form=="T" & xft$Group=="Gymnosperm" & xft$Phenology!="W" & xft$Plant.organ =="S")], na.rm=T) 
#ever needle = biom=TMR/TMS, growth=T, group=Gymnosperm, Phenology = E


########. Deep LFTs #############3
Deep.LFT.params <- traits[-which(is.na(traits$Deep.LFT)),] %>% group_by(Deep.LFT) %>% summarise (LMA = mean(LMA, na.rm=T), LeafLife=mean(LEAF_LIFE, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T))
Deep.LFT.params$P50 <- NA
Deep.LFT.params$P50[which(Deep.LFT.params$Deep.LFT=="Angio")] <- mean(xft$P50[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form %in% c("T","S") & xft$Group=="Angiosperm" & xft$Plant.organ =="S")], na.rm=T) 
#angio = biom=TMR/TMS, growth form=T/S group=Angiosperm
Deep.LFT.params$P50[which(Deep.LFT.params$Deep.LFT=="Cupressaceae")] <- mean(xft$P50[which(xft$Cleaned.family =="Cupressaceae")], na.rm=T)
# Cuppressaceae 
Deep.LFT.params$P50[which(Deep.LFT.params$Deep.LFT=="Pinaceae")] <- mean(xft$P50[which(xft$Cleaned.family=="Pinaceae")], na.rm=T)
# Pinaceae

Deep.LFT.params$Ks <- NA
Deep.LFT.params$Ks[which(Deep.LFT.params$Deep.LFT=="Angio")] <- mean(xft$Ks[which(xft$Biome %in% c("TMR","TMS") & xft$Growth.form %in% c("T","S") & xft$Group=="Angiosperm" & xft$Plant.organ =="S")], na.rm=T) 
#angio = biom=TMR/TMS, growth form=T/S group=Angiosperm
Deep.LFT.params$Ks[which(Deep.LFT.params$Deep.LFT=="Cupressaceae")] <- mean(xft$Ks[which(xft$Cleaned.family =="Cupressaceae")], na.rm=T)
# Cuppressaceae 
Deep.LFT.params$Ks[which(Deep.LFT.params$Deep.LFT=="Pinaceae")] <- mean(xft$Ks[which(xft$Cleaned.family=="Pinaceae")], na.rm=T)
# Pinaceae

###############. Shallow LFTs ################

Shallow.LFT.params <- traits[-which(is.na(traits$Shallow.LFT.larix)),] %>% group_by(Shallow.LFT.larix) %>% summarise (LMA = mean(LMA, na.rm=T), LeafLife=mean(LEAF_LIFE, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T))



