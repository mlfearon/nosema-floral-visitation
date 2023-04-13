# Code for the analyses for "Honeybee visitation behavior on shared floral resources 
# increases Vairimorpha ceranae prevalence in bumblebees"

# Manuscript submitted to: 


# Vairimorpha (=Nosema) ceranae prevalence in Apis and Bombus analysis

# This script includes the full analyses and figures for how pollinator visitation number and duration of visits to flowers
# impact V. ceranae prevalence in Apis mellifera and Bombus spp. 


# Written by: Michelle Fearon and Maryellen Zbrozek
# Last updated: 13 April 2023




# set the working directory
setwd("C:/Users/Maryellen/Documents/School/Honors Thesis/Manuscript/Analysis")


# Import libraries needed 

library(car)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(dplyr)
library(ggplot2)
library(ape)
library(epiR)
library(here)


# set the path to the script relative to the project root directory
here::i_am("scripts/NosemaAnalysis_2.20.2023.R")

###############################
### LOAD THE DATA FOR ANALYSIS
###############################
data <- read.csv(here("data/NosemaAnalysis_18Sremoved.csv"), stringsAsFactors = F)

head(data)
summary(data)

# convert variables into characters that were misread in by R
data$Year <- as.character(data$Year)
data$Visit <- as.character(data$Visit)



#####################################
# Calculate total Nosema prevalence in Apis and Bombus
#####################################
spp_prev <- matrix(NA, ncol = 5, nrow = 2)
rownames(spp_prev) <- c("Apis", "Bombus")
colnames(spp_prev) <- c("N_tested", "N_Nosema", "Nosema_prev", "Nosema_Lower", "Nosema_Upper")

GENUS <- unique(data$Genus)

for (i in GENUS){
  tempxx <- data[ data$Genus == i, ]
  
  ### Calculate Nosema prevalence for both species at each site
  # calculate the total number of infected ind for nosema (store in site_summary matrix) [BOTH SPECIES]
  infectedNosema <- sum(tempxx$Nosema)
  spp_prev[ i , "N_Nosema"] <- infectedNosema
  
  # calculate the total number tested
  tested <- length(tempxx$Genus)
  spp_prev[ i , "N_tested"] <- tested
  
  # unlist prevalence estimates to make easier to work with
  Nos <- unlist(epi.prev(pos = infectedNosema, tested = tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))
  
  
  # store total nosema prevalence and lower and upper confidence interval in site_summary matrix
  spp_prev[ i, "Nosema_prev"] <- Nos[1]
  spp_prev[ i, "Nosema_Lower"] <- Nos[2]
  spp_prev[ i, "Nosema_Upper"] <- Nos[3]
  
  
}
spp_prev



#####################################
####### Test of two proportions between total apis and bombus nosema prevalence ########
######################################

######## Nosema Prevalence

#Apis vs Bombus Nosema
infected <- c(spp_prev[ "Apis", "N_Nosema"], spp_prev[ "Bombus", "N_Nosema"])
total <- c(spp_prev[ "Apis", "N_tested"], spp_prev[ "Bombus", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)
# not significant p = 0.84


# subset data by species
data_Apis <- data[ data$Genus == "Apis", ]
data_Bombus <- data[ data$Genus == "Bombus", ]

unique(data_Bombus$Site)
unique(data_Apis$Site)
unique(data$Site)


#############################################
### Import functions needed for analyses
#############################################

# vif.mer function developed by Austin F. Frank https://raw.github.com/aufrank/R-hacks/master/mer-utils.R

vif.mer<-function(fit){   ## adapted from rms:vif
  v<-vcov(fit)
  nam<-names(fixef(fit)) ## exclude intercepts
  ns<-sum(1*(nam=="Intercept"|nam=="Intercept"))
  if(ns>0){
    v<-v[-(1:ns),-(1:ns),drop=FALSE]
    nam<-nam[-(1:ns)]
  }
  d<-diag(v)^(0.5)
  v<-diag(solve(v/(d%o%d)))
  names(v)<-nam
  v
}

# overdispersion test for glmer models (from Mark Hill 5/21/18)
# if overdispersed in a glm --> use quasipoisson and an F test to compare models
# if overdispersed in a glmer --> add ID as a random effect

# Overdispersion parameter estimation function from Ben Bolker: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
# Comment additions from http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
## updated overdispersion function from Ben Bolker
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
} # Generates a p-value. If less than 0.05, the data are overdispersed.



##############################################
### GLMMs for Nosema presence/absence
##############################################

# test of All nosema across sites and species
fit_cat <- glm(Nosema ~ Genus * Site, family = binomial, data = data)
summary(fit_cat)
Anova(fit_cat)



####Redone model for spatial autocorrelation

#moran's I
#Test for spatial autocorrelation
mod1_resid <- simulateResiduals(fit_cat)
mod1_resid2 <- recalculateResiduals(mod1_resid, group = data$Site) # group residuals by site
testSpatialAutocorrelation(mod1_resid2, unique(data$Long), unique(data$Lat))
# is significant! 

# Redo of model because of significant spatial autocorrelation

# first we need to create a numeric factor recording the coordinates of the sampled locations
data$pos <- numFactor(scale(data$Lat), scale(data$Long))
# then create a dummy group factor to be used as a random term
data$grp <- factor(rep(1, nrow(data)))


# fit the updated model including spatial random effects 
mod1_redo <- glmmTMB(Nosema ~ Genus * Site + mat(pos + 0|grp), family = binomial, data)
# model summary of fixed effects
summary(mod1_redo)
Anova(mod1_redo)


#recheck the spatial autocorrelation
mod1_resid_redo <- simulateResiduals(mod1_redo)
mod1_resid2_redo <- recalculateResiduals(mod1_resid_redo, group = data$Site) # group residuals by site
testSpatialAutocorrelation(mod1_resid2_redo, unique(data$Long), unique(data$Lat))
##no longer significant


### Pairwise tests of Nosema prevalence across host species and sites

# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Nosema prev between each host spp within each site)
a <- emmeans(mod1_redo, specs = pairwise ~ Genus|Site, type = "response")
a

# calculates the pairwise tests among each pair of sites for both host species 
#(determines if there are sig differences in Nosema prev between each site for each host)
b <- emmeans(mod1_redo, specs = pairwise ~ Site|Genus, type = "response")
b

##########################
# Appendix Table S7
##########################
PrevByHost <- a$emmeans %>%
  as.data.frame()
PrevByHost
write.csv(PrevByHost, file = "tables/AppendixTableS7_Nosema_Predicted_Prevalence_By_Host_and_Site.csv", quote = F)

##########################
# Appendix Table S8
##########################
# data frame of pairwise tests for each site for each species
pairwise_sites <- b$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()
pairwise_sites
write.csv(pairwise_sites, file = "tables/AppendixTableS8_Nosema_Pairwise_Tests__BySiteForEachHost.csv", quote = F)

##########################
# Appendix Table S9
##########################
# data frame of pairwise tests for each species at each site
pairwise_hosts <- a$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()
pairwise_hosts
write.csv(pairwise_hosts, file = "tables/AppendixTableS9_Nosema_Pairwise_Tests__ByHostWithinEachSite.csv", quote = F)



##########################
## Figure 1
##########################

##### plot of nesema prevalence by site for each host species
# model predictions of nosema prevalence for each site and host species
me <- ggpredict(mod1_redo, c("Site", "Genus"))

# color scheme for Apis and Bombus
Tol_muted2 <- c('#882255','#44AA99')

prev_site <- plot(me) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  labs(x = "Site", y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence (%)"), title = NULL) +
  annotate(geom= "text", x = 5, y = 0.05, label = expression("Genus" %*% "Site: p = 0.4"), color = "black", size = 3.5) +
  coord_cartesian(ylim = c(0,1)) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title = element_text(size=10, color="black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9))
print(prev_site)
ggsave(here("figures/Fig1_NosemaPrev_BY_Site_BOTH_Hosts_predict.tiff"), plot = prev_site, dpi = 300, width = 12, height = 8, units = "cm", compression="lzw")






##########################################
##########GLMM for Apis only: evaluating effects of visitation to flowers on Nosema in Apis
##########################################

#pearson's for apis

cor(data_Apis$VisitRichnessPerFlower, data_Apis$APIS_visits)
cor(data_Apis$VisitRichnessPerFlower, data_Apis$BOMB_visits)
cor(data_Apis$VisitRichnessPerFlower, data_Apis$Other_visits)
cor(data_Apis$APIS_visits, data_Apis$BOMB_visits)


# Scale transform all main factors
data_Apis$VisitRichnessPerFlower_z <- as.numeric(scale(data_Apis$VisitRichnessPerFlower))
data_Apis$VisitDur_z <- as.numeric(scale(log(data_Apis$VisitDur+1)))
data_Apis$VisitNum_z <- as.numeric(scale(log(data_Apis$VisitNum+1)))
data_Apis$APIS_visits_z <- as.numeric(scale(log(data_Apis$APIS_visits+1)))
data_Apis$BOMB_visits_z <- as.numeric(scale(log(data_Apis$BOMB_visits+1)))
data_Apis$APIS_freq_z <- as.numeric(scale(data_Apis$APIS_freq))
data_Apis$BOMB_freq_z <- as.numeric(scale(data_Apis$BOMB_freq))
data_Apis$Native_visits_z <- as.numeric(scale(log(data_Apis$Native_visits+1)))
data_Apis$Other_visits_z <- as.numeric(scale(log(data_Apis$Other_visits+1)))
data_Apis$Native_freq_z <- as.numeric(scale(data_Apis$Native_freq))
data_Apis$Other_freq_z <- as.numeric(scale(data_Apis$Other_freq))
data_Apis$APIS_rate_z <- as.numeric(scale(log(data_Apis$APIS_rate+1)))
data_Apis$BOMB_rate_z <- as.numeric(scale(log(data_Apis$BOMB_rate+1)))
data_Apis$Other_rate_z <- as.numeric(scale(log(data_Apis$Other_rate+1)))
data_Apis$APIS_visitdur_z <- as.numeric(scale(log(data_Apis$APIS_visitdur+1)))
data_Apis$BOMB_visitdur_z <- as.numeric(scale(log(data_Apis$BOMB_visitdur+1)))
data_Apis$Other_visitdur_z <- as.numeric(scale(log(data_Apis$Other_visitdur+1)))
data_Apis$APIS_visitdur2_z <- as.numeric(scale(log(data_Apis$APIS_visitdur2+1)))
data_Apis$APIS_visitdur3_z <- as.numeric(scale(log(data_Apis$APIS_visitdur3+1)))
data_Apis$APIS_visitdur4_z <- as.numeric(scale(log(data_Apis$APIS_visitdur4+1)))
data_Apis$APIS_visitdur5_z <- as.numeric(scale(log(data_Apis$APIS_visitdur5+1)))
data_Apis$BOMB_visitdur2_z <- as.numeric(scale(log(data_Apis$BOMB_visitdur2+1)))
data_Apis$BOMB_visitdur3_z <- as.numeric(scale(log(data_Apis$BOMB_visitdur3+1)))
data_Apis$BOMB_visitdur4_z <- as.numeric(scale(log(data_Apis$BOMB_visitdur4+1)))
data_Apis$BOMB_visitdur5_z <- as.numeric(scale(log(data_Apis$BOMB_visitdur5+1)))
data_Apis$Other_visitdur2_z <- as.numeric(scale(log(data_Apis$Other_visitdur2+1)))
data_Apis$Other_visitdur3_z <- as.numeric(scale(log(data_Apis$Other_visitdur3+1)))
data_Apis$Other_visitdur4_z <- as.numeric(scale(log(data_Apis$Other_visitdur4+1)))
data_Apis$Other_visitdur5_z <- as.numeric(scale(log(data_Apis$Other_visitdur5+1)))
data_Apis$APIS_dur_z <- as.numeric(scale(log(data_Apis$APIS_dur+1)))
data_Apis$BOMB_dur_z <- as.numeric(scale(log(data_Apis$BOMB_dur+1)))
data_Apis$Other_dur_z <- as.numeric(scale(log(data_Apis$Other_dur+1)))
data_Apis$APIS_dur2_z <- as.numeric(scale(log(data_Apis$APIS_dur2+1)))
data_Apis$APIS_dur3_z <- as.numeric(scale(log(data_Apis$APIS_dur3+1)))
data_Apis$APIS_dur4_z <- as.numeric(scale(log(data_Apis$APIS_dur4+1)))
data_Apis$APIS_dur5_z <- as.numeric(scale(log(data_Apis$APIS_dur5+1)))
data_Apis$BOMB_dur2_z <- as.numeric(scale(log(data_Apis$BOMB_dur2+1)))
data_Apis$BOMB_dur3_z <- as.numeric(scale(log(data_Apis$BOMB_dur3+1)))
data_Apis$BOMB_dur4_z <- as.numeric(scale(log(data_Apis$BOMB_dur4+1)))
data_Apis$BOMB_dur5_z <- as.numeric(scale(log(data_Apis$BOMB_dur5+1)))
data_Apis$Other_dur2_z <- as.numeric(scale(log(data_Apis$Other_dur2+1)))
data_Apis$Other_dur3_z <- as.numeric(scale(log(data_Apis$Other_dur3+1)))
data_Apis$Other_dur4_z <- as.numeric(scale(log(data_Apis$Other_dur4+1)))
data_Apis$Other_dur5_z <- as.numeric(scale(log(data_Apis$Other_dur5+1)))



# Model of visit number for Nosema in Apis
fit24 <- glmer(Nosema ~ APIS_visits_z + BOMB_visits_z + Other_visits_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit24)
vif(fit24)
vif.mer(fit24)
plot(fit24)
qqnorm(resid(fit24))
qqline(resid(fit24))
overdisp_fun(fit24)

# Model of visit frequency for Nosema in Apis  (not shown in manuscript)
fit26 <- glmer(Nosema ~ APIS_freq_z + BOMB_freq_z + Other_freq_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit26)
vif(fit26)
vif.mer(fit26)
plot(fit26)
qqnorm(resid(fit26))
qqline(resid(fit26))
overdisp_fun(fit26)



# Model of visit duration for Nosema in Apis
fit28 <- glmer(Nosema ~ APIS_visitdur_z + BOMB_visitdur_z + Other_visitdur_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit28)
vif(fit28)
plot(fit28)
qqnorm(resid(fit28))
qqline(resid(fit28))
overdisp_fun(fit28)

cor(data$APIS_dur, data$Other_dur)


fit28a <- glmer(Nosema ~ APIS_visitdur2_z + BOMB_visitdur2_z + Other_visitdur2_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit28a)
vif(fit28a)
overdisp_fun(fit28a)

fit28b <- glmer(Nosema ~ APIS_visitdur3_z + BOMB_visitdur3_z + Other_visitdur3_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit28b)
vif(fit28b)
overdisp_fun(fit28b)

fit28c <- glmer(Nosema ~ APIS_visitdur4_z + BOMB_visitdur4_z + Other_visitdur4_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit28c)
vif(fit28c)
overdisp_fun(fit28c)

fit28d <- glmer(Nosema ~ APIS_visitdur5_z + BOMB_visitdur5_z + Other_visitdur5_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit28d)
vif(fit28d)
overdisp_fun(fit28d)


##################################
## Appendix Table S4: Apis visitation model outputs
##################################
summary(fit24)  # visit number
summary(fit28)  # total visit duration
summary(fit28a) # duration of petal only visits
summary(fit28b) # duration of nectar only visits
summary(fit28c) # duration of pollen only visits
summary(fit28d) # duration of pollen+nectar visits







##########################################
##########GLMM for Bombus only: evaluating effects of visitation to flowers on Nosema in Bombus
##########################################

head(data_Bombus)

data_Bombus$VisitRichnessPerFlower_z <- as.numeric(scale(data_Bombus$VisitRichnessPerFlower))
data_Bombus$VisitDur_z <- as.numeric(scale(log(data_Bombus$VisitDur+1)))
data_Bombus$VisitNum_z <- as.numeric(scale(log(data_Bombus$VisitNum+1)))
data_Bombus$APIS_visits_z <- as.numeric(scale(log(data_Bombus$APIS_visits+1)))
data_Bombus$BOMB_visits_z <- as.numeric(scale(log(data_Bombus$BOMB_visits+1)))
data_Bombus$APIS_freq_z <- as.numeric(scale(data_Bombus$APIS_freq))
data_Bombus$BOMB_freq_z <- as.numeric(scale(data_Bombus$BOMB_freq))
data_Bombus$Native_visits_z <- as.numeric(scale(log(data_Bombus$Native_visits+1)))
data_Bombus$Other_visits_z <- as.numeric(scale(log(data_Bombus$Other_visits+1)))
data_Bombus$Native_freq_z <- as.numeric(scale(data_Bombus$Native_freq))
data_Bombus$Other_freq_z <- as.numeric(scale(data_Bombus$Other_freq))
data_Bombus$APIS_visitdur_z <- as.numeric(scale(log(data_Bombus$APIS_visitdur+1)))
data_Bombus$BOMB_visitdur_z <- as.numeric(scale(log(data_Bombus$BOMB_visitdur+1)))
data_Bombus$Other_visitdur_z <- as.numeric(scale(log(data_Bombus$Other_visitdur+1)))
data_Bombus$APIS_visitdur2_z <- as.numeric(scale(log(data_Bombus$APIS_visitdur2+1)))
data_Bombus$APIS_visitdur3_z <- as.numeric(scale(log(data_Bombus$APIS_visitdur3+1)))
data_Bombus$APIS_visitdur4_z <- as.numeric(scale(log(data_Bombus$APIS_visitdur4+1)))
data_Bombus$APIS_visitdur5_z <- as.numeric(scale(log(data_Bombus$APIS_visitdur5+1)))
data_Bombus$BOMB_visitdur2_z <- as.numeric(scale(log(data_Bombus$BOMB_visitdur2+1)))
data_Bombus$BOMB_visitdur3_z <- as.numeric(scale(log(data_Bombus$BOMB_visitdur3+1)))
data_Bombus$BOMB_visitdur4_z <- as.numeric(scale(log(data_Bombus$BOMB_visitdur4+1)))
data_Bombus$BOMB_visitdur5_z <- as.numeric(scale(log(data_Bombus$BOMB_visitdur5+1)))
data_Bombus$Other_visitdur2_z <- as.numeric(scale(log(data_Bombus$Other_visitdur2+1)))
data_Bombus$Other_visitdur3_z <- as.numeric(scale(log(data_Bombus$Other_visitdur3+1)))
data_Bombus$Other_visitdur4_z <- as.numeric(scale(log(data_Bombus$Other_visitdur4+1)))
data_Bombus$Other_visitdur5_z <- as.numeric(scale(log(data_Bombus$Other_visitdur5+1)))


# pearson's correlations
cor(data_Bombus$VisitRichnessPerFlower, data_Bombus$APIS_visits)
cor(data_Bombus$VisitRichnessPerFlower, data_Bombus$BOMB_visits)
cor(data_Bombus$VisitRichnessPerFlower, data_Bombus$Other_visits)

cor(data_Bombus$APIS_visits, data_Bombus$BOMB_visits)


# Model of visit number for Nosema in Bombus 
fit36 <- glmer(Nosema ~ APIS_visits_z + BOMB_visits_z + Other_visits_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
vif(fit36)
summary(fit36)
plot(fit36)
qqnorm(resid(fit36))
qqline(resid(fit36))
overdisp_fun(fit36)


# Model of visit frequency for Nosema in Bombus (Not shown in manuscript)
fit38 <- glmer(Nosema ~ APIS_freq_z + BOMB_freq_z + Other_freq_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
vif(fit38)
summary(fit38)
plot(fit38)
qqnorm(resid(fit38))
qqline(resid(fit38))
overdisp_fun(fit38)


# Model of visit duration for Nosema in Bombus
fit40 <- glmer(Nosema ~ APIS_visitdur_z + BOMB_visitdur_z + Other_visitdur_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit40)
vif(fit40)
plot(fit40)
qqnorm(resid(fit40))
qqline(resid(fit40))
overdisp_fun(fit40)



fit40a <- glmer(Nosema ~ APIS_visitdur2_z + BOMB_visitdur2_z + Other_visitdur2_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit40a)
vif(fit40a)
overdisp_fun(fit40a)

fit40b <- glmer(Nosema ~ APIS_visitdur3_z + BOMB_visitdur3_z + Other_visitdur3_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit40b)
vif(fit40b)
overdisp_fun(fit40b)

fit40c <- glmer(Nosema ~ APIS_visitdur4_z + BOMB_visitdur4_z + Other_visitdur4_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit40c)
vif(fit40c)
overdisp_fun(fit40c)

fit40d <- glmer(Nosema ~ APIS_visitdur5_z + BOMB_visitdur5_z + Other_visitdur5_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit40d)
vif(fit40d)
overdisp_fun(fit40d)


##################################
## Appendix Table S5: Bombus visitation model outputs
##################################
summary(fit36)  # visit number
summary(fit40)  # total visit duration
summary(fit40a) # duration of petal only visits
summary(fit40b) # duration of nectar only visits
summary(fit40c) # duration of pollen only visits
summary(fit40d) # duration of pollen+nectar visits



##################################
## Appendix Table S6: Moran's I
##################################
#Test for spatial autocorrelation for all Nosema in Apis ande Bombus models

# Nosema in Apis vs visitaiton number
fit24_resid <- simulateResiduals(fit24)
fit24_resid2 <- recalculateResiduals(fit24_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit24_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs total visit duration
fit28_resid <- simulateResiduals(fit28)
fit28_resid2 <- recalculateResiduals(fit28_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit28_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs petal only visit duration
fit28a_resid <- simulateResiduals(fit28a)
fit28a_resid2 <- recalculateResiduals(fit28a_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit28a_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs nectar only visit duration
fit28b_resid <- simulateResiduals(fit28b)
fit28b_resid2 <- recalculateResiduals(fit28b_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit28b_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs pollen only visit duration
fit28c_resid <- simulateResiduals(fit28c)
fit28c_resid2 <- recalculateResiduals(fit28c_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit28c_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs pollen+nectar visit duration
fit28d_resid <- simulateResiduals(fit28d)
fit28d_resid2 <- recalculateResiduals(fit28d_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit28d_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 


# Nosema in Bombus vs visitaiton number
fit36_resid <- simulateResiduals(fit36)
fit36_resid2 <- recalculateResiduals(fit36_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit36_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs total visit duration
fit40_resid <- simulateResiduals(fit40)
fit40_resid2 <- recalculateResiduals(fit40_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit40_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs petal only visit duration
fit40a_resid <- simulateResiduals(fit40a)
fit40a_resid2 <- recalculateResiduals(fit40a_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit40a_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs nectar only visit duration
fit40b_resid <- simulateResiduals(fit40b)
fit40b_resid2 <- recalculateResiduals(fit40b_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit40b_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs pollen only visit duration
fit40c_resid <- simulateResiduals(fit40c)
fit40c_resid2 <- recalculateResiduals(fit40c_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit40c_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs pollen+nectar visit duration
fit40d_resid <- simulateResiduals(fit40d)
fit40d_resid2 <- recalculateResiduals(fit40d_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit40d_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 



#### Figures of Nosema prevalence vs visit number and duration

# Load libraries
library(grid)
library(gridExtra)
library(reshape)
library(ggiraph)
library(ggiraphExtra)
library(dplyr)
library(moonBook)
library(ggeffects)
library(colorBlindness)
library(emmeans)


# function to make a grid plot with a shared legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

# color palette
Tol_muted <- c('#882255','#44AA99', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#DDDDDD', '#999933')
Tol_muted2 <- c('#882255','#44AA99')


# extract model marginal means for all models

#################################
#### Apis visit number figure
#################################
# Nosema in Apis
me_ApisPrev_ApisNum <- ggpredict(fit24, c("APIS_visits_z"))
plot(me_ApisPrev_ApisNum, add.data = T)
me_ApisPrev_ApisNum$Host_Species <- "Apis"

# recalculate original Apis visit number 
me_ApisPrev_ApisNum$x # scaled Apis visit number
ApisNumA_mean <- mean(log(data_Apis$APIS_visits+1)) # mean of original richness from disease data set
ApisNumA_sd <- sd(log(data_Apis$APIS_visits+1)) # sd of original richness from disease data set

me_ApisPrev_ApisNum$ApisVisits <- t((t(me_ApisPrev_ApisNum$x) * ApisNumA_sd) + ApisNumA_mean)

# Nosema in Bombus
me_BombPrev_ApisNum <- ggpredict(fit36, c("APIS_visits_z"))
plot(me_BombPrev_ApisNum, add.data = T)
me_BombPrev_ApisNum$Host_Species <- "Bombus"

# recalculate original Apis visit number 
me_BombPrev_ApisNum$x # scaled Apis visit number
ApisNumB_mean <- mean(log(data_Bombus$APIS_visits+1)) # mean of original richness from disease data set
ApisNumB_sd <- sd(log(data_Bombus$APIS_visits+1)) # sd of original richness from disease data set

me_BombPrev_ApisNum$ApisVisits <- t((t(me_BombPrev_ApisNum$x) * ApisNumB_sd) + ApisNumB_mean)


# combine by rows
me_ApisNum <- rbind(me_ApisPrev_ApisNum, me_BombPrev_ApisNum)


# plot of Apis visit number vs nosema prevalence
prev_ApisNum <- ggplot(me_ApisNum, aes(x = ApisVisits, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2, outline.type = NULL) +
  scale_linetype_manual(values=c("dashed", "solid"), guide = "none") +
  labs(color = "Species", tag = "a", x = bquote(atop("Avg. Number of", "Honeybee Visits (per 30 min)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisNum)
#ggsave("VirusPrev_Richness_predict.png", plot = virus_rich, dpi = 300, width = 18, height = 8, units = "cm")



######################################
## Apis visit duration (total) figure
#####################################
# Nosema in Apis
me_ApisPrev_ApisDur <- ggpredict(fit28, c("APIS_visitdur_z"))
plot(me_ApisPrev_ApisDur, add.data = T)
me_ApisPrev_ApisDur$Host_Species <- "Apis"

# recalculate original Apis visit duration 
me_ApisPrev_ApisDur$x # scaled Apis visit duration
ApisDurA_mean <- mean(log(data_Apis$APIS_visitdur+1)) # mean of original apis visit duration
ApisDurA_sd <- sd(log(data_Apis$APIS_visitdur+1)) # sd of original apis visit duration

me_ApisPrev_ApisDur$ApisDur <- t((t(me_ApisPrev_ApisDur$x) * ApisDurA_sd) + ApisDurA_mean)

# Nosema in Bombus
me_BombPrev_ApisDur <- ggpredict(fit40, c("APIS_visitdur_z"))
plot(me_BombPrev_ApisDur, add.data = T)
me_BombPrev_ApisDur$Host_Species <- "Bombus"

# recalculate original Apis visit duration 
me_BombPrev_ApisDur$x # scaled Apis visit duration
ApisDurB_mean <- mean(log(data_Bombus$APIS_visitdur+1)) # mean of original apis visit duration 
ApisDurB_sd <- sd(log(data_Bombus$APIS_visitdur+1)) # sd of original apis visit duration

me_BombPrev_ApisDur$ApisDur <- t((t(me_BombPrev_ApisDur$x) * ApisDurB_sd) + ApisDurB_mean)


# combine by rows
me_ApisDur <- rbind(me_ApisPrev_ApisDur, me_BombPrev_ApisDur)


# plot of Apis visit duration total vs nosema prevalence
prev_ApisDur <- ggplot(me_ApisDur, aes(x = ApisDur, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2, outline.type = NULL) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "b", x = bquote(atop("Avg. Duration of", "Honeybee Visits (seconds)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
                                 legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisDur)


###############################################
## Apis visit duration on pollen+nectar figure
###############################################
# Nosema in Apis
me_ApisPrev_ApisDur_PN <- ggpredict(fit28d, c("APIS_visitdur5_z"))
plot(me_ApisPrev_ApisDur_PN, add.data = T)
me_ApisPrev_ApisDur_PN$Host_Species <- "Apis"

# recalculate original Apis visit duration on pollen+nectar 
me_ApisPrev_ApisDur_PN$x # scaled Apis visit duration on pollen+nectar 
ApisDurA_mean <- mean(log(data_Apis$APIS_visitdur5+1)) # mean of original apis visit duration
ApisDurA_sd <- sd(log(data_Apis$APIS_visitdur5+1)) # sd of original apis visit duration

me_ApisPrev_ApisDur_PN$ApisDurPN <- t((t(me_ApisPrev_ApisDur_PN$x) * ApisDurA_sd) + ApisDurA_mean)

# Nosema in Bombus
me_BombPrev_ApisDur_PN <- ggpredict(fit40d, c("APIS_visitdur5_z"))
plot(me_BombPrev_ApisDur_PN, add.data = T)
me_BombPrev_ApisDur_PN$Host_Species <- "Bombus"

# recalculate original Apis visit duration on pollen+nectar 
me_BombPrev_ApisDur_PN$x # scaled Apis visit duration on pollen+nectar 
ApisDurB_mean <- mean(log(data_Bombus$APIS_visitdur5+1)) # mean of original apis visit duration 
ApisDurB_sd <- sd(log(data_Bombus$APIS_visitdur5+1)) # sd of original apis visit duration

me_BombPrev_ApisDur_PN$ApisDurPN <- t((t(me_BombPrev_ApisDur_PN$x) * ApisDurB_sd) + ApisDurB_mean)


# combine by rows
me_ApisDur_PN <- rbind(me_ApisPrev_ApisDur_PN, me_BombPrev_ApisDur_PN)


# plot of Apis visit duration Pollen+Nectar vs nosema prevalence
prev_ApisDur_PN <- ggplot(me_ApisDur_PN, aes(x = ApisDurPN, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2, outline.type = NULL) +
  scale_linetype_manual(values=c("dashed", "solid"), guide = "none") +
  labs(color = "Species", tag = "c", x = bquote(atop("Avg. Duration of Honeybee", "Pollen + Nectar Visits (seconds)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisDur_PN)




###############################
### FULL FIGURE 2 in Manuscript
###############################
# Grid of both Apis visit number and duration plots
prev_fig2 <- grid_arrange_shared_legend(prev_ApisNum, prev_ApisDur, prev_ApisDur_PN, nrow=1, ncol = 3, position = "right")
ggsave(here("figures/Fig2.tiff"), plot = prev_fig2, dpi = 600, width = 8, height = 3.5, units = "in", compression="lzw")

# color blind check
cvdPlot(prev_fig2)




##### Additional Figures not included in the manuscript

###########################################
#### Apis visit duration on PETALS figure
###########################################
# Nosema in Apis
me_ApisPrev_ApisDur2 <- ggpredict(fit28a, c("APIS_visitdur2_z"))
plot(me_ApisPrev_ApisDur2, add.data = T)
me_ApisPrev_ApisDur2$Host_Species <- "Apis"

# recalculate original Apis visit duration of behavior 2 
me_ApisPrev_ApisDur2$x # scaled Apis duration of behavior 2 
ApisDur2A_mean <- mean(log(data_Apis$APIS_visitdur2+1)) # mean of original Apis duration of behavior 2 
ApisDur2A_sd <- sd(log(data_Apis$APIS_visitdur2+1)) # sd of original Apis duration of behavior 2 

me_ApisPrev_ApisDur2$ApisDur2 <- t((t(me_ApisPrev_ApisDur2$x) * ApisDur2A_sd) + ApisDur2A_mean)

# Nosema in Bombus
me_BombPrev_ApisDur2 <- ggpredict(fit40a, c("APIS_visitdur2_z"))
plot(me_BombPrev_ApisDur2, add.data = T)
me_BombPrev_ApisDur2$Host_Species <- "Bombus"

# recalculate original Apis visit duration of behavior 2  
me_BombPrev_ApisDur2$x # scaled Apis visit duration of behavior 2 
ApisDur2B_mean <- mean(log(data_Bombus$APIS_visitdur2+1)) # mean of original Apis duration of behavior 2 
ApisDur2B_sd <- sd(log(data_Bombus$APIS_visitdur2+1)) # sd of original Apis duration of behavior 2 

me_BombPrev_ApisDur2$ApisDur2 <- t((t(me_BombPrev_ApisDur2$x) * ApisDur2B_sd) + ApisDur2B_mean)


# combine by rows
me_ApisDur2 <- rbind(me_ApisPrev_ApisDur2, me_BombPrev_ApisDur2)


# plot of Apis visit duration total vs nosema prevalence
prev_ApisDur2 <- ggplot(me_ApisDur2, aes(x = ApisDur2, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  geom_line(aes(color = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2, outline.type = NULL) +
  labs(color = "Species", tag = "a", x = bquote("Avg. Duration of" ~ italic(Apis) ~ "on petals (seconds)"), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisDur2)

###########################################
## Bombus visit duration on PETALS figure
###########################################
# Nosema in Apis
me_ApisPrev_BombusDur2 <- ggpredict(fit28a, c("BOMB_visitdur2_z"))
plot(me_ApisPrev_BombusDur2, add.data = T)
me_ApisPrev_BombusDur2$Host_Species <- "Apis"

# recalculate original Bombus visit duration of behavior 2  
me_ApisPrev_BombusDur2$x # scaled Bombus visit duration of behavior 2 
BombusDurA_mean <- mean(log(data_Apis$BOMB_visitdur2+1)) # mean of original Bombus visit duration of behavior 2  
BombusDurA_sd <- sd(log(data_Apis$BOMB_visitdur2+1)) # sd of original Bombus visit duration of behavior 2  

me_ApisPrev_BombusDur2$BombusDur2 <- t((t(me_ApisPrev_BombusDur2$x) * BombusDurA_sd) + BombusDurA_mean)

# Nosema in Bombus
me_BombPrev_BombusDur2 <- ggpredict(fit40a, c("BOMB_visitdur2_z"))
plot(me_BombPrev_BombusDur2, add.data = T)
me_BombPrev_BombusDur2$Host_Species <- "Bombus"

# recalculate original Apis visit number 
me_BombPrev_BombusDur2$x # scaled Apis visit number
BombusDurB_mean <- mean(log(data_Bombus$BOMB_visitdur2+1)) # mean of original richness from disease data set
BombusDurB_sd <- sd(log(data_Bombus$BOMB_visitdur2+1)) # sd of original richness from disease data set

me_BombPrev_BombusDur2$BombusDur2 <- t((t(me_BombPrev_BombusDur2$x) * BombusDurB_sd) + BombusDurB_mean)


# combine by rows
me_BombusDur2 <- rbind(me_ApisPrev_BombusDur2, me_BombPrev_BombusDur2)


# plot of Apis visit duration total vs nosema prevalence
prev_BombusDur2 <- ggplot(me_BombusDur2, aes(x = BombusDur2, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  geom_line(aes(color = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2, outline.type = NULL) +
  labs(color = "Species", tag = "b", x = bquote("Avg. Duration of" ~ italic(Bombus) ~ "on petals (seconds)"), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  #facet_wrap(~facet) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_BombusDur2)



# Grid of both Apis and bombus duration on PETALS
prev_fig3 <- grid_arrange_shared_legend(prev_ApisDur2, prev_BombusDur2, nrow=1, ncol = 2, position = "right")
ggsave("Fig3.tiff", plot = prev_fig3, dpi = 600, width = 6.5, height = 3, units = "in", compression="lzw")

# color blind check
cvdPlot(prev_fig3)


