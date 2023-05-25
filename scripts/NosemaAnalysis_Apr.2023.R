# Code for the analyses for "Honeybee visitation behavior on shared flowers 
# increases Vairimorpha ceranae prevalence in bumblebees"

# Manuscript submitted to: Ecology and Evolution


# Vairimorpha (=Nosema) ceranae prevalence in Apis and Bombus analysis

# This script includes the full analyses and figures for how pollinator visitation number and duration of visits to flowers
# impact V. ceranae prevalence in Apis mellifera and Bombus spp. 


# Written by: Michelle Fearon and Maryellen Zbrozek
# Last updated: 25 May 2023



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
here::i_am("scripts/NosemaAnalysis_Apr.2023.R")

###############################
### LOAD THE DATA FOR ANALYSIS
###############################
data <- read.csv(here("data/NosemaAnalysis.csv"), stringsAsFactors = F)

head(data)
summary(data)

# convert variables into characters that were misread in by R
data$Year <- as.character(data$Year)
data$Visit <- as.character(data$Visit)


#########################################
#  Appendix S1, Table S4: Count of Number tested and Nosema positive per genus per visit per site
#########################################
sum.data <- data %>%
  dplyr::group_by(Genus, Site, Visit) %>%
  summarise(N_tested = length(Nosema),
          N_postive = sum(Nosema))

write.csv(sum.data, file = "tables/AppendixTableS4_Counts_of_tested+pos_perSiteVisitGenus.csv", quote = F)


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


#####################################
# Calculate total Nosema prevalence in Apis and Bombus PER SITE
#####################################
Apis_site_prev <- matrix(NA, ncol = 7, nrow = 6)
rownames(Apis_site_prev) <- c(unique(data$Site))
colnames(Apis_site_prev) <- c("Site", "Genus", "N_tested", "N_Nosema", "Nosema_prev", "Nosema_Lower", "Nosema_Upper")

SITE <- unique(data$Site)


for (i in SITE){
  temp_site <- data[ data$Site == i & data$Genus == "Apis", ]
  Apis_site_prev[i , "Site"] <- i
  Apis_site_prev[i , "Genus"] <- "Apis"
  
  ### Calculate Nosema prevalence for both species at each site
  # calculate the total number of infected ind for nosema (store in site_summary matrix) [BOTH SPECIES]
  infectedNosema <- sum(temp_site$Nosema)
  Apis_site_prev[ i , "N_Nosema"] <- infectedNosema
  
  # calculate the total number tested
  tested <- length(temp_site$Genus)
  Apis_site_prev[ i , "N_tested"] <- tested
  
  # unlist prevalence estimates to make easier to work with
  Nos <- unlist(epi.prev(pos = infectedNosema, tested = tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))
  
  
  # store total nosema prevalence and lower and upper confidence interval in site_summary matrix
  Apis_site_prev[ i, "Nosema_prev"] <- Nos[1]
  Apis_site_prev[ i, "Nosema_Lower"] <- Nos[2]
  Apis_site_prev[ i, "Nosema_Upper"] <- Nos[3]
}
Apis_site_prev


Bombus_site_prev <- matrix(NA, ncol = 7, nrow = 6)
rownames(Bombus_site_prev) <- c(unique(data$Site))
colnames(Bombus_site_prev) <- c("Site", "Genus", "N_tested", "N_Nosema", "Nosema_prev", "Nosema_Lower", "Nosema_Upper")

for (i in SITE){
  temp_site <- data[ data$Site == i & data$Genus == "Bombus", ]
  Bombus_site_prev[i , "Site"] <- i
  Bombus_site_prev[i , "Genus"] <- "Bombus"
  
  ### Calculate Nosema prevalence for both species at each site
  # calculate the total number of infected ind for nosema (store in site_summary matrix) [BOTH SPECIES]
  infectedNosema <- sum(temp_site$Nosema)
  Bombus_site_prev[ i , "N_Nosema"] <- infectedNosema
  
  # calculate the total number tested
  tested <- length(temp_site$Genus)
  Bombus_site_prev[ i , "N_tested"] <- tested
  
  # unlist prevalence estimates to make easier to work with
  Nos <- unlist(epi.prev(pos = infectedNosema, tested = tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))
  
  
  # store total nosema prevalence and lower and upper confidence interval in site_summary matrix
  Bombus_site_prev[ i, "Nosema_prev"] <- Nos[1]
  Bombus_site_prev[ i, "Nosema_Lower"] <- Nos[2]
  Bombus_site_prev[ i, "Nosema_Upper"] <- Nos[3]
}
Bombus_site_prev


##########################
# Appendix S1, Table S9
##########################
# combine matrices into a data frame and convert data to numeric classes
site_prev <- as.data.frame(rbind(Apis_site_prev, Bombus_site_prev))
rownames(site_prev) <- NULL
site_prev$N_tested <- as.numeric(site_prev$N_tested)
site_prev$Nosema_prev <- as.numeric(site_prev$Nosema_prev)
site_prev$Nosema_Lower <- as.numeric(site_prev$Nosema_Lower)
site_prev$Nosema_Upper <- as.numeric(site_prev$Nosema_Upper)


write.csv(site_prev, file = "tables/AppendixTableS10_Nosema_Predicted_Prevalence_By_Host_and_Site.csv", quote = F)









####################
# Appendix S1: Figure S2
#################### 

###### Figure of Nosema prevalence per site for each species

# color scheme for Apis and Bombus
Tol_muted2 <- c('#882255','#44AA99')

prev_site_plot <- ggplot(site_prev, aes(x = Site, y = Nosema_prev, color = Genus)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_errorbar(aes(ymin = Nosema_Lower, ymax = Nosema_Upper), position = position_dodge(width = 0.4), width = 0.2) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  labs(x = "Site", y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence (%)"), title = NULL) +
  coord_cartesian(ylim = c(0,100)) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title = element_text(size=10, color="black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9))
print(prev_site_plot)
ggsave(here("figures/FigS2_NosemaPrev_BY_Site_BOTH_Hosts.tiff"), plot = prev_site_plot, dpi = 300, width = 5, height = 4, units = "in", compression="lzw")




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
fit_Apis_visits <- glmer(Nosema ~ APIS_visits_z + BOMB_visits_z + Other_visits_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
#fit_Apis_visits <- glmer(Nosema ~ APIS_visits_z + BOMB_visits_z + Other_visits_z + Visit + (1|Site), family = binomial, data = data_Apis)
summary(fit_Apis_visits)
vif(fit_Apis_visits)
vif.mer(fit_Apis_visits)
plot(fit_Apis_visits)
qqnorm(resid(fit_Apis_visits))
qqline(resid(fit_Apis_visits))
overdisp_fun(fit_Apis_visits)
testDispersion(fit_Apis_visits)
testZeroInflation(fit_Apis_visits)
ApisVisit_simResid <- simulateResiduals(fittedModel = fit_Apis_visits)
plot(ApisVisit_simResid) 

# Model of visit rate for Nosema in Apis  (not shown in manuscript)
fit_Apis_rate <- glmer(Nosema ~ APIS_rate_z + BOMB_rate_z + Other_rate_z + Visit + (1|Site), family = binomial, data = data_Apis)
summary(fit_Apis_rate)
vif(fit_Apis_rate)
vif.mer(fit_Apis_rate)
plot(fit_Apis_rate)
qqnorm(resid(fit_Apis_rate))
qqline(resid(fit_Apis_rate))
overdisp_fun(fit_Apis_rate)
testDispersion(fit_Apis_rate)
testZeroInflation(fit_Apis_rate)
ApisRate_simResid <- simulateResiduals(fittedModel = fit_Apis_rate)
plot(ApisRate_simResid) 
plot(fit_Apis_rate_resid2)

# Nosema in Apis vs visitation rate
fit_Apis_rate_resid <- simulateResiduals(fit_Apis_rate)
fit_Apis_rate_resid2 <- recalculateResiduals(fit_Apis_rate_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_rate_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
testSpatialAutocorrelation(fit_Apis_rate_resid, data_Apis$Long, data_Apis$Lat)
# significant spatial autocorrelation! 



# Spatial autocorrelation test below indicates that there is significant spatial autocorrelation
# Update visit rate model for Nosema in Apis
library(spaMM)
fit_Apis_rate_update <- fitme(Nosema ~ APIS_rate_z + BOMB_rate_z + Other_rate_z + (1|Site) + (1|Site:Visit) + Matern(1 | Lat + Long), family = binomial, data_Apis)
summary(fit_Apis_rate_update)
AIC(fit_Apis_rate_update)
overdisp_fun(fit_Apis_rate_update)

dd <- dist(data_Apis[,c("Lat","Long")])*111.139  # Multiply the degrees of separation of longitude and latitude by 111.139 to get the corresponding linear distances in kilometers.
mm <- MaternCorr(dd, nu = 16.6, rho = 16.87355)
plot(as.numeric(dd), as.numeric(mm), xlab = "Distance between pairs of location [in km]", ylab = "Estimated correlation")
# Sites less than 10 km away from each other are very correlated... looks like mostly within a site there is a lot of correlation, but as you move away, it is essentially zero


# UPDATED Nosema in Apis vs visitation rate that corrected for spatial autocorrelation in the model
fit_Apis_rate2_resid <- simulateResiduals(fit_Apis_rate_update)
fit_Apis_rate2_resid2 <- recalculateResiduals(fit_Apis_rate2_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_rate2_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# significant spatial autocorrelation!




# Model of visit frequency for Nosema in Apis  (not shown in manuscript)
fit_Apis_freq <- glmer(Nosema ~ APIS_freq_z + BOMB_freq_z + Other_freq_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_freq)
vif(fit_Apis_freq)
vif.mer(fit_Apis_freq) # multi-collinearity
plot(fit_Apis_freq)
qqnorm(resid(fit_Apis_freq))
qqline(resid(fit_Apis_freq))
overdisp_fun(fit_Apis_freq)



# Model of visit duration (i.e. behavior duration per each bee visit) for Nosema in Apis
fit_Apis_visitdur <- glmer(Nosema ~ APIS_visitdur_z + BOMB_visitdur_z + Other_visitdur_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_visitdur)
vif(fit_Apis_visitdur)
plot(fit_Apis_visitdur)
qqnorm(resid(fit_Apis_visitdur))
qqline(resid(fit_Apis_visitdur))
overdisp_fun(fit_Apis_visitdur)

# visit duration to petals
fit_Apis_visitdur2 <- glmer(Nosema ~ APIS_visitdur2_z + BOMB_visitdur2_z + Other_visitdur2_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_visitdur2)
vif(fit_Apis_visitdur2)
overdisp_fun(fit_Apis_visitdur2)

# visit duration to nectar
fit_Apis_visitdur3 <- glmer(Nosema ~ APIS_visitdur3_z + BOMB_visitdur3_z + Other_visitdur3_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_visitdur3)
vif(fit_Apis_visitdur3)
overdisp_fun(fit_Apis_visitdur3)

# visit duration to pollen
fit_Apis_visitdur4 <- glmer(Nosema ~ APIS_visitdur4_z + BOMB_visitdur4_z + Other_visitdur4_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_visitdur4)
vif(fit_Apis_visitdur4)
overdisp_fun(fit_Apis_visitdur4)

# visit duration to pollen+nectar
fit_Apis_visitdur5 <- glmer(Nosema ~ APIS_visitdur5_z + BOMB_visitdur5_z + Other_visitdur5_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_visitdur5)
vif(fit_Apis_visitdur5)
overdisp_fun(fit_Apis_visitdur5)




# Model of sum duration (i.e. sum of all behavior duration that occurred within the 30 min video) for Nosema in Apis
fit_Apis_dur <- glmer(Nosema ~ APIS_dur_z + BOMB_dur_z + Other_dur_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_dur)
vif(fit_Apis_dur)
plot(fit_Apis_dur)
qqnorm(resid(fit_Apis_dur))
qqline(resid(fit_Apis_dur))
overdisp_fun(fit_Apis_dur)

# sum duration to petals
fit_Apis_dur2 <- glmer(Nosema ~ APIS_dur2_z + BOMB_dur2_z + Other_dur2_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_dur2)
vif(fit_Apis_dur2)
overdisp_fun(fit_Apis_dur2)

# Spatial autocorrelation test below indicates that there is significant spatial autocorrelation
# Update sum duration to petals model for Nosema in Apis
fit_Apis_dur2_update <- fitme(Nosema ~ APIS_dur2_z + BOMB_dur2_z + Other_dur2_z +  Matern(1 | Lat + Long), family = binomial, data_Apis)
summary(fit_Apis_dur2_update)
AIC(fit_Apis_dur2_update)
overdisp_fun(fit_Apis_dur2_update)

dd <- dist(data_Apis[,c("Lat","Long")])*111.139  # Multiply the degrees of separation of longitude and latitude by 111.139 to get the corresponding linear distances in kilometers.
mm <- MaternCorr(dd, nu = 16.6, rho = 20.097)
plot(as.numeric(dd), as.numeric(mm), xlab = "Distance between pairs of location [in km]", ylab = "Estimated correlation")
# Sites less than 10 km away from each other are very correlated... looks like mostly within a site there is a lot of correlation, but as you move away, it is essentially zero

# sum duration to nectar
fit_Apis_dur3 <- glmer(Nosema ~ APIS_dur3_z + BOMB_dur3_z + Other_dur3_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_dur3)
vif(fit_Apis_dur3)
overdisp_fun(fit_Apis_dur3)

# sum duration to pollen
fit_Apis_dur4 <- glmer(Nosema ~ APIS_dur4_z + BOMB_dur4_z + Other_dur4_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_dur4)
vif(fit_Apis_dur4)
overdisp_fun(fit_Apis_dur4)

# sum duration to pollen+nectar
fit_Apis_dur5 <- glmer(Nosema ~ APIS_dur5_z + BOMB_dur5_z + Other_dur5_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Apis)
summary(fit_Apis_dur5)
vif(fit_Apis_dur5)
overdisp_fun(fit_Apis_dur5)


##################################
## Appendix S1, Table S7: Honeybee visitation model outputs
##################################
summary(fit_Apis_visits)    # visit number
summary(fit_Apis_rate)      # visit rate (visit number / total time)
summary(fit_Apis_rate_update)# visit rate (visit number / total time) - updated to account for spatial autocorrelation
summary(fit_Apis_visitdur)  # total duration per visit
summary(fit_Apis_visitdur2) # duration per visit of petal only visits
summary(fit_Apis_visitdur3) # duration per visit of nectar only visits
summary(fit_Apis_visitdur4) # duration per visit of pollen only visits
summary(fit_Apis_visitdur5) # duration per visit of pollen+nectar visits
summary(fit_Apis_dur)  # total sum duration per 30 min
summary(fit_Apis_dur2) # sum duration per 30 min of petal only visits
summary(fit_Apis_dur3) # sum duration per 30 min of nectar only visits
summary(fit_Apis_dur4) # sum duration per 30 min of pollen only visits
summary(fit_Apis_dur5) # sum duration per 30 min of pollen+nectar visits






##########################################
##########GLMM for Bombus only: evaluating effects of visitation to flowers on Nosema in Bombus
##########################################

head(data_Bombus)

data_Bombus$VisitRichnessPerFlower_z <- as.numeric(scale(data_Bombus$VisitRichnessPerFlower))
data_Bombus$VisitDur_z <- as.numeric(scale(log(data_Bombus$VisitDur+1)))
data_Bombus$VisitNum_z <- as.numeric(scale(log(data_Bombus$VisitNum+1)))
data_Bombus$APIS_visits_z <- as.numeric(scale(log(data_Bombus$APIS_visits+1)))
data_Bombus$BOMB_visits_z <- as.numeric(scale(log(data_Bombus$BOMB_visits+1)))
data_Bombus$Native_visits_z <- as.numeric(scale(log(data_Bombus$Native_visits+1)))
data_Bombus$Other_visits_z <- as.numeric(scale(log(data_Bombus$Other_visits+1)))
data_Bombus$APIS_freq_z <- as.numeric(scale(data_Bombus$APIS_freq))
data_Bombus$BOMB_freq_z <- as.numeric(scale(data_Bombus$BOMB_freq))
data_Bombus$Native_freq_z <- as.numeric(scale(data_Bombus$Native_freq))
data_Bombus$Other_freq_z <- as.numeric(scale(data_Bombus$Other_freq))
data_Bombus$APIS_rate_z <- as.numeric(scale(log(data_Bombus$APIS_rate+1)))
data_Bombus$BOMB_rate_z <- as.numeric(scale(log(data_Bombus$BOMB_rate+1)))
data_Bombus$Other_rate_z <- as.numeric(scale(log(data_Bombus$Other_rate+1)))
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
data_Bombus$APIS_dur_z <- as.numeric(scale(log(data_Bombus$APIS_dur+1)))
data_Bombus$BOMB_dur_z <- as.numeric(scale(log(data_Bombus$BOMB_dur+1)))
data_Bombus$Other_dur_z <- as.numeric(scale(log(data_Bombus$Other_dur+1)))
data_Bombus$APIS_dur2_z <- as.numeric(scale(log(data_Bombus$APIS_dur2+1)))
data_Bombus$APIS_dur3_z <- as.numeric(scale(log(data_Bombus$APIS_dur3+1)))
data_Bombus$APIS_dur4_z <- as.numeric(scale(log(data_Bombus$APIS_dur4+1)))
data_Bombus$APIS_dur5_z <- as.numeric(scale(log(data_Bombus$APIS_dur5+1)))
data_Bombus$BOMB_dur2_z <- as.numeric(scale(log(data_Bombus$BOMB_dur2+1)))
data_Bombus$BOMB_dur3_z <- as.numeric(scale(log(data_Bombus$BOMB_dur3+1)))
data_Bombus$BOMB_dur4_z <- as.numeric(scale(log(data_Bombus$BOMB_dur4+1)))
data_Bombus$BOMB_dur5_z <- as.numeric(scale(log(data_Bombus$BOMB_dur5+1)))
data_Bombus$Other_dur2_z <- as.numeric(scale(log(data_Bombus$Other_dur2+1)))
data_Bombus$Other_dur3_z <- as.numeric(scale(log(data_Bombus$Other_dur3+1)))
data_Bombus$Other_dur4_z <- as.numeric(scale(log(data_Bombus$Other_dur4+1)))
data_Bombus$Other_dur5_z <- as.numeric(scale(log(data_Bombus$Other_dur5+1)))


# pearson's correlations
cor(data_Bombus$VisitRichnessPerFlower, data_Bombus$APIS_visits)
cor(data_Bombus$VisitRichnessPerFlower, data_Bombus$BOMB_visits)
cor(data_Bombus$VisitRichnessPerFlower, data_Bombus$Other_visits)

cor(data_Bombus$APIS_visits, data_Bombus$BOMB_visits)


# Model of visit number for Nosema in Bombus 
fit_Bombus_visits <- glmer(Nosema ~ APIS_visits_z + BOMB_visits_z + Other_visits_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
vif(fit_Bombus_visits)
summary(fit_Bombus_visits)
plot(fit_Bombus_visits)
qqnorm(resid(fit_Bombus_visits))
qqline(resid(fit_Bombus_visits))
overdisp_fun(fit_Bombus_visits)


# Model of visit rate for Nosema in Bombus  (not shown in manuscript)
fit_Bombus_rate <- glmer(Nosema ~ APIS_rate_z + BOMB_rate_z + Other_rate_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_rate)
vif(fit_Bombus_rate)
vif.mer(fit_Bombus_rate)
plot(fit_Bombus_rate)
qqnorm(resid(fit_Bombus_rate))
qqline(resid(fit_Bombus_rate))
overdisp_fun(fit_Bombus_rate)


# Model of visit frequency for Nosema in Bombus (Not shown in manuscript)
fit_Bombus_freq <- glmer(Nosema ~ APIS_freq_z + BOMB_freq_z + Other_freq_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
vif(fit_Bombus_freq) # highest VIF is 4.05
summary(fit_Bombus_freq)
plot(fit_Bombus_freq)
qqnorm(resid(fit_Bombus_freq))
qqline(resid(fit_Bombus_freq))
overdisp_fun(fit_Bombus_freq)


# Model of visit duration (i.e. behavior duration per each bee visit) for Nosema in Bombus
fit_Bombus_visitdur <- glmer(Nosema ~ APIS_visitdur_z + BOMB_visitdur_z + Other_visitdur_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_visitdur)
vif(fit_Bombus_visitdur)
plot(fit_Bombus_visitdur)
qqnorm(resid(fit_Bombus_visitdur))
qqline(resid(fit_Bombus_visitdur))
overdisp_fun(fit_Bombus_visitdur)


# visit duration to petals
fit_Bombus_visitdur2 <- glmer(Nosema ~ APIS_visitdur2_z + BOMB_visitdur2_z + Other_visitdur2_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_visitdur2)
vif(fit_Bombus_visitdur2)
overdisp_fun(fit_Bombus_visitdur2)

# visit duration to nectar
fit_Bombus_visitdur3 <- glmer(Nosema ~ APIS_visitdur3_z + BOMB_visitdur3_z + Other_visitdur3_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_visitdur3)
vif(fit_Bombus_visitdur3)
overdisp_fun(fit_Bombus_visitdur3)

# visit duration to pollen
fit_Bombus_visitdur4 <- glmer(Nosema ~ APIS_visitdur4_z + BOMB_visitdur4_z + Other_visitdur4_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_visitdur4)
vif(fit_Bombus_visitdur4)
overdisp_fun(fit_Bombus_visitdur4)

# visit duration to pollen+nectar
fit_Bombus_visitdur5 <- glmer(Nosema ~ APIS_visitdur5_z + BOMB_visitdur5_z + Other_visitdur5_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_visitdur5)
vif(fit_Bombus_visitdur5)
overdisp_fun(fit_Bombus_visitdur5)




# Model of sum duration (i.e. sum of all behavior duration that occurred within the 30 min video) for Nosema in Bombus
fit_Bombus_dur <- glmer(Nosema ~ APIS_dur_z + BOMB_dur_z + Other_dur_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_dur)
vif(fit_Bombus_dur)
plot(fit_Bombus_dur)
qqnorm(resid(fit_Bombus_dur))
qqline(resid(fit_Bombus_dur))
overdisp_fun(fit_Bombus_dur)


# sum duration to petals
fit_Bombus_dur2 <- glmer(Nosema ~ APIS_dur2_z + BOMB_dur2_z + Other_dur2_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_dur2)
vif(fit_Bombus_dur2)
overdisp_fun(fit_Bombus_dur2)

# sum duration to nectar
fit_Bombus_dur3 <- glmer(Nosema ~ APIS_dur3_z + BOMB_dur3_z + Other_dur3_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_dur3)
vif(fit_Bombus_dur3)
overdisp_fun(fit_Bombus_dur3)

# sum duration to pollen
fit_Bombus_dur4 <- glmer(Nosema ~ APIS_dur4_z + BOMB_dur4_z + Other_dur4_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_dur4)
vif(fit_Bombus_dur4)
overdisp_fun(fit_Bombus_dur4)

# sum duration to pollen+nectar
fit_Bombus_dur5 <- glmer(Nosema ~ APIS_dur5_z + BOMB_dur5_z + Other_dur5_z + (1|Site) + (1|Site:Visit), family = binomial, data = data_Bombus)
summary(fit_Bombus_dur5)
vif(fit_Bombus_dur5)
overdisp_fun(fit_Bombus_dur5)




##################################
## Appendix S1, Table S8: Bumblebee visitation model outputs
##################################
summary(fit_Bombus_visits)    # visit number
summary(fit_Bombus_rate)      # visit number / total video time
summary(fit_Bombus_visitdur)  # total duration per visit
summary(fit_Bombus_visitdur2) # duration per visit of petal only visits
summary(fit_Bombus_visitdur3) # duration per visit of nectar only visits
summary(fit_Bombus_visitdur4) # duration per visit of pollen only visits
summary(fit_Bombus_visitdur5) # duration per visit of pollen+nectar visits
summary(fit_Bombus_dur)  # total sum duration per 30 min
summary(fit_Bombus_dur2) # sum duration per 30 min of petal only visits
summary(fit_Bombus_dur3) # sum duration per 30 min of nectar only visits
summary(fit_Bombus_dur4) # sum duration per 30 min of pollen only visits
summary(fit_Bombus_dur5) # sum duration per 30 min of pollen+nectar visits




##################################
## Appendix S1, Table S9: Moran's I
##################################
#Test for spatial autocorrelation for all Nosema in Apis and Bombus models

# Nosema in Apis vs visitation number
fit_Apis_visits_resid <- simulateResiduals(fit_Apis_visits)
fit_Apis_visits_resid2 <- recalculateResiduals(fit_Apis_visits_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_visits_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs visitation rate
fit_Apis_rate_resid <- simulateResiduals(fit_Apis_rate)
fit_Apis_rate_resid2 <- recalculateResiduals(fit_Apis_rate_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_rate_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# significant spatial autocorrelation! 

# UPDATED Nosema in Apis vs visitation rate that corrected for spatial autocorrelation in the model
fit_Apis_rate2_resid <- simulateResiduals(fit_Apis_rate_update)
fit_Apis_rate2_resid2 <- recalculateResiduals(fit_Apis_rate2_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_rate2_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# significant spatial autocorrelation!

# Nosema in Apis vs total duration per visit
fit_Apis_visitdur_resid <- simulateResiduals(fit_Apis_visitdur)
fit_Apis_visitdur_resid2 <- recalculateResiduals(fit_Apis_visitdur_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_visitdur_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs petal only duration per visit
fit_Apis_visitdur2_resid <- simulateResiduals(fit_Apis_visitdur2)
fit_Apis_visitdur2_resid2 <- recalculateResiduals(fit_Apis_visitdur2_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_visitdur2_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs nectar only duration per visit
fit_Apis_visitdur3_resid <- simulateResiduals(fit_Apis_visitdur3)
fit_Apis_visitdur3_resid2 <- recalculateResiduals(fit_Apis_visitdur3_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_visitdur3_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs pollen only duration per visit
fit_Apis_visitdur4_resid <- simulateResiduals(fit_Apis_visitdur4)
fit_Apis_visitdur4_resid2 <- recalculateResiduals(fit_Apis_visitdur4_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_visitdur4_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs pollen+nectar duration per visit
fit_Apis_visitdur5_resid <- simulateResiduals(fit_Apis_visitdur5)
fit_Apis_visitdur5_resid2 <- recalculateResiduals(fit_Apis_visitdur5_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_visitdur5_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs total sum duration per 30 min
fit_Apis_dur_resid <- simulateResiduals(fit_Apis_dur)
fit_Apis_dur_resid2 <- recalculateResiduals(fit_Apis_dur_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_dur_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs petal only sum duration per 30 min
fit_Apis_dur2_resid <- simulateResiduals(fit_Apis_dur2)
fit_Apis_dur2_resid2 <- recalculateResiduals(fit_Apis_dur2_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_dur2_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# significant spatial autocorrelation!  

# UPDATED Nosema in Apis vs petal only sum duration per 30 min
fit_Apis_dur2_resid <- simulateResiduals(fit_Apis_dur2_update)
fit_Apis_dur2_resid2 <- recalculateResiduals(fit_Apis_dur2_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_dur2_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# significant spatial autocorrelation! 

# Nosema in Apis vs nectar only sum duration per 30 min
fit_Apis_dur3_resid <- simulateResiduals(fit_Apis_dur3)
fit_Apis_dur3_resid2 <- recalculateResiduals(fit_Apis_dur3_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_dur3_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs pollen only sum duration per 30 min
fit_Apis_dur4_resid <- simulateResiduals(fit_Apis_dur4)
fit_Apis_dur4_resid2 <- recalculateResiduals(fit_Apis_dur4_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_dur4_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 

# Nosema in Apis vs pollen+nectar sum duration per 30 min
fit_Apis_dur5_resid <- simulateResiduals(fit_Apis_dur5)
fit_Apis_dur5_resid2 <- recalculateResiduals(fit_Apis_dur5_resid, group = data_Apis$Site) # group residuals by site
testSpatialAutocorrelation(fit_Apis_dur5_resid2, unique(data_Apis$Long), unique(data_Apis$Lat))
# not significant! 


# Nosema in Bombus vs visitation number
fit_Bombus_visits_resid <- simulateResiduals(fit_Bombus_visits)
fit_Bombus_visits_resid2 <- recalculateResiduals(fit_Bombus_visits_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_visits_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs visitation rate
fit_Bombus_rate_resid <- simulateResiduals(fit_Bombus_rate)
fit_Bombus_rate_resid2 <- recalculateResiduals(fit_Bombus_rate_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_rate_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs total duration per visit
fit_Bombus_visitdur_resid <- simulateResiduals(fit_Bombus_visitdur)
fit_Bombus_visitdur_resid2 <- recalculateResiduals(fit_Bombus_visitdur_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_visitdur_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs petal only duration per visit
fit_Bombus_visitdur2_resid <- simulateResiduals(fit_Bombus_visitdur2)
fit_Bombus_visitdur2_resid2 <- recalculateResiduals(fit_Bombus_visitdur2_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_visitdur2_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs nectar only duration per visit
fit_Bombus_visitdur3_resid <- simulateResiduals(fit_Bombus_visitdur3)
fit_Bombus_visitdur3_resid2 <- recalculateResiduals(fit_Bombus_visitdur3_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_visitdur3_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs pollen only duration per visit
fit_Bombus_visitdur4_resid <- simulateResiduals(fit_Bombus_visitdur4)
fit_Bombus_visitdur4_resid2 <- recalculateResiduals(fit_Bombus_visitdur4_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_visitdur4_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs pollen+nectar duration per visit
fit_Bombus_visitdur5_resid <- simulateResiduals(fit_Bombus_visitdur5)
fit_Bombus_visitdur5_resid2 <- recalculateResiduals(fit_Bombus_visitdur5_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_visitdur5_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 


# Nosema in Bombus vs total sum duration per 30 min
fit_Bombus_dur_resid <- simulateResiduals(fit_Bombus_dur)
fit_Bombus_dur_resid2 <- recalculateResiduals(fit_Bombus_dur_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_dur_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs petal only sum duration per 30 min
fit_Bombus_dur2_resid <- simulateResiduals(fit_Bombus_dur2)
fit_Bombus_dur2_resid2 <- recalculateResiduals(fit_Bombus_dur2_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_dur2_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs nectar only sum duration per 30 min
fit_Bombus_dur3_resid <- simulateResiduals(fit_Bombus_dur3)
fit_Bombus_dur3_resid2 <- recalculateResiduals(fit_Bombus_dur3_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_dur3_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs pollen only sum duration per 30 min
fit_Bombus_dur4_resid <- simulateResiduals(fit_Bombus_dur4)
fit_Bombus_dur4_resid2 <- recalculateResiduals(fit_Bombus_dur4_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_dur4_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
# not significant! 

# Nosema in Bombus vs pollen+nectar sum duration per 30 min
fit_Bombus_dur5_resid <- simulateResiduals(fit_Bombus_dur5)
fit_Bombus_dur5_resid2 <- recalculateResiduals(fit_Bombus_dur5_resid, group = data_Bombus$Site) # group residuals by site
testSpatialAutocorrelation(fit_Bombus_dur5_resid2, unique(data_Bombus$Long), unique(data_Bombus$Lat))
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
#### Apis visit number figure  (Figure 2 A)
#################################
# Nosema in Apis
me_ApisPrev_ApisNum <- ggpredict(fit_Apis_visits, c("APIS_visits_z"))
plot(me_ApisPrev_ApisNum, add.data = T)
me_ApisPrev_ApisNum$Host_Species <- "Apis"

# recalculate original Apis visit number 
me_ApisPrev_ApisNum$x # scaled Apis visit number
ApisNumA_mean <- mean(log(data_Apis$APIS_visits+1)) # mean of original richness from disease data set
ApisNumA_sd <- sd(log(data_Apis$APIS_visits+1)) # sd of original richness from disease data set

me_ApisPrev_ApisNum$ApisVisits <- t((t(me_ApisPrev_ApisNum$x) * ApisNumA_sd) + ApisNumA_mean)

# Nosema in Bombus
me_BombPrev_ApisNum <- ggpredict(fit_Bombus_visits, c("APIS_visits_z"))
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
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "solid"), guide = "none") +
  labs(color = "Species", tag = "a", x = bquote(atop("Avg. Number of", "Honeybee Visits (per 30 min)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisNum)
#ggsave("VirusPrev_Richness_predict.png", plot = virus_rich, dpi = 300, width = 18, height = 8, units = "cm")



######################################
## Apis visit duration (total) figure (Figure 2 B)
#####################################
# Nosema in Apis
me_ApisPrev_ApisDur <- ggpredict(fit_Apis_visitdur, c("APIS_visitdur_z"))
plot(me_ApisPrev_ApisDur, add.data = T)
me_ApisPrev_ApisDur$Host_Species <- "Apis"

# recalculate original Apis visit duration 
me_ApisPrev_ApisDur$x # scaled Apis visit duration
ApisDurA_mean <- mean(log(data_Apis$APIS_visitdur+1)) # mean of original apis visit duration
ApisDurA_sd <- sd(log(data_Apis$APIS_visitdur+1)) # sd of original apis visit duration

me_ApisPrev_ApisDur$ApisDur <- t((t(me_ApisPrev_ApisDur$x) * ApisDurA_sd) + ApisDurA_mean)

# Nosema in Bombus
me_BombPrev_ApisDur <- ggpredict(fit_Bombus_visitdur, c("APIS_visitdur_z"))
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
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "b", x = bquote(atop("Avg. Duration of Honeybee", " Visits (seconds/visit)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
                                 legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisDur)


###############################################
## Apis visit duration on pollen+nectar figure  (Figure 2 C)
###############################################
# Nosema in Apis
me_ApisPrev_ApisDur_PN <- ggpredict(fit_Apis_visitdur5, c("APIS_visitdur5_z"))
plot(me_ApisPrev_ApisDur_PN, add.data = T)
me_ApisPrev_ApisDur_PN$Host_Species <- "Apis"

# recalculate original Apis visit duration on pollen+nectar 
me_ApisPrev_ApisDur_PN$x # scaled Apis visit duration on pollen+nectar 
ApisDurA_mean <- mean(log(data_Apis$APIS_visitdur5+1)) # mean of original apis visit duration
ApisDurA_sd <- sd(log(data_Apis$APIS_visitdur5+1)) # sd of original apis visit duration

me_ApisPrev_ApisDur_PN$ApisDurPN <- t((t(me_ApisPrev_ApisDur_PN$x) * ApisDurA_sd) + ApisDurA_mean)

# Nosema in Bombus
me_BombPrev_ApisDur_PN <- ggpredict(fit_Bombus_visitdur5, c("APIS_visitdur5_z"))
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
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "solid"), guide = "none") +
  labs(color = "Species", tag = "c", x = bquote(atop("Avg. Duration of Honeybee", "Pollen + Nectar Visits (seconds/visit)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
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





#################################
#### Bombus visit number figure   (Figure 3 A)
#################################
# Nosema in Apis
me_ApisPrev_BombusNum <- ggpredict(fit_Apis_visits, c("BOMB_visits_z"))
plot(me_ApisPrev_BombusNum, add.data = T)
me_ApisPrev_BombusNum$Host_Species <- "Apis"

# recalculate original Bombus visit number 
me_ApisPrev_BombusNum$x # scaled Bombus visit number 
BombusNumA_mean <- mean(log(data_Apis$BOMB_visits+1)) # mean of original Bombus visit number from disease data set
BombusNumA_sd <- sd(log(data_Apis$BOMB_visits+1)) # sd of original Bombus visit number from disease data set

me_ApisPrev_BombusNum$BombusVisits <- t((t(me_ApisPrev_BombusNum$x) * BombusNumA_sd) + BombusNumA_mean)

# Nosema in Bombus
me_BombPrev_BombusNum <- ggpredict(fit_Bombus_visits, c("BOMB_visits_z"))
plot(me_BombPrev_BombusNum, add.data = T)
me_BombPrev_BombusNum$Host_Species <- "Bombus"

# recalculate original Bombus visit number 
me_BombPrev_BombusNum$x # scaled Bombus visit number 
BombusNumB_mean <- mean(log(data_Bombus$BOMB_visits+1)) # mean of original Bombus visit number from disease data set
BombusNumB_sd <- sd(log(data_Bombus$BOMB_visits+1)) # sd of original Bombus visit number from disease data set

me_BombPrev_BombusNum$BombusVisits <- t((t(me_BombPrev_BombusNum$x) * BombusNumB_sd) + BombusNumB_mean)


# combine by rows
me_BombusNum <- rbind(me_ApisPrev_BombusNum, me_BombPrev_BombusNum)


# plot of Bombus visit number vs nosema prevalence
prev_BombusNum <- ggplot(me_BombusNum, aes(x = BombusVisits, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "a", x = bquote(atop("Avg. Number of", "Bumblebee Visits (per 30 min)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  xlim(0,3) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_BombusNum)

#################################
#### Other visit number figure     (Figure 3 B)
#################################
# Nosema in Apis
me_ApisPrev_OtherNum <- ggpredict(fit_Apis_visits, c("Other_visits_z"))
plot(me_ApisPrev_OtherNum, add.data = T)
me_ApisPrev_OtherNum$Host_Species <- "Apis"

# recalculate original Other visit number 
me_ApisPrev_OtherNum$x # scaled Other visit number 
OtherNumA_mean <- mean(log(data_Apis$Other_visits+1)) # mean of original Other visit number from disease data set
OtherNumA_sd <- sd(log(data_Apis$Other_visits+1)) # sd of original Other visit number from disease data set

me_ApisPrev_OtherNum$OtherVisits <- t((t(me_ApisPrev_OtherNum$x) * OtherNumA_sd) + OtherNumA_mean)

# Nosema in Other
me_BombPrev_OtherNum <- ggpredict(fit_Bombus_visits, c("Other_visits_z"))
plot(me_BombPrev_OtherNum, add.data = T)
me_BombPrev_OtherNum$Host_Species <- "Bombus"

# recalculate original Other visit number 
me_BombPrev_OtherNum$x # scaled Other visit number 
OtherNumB_mean <- mean(log(data_Bombus$Other_visits+1)) # mean of original Other visit number from disease data set
OtherNumB_sd <- sd(log(data_Bombus$Other_visits+1)) # sd of original Other visit number from disease data set

me_BombPrev_OtherNum$OtherVisits <- t((t(me_BombPrev_OtherNum$x) * OtherNumB_sd) + OtherNumB_mean)


# combine by rows
me_OtherNum <- rbind(me_ApisPrev_OtherNum, me_BombPrev_OtherNum)


# plot of Other visit number vs nosema prevalence
prev_OtherNum <- ggplot(me_OtherNum, aes(x = OtherVisits, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "b", x = bquote(atop("Avg. Number of Other", "Pollinator Visits (per 30 min)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_OtherNum)



#################################
#### Bombus duration per visit pollen + nectar figure  (Figure 3 C)
#################################
# Nosema in Apis
me_ApisPrev_BombusDur_PN <- ggpredict(fit_Apis_visitdur5, c("BOMB_visitdur5_z"))
plot(me_ApisPrev_BombusDur_PN, add.data = T)
me_ApisPrev_BombusDur_PN$Host_Species <- "Apis"

# recalculate original Bombus duration per visit pollen + nectar 
me_ApisPrev_BombusDur_PN$x # scaled Bombus duration per visit pollen + nectar 
BombusDur5A_mean <- mean(log(data_Apis$BOMB_visitdur5+1)) # mean of original Bombus duration per visit pollen + nectar from disease data set
BombusDur5A_sd <- sd(log(data_Apis$BOMB_visitdur5+1)) # sd of original Bombus duration per visit pollen + nectar from disease data set

me_ApisPrev_BombusDur_PN$BombusVisitDur5 <- t((t(me_ApisPrev_BombusDur_PN$x) * BombusDur5A_sd) + BombusDur5A_mean)

# Nosema in Bombus
me_BombPrev_BombusDur_PN <- ggpredict(fit_Bombus_visitdur5, c("BOMB_visitdur5_z"))
plot(me_BombPrev_BombusDur_PN, add.data = T)
me_BombPrev_BombusDur_PN$Host_Species <- "Bombus"

# recalculate original Bombus duration per visit pollen + nectar 
me_BombPrev_BombusDur_PN$x # scaled Bombus duration per visit pollen + nectar 
BombusDur5B_mean <- mean(log(data_Bombus$BOMB_visitdur5+1)) # mean of original Bombus duration per visit pollen + nectar from disease data set
BombusDur5B_sd <- sd(log(data_Bombus$BOMB_visitdur5+1)) # sd of original Bombus duration per visit pollen + nectar from disease data set

me_BombPrev_BombusDur_PN$BombusVisitDur5 <- t((t(me_BombPrev_BombusDur_PN$x) * BombusDur5B_sd) + BombusDur5B_mean)


# combine by rows
me_BombusDur5 <- rbind(me_ApisPrev_BombusDur_PN, me_BombPrev_BombusDur_PN)


# plot of Bombus duration per visit pollen + nectar vs nosema prevalence
prev_BombusDur5 <- ggplot(me_BombusDur5, aes(x = BombusVisitDur5, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "c", x = bquote(atop("Avg. Duration of Bumblebee", "Pollen+Nectar Visits (second/visit)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  xlim(0, 4) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_BombusDur5)

#################################
#### Other duration per visit pollen + nectar figure  (Figure 3 D)
#################################
# Nosema in Apis
me_ApisPrev_OtherDur_PN <- ggpredict(fit_Apis_visitdur5, c("Other_visitdur5_z"))
plot(me_ApisPrev_OtherDur_PN, add.data = T)
me_ApisPrev_OtherDur_PN$Host_Species <- "Apis"

# recalculate original Other duration per visit pollen + nectar 
me_ApisPrev_OtherDur_PN$x # scaled Other duration per visit pollen + nectar 
OtherDur5A_mean <- mean(log(data_Apis$Other_visitdur5+1)) # mean of original Other duration per visit pollen + nectar from disease data set
OtherDur5A_sd <- sd(log(data_Apis$Other_visitdur5+1)) # sd of original Other duration per visit pollen + nectar from disease data set

me_ApisPrev_OtherDur_PN$OtherVisitDur5 <- t((t(me_ApisPrev_OtherDur_PN$x) * OtherDur5A_sd) + OtherDur5A_mean)

# Nosema in Other
me_BombPrev_OtherDur_PN <- ggpredict(fit_Bombus_visitdur5, c("Other_visitdur5_z"))
plot(me_BombPrev_OtherDur_PN, add.data = T)
me_BombPrev_OtherDur_PN$Host_Species <- "Bombus"

# recalculate original Other duration per visit pollen + nectar 
me_BombPrev_OtherDur_PN$x # scaled Other duration per visit pollen + nectar 
OtherDur5B_mean <- mean(log(data_Bombus$Other_visitdur5+1)) # mean of original Other duration per visit pollen + nectar from disease data set
OtherDur5B_sd <- sd(log(data_Bombus$Other_visitdur5+1)) # sd of original Other duration per visit pollen + nectar from disease data set

me_BombPrev_OtherDur_PN$OtherVisitDur5 <- t((t(me_BombPrev_OtherDur_PN$x) * OtherDur5B_sd) + OtherDur5B_mean)


# combine by rows
me_OtherDur5 <- rbind(me_ApisPrev_OtherDur_PN, me_BombPrev_OtherDur_PN)


# plot of Other duration per visit pollen + nectar vs nosema prevalence
prev_OtherDur5 <- ggplot(me_OtherDur5, aes(x = OtherVisitDur5, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "d", x = bquote(atop("Avg. Duration of Other Pollinator", "Pollen+Nectar Visits (second/visit)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  xlim(-0.1, 4) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_OtherDur5)


###############################
### Full figure 3 of Bombus and Other pollinator visit number and duration per vistit on pollen+nectar
###############################
# Grid of both isit rate and sum duration plots
prev_fig3 <- grid_arrange_shared_legend(prev_BombusNum, prev_OtherNum, prev_BombusDur5, prev_OtherDur5, nrow=1, ncol = 4, position = "right")
ggsave(here("figures/Fig3.tiff"), plot = prev_fig3, dpi = 600, width = 10, height = 3, units = "in", compression="lzw")

# color blind check
cvdPlot(prev_fig3)







##################################################################
##### Additional Figures not included in the manuscript




###########################################
#### Apis visit duration on PETALS figure
###########################################
# Nosema in Apis
me_ApisPrev_ApisDur2 <- ggpredict(fit_Apis_visitdur2, c("APIS_visitdur2_z"))
plot(me_ApisPrev_ApisDur2, add.data = T)
me_ApisPrev_ApisDur2$Host_Species <- "Apis"

# recalculate original Apis visit duration of behavior 2 
me_ApisPrev_ApisDur2$x # scaled Apis duration of behavior 2 
ApisDur2A_mean <- mean(log(data_Apis$APIS_visitdur2+1)) # mean of original Apis duration of behavior 2 
ApisDur2A_sd <- sd(log(data_Apis$APIS_visitdur2+1)) # sd of original Apis duration of behavior 2 

me_ApisPrev_ApisDur2$ApisDur2 <- t((t(me_ApisPrev_ApisDur2$x) * ApisDur2A_sd) + ApisDur2A_mean)

# Nosema in Bombus
me_BombPrev_ApisDur2 <- ggpredict(fit_Bombus_visitdur2, c("APIS_visitdur2_z"))
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
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "a", x = bquote(atop("Avg. Duration of Honeybee", "Petal visits (seconds/visit)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisDur2)

###########################################
## Bombus visit duration on PETALS figure
###########################################
# Nosema in Apis
me_ApisPrev_BombusDur2 <- ggpredict(fit_Apis_visitdur2, c("BOMB_visitdur2_z"))
plot(me_ApisPrev_BombusDur2, add.data = T)
me_ApisPrev_BombusDur2$Host_Species <- "Apis"

# recalculate original Bombus visit duration of behavior 2  
me_ApisPrev_BombusDur2$x # scaled Bombus visit duration of behavior 2 
BombusDurA_mean <- mean(log(data_Apis$BOMB_visitdur2+1)) # mean of original Bombus visit duration of behavior 2  
BombusDurA_sd <- sd(log(data_Apis$BOMB_visitdur2+1)) # sd of original Bombus visit duration of behavior 2  

me_ApisPrev_BombusDur2$BombusDur2 <- t((t(me_ApisPrev_BombusDur2$x) * BombusDurA_sd) + BombusDurA_mean)

# Nosema in Bombus
me_BombPrev_BombusDur2 <- ggpredict(fit_Bombus_visitdur2, c("BOMB_visitdur2_z"))
plot(me_BombPrev_BombusDur2, add.data = T)
me_BombPrev_BombusDur2$Host_Species <- "Bombus"

# recalculate original Bombus visit duration of behavior 2
me_BombPrev_BombusDur2$x # scaled Bombus visit duration of behavior 2
BombusDurB_mean <- mean(log(data_Bombus$BOMB_visitdur2+1)) # mean of original Bombus visit duration of behavior 2 from data set
BombusDurB_sd <- sd(log(data_Bombus$BOMB_visitdur2+1)) # sd of original Bombus visit duration of behavior 2 from data set

me_BombPrev_BombusDur2$BombusDur2 <- t((t(me_BombPrev_BombusDur2$x) * BombusDurB_sd) + BombusDurB_mean)


# combine by rows
me_BombusDur2 <- rbind(me_ApisPrev_BombusDur2, me_BombPrev_BombusDur2)


# plot of Bombus visit duration total vs nosema prevalence
prev_BombusDur2 <- ggplot(me_BombusDur2, aes(x = BombusDur2, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "b", x = bquote(atop("Avg. Duration of Bumblebee", "Petal Visits (seconds/visit)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_BombusDur2)



# Grid of both Apis and bombus duration on PETALS
prev_figS4 <- grid_arrange_shared_legend(prev_ApisDur2, prev_BombusDur2, nrow=1, ncol = 2, position = "right")
ggsave("figures/FigS4.tiff", plot = prev_figS4, dpi = 600, width = 6.5, height = 3, units = "in", compression="lzw")

# color blind check
cvdPlot(prev_figS2)



#################################
#### Apis visit rate figure
#################################
# Nosema in Apis
me_ApisPrev_ApisRate <- ggpredict(fit_Apis_rate, c("APIS_rate_z"))
plot(me_ApisPrev_ApisRate, add.data = T)
me_ApisPrev_ApisRate$Host_Species <- "Apis"

# recalculate original Apis visit rate 
me_ApisPrev_ApisRate$x # scaled Apis visit rate
ApisRateA_mean <- mean(log(data_Apis$APIS_rate+1)) # mean of original Apis visit rate from data set
ApisRateA_sd <- sd(log(data_Apis$APIS_rate+1)) # sd of original Apis visit rate from data set

me_ApisPrev_ApisRate$ApisRate <- t((t(me_ApisPrev_ApisRate$x) * ApisRateA_sd) + ApisRateA_mean)

# Nosema in Bombus
me_BombPrev_ApisRate <- ggpredict(fit_Bombus_rate, c("APIS_rate_z"))
plot(me_BombPrev_ApisRate, add.data = T)
me_BombPrev_ApisRate$Host_Species <- "Bombus"

# recalculate original Apis visit rate 
me_BombPrev_ApisRate$x # scaled Apis visit rate
ApisRateB_mean <- mean(log(data_Bombus$APIS_rate+1)) # mean of original Apis visit rate from data set
ApisRateB_sd <- sd(log(data_Bombus$APIS_rate+1)) # sd of original Apis visit rate from data set

me_BombPrev_ApisRate$ApisRate <- t((t(me_BombPrev_ApisRate$x) * ApisRateB_sd) + ApisRateB_mean)


# combine by rows
me_ApisRate <- rbind(me_ApisPrev_ApisRate, me_BombPrev_ApisRate)


# plot of Apis visit rate vs nosema prevalence
prev_ApisRate <- ggplot(me_ApisRate, aes(x = ApisRate, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "solid"), guide = "none") +
  labs(color = "Species", tag = "a", x = bquote(atop("Avg. Rate of", "Honeybee Visits (per min)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisRate)
#ggsave("VirusPrev_Richness_predict.png", plot = virus_rich, dpi = 300, width = 18, height = 8, units = "cm")



######################################
## Apis sum duration (total) per 30 min figure
#####################################
# Nosema in Apis
me_ApisPrev_ApisSumDur <- ggpredict(fit_Apis_dur, c("APIS_dur_z"))
plot(me_ApisPrev_ApisSumDur, add.data = T)
me_ApisPrev_ApisSumDur$Host_Species <- "Apis"

# recalculate original Apis sum duration per 30 min
me_ApisPrev_ApisSumDur$x # scaled Apis sum duration
ApisSumDurA_mean <- mean(log(data_Apis$APIS_dur+1)) # mean of original apis sum duration
ApisSumDurA_sd <- sd(log(data_Apis$APIS_dur+1)) # sd of original apis sum duration

me_ApisPrev_ApisSumDur$ApisSumDur <- t((t(me_ApisPrev_ApisSumDur$x) * ApisSumDurA_sd) + ApisSumDurA_mean)

# Nosema in Bombus
me_BombPrev_ApisSumDur <- ggpredict(fit_Bombus_dur, c("APIS_dur_z"))
plot(me_BombPrev_ApisSumDur, add.data = T)
me_BombPrev_ApisSumDur$Host_Species <- "Bombus"

# recalculate original Apis sum duration 
me_BombPrev_ApisSumDur$x # scaled Apis sum duration
ApisSumDurB_mean <- mean(log(data_Bombus$APIS_dur+1)) # mean of original apis sum duration 
ApisSumDurB_sd <- sd(log(data_Bombus$APIS_dur+1)) # sd of original apis sum duration

me_BombPrev_ApisSumDur$ApisSumDur <- t((t(me_BombPrev_ApisSumDur$x) * ApisSumDurB_sd) + ApisSumDurB_mean)


# combine by rows
me_ApisSumDur <- rbind(me_ApisPrev_ApisSumDur, me_BombPrev_ApisSumDur)


# plot of Apis visit duration total vs nosema prevalence
prev_ApisSumDur <- ggplot(me_ApisSumDur, aes(x = ApisSumDur, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "dashed"), guide = "none") +
  labs(color = "Species", tag = "b", x = bquote(atop("Avg. Sum Duration of", "Honeybee Visits (seconds)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisSumDur)



###########################################
## Bombus Sum duration on PETALS per 30 min figure
###########################################
# Nosema in Apis
me_ApisPrev_BombusSumDur2 <- ggpredict(fit_Apis_dur2, c("BOMB_dur2_z"))  # generated figure based on initial, non-spatial model because ggpredict does not work on spaMM models, estimates are the same though
plot(me_ApisPrev_BombusSumDur2, add.data = T)
me_ApisPrev_BombusSumDur2$Host_Species <- "Apis"

# recalculate original Bombus sum duration of behavior 2  
me_ApisPrev_BombusSumDur2$x # scaled Bombus sum duration of behavior 2 
BombusSumDurA_mean <- mean(log(data_Apis$BOMB_dur2+1)) # mean of original Bombus sum duration of behavior 2  
BombusSumDurA_sd <- sd(log(data_Apis$BOMB_dur2+1)) # sd of original Bombus sum duration of behavior 2  

me_ApisPrev_BombusSumDur2$BombusSumDur2 <- t((t(me_ApisPrev_BombusSumDur2$x) * BombusSumDurA_sd) + BombusSumDurA_mean)

# Nosema in Bombus
me_BombPrev_BombusSumDur2 <- ggpredict(fit_Bombus_dur2, c("BOMB_dur2_z"))
plot(me_BombPrev_BombusSumDur2, add.data = T)
me_BombPrev_BombusSumDur2$Host_Species <- "Bombus"

# recalculate original Bombus sum duration on petals
me_BombPrev_BombusSumDur2$x # scaled Bombus sum duration on petals
BombusSumDurB_mean <- mean(log(data_Bombus$BOMB_dur2+1)) # mean of original Bombus sum duration on petals from disease data set
BombusSumDurB_sd <- sd(log(data_Bombus$BOMB_dur2+1)) # sd of original Bombus sum duration on petals from disease data set

me_BombPrev_BombusSumDur2$BombusSumDur2 <- t((t(me_BombPrev_BombusSumDur2$x) * BombusSumDurB_sd) + BombusSumDurB_mean)


# combine by rows
me_BombusSumDur2 <- rbind(me_ApisPrev_BombusSumDur2, me_BombPrev_BombusSumDur2)


# plot of Bombus sum duration total vs nosema prevalence
prev_BombusSumDur2 <- ggplot(me_BombusSumDur2, aes(x = BombusSumDur2, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote(italic("Apis")), bquote(italic("Bombus")))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("solid", "solid"), guide = "none") +
  labs(color = "Species", tag = "c", x = bquote(atop("Avg. Sum Duration of", "Bumblebee Petal Visits (seconds)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  xlim(0, 3) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_BombusSumDur2)


###############################################
## Apis Sum duration on pollen+nectar figure
###############################################
# Nosema in Apis
me_ApisPrev_ApisSumDur_PN <- ggpredict(fit_Apis_dur5, c("APIS_dur5_z"))
plot(me_ApisPrev_ApisSumDur_PN, add.data = T)
me_ApisPrev_ApisSumDur_PN$Host_Species <- "Apis"

# recalculate original Apis sum duration on pollen+nectar 
me_ApisPrev_ApisSumDur_PN$x # scaled Apis sum duration on pollen+nectar 
ApisSumDurA_mean <- mean(log(data_Apis$APIS_dur5+1)) # mean of original Apis sum duration on pollen+nectar 
ApisSumDurA_sd <- sd(log(data_Apis$APIS_dur5+1)) # sd of original Apis sum duration on pollen+nectar 

me_ApisPrev_ApisSumDur_PN$ApisSumDurPN <- t((t(me_ApisPrev_ApisSumDur_PN$x) * ApisSumDurA_sd) + ApisSumDurA_mean)

# Nosema in Bombus
me_BombPrev_ApisSumDur_PN <- ggpredict(fit_Bombus_dur5, c("APIS_dur5_z"))
plot(me_BombPrev_ApisSumDur_PN, add.data = T)
me_BombPrev_ApisSumDur_PN$Host_Species <- "Bombus"

# recalculate original Apis visit duration on pollen+nectar 
me_BombPrev_ApisSumDur_PN$x # scaled Apis visit duration on pollen+nectar 
ApisSumDurB_mean <- mean(log(data_Bombus$APIS_dur5+1)) # mean of original Apis sum duration on pollen+nectar 
ApisSumDurB_sd <- sd(log(data_Bombus$APIS_dur5+1)) # sd of original Apis sum duration on pollen+nectar 

me_BombPrev_ApisSumDur_PN$ApisSumDurPN <- t((t(me_BombPrev_ApisSumDur_PN$x) * ApisSumDurB_sd) + ApisSumDurB_mean)


# combine by rows
me_ApisSumDur_PN <- rbind(me_ApisPrev_ApisSumDur_PN, me_BombPrev_ApisSumDur_PN)


# plot of Apis sum duration Pollen+Nectar vs nosema prevalence
prev_ApisSumDur_PN <- ggplot(me_ApisSumDur_PN, aes(x = ApisSumDurPN, y = predicted)) +
  scale_fill_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  scale_color_manual(values = Tol_muted2, name = "Species", labels = c(bquote("Honeybees"), bquote("Bumblebees"))) +
  geom_line(aes(color = Host_Species, linetype = Host_Species), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Host_Species), alpha = 0.2) +
  scale_linetype_manual(values=c("dashed", "solid"), guide = "none") +
  labs(color = "Species", tag = "d", x = bquote(atop("Avg. Sum Duration of Honeybee", "Pollen + Nectar Visits (seconds)")), y = bquote(italic(V.) ~ italic(ceranae) ~ " Prevalence")) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8))
print(prev_ApisSumDur_PN)




###############################
### Full Figure of Visit rate and sum duration plots [NOT IN MANUSCRIPT]
###############################
# Grid of both isit rate and sum duration plots
prev_fig5 <- grid_arrange_shared_legend(prev_ApisRate, prev_ApisSumDur, prev_BombusSumDur2, prev_ApisSumDur_PN, nrow=1, ncol = 4, position = "right")
ggsave(here("figures/Fig5.tiff"), plot = prev_fig5, dpi = 600, width = 10, height = 3.5, units = "in", compression="lzw")

# color blind check
cvdPlot(prev_fig4)




