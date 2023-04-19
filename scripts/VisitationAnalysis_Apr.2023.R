# Code for the analyses for "Honeybee visitation behavior on shared floral resources 
# increases Vairimorpha ceranae prevalence in bumblebees"

# Manuscript submitted to: 


# Vairimorpha (=Nosema) ceranae prevalence in Apis and Bombus analysis

# This script includes the full analyses and figures for how pollinator visitation number and duration of visits to flowers
# impact V. ceranae prevalence in Apis mellifera and Bombus spp. 


# Written by: Michelle Fearon
# Last updated: 19 April 2023


# Import libraries needed 

library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(EnvStats)
library(here)



# set the path to the script relative to the project root directory
here::i_am("scripts/VisitationAnalysis_Apr.2023.R")


###############################
### LOAD THE DATA FOR ANALYSIS
###############################
visitation_bySpp <- read.csv(here("data/Visitation_bySpp.csv"), stringsAsFactors = F)

head(visitation_bySpp)
summary(visitation_bySpp)


avgs_bySpp <- read.csv(here("data/VisitationAvgs_bySpp.csv"), stringsAsFactors = F)


## histograms
hist(visitation_bySpp$visits)
hist(visitation_bySpp$rate)
hist(visitation_bySpp$dur)

###############################
# Tests of differences in visitation metrics between species (Apis, Bombus, Other)
# For each type of visitation metric, compare whether we get the same results from the per flower (as collected) and averaged per site/visit and per site
##############################


# color palette for figures
# color palette
Tol_muted <- c('#882255','#44AA99', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#DDDDDD', '#999933')
Tol_muted3 <- c('#882255','#44AA99','#88CCEE')


#### Visit number per 30 min

# remove one outlier from the Other pollinators that has 136 visits to one flower (many combined species)
visitation_bySpp_noVisitNumOutlier <- filter(visitation_bySpp, visits < 100)

visitnum_mod1 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp_noVisitNumOutlier)
visitnum_mod2 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp_noVisitNumOutlier)
AIC(visitnum_mod1, visitnum_mod2)  # negbinomial model has much lower AIC
summary(visitnum_mod2)
Anova(visitnum_mod2)
testDispersion(visitnum_mod2)
testZeroInflation(visitnum_mod2)
visitnum_simResid1 <- simulateResiduals(fittedModel = visitnum_mod2)
plot(visitnum_simResid1) 
visitnum_contrasts <- emmeans(visitnum_mod2, spec = pairwise ~ Genus, type = "response")
visitnum_contrasts    #  Apis is sig lower than Bombus, but not with Other. Bombus and Other not sig different
me_visitnum <- ggpredict(visitnum_mod2, "Genus")
plot(me_visitnum, add.data = T) +
  scale_y_log10()

#  Visit number average by visit to each site produces similar results
# response variable is challenging b/c it has a count-link distribution but it is an average and therefore not an integer. Better to use the models above.
visitnum_mod <- glmmTMB(visits ~ Genus + (1|Site/Visit), family = nbinom2(), data = avgs_bySpp)
summary(visitnum_mod2)


#### Total sum duration per 30 min

# remove one outlier from the Other pollinators that has 2800+ seconds of visits to one flower (many combined species)
visitation_bySpp_noDurOutlier <- filter(visitation_bySpp, dur < 2500)

totdur_mod1 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp_noDurOutlier)
totdur_mod2 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp_noDurOutlier)
AIC(totdur_mod1, totdur_mod2)  # negbinomial model has much lower AIC
summary(totdur_mod2)
Anova(totdur_mod2)
testDispersion(totdur_mod2)
testZeroInflation(totdur_mod2)
totdur_simResid1 <- simulateResiduals(fittedModel = totdur_mod2)
plot(totdur_simResid1) 
totdur_contrasts <- emmeans(totdur_mod2, spec = pairwise ~ Genus, type = "response")
totdur_contrasts     # no sig diff based on genus
me_totdur <- ggpredict(totdur_mod2, "Genus")
plot(me_totdur, add.data = T) +
  scale_y_log10()



#### Total sum duration on petals per 30 min

totdur2_mod1 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp)
totdur2_mod2 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp)
AIC(totdur2_mod1, totdur2_mod2)  # poisson model doesn't converge, but the negbinomial model does
summary(totdur2_mod2)
Anova(totdur2_mod2)
testDispersion(totdur2_mod2)
testZeroInflation(totdur2_mod2)
totdur2_simResid1 <- simulateResiduals(fittedModel = totdur2_mod2)
plot(totdur2_simResid1) 
totdur2_contrasts <- emmeans(totdur2_mod2, spec = pairwise ~ Genus, type = "response")
totdur2_contrasts   # Bombus is sig lower than Apis or Other, Apis and Other not sig different
me_totdur2 <- ggpredict(totdur2_mod2, "Genus")
plot(me_totdur2, add.data = T) +
  scale_y_log10()



#### Total sum duration on pollen+nectar per 30 min

totdur5_mod1 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp)
totdur5_mod2 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp)
AIC(totdur5_mod1, totdur5_mod2)
summary(totdur5_mod2)
Anova(totdur5_mod2)
testDispersion(totdur5_mod2)
testZeroInflation(totdur5_mod2)
totdur5_simResid1 <- simulateResiduals(fittedModel = totdur5_mod2)
plot(totdur5_simResid1) 
totdur5_contrasts <- emmeans(totdur5_mod2, spec = pairwise ~ Genus, type = "response")
totdur5_contrasts    # No sig diff based on genus
me_totdur5 <- ggpredict(totdur5_mod2, "Genus")
plot(me_totdur5, add.data = T) +
  scale_y_log10()



#### Visit rate (visit per min)

# remove one outlier from the Other pollinators that has 136 visits to one flower (many combined species)
visitation_bySpp_noVisitNumOutlier <- filter(visitation_bySpp, visits < 100)

visitrate_mod1 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(totdur_min), family = poisson, data = visitation_bySpp_noVisitNumOutlier)
visitrate_mod2 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(totdur_min), family = nbinom2(), data = visitation_bySpp_noVisitNumOutlier)
AIC(visitrate_mod1, visitrate_mod2)  # negbinomial model has much lower AIC
summary(visitrate_mod2)
Anova(visitrate_mod2)
testDispersion(visitrate_mod2)
testZeroInflation(visitrate_mod2)
visitrate_simResid1 <- simulateResiduals(fittedModel = visitrate_mod2)
plot(visitrate_simResid1) 
visitrate_contrasts <- emmeans(visitrate_mod2, spec = pairwise ~ Genus, type = "response", offset = T)
visitrate_contrasts     #  Apis is sig lower than Bombus, but not with Other. Bombus and Other not sig different
me_visitrate <- ggpredict(visitrate_mod2, "Genus")
plot(me_visitrate, add.data = T) +  # need to add in rate raw data points instead of visits data points
  scale_y_log10()



#### Total duration per visit

# remove one outlier from the Other pollinators that has 2800+ seconds of visits to one flower (many combined species)
visitation_bySpp_noDurOutlier <- filter(visitation_bySpp, dur < 2500)

totvisitdur_mod1 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = poisson, data = visitation_bySpp)
totvisitdur_mod2 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = nbinom2(), data = visitation_bySpp)
AIC(totvisitdur_mod1, totvisitdur_mod2)  # negbinomial model has much lower AIC
summary(totvisitdur_mod2)
Anova(totvisitdur_mod2)
testDispersion(totvisitdur_mod2)
testZeroInflation(totvisitdur_mod2)
totvisitdur_simResid1 <- simulateResiduals(fittedModel = totvisitdur_mod2)
plot(totvisitdur_simResid1)   # residual v predicted plot has some problems; quantile deviations detected
totvisitdur_contrasts <- emmeans(totvisitdur_mod2, spec = pairwise ~ Genus, type = "response", offset = T)
totvisitdur_contrasts    # No sig dif based on genus
me_totvisitdur <- ggpredict(totvisitdur_mod2, "Genus")
plot(me_totvisitdur, add.data = T) +  # need to add in duration per visit raw data points instead of total dur data points
  scale_y_log10()



#### Total duration per visit on petals per 30 min

totvisitdur2_mod1 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = poisson, data = visitation_bySpp)
totvisitdur2_mod2 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = nbinom2(), data = visitation_bySpp)
AIC(totvisitdur2_mod2)  # poisson model doesn't converge, but the negbinomial model does
summary(totvisitdur2_mod2)
Anova(totvisitdur2_mod2)
testDispersion(totvisitdur2_mod2)
testZeroInflation(totvisitdur2_mod2)
totvisitdur2_simResid1 <- simulateResiduals(fittedModel = totvisitdur2_mod2)
plot(totvisitdur2_simResid1)  # residual v predicted plot has some problems
totvisitdur2_contrasts <- emmeans(totvisitdur2_mod2, spec = pairwise ~ Genus, type = "response", offset = T)
totvisitdur2_contrasts     # Bombus is sig lower than Apis or Other, Apis and Other not sig different
me_totvisitdur2 <- ggpredict(totvisitdur2_mod2, "Genus")
plot(me_totvisitdur2) +    # need to add in duration per visit to petals raw data points instead of total dur data points
  scale_y_log10()



#### Total duration per visit on pollen+nectar per 30 min

totvisitdur5_mod1 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = poisson, data = visitation_bySpp)
totvisitdur5_mod2 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = nbinom2(), data = visitation_bySpp)
AIC(totvisitdur5_mod1, totvisitdur5_mod2)
summary(totvisitdur5_mod2)
Anova(totvisitdur5_mod2)
testDispersion(totvisitdur5_mod2)
testZeroInflation(totvisitdur5_mod2)
totvisitdur5_simResid1 <- simulateResiduals(fittedModel = totvisitdur5_mod2)
plot(totvisitdur5_simResid1)  # residual v predicted plot has some problems; quantile deviations detected
totvisitdur5_contrasts <- emmeans(totvisitdur5_mod2, spec = pairwise ~ Genus, type = "response", offset  = T)
totvisitdur5_contrasts    # No sig differences between genus
me_totvisitdur5 <- ggpredict(totvisitdur5_mod2, "Genus")
plot(me_totvisitdur5, add.data = T) +
  scale_y_log10()
