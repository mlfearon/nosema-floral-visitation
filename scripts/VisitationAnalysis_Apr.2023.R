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
library(ggeffects)
library(ggplot2)
library(ggpubr)
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
head(avgs_bySpp)

## histograms
hist(visitation_bySpp$visits)
hist(visitation_bySpp$rate)
hist(visitation_bySpp$dur)


# color palette for figures
# color palette
Tol_muted <- c('#882255','#44AA99', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#DDDDDD', '#999933')
Tol_muted3 <- c('#882255','#44AA99','#88CCEE')


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


###############################
# Tests of differences in visitation metrics between species (Apis, Bombus, Other)
# For each type of visitation metric, compare whether we get the same results from the per flower (as collected) and averaged per site/visit and per site
##############################

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
visitnum_contrasts$contrasts
visitnum_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0003, "***",
  "APIS",     "Other", 0.1467, "NS",
  "BOMB",     "Other", 0.6152, "NS")

me_visitnum <- ggpredict(visitnum_mod2, "Genus")
plot(me_visitnum, add.data = T) +
  scale_y_log10()


# Visit number per 30 min by genus figure
visitnum_plot <- ggplot(me_visitnum, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp_noVisitNumOutlier, aes(x= Genus, y = visits, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators" )) +
  scale_y_log10(breaks = c(1, 3, 10, 30, 100)) +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  stat_pvalue_manual(visitnum_sig, y.position = 1.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "a", x = "Species", y = "Number of Visits (per 30 min)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
visitnum_plot


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
totdur_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.4282, "NS",
  "APIS",     "Other", 0.2864, "NS",
  "BOMB",     "Other", 0.8332, "NS")
me_totdur <- ggpredict(totdur_mod2, "Genus")
plot(me_totdur, add.data = T) +
  scale_y_log10()

# Total sum duration per 30 min by genus figure
totdur_plot <- ggplot(me_totdur, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp_noDurOutlier, aes(x= Genus, y = dur, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  annotate(geom = "text", x = 1.8, y = 1500, label = bquote(paste("Species, ", chi^2, " = 2.3, p = NS")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  #stat_pvalue_manual(totdur_sig, y.position = 3.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 3) +
  labs(color = "Species", tag = "f", x = "Species", y = "Sum Duration of Visits (seconds)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totdur_plot


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
totdur2_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0126, "*",
  "APIS",     "Other", 0.2748, "NS",
  "BOMB",     "Other", 0.0001, "***")
me_totdur2 <- ggpredict(totdur2_mod2, "Genus")
plot(me_totdur2, add.data = T) +
  scale_y_log10()

# Total sum duration on petals per 30 min by genus figure
totdur2_plot <- ggplot(me_totdur2, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = dur2, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  #annotate(geom = "text", x = 1, y = 1400, label = bquote(paste("Species, ", chi^2, " = 2.3, p = NS")), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  stat_pvalue_manual(totdur2_sig, y.position = 2.8, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "g", x = "Species", y = "Sum Duration of Petal Visits (seconds)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totdur2_plot


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
totdur5_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.3249, "NS",
  "APIS",     "Other", 0.9848, "NS",
  "BOMB",     "Other", 0.1401, "NS")
me_totdur5 <- ggpredict(totdur5_mod2, "Genus")
plot(me_totdur5, add.data = T) +
  scale_y_log10()


# Total sum duration on pollen+nectar per 30 min by genus figure
totdur5_plot <- ggplot(me_totdur5, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = dur5, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  annotate(geom = "text", x = 1.8, y = 1500, label = bquote(paste("Species, ", chi^2, " = 4.9, p = 0.08")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  #stat_pvalue_manual(totdur5_sig, y.position = 3.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 3) +
  labs(color = "Species", tag = "h", x = "Species", y = "Sum Duration of Visits (seconds)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totdur5_plot



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
visitrate_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0003, "***",
  "APIS",     "Other", 0.1242, "NS",
  "BOMB",     "Other", 0.9465, "NS")
me_visitrate <- ggpredict(visitrate_mod2, "Genus")
plot(me_visitrate) +  # need to add in rate raw data points instead of visits data points
  scale_y_log10()


# Visit rate (visits per min) by genus figure
visitrate_plot <- ggplot(me_visitrate, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = rate, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators" )) +
  scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10)) +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  stat_pvalue_manual(visitnum_sig, y.position = 0.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "e", x = "Species", y = "Rate of Visits (per min)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
visitrate_plot




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
totvisitdur_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.5539, "NS",
  "APIS",     "Other", 0.8760, "NS",
  "BOMB",     "Other", 0.0738, "NS")
me_totvisitdur <- ggpredict(totvisitdur_mod2, "Genus")
plot(me_totvisitdur, add.data = T) +  # need to add in duration per visit raw data points instead of total dur data points
  scale_y_log10()


# Total duration per visit by genus figure
totvisitdur_plot <- ggplot(me_totvisitdur, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = visitdur, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  annotate(geom = "text", x = 1.8, y = 300, label = bquote(paste("Species, ", chi^2, " = 5.0, p = 0.08")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  #stat_pvalue_manual(totvisitdur_sig, y.position = 3.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 3) +
  labs(color = "Species", tag = "b", x = "Species", y = "Duration per Visit (seconds/visit)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totvisitdur_plot


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
totvisitdur2_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0001, "***",
  "APIS",     "Other", 0.8722, "NS",
  "BOMB",     "Other", 0.0001, "***")
me_totvisitdur2 <- ggpredict(totvisitdur2_mod2, "Genus")
plot(me_totvisitdur2) +    # need to add in duration per visit to petals raw data points instead of total dur data points
  scale_y_log10()

# Total duration per visit on petals by genus figure
totvisitdur2_plot <- ggplot(me_totvisitdur2, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = visitdur2, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  #annotate(geom = "text", x = 1.8, y = 300, label = bquote(paste("Species, ", chi^2, " = 5.0, p = 0.08")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  stat_pvalue_manual(totvisitdur2_sig, y.position = 1.9, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "c", x = "Species", y = "Duration on Petals per Visit (seconds/visit)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totvisitdur2_plot


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
totvisitdur5_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.9722, "NS",
  "APIS",     "Other", 0.4824, "NS",
  "BOMB",     "Other", 0.1458, "NS")
me_totvisitdur5 <- ggpredict(totvisitdur5_mod2, "Genus")
plot(me_totvisitdur5, add.data = T) +   # need to add in duration per visit to pollen + nectar raw data points instead of total dur data points
  scale_y_log10()

# Total duration per visit on pollen + nectar by genus figure
totvisitdur5_plot <- ggplot(me_totvisitdur5, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = visitdur5, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  annotate(geom = "text", x = 1.8, y = 300, label = bquote(paste("Species, ", chi^2, " = 3.6, p = NS")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Other")) +
  #stat_pvalue_manual(totvisitdur5_sig, y.position = 1.9, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "d", x = "Species", y = "Duration on Pollen + Nectar per Visit (seconds/visit)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totvisitdur5_plot



#############################
# Combined figure of all visitation behaviors by genus
############################

# Grid of both Apis visit number and duration plots
visitation_fig <- grid_arrange_shared_legend(visitnum_plot, totvisitdur_plot, totvisitdur2_plot, totvisitdur5_plot, visitrate_plot, totdur_plot, totdur2_plot, totdur5_plot, nrow=2, ncol = 4, position = "right")
ggsave(here("figures/Fig6.tiff"), plot = visitation_fig, dpi = 600, width = 11, height = 7, units = "in", compression="lzw")


