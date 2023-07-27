# Code for the analyses for "Honeybee visitation behavior on shared flowers
# increases Vairimorpha ceranae prevalence in bumblebees"

# Manuscript submitted to: Ecology and Evolution


# Vairimorpha (=Nosema) ceranae prevalence in Apis and Bombus analysis

# This script includes the full analyses and figures for how pollinator visitation number and duration per visit to flowers
# vary based on pollinator group (honeybees, bumblebees, and other pollinators)


# Written by: Michelle Fearon
# Last updated: 27 July 2023


# Import libraries needed 

library(lme4)
library(lmerTest)
library(glmmTMB)
library(car)
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
visitation_bySpp$Genus <- factor(visitation_bySpp$Genus, levels = c("APIS", "BOMB", "PEPO", "Other"))

## histograms
hist(visitation_bySpp$visits)
hist(visitation_bySpp$rate)
hist(visitation_bySpp$dur)


# color palette for figures
# color palette
Tol_muted <- c('#882255','#44AA99', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#DDDDDD', '#999933')
Tol_muted3 <- c('#882255','#44AA99','#332288', '#88CCEE')


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

#### Visit number per 30 min [INCLUDED IN THE MANUSCRIPT]

# remove one outlier from the Other pollinators that has 136 visits to one flower (many combined species)
visitation_bySpp_noVisitNumOutlier <- filter(visitation_bySpp, visits < 100)

range(visitation_bySpp_noVisitNumOutlier$visits)
range(visitation_bySpp$visits)
visitnum_mod1 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp_noVisitNumOutlier)
visitnum_mod2 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp_noVisitNumOutlier)  # use this model
AIC(visitnum_mod1, visitnum_mod2)  # negbinomial model has much lower AIC
summary(visitnum_mod2)
Anova(visitnum_mod2)

# tests of model assumptions and diagnostics
testDispersion(visitnum_mod2)
testZeroInflation(visitnum_mod2)

# spatial autocorrelation test
visitnum_resid <- simulateResiduals(visitnum_mod2)
visitnum_resid2 <- recalculateResiduals(visitnum_resid, group = visitation_bySpp_noVisitNumOutlier$Site) # group residuals by site
testSpatialAutocorrelation(visitnum_resid2, unique(visitation_bySpp_noVisitNumOutlier$Long), unique(visitation_bySpp_noVisitNumOutlier$Lat)) # not sig!

# post hoc test
visitnum_contrasts <- emmeans(visitnum_mod2, spec = pairwise ~ Genus, type = "response")
visitnum_contrasts    #  Apis is sig lower than Bombus, but not with Other. Bombus and Other not sig different
visitnum_contrasts$contrasts
visitnum_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0008, "***",
  "APIS",     "Other", 0.9906, "NS",
  "APIS",     "PEPO", 0.8316, "NS",
  "BOMB",     "Other", 0.0016, "**",
  "BOMB",     "PEPO", 0.1843, "NS",
  "Other",     "PEPO", 0.8667, "NS")

visitnum_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0008, "***",
  "BOMB",     "Other", 0.0016, "**")

me_visitnum <- ggpredict(visitnum_mod2, "Genus")
plot(me_visitnum, add.data = T) +
  scale_y_log10()



# Visit number per 30 min by genus figure
visitnum_plot <- ggplot(me_visitnum, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp_noVisitNumOutlier, aes(x= Genus, y = visits, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10(breaks = c(1, 3, 10, 30, 100)) +
  scale_x_discrete(labels = NULL) +
  stat_pvalue_manual(visitnum_sig, y.position = 1.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.8) +
  labs(color = "Species", tag = "a", x = NULL, y = "Number of Visits (per 30 min)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
visitnum_plot


#  Visit number average by visit to each site produces similar results
# response variable is challenging b/c it has a count-link distribution but it is an average and therefore not an integer. Better to use the models above.
visitnum_mod <- glmmTMB(visits ~ Genus + (1|Site/Visit), family = nbinom2(), data = avgs_bySpp)
summary(visitnum_mod2)




#### Total duration per visit [INCLUDED IN THE MANUSCRIPT]

totvisitdur_mod1 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = poisson, data = visitation_bySpp)
totvisitdur_mod2 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = nbinom2(), data = visitation_bySpp)
AIC(totvisitdur_mod1, totvisitdur_mod2)  # negbinomial model has much lower AIC
summary(totvisitdur_mod2)
Anova(totvisitdur_mod2)
# test of model diagnostics
testDispersion(totvisitdur_mod2)
testZeroInflation(totvisitdur_mod2)

# spatial autocorrelation test
totvisitdur_resid <- simulateResiduals(totvisitdur_mod2)
totvisitdur_resid2 <- recalculateResiduals(totvisitdur_resid, group = visitation_bySpp$Site) # group residuals by site
testSpatialAutocorrelation(totvisitdur_resid2, unique(visitation_bySpp$Long), unique(visitation_bySpp$Lat)) # not sig!

# post hoc tests
totvisitdur_contrasts <- emmeans(totvisitdur_mod2, spec = pairwise ~ Genus, type = "response", offset = T)
totvisitdur_contrasts    # No sig dif based on genus
totvisitdur_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.8311, "NS",
  "APIS",     "PEPO", 0.8321, "NS",
  "APIS",     "Other", 0.7975, "NS",
  "BOMB",     "PEPO", 0.9984, "NS",
  "BOMB",     "Other", 0.0721, "NS",
  "PEPO",     "Other", 0.1815, "NS")
me_totvisitdur <- ggpredict(totvisitdur_mod2, "Genus")
plot(me_totvisitdur, add.data = T) +  # need to add in duration per visit raw data points instead of total dur data points
  scale_y_log10()


# Total duration per visit by genus figure
totvisitdur_plot <- ggplot(me_totvisitdur, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = visitdur, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  annotate(geom = "text", x = 1.8, y = 300, label = bquote(paste("Species, ", chi^2, " = 7.4, p = 0.06")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = NULL) +
  #stat_pvalue_manual(totvisitdur_sig, y.position = 3.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 3) +
  labs(color = "Species", tag = "b", x = NULL, y = "Duration per Visit (seconds/visit)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totvisitdur_plot




#### Total duration per visit on petals per 30 min [INCLUDED IN THE MANUSCRIPT]

#totvisitdur2_mod1 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = poisson, data = visitation_bySpp)
totvisitdur2_mod2 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = nbinom2(), data = visitation_bySpp)
AIC(totvisitdur2_mod2)  # poisson model doesn't converge, but the negbinomial model does
summary(totvisitdur2_mod2)
Anova(totvisitdur2_mod2)
# tests of model diagnostics
testDispersion(totvisitdur2_mod2)
testZeroInflation(totvisitdur2_mod2)

# spatial autocorrelation test
totvisitdur2_resid <- simulateResiduals(totvisitdur2_mod2)
totvisitdur2_resid2 <- recalculateResiduals(totvisitdur2_resid, group = visitation_bySpp$Site) # group residuals by site
testSpatialAutocorrelation(totvisitdur2_resid2, unique(visitation_bySpp$Long), unique(visitation_bySpp$Lat)) # not sig!

# post hoc tests
totvisitdur2_contrasts <- emmeans(totvisitdur2_mod2, spec = pairwise ~ Genus, type = "response", offset = T)
totvisitdur2_contrasts     # Bombus is sig lower than Apis or Other, Apis and Other not sig different
totvisitdur2_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0001, "***",
  "APIS",     "PEPO", 0.0001, "***",
  "APIS",     "Other", 0.3497, "NS",
  "BOMB",     "PEPO", 0.9694, "NS",
  "BOMB",     "Other", 0.0001, "***",
  "PEPO",     "Other", 0.0001, "***")

totvisitdur2_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0001, "***",
  "APIS",     "PEPO", 0.0001, "***",
  "BOMB",     "Other", 0.0001, "***",
  "PEPO",     "Other", 0.0001, "***")
me_totvisitdur2 <- ggpredict(totvisitdur2_mod2, "Genus")
plot(me_totvisitdur2) +    # need to add in duration per visit to petals raw data points instead of total dur data points
  scale_y_log10()

# Total duration per visit on petals by genus figure
totvisitdur2_plot <- ggplot(me_totvisitdur2, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = visitdur2, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  #annotate(geom = "text", x = 1.8, y = 300, label = bquote(paste("Species, ", chi^2, " = 5.0, p = 0.08")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Squash bees", "Other")) +
  stat_pvalue_manual(totvisitdur2_sig, y.position = 1.9, step.increase = 0.09, label = "p.adj.signif", hide.ns = T, label.size = 2.8) +
  labs(color = "Species", tag = "c", x = "Species", y = "Duration on Petals per Visit (seconds/visit)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black", angle = 20, vjust = 1, hjust = 1), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totvisitdur2_plot



#### Total duration per visit on pollen per 30 min [INCLUDED IN THE MANUSCRIPT]

totvisitdur4_mod1 <- glmmTMB(dur4 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = poisson, data = visitation_bySpp)
totvisitdur4_mod2 <- glmmTMB(dur4 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = nbinom2(), data = visitation_bySpp)
AIC(totvisitdur4_mod1, totvisitdur4_mod2) # neg binom model is better
summary(totvisitdur4_mod2)
Anova(totvisitdur4_mod2)
# test model diagnostics
testDispersion(totvisitdur4_mod2)
testZeroInflation(totvisitdur4_mod2)  

# spatial autocorrelation test
totvisitdur4_resid <- simulateResiduals(totvisitdur4_mod2)
totvisitdur4_resid2 <- recalculateResiduals(totvisitdur4_resid, group = visitation_bySpp$Site) # group residuals by site
testSpatialAutocorrelation(totvisitdur4_resid2, unique(visitation_bySpp$Long), unique(visitation_bySpp$Lat)) # not sig!

# post hoc tests
totvisitdur4_contrasts <- emmeans(totvisitdur4_mod2, spec = pairwise ~ Genus, type = "response", offset  = T)
totvisitdur4_contrasts    
totvisitdur4_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.4670, "NS",
  "APIS",     "PEPO", 0.1657, "NS",
  "APIS",     "Other", 0.0001, "***",
  "BOMB",     "PEPO", 0.6151, "NS",
  "BOMB",     "Other", 0.0001, "***",
  "PEPO",     "Other", 0.0005, "***")

totvisitdur4_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "Other", 0.0001, "***",
  "BOMB",     "Other", 0.0001, "***",
  "PEPO",     "Other", 0.0005, "***")
me_totvisitdur4 <- ggpredict(totvisitdur4_mod2, "Genus")
plot(me_totvisitdur4, add.data = T) +   # need to add in duration per visit to pollen raw data points instead of total dur data points
  scale_y_log10()

# Total duration per visit on pollen by genus figure
totvisitdur4_plot <- ggplot(me_totvisitdur4, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = visitdur4, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  #annotate(geom = "text", x = 1.8, y = 300, label = bquote(paste("Species, ", chi^2, " = 3.6, p = NS")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Squash bees", "Other")) +
  stat_pvalue_manual(totvisitdur4_sig, y.position = 1.7, step.increase = 0.07, label = "p.adj.signif", hide.ns = T, label.size = 2.8) +
  labs(color = "Species", x = "Species", y = "Duration on Pollen per Visit (seconds/visit)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black", angle = 30, vjust = 1, hjust = 1), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=9, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totvisitdur4_plot
ggsave(here("figures/FigS1.tiff"), plot = totvisitdur4_plot, dpi = 600, width = 4, height = 3.5, units = "in", compression="lzw")




#### Total duration per visit on pollen+nectar per 30 min [INCLUDED IN THE MANUSCRIPT]

totvisitdur5_mod1 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = poisson, data = visitation_bySpp)
totvisitdur5_mod2 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(visits+1), family = nbinom2(), data = visitation_bySpp)
AIC(totvisitdur5_mod1, totvisitdur5_mod2)  # neg binom model is better
summary(totvisitdur5_mod2)
Anova(totvisitdur5_mod2)
# tests of model diagnostics
testDispersion(totvisitdur5_mod2)
testZeroInflation(totvisitdur5_mod2)

# spatial autocorrelation test
totvisitdur5_resid <- simulateResiduals(totvisitdur5_mod2)
totvisitdur5_resid2 <- recalculateResiduals(totvisitdur5_resid, group = visitation_bySpp$Site) # group residuals by site
testSpatialAutocorrelation(totvisitdur5_resid2, unique(visitation_bySpp$Long), unique(visitation_bySpp$Lat)) # not sig!

# post hoc tests
totvisitdur5_contrasts <- emmeans(totvisitdur5_mod2, spec = pairwise ~ Genus, type = "response", offset  = T)
totvisitdur5_contrasts    # No sig differences between genus
totvisitdur5_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.9464, "NS",
  "APIS",     "PEPO", 0.5738, "NS",
  "APIS",     "Other", 0.9778, "NS",
  "BOMB",     "PEPO", 0.1814, "NS",
  "BOMB",     "Other", 0.9995, "NS",
  "PEPO",     "Other", 0.0930, "NS")
me_totvisitdur5 <- ggpredict(totvisitdur5_mod2, "Genus")
plot(me_totvisitdur5, add.data = T) +   # need to add in duration per visit to pollen + nectar raw data points instead of total dur data points
  scale_y_log10()

# Total duration per visit on pollen + nectar by genus figure
totvisitdur5_plot <- ggplot(me_totvisitdur5, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = visitdur5, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  annotate(geom = "text", x = 1.8, y = 300, label = bquote(paste("Species, ", chi^2, " = 6.75, p = 0.08")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  #stat_pvalue_manual(totvisitdur5_sig, y.position = 1.9, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "d", x = "Species", y = "Duration on Pollen + Nectar per Visit (seconds/visit)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black", angle = 20, vjust = 1, hjust = 1), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totvisitdur5_plot




#############################
# Figure 2 of all visitation behaviors by genus
############################

# Grid of both Apis visit number and duration plots
visitation_fig2 <- ggarrange(visitnum_plot, totvisitdur_plot, totvisitdur2_plot, totvisitdur5_plot, nrow=2, ncol = 2, heights = c(2.6,3), legend = "none", common.legend = T)
ggsave(here("figures/Fig2.tiff"), plot = visitation_fig2, dpi = 600, width = 6, height = 6, units = "in", compression="lzw")




###########################
# Additional analyses not in the manuscript
###########################


#### Total sum duration per 30 min [NOT IN MANUSCRIPT]
        # Note: We ended up using duration per visit instead of the total sum duration per 30 min in the manuscript (duration per visit models are below starting at line 344)

# remove one outlier from the Other pollinators that has 2800+ seconds of visits to one flower (many combined species)
visitation_bySpp_noDurOutlier <- filter(visitation_bySpp, dur < 2500)

totdur_mod1 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp_noDurOutlier)
totdur_mod2 <- glmmTMB(dur ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp_noDurOutlier)
AIC(totdur_mod1, totdur_mod2)  # negbinomial model has much lower AIC
summary(totdur_mod2)
Anova(totdur_mod2)
# test of model diagnostics
testDispersion(totdur_mod2)
testZeroInflation(totdur_mod2)

# spatial autocorrelation test
totdur_resid <- simulateResiduals(totdur_mod2)
totdur_resid2 <- recalculateResiduals(totdur_resid, group = visitation_bySpp_noDurOutlier$Site) # group residuals by site
testSpatialAutocorrelation(totdur_resid2, unique(visitation_bySpp_noDurOutlier$Long), unique(visitation_bySpp_noDurOutlier$Lat)) # not sig!

# post hoc tests
totdur_contrasts <- emmeans(totdur_mod2, spec = pairwise ~ Genus, type = "response")
totdur_contrasts     # no sig diff based on genus
totdur_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.6704, "NS",
  "APIS",     "Other", 0.6504, "NS",
  "APIS",     "PEPO", 0.9080, "NS",
  "BOMB",     "Other", 0.9972, "NS",
  "BOMB",     "PEPO", 0.2387, "NS",
  "Other",     "PEPO", 0.2021, "NS")
me_totdur <- ggpredict(totdur_mod2, "Genus")
plot(me_totdur, add.data = T) +
  scale_y_log10()

# Total sum duration per 30 min by genus figure
totdur_plot <- ggplot(me_totdur, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp_noDurOutlier, aes(x= Genus, y = dur, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  annotate(geom = "text", x = 1.8, y = 1500, label = bquote(paste("Species, ", chi^2, " = 5.3, p = NS")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  #stat_pvalue_manual(totdur_sig, y.position = 3.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 3) +
  labs(color = "Species", tag = "f", x = "Species", y = "Sum Duration of Visits (seconds)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totdur_plot


#### Total sum duration on petals per 30 min [NOT IN MANUSCRIPT]

totdur2_mod1 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp)
totdur2_mod2 <- glmmTMB(dur2 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp)
AIC(totdur2_mod1, totdur2_mod2)  # poisson model doesn't converge, but the negbinomial model does
summary(totdur2_mod2)
Anova(totdur2_mod2)
# tests of model diagnostics
testDispersion(totdur2_mod2)
testZeroInflation(totdur2_mod2)

# spatial autocorrelation test
totdur2_resid <- simulateResiduals(totdur2_mod2)
totdur2_resid2 <- recalculateResiduals(totdur2_resid, group = visitation_bySpp$Site) # group residuals by site
testSpatialAutocorrelation(totdur2_resid2, unique(visitation_bySpp$Long), unique(visitation_bySpp$Lat)) # not sig!

# post hoc tests
totdur2_contrasts <- emmeans(totdur2_mod2, spec = pairwise ~ Genus, type = "response")
totdur2_contrasts   # Bombus is sig lower than Apis or Other, Apis and Other not sig different
totdur2_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0109, "*",
  "APIS",     "PEPO", 0.0034, "**",
  "APIS",     "Other", 0.4754, "NS",
  "BOMB",     "PEPO", 0.4714, "NS",
  "BOMB",     "Other", 0.0001, "***",
  "PEPO",     "Other", 0.0001, "***")
totdur2_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0109, "*",
  "APIS",     "PEPO", 0.0034, "**",
  "BOMB",     "Other", 0.0001, "***",
  "PEPO",     "Other", 0.0001, "***")
me_totdur2 <- ggpredict(totdur2_mod2, "Genus")
plot(me_totdur2, add.data = T) +
  scale_y_log10()

# Total sum duration on petals per 30 min by genus figure
totdur2_plot <- ggplot(me_totdur2, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = dur2, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  #annotate(geom = "text", x = 1, y = 1400, label = bquote(paste("Species, ", chi^2, " = 2.3, p = NS")), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  stat_pvalue_manual(totdur2_sig, y.position = 2.8, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "g", x = "Species", y = "Sum Duration of Petal Visits (seconds)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totdur2_plot


#### Total sum duration on pollen+nectar per 30 min [NOT IN MANUSCRIPT]

totdur5_mod1 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = poisson, data = visitation_bySpp)
totdur5_mod2 <- glmmTMB(dur5 ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, family = nbinom2(), data = visitation_bySpp)
AIC(totdur5_mod1, totdur5_mod2)
summary(totdur5_mod2)
Anova(totdur5_mod2)
# test of model diagnostics
testDispersion(totdur5_mod2)
testZeroInflation(totdur5_mod2)

# spatial autocorrelation test
totdur5_resid <- simulateResiduals(totdur5_mod2)
totdur5_resid2 <- recalculateResiduals(totdur5_resid, group = visitation_bySpp$Site) # group residuals by site
testSpatialAutocorrelation(totdur5_resid2, unique(visitation_bySpp$Long), unique(visitation_bySpp$Lat)) # not sig!

# post hoc tests
totdur5_contrasts <- emmeans(totdur5_mod2, spec = pairwise ~ Genus, type = "response")
totdur5_contrasts    # No sig diff based on genus
totdur5_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.2211, "NS",
  "APIS",     "PEPO", 0.7360, "NS",
  "APIS",     "Other", 0.8940, "NS",
  "BOMB",     "PEPO", 0.0079, "**",
  "BOMB",     "Other", 0.0148, "*",
  "PEPO",     "Other", 0.9524, "NS")
totdur5_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "BOMB",     "PEPO", 0.0079, "**",
  "BOMB",     "Other", 0.0148, "*")
me_totdur5 <- ggpredict(totdur5_mod2, "Genus")
plot(me_totdur5, add.data = T) +
  scale_y_log10()


# Total sum duration on pollen+nectar per 30 min by genus figure
totdur5_plot <- ggplot(me_totdur5, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = dur5, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  #annotate(geom = "text", x = 1.8, y = 1500, label = bquote(paste("Species, ", chi^2, " = 4.9, p = 0.08")), size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  stat_pvalue_manual(totdur5_sig, y.position = 3.5, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 3) +
  labs(color = "Species", tag = "h", x = "Species", y = "Sum Duration of Pollen + Nectar Visits (seconds)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
totdur5_plot



#### Visit rate (visit per min) [NOT IN MANUSCRIPT]

# remove one outlier from the Other pollinators that has 136 visits to one flower (many combined species)
visitation_bySpp_noVisitNumOutlier <- filter(visitation_bySpp, visits < 100)

visitrate_mod1 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(totdur_min), family = poisson, data = visitation_bySpp_noVisitNumOutlier)
visitrate_mod2 <- glmmTMB(visits ~ Genus + (1|Site/Visit/FlowerID), ziformula = ~1 + Genus, offset = log(totdur_min), family = nbinom2(), data = visitation_bySpp_noVisitNumOutlier)
AIC(visitrate_mod1, visitrate_mod2)  # negbinomial model has much lower AIC
summary(visitrate_mod2)
Anova(visitrate_mod2)
# test of model diagnostics
testDispersion(visitrate_mod2)
testZeroInflation(visitrate_mod2)

# spatial autocorrelation test
visitrate_resid <- simulateResiduals(visitrate_mod2)
visitrate_resid2 <- recalculateResiduals(visitrate_resid, group = visitation_bySpp$Site) # group residuals by site
testSpatialAutocorrelation(visitrate_resid2, unique(visitation_bySpp$Long), unique(visitation_bySpp$Lat)) # not sig!

# post hoc tests
visitrate_contrasts <- emmeans(visitrate_mod2, spec = pairwise ~ Genus, type = "response", offset = T)
visitrate_contrasts     #  Apis is sig lower than Bombus, but not with Other. Bombus and Other not sig different
visitrate_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0008, "***",
  "APIS",     "PEPO", 0.7545, "NS",
  "APIS",     "Other", 0.9639, "NS",
  "BOMB",     "PEPO", 0.4563, "NS",
  "BOMB",     "Other", 0.0140, "*",
  "PEPO",     "Other", 0.8258, "NS")
visitrate_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~p.adj.signif,
  "APIS",     "BOMB", 0.0008, "***",
  "BOMB",     "Other", 0.0140, "*")
me_visitrate <- ggpredict(visitrate_mod2, "Genus")
plot(me_visitrate) +  # need to add in rate raw data points instead of visits data points
  scale_y_log10()


# Visit rate (visits per min) by genus figure
visitrate_plot <- ggplot(me_visitrate, aes(x = x, y= predicted)) +
  geom_jitter(data = visitation_bySpp, aes(x= Genus, y = rate, color = Genus), alpha = 0.6, width = 0.2, size = 1.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.05) +
  scale_color_manual(values = Tol_muted3, name = "Species", labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10)) +
  scale_x_discrete(labels = c("Honeybees", "Bumblebees", "Squash bees", "Other Pollinators")) +
  stat_pvalue_manual(visitnum_sig, y.position = 0.45, step.increase = 0.1, label = "p.adj.signif", hide.ns = T, label.size = 2.5) +
  labs(color = "Species", tag = "e", x = "Species", y = "Rate of Visits (per min)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=7, color = "black"), axis.text.y = element_text(size=8, color = "black"), axis.title = element_text(size=8, color="black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))
visitrate_plot









#############################
# Combined figure of all visitation behaviors by genus
############################

# Grid of both Apis visit number and duration plots
visitation_all_fig <- grid_arrange_shared_legend(visitnum_plot, totvisitdur_plot, totvisitdur2_plot, totvisitdur5_plot, visitrate_plot, totdur_plot, totdur2_plot, totdur5_plot, nrow=2, ncol = 4, position = "right")
ggsave(here("figures/All visitation models.tiff"), plot = visitation_all_fig, dpi = 600, width = 11, height = 7, units = "in", compression="lzw")


