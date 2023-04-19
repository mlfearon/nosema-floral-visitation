# Calculations and data processing of Vairimorpha (=Nosema) ceranae prevalence in Apis and Bombus
# and pollinator visitation behavior to flowers

# This script produces the full compiled data for "Honeybee visitation behavior on shared floral resources 
# increases Vairimorpha ceranae prevalence in bumblebees". The resulting data set is analyzed using the code
# titled NosemaAnalysis_10.28.2022.R


# Written by: Michelle Fearon and Maryellen Zbrozek
# Last updated: 13 April 2023




# set the working directory
setwd("C:/Users/Maryellen/Documents/School/Honors Thesis/Manuscript/Analysis")


# Import libraries needed 
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(epiR)
library(here)

# set the path to the script relative to the project root directory
here::i_am("scripts/NosemaAnalysis_2.20.2023.R")


# load binary nosema data
nosema <- read.csv("data/Nosema_pos_update_missing18Sremoved.csv", stringsAsFactors = F)
head(nosema)
summary(nosema)
dim(nosema)
nosema

# check for any samples that did not produce an 18S band
No18S <- filter(nosema, X18S == 0 | is.na(X18S))
   #found 13 examples, correspond to missing data
check18S <- filter(nosema, is.na(X18S)) #all are n/a
check18S2 <- filter(nosema, X18S == 0)

#filter nosema dataset to only include samples with 18S data
nosema <- filter(nosema, X18S == 1)



# this site does not have any Bombus samples, so it was removed from the analysis
nosema <- filter(nosema, Site != "SP")

unique(nosema$Site)


###### add in visitation covariates to data set 
visitation <- read.csv("data/Video_Summary_2016_forNosemaAnalysis.csv", stringsAsFactors = F)
summary(visitation)


# calculate proportion of apis visits, prop bombus visits, total non-apis or bombus visits, and duration of time per visit doing each behavior
visitation <- visitation %>%
  select(FlowerID, Site, Year, Visit, Date, StartHour, totdur_min:Other_dur, APIS_dur2:Other_dur5) %>%
  mutate(APIS_freq = if_else(VisitNum == 0, 0, APIS_visits/VisitNum), BOMB_freq = if_else(VisitNum == 0, 0, BOMB_visits/VisitNum),
         Native_freq = if_else(VisitNum == 0, 0, Native_visits/VisitNum), Other_freq = if_else(VisitNum == 0, 0, Other_visits/VisitNum), Other_rate = Other_visits/totdur_min,
         APIS_visitdur = if_else(APIS_visits == 0, 0, APIS_dur/APIS_visits), APIS_visitdur2 = if_else(APIS_visits == 0, 0, APIS_dur2/APIS_visits),
         APIS_visitdur3 = if_else(APIS_visits == 0, 0, APIS_dur3/APIS_visits), APIS_visitdur4 = if_else(APIS_visits == 0, 0, APIS_dur4/APIS_visits),
         APIS_visitdur5 = if_else(APIS_visits == 0, 0, APIS_dur5/APIS_visits), BOMB_visitdur = if_else(BOMB_visits == 0, 0, BOMB_dur/BOMB_visits),
         BOMB_visitdur2 = if_else(BOMB_visits == 0, 0, BOMB_dur2/BOMB_visits), BOMB_visitdur3 = if_else(BOMB_visits == 0, 0, BOMB_dur3/BOMB_visits),
         BOMB_visitdur4 = if_else(BOMB_visits == 0, 0, BOMB_dur4/BOMB_visits), BOMB_visitdur5 = if_else(BOMB_visits == 0, 0, BOMB_dur5/BOMB_visits),
         Other_visitdur = if_else(Other_visits == 0, 0, Other_dur/Other_visits), Other_visitdur2 = if_else(Other_visits == 0, 0, Other_dur2/Other_visits),
         Other_visitdur3 = if_else(Other_visits == 0, 0, Other_dur3/Other_visits), Other_visitdur4 = if_else(Other_visits == 0, 0, Other_dur4/Other_visits), 
         Other_visitdur5 = if_else(Other_visits == 0, 0, Other_dur5/Other_visits))

## Note (these durations are all in seconds per visit and calculated for Apis, Bombus, and all other species)
# visitdur = total duration
# visitdur2 = duration on petals only
# visitdur3 = duration on nectar only
# visitdur4 = duration on pollen only
# visitdur5 = duration on pollen and nectar simultaneously

# total time in minutes of each video recording
mean(visitation$totdur_min)
median(visitation$totdur_min)
sd(visitation$totdur_min)
hist(visitation$totdur_min)

# calculate shannon index per flower
spp_visits <- visitation %>%
  select(APIS_visits:TRIE_visits)
visitation$VisitShannon = diversity(spp_visits, "shannon")


View(visitation)



cor.test(visitation$APIS_visits, visitation$BOMB_visits)
cor.test(visitation$APIS_visits, visitation$Other_visits)
cor.test(visitation$Other_visits, visitation$BOMB_visits)

# average based on visit (1 or 2) to each site
avgs <- visitation %>%
  group_by(Site, Visit) %>%
  summarize(VisitDur = mean(VisitDur), VisitNum = mean(VisitNum), VisitRichnessPerFlower = mean(VisitRichness),
            APIS_visits = mean(APIS_visits), BOMB_visits = mean(BOMB_visits), Other_visits = mean(Other_visits), Native_visits = mean(Native_visits),
            APIS_freq = mean(APIS_freq), BOMB_freq = mean(BOMB_freq), Other_freq = mean(Other_freq), Native_freq = mean(Native_freq),
            APIS_rate = mean(APIS_rate), BOMB_rate = mean(BOMB_rate), Other_rate = mean(Other_rate),
            APIS_dur = mean(APIS_dur), APIS_dur2 = mean(APIS_dur2), APIS_dur3 = mean(APIS_dur3), APIS_dur4 = mean(APIS_dur4),
            APIS_dur5 = mean(APIS_dur5), APIS_visitdur = mean(APIS_visitdur), APIS_visitdur2 = mean(APIS_visitdur2), 
            APIS_visitdur3 = mean(APIS_visitdur3), APIS_visitdur4 = mean(APIS_visitdur4), APIS_visitdur5 = mean(APIS_visitdur5),
            BOMB_dur = mean(BOMB_dur), BOMB_dur2 = mean(BOMB_dur2), BOMB_dur3 = mean(BOMB_dur3), BOMB_dur4 = mean(BOMB_dur4), 
            BOMB_dur5 = mean(BOMB_dur5), BOMB_visitdur = mean(BOMB_visitdur), BOMB_visitdur2 = mean(BOMB_visitdur2), 
            BOMB_visitdur3 = mean(BOMB_visitdur3), BOMB_visitdur4 = mean(BOMB_visitdur4), BOMB_visitdur5 = mean(BOMB_visitdur5),
            Other_dur = mean(Other_dur), Other_dur2 = mean(Other_dur2), Other_dur3 = mean(Other_dur3), Other_dur4 = mean(Other_dur4),
            Other_dur5 = mean(Other_dur5), Other_visitdur = mean(Other_visitdur), Other_visitdur2 = mean(Other_visitdur2), 
            Other_visitdur3 = mean(Other_visitdur3), Other_visitdur4 = mean(Other_visitdur4), Other_visitdur5 = mean(Other_visitdur5),
            VisitShannon = mean(VisitShannon))
View(avgs)

hist(avgs$VisitNum)

# add additional visitation variables to data set
data <- left_join(nosema, avgs, by = c("Site", "Visit"))
head(data)

# convert variables into characters that were misread in by R
data$Year <- as.character(data$Year)
data$Visit <- as.character(data$Visit)
data$Sex <- as.character(data$Sex)


#Check that the variables updated correctly
summary(data)

unique(data$Site)


# Data file to be used in analyses
write.csv(data, "data/NosemaAnalysis_18Sremoved.csv", quote = F)




# Table 2 in manuscript
# create summary table of visitation factors for each field site
avgs_perSite <- visitation %>%
  group_by(Site) %>%
  summarize(VisitRichnessPerFlower = mean(VisitRichness), VisitRichnessPerFlower_min = min(VisitRichness), VisitRichnessPerFlower_max = max(VisitRichness),
            APIS_visits_mean = mean(APIS_visits), APIS_visits_min = min(APIS_visits), APIS_visits_max = max(APIS_visits), 
            BOMB_visits_mean = mean(BOMB_visits), BOMB_visits_min = min(BOMB_visits), BOMB_visits_max = max(BOMB_visits),
            Other_visits_mean = mean(Other_visits), Other_visits_min = min(Other_visits), Other_visits_max = max(Other_visits),
            APIS_rate_mean = mean(APIS_rate), APIS_rate_min = min(APIS_rate), APIS_rate_max = max(APIS_rate), 
            BOMB_rate_mean = mean(BOMB_rate), BOMB_rate_min = min(BOMB_rate), BOMB_rate_max = max(BOMB_rate),
            Other_rate_mean = mean(Other_rate), Other_rate_min = min(Other_rate), Other_rate_max = max(Other_rate),
            APIS_visitdur_mean = mean(APIS_visitdur), APIS_visitdur_min = min(APIS_visitdur), APIS_visitdur_max = max(APIS_visitdur),
            BOMB_visitdur_mean = mean(BOMB_visitdur), BOMB_visitdur_min = min(BOMB_visitdur), BOMB_visitdur_max = max(BOMB_visitdur),
            Other_visitdur_mean = mean(Other_visitdur), Other_visitdur_min = min(Other_visitdur), Other_visitdur_max = max(Other_visitdur),
            APIS_visitdur2_mean = mean(APIS_visitdur2), APIS_visitdur2_min = min(APIS_visitdur2), APIS_visitdur2_max = max(APIS_visitdur2),
            APIS_visitdur3_mean = mean(APIS_visitdur3), APIS_visitdur3_min = min(APIS_visitdur3), APIS_visitdur3_max = max(APIS_visitdur3),
            APIS_visitdur4_mean = mean(APIS_visitdur4), APIS_visitdur4_min = min(APIS_visitdur4), APIS_visitdur4_max = max(APIS_visitdur4),
            APIS_visitdur5_mean = mean(APIS_visitdur5), APIS_visitdur5_min = min(APIS_visitdur5), APIS_visitdur5_max = max(APIS_visitdur5),
            BOMB_visitdur2_mean = mean(BOMB_visitdur2), BOMB_visitdur2_min = min(BOMB_visitdur2), BOMB_visitdur2_max = max(BOMB_visitdur2), 
            BOMB_visitdur3_mean = mean(BOMB_visitdur3), BOMB_visitdur3_min = min(BOMB_visitdur3), BOMB_visitdur3_max = max(BOMB_visitdur3),
            BOMB_visitdur4_mean = mean(BOMB_visitdur4), BOMB_visitdur4_min = min(BOMB_visitdur4), BOMB_visitdur4_max = max(BOMB_visitdur4),
            BOMB_visitdur5_mean = mean(BOMB_visitdur5),BOMB_visitdur5_min = min(BOMB_visitdur5), BOMB_visitdur5_max = max(BOMB_visitdur5),
            Other_visitdur2_mean = mean(Other_visitdur2), Other_visitdur2_min = min(Other_visitdur2), Other_visitdur2_max = max(Other_visitdur2),
            Other_visitdur3_mean = mean(Other_visitdur3), Other_visitdur3_min = min(Other_visitdur3), Other_visitdur3_max = max(Other_visitdur3),
            Other_visitdur4_mean = mean(Other_visitdur4), Other_visitdur4_min = min(Other_visitdur4), Other_visitdur4_max = max(Other_visitdur4),
            Other_visitdur5_mean = mean(Other_visitdur5), Other_visitdur5_min = min(Other_visitdur5), Other_visitdur5_max = max(Other_visitdur5),)
View(avgs_perSite)


# export file to make the table
write.csv(avgs_perSite, "tables/MeanBeeVisitationPerSite.csv", quote = F)



######
# Figures of Apis vs Bombus visitation metrics 
#####

plot <- ggplot(visitation, aes(x = BOMB_visits, y = APIS_visits)) +
  geom_jitter(width = 0.3, height = 0.3, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot)

plot2 <- ggplot(avgs, aes(x = BOMB_visits, y = APIS_visits)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot2)

plot3 <- ggplot(visitation, aes(x = BOMB_dur, y = APIS_dur)) +
  geom_jitter(width = 5, height = 3, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot3)

plot4 <- ggplot(avgs, aes(x = BOMB_dur, y = APIS_dur)) +
  geom_jitter(width = 5, height = 3, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot4)

plot5 <- ggplot(avgs, aes(x = BOMB_visitdur, y = APIS_visitdur)) +
  geom_jitter(width = 1, height = 1, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot5)

plot6 <- ggplot(visitation, aes(x = BOMB_dur2, y = APIS_dur2)) +
  geom_jitter(width = 1, height = 2, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot6)

plot7 <- ggplot(avgs, aes(x = BOMB_dur2, y = APIS_dur2)) +
  geom_jitter(width = 1, height = 2, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot7)

plot8 <- ggplot(avgs, aes(x = BOMB_visitdur2, y = APIS_visitdur2)) +
  geom_jitter(width = 0, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot8)


plot9 <- ggplot(visitation, aes(x = BOMB_dur5, y = APIS_dur5)) +
  geom_jitter(width = 1, height = 2, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot9)

plot10 <- ggplot(avgs, aes(x = BOMB_dur5, y = APIS_dur5)) +
  geom_jitter(width = 1, height = 2, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot10)

plot11 <- ggplot(avgs, aes(x = BOMB_visitdur5, y = APIS_visitdur5)) +
  geom_jitter(width = 0, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, linetype = "dashed")
print(plot11)

##############################
# Manipulate visitation data set to analyze differences in visitation by species
#############################

# Data collected per flower
dim(visitation)
vistation_meta <- visitation %>%
  select(FlowerID:behavdur_sec)

visitation_bySpp <- visitation %>%
  select(FlowerID:totdur_sec, matches("APIS"), matches("BOMB"), matches("Other")) %>%
  pivot_longer(APIS_visits:Other_visitdur5, names_to = c("Genus", "visit_metric"), names_sep = "_", values_to = "value") %>%
  filter(Genus != "APISBOMB", visit_metric != "freq") %>%
  pivot_wider(names_from = "visit_metric", values_from = "value")
dim(visitation_bySpp)
View(visitation_bySpp)

# export file for visitation analyses
write.csv(visitation_bySpp, "data/Visitation_bySpp.csv", quote = F)

avgs_bySpp <- avgs %>%
  select(Site:Visit, matches("APIS"), matches("BOMB"), matches("Other")) %>%
  pivot_longer(APIS_visits:Other_visitdur5, names_to = c("Genus", "visit_metric"), names_sep = "_", values_to = "value") %>%
  filter(visit_metric != "freq") %>%
  pivot_wider(names_from = "visit_metric", values_from = "value")
View(avgs_bySpp)

# export file for visitation analyses
write.csv(avgs_bySpp, "data/VisitationAvgs_bySpp.csv", quote = F)

 



