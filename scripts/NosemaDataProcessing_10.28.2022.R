# Calculations and data processing of Vairimorpha (=Nosema) ceranae prevalence in Apis and Bombus
# and pollinator visitation behavior to flowers

# This script produces the full compiled data for "Honeybee visitation behavior on shared floral resources 
# increases Vairimorpha ceranae prevalence in bumblebees". The resulting data set is analyzed using the code
# titled NosemaAnalysis_10.28.2022.R


# Written by: Michelle Fearon and Maryellen Zbrozek
# Last updated: 28 October 2022




rm(list = ls())
# set the working directory
setwd("C:/Users/mlfearon/OneDrive - Umich/PROJECTS/Manuscripts/Nosema and visitation patterns/Data and Codes")
setwd("C:/Users/Maryellen/Documents/School/Honors Thesis/Manuscript/Analysis")


# Import libraries needed 
library(dplyr)
library(ggplot2)
library(vegan)
library(epiR)


# load binary nosema data
nosema <- read.csv("Nosema_pos_update_missing18Sremoved.csv", stringsAsFactors = F)
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
visitation <- read.csv("Video_Summary_2016_forNosemaAnalysis.csv", stringsAsFactors = F)
summary(visitation)


# calculate proportion of apis visits, prop bombus visits, total non-apis or bombus visits, and duration of time per visit doing each behavior
visitation <- visitation %>%
  select(FlowerID, Site, Year, Visit, Date, StartHour, totdur_min:Other_dur, APIS_dur2:Other_dur5) %>%
  mutate(APIS_freq = if_else(VisitNum == 0, 0, APIS_visits/VisitNum), BOMB_freq = if_else(VisitNum == 0, 0, BOMB_visits/VisitNum),
         Native_freq = if_else(VisitNum == 0, 0, Native_visits/VisitNum), Other_freq = if_else(VisitNum == 0, 0, Other_visits/VisitNum),
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



# calculate shannon index per flower
spp_visits <- visitation %>%
  select(APIS_visits:TRIE_visits)
visitation$VisitShannon = diversity(spp_visits, "shannon")


View(visitation)



# average based on visit (1 or 2) to each site
avgs <- visitation %>%
  group_by(Site, Visit) %>%
  summarize(VisitDur = mean(VisitDur), VisitNum = mean(VisitNum), VisitRichnessPerFlower = mean(VisitRichness),
            APIS_visits = mean(APIS_visits), BOMB_visits = mean(BOMB_visits),
            APIS_freq = mean(APIS_freq), BOMB_freq = mean(BOMB_freq), Native_visits = mean(Native_visits), 
            Other_visits = mean(Other_visits), Native_freq = mean(Native_freq), Other_freq = mean(Other_freq), 
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
write.csv(data, "NosemaAnalysis_17Oct2022_18Sremoved.csv", quote = F)




# Table 2 in manuscript
# create summary table of visitation factors for each field site
avgs_perSite <- visitation %>%
  group_by(Site) %>%
  summarize(VisitRichnessPerFlower = mean(VisitRichness), VisitRichnessPerFlower_min = min(VisitRichness), VisitRichnessPerFlower_max = max(VisitRichness),
            APIS_visits_mean = mean(APIS_visits), APIS_visits_min = min(APIS_visits), APIS_visits_max = max(APIS_visits), 
            BOMB_visits_mean = mean(BOMB_visits), BOMB_visits_min = min(BOMB_visits), BOMB_visits_max = max(BOMB_visits),
            Other_visits_mean = mean(Other_visits), Other_visits_min = min(Other_visits), Other_visits_max = max(Other_visits),
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
write.csv(avgs_perSite, "MeanBeeVisitationPerSite_17Oct2022.csv", quote = F)





