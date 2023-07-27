# Calculations and data processing of Vairimorpha (=Nosema) ceranae prevalence in Apis and Bombus
# and pollinator visitation behavior to flowers

# This script produces the full compiled data for "Honeybee visitation behavior on shared flowers 
# increases Vairimorpha ceranae prevalence in bumblebees". The resulting data set is analyzed using the code
# titled NosemaAnalysis_Apr.2023.R and VisitationAnalysis_Apr.2023.R


# Written by: Michelle Fearon and Maryellen Zbrozek
# Last updated: 27 July 2023



# Import libraries needed 
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(epiR)
library(here)

# set the path to the script relative to the project root directory
here::i_am("scripts/NosemaDataProcessing.R")


## Reviewer had a question about Bombus species sampled from each site and visit, requested a table for the Appendix
# load pollinator community data previously published in Fearon and Tibbetts 2021, https://doi.org/10.5061/dryad.zpc866t7g
poll_comm <- read.csv("data/PollinatorComm2015_2016_publish.csv", stringsAsFactors = F)
head(poll_comm)
poll_comm_bombus <- filter(poll_comm, Genus == "Bombus", Year == 2016, Site != "T", Site != "SP") # only include sampling visits included in this study

#################################
## Appendix S1: Table S2 - Bombus sampling at each site visit
#################################
bombus_sampling <- table(poll_comm_bombus$SiteVisit, poll_comm_bombus$Code)
write.csv(bombus_sampling, "tables/AppendixTableS2_BombusSampling.csv", quote = F)


# Squash bee samping at each site
poll_comm_pepo <- filter(poll_comm, Genus == "Peponapis", Year == 2016, Site != "T", Site != "SP") # only include sampling visits included in this study
squashbee_sampling <- table(poll_comm_pepo$Site, poll_comm_pepo$Code)



#### load binary nosema data
nosema <- read.csv("data/Nosema_pos.csv", stringsAsFactors = F)
head(nosema)
summary(nosema)
dim(nosema)
nosema

# check for any samples that did not produce an 18S band
No18S <- filter(nosema, X18S == 0 | is.na(X18S))
No18S
   #None, all samples have a good 18S band


unique(nosema$Site)

# remove unnecessary columns
nosema <- select(nosema, Sample:Sex)



# get unique Latitude and Longitude per site
LatLong <- nosema %>% select(Site, Lat, Long) %>%
  group_by(Site) %>%
  summarize(Site = unique(Site),
            Lat = unique(Lat),
            Long = unique(Long))
LatLong



###### add in visitation covariates to data set 
visitation <- read.csv("data/Video_Summary_2016_forNosemaAnalysis.csv", stringsAsFactors = F)
summary(visitation)

# Note:
# Other1 = All pollinator taxa excluding Apis and Bombus
# Other 2 = App pollinator taxa excluding Apis, Bombus, and Eucera (squash bees, PEPO)


# calculate frequency of apis visits, frequency of bombus visits, total non-apis and non-bombus visits, and duration of time per visit doing each behavior
visitation <- visitation %>%
  select(FlowerID, Site, Year, Visit, Date, StartHour, totdur_min:Other2_dur, APIS_dur2:Other2_dur5) %>%
  mutate(APIS_freq = if_else(VisitNum == 0, 0, APIS_visits/VisitNum), BOMB_freq = if_else(VisitNum == 0, 0, BOMB_visits/VisitNum), PEPO_freq = if_else(VisitNum == 0, 0, PEPO_visits/VisitNum),
         Native_freq = if_else(VisitNum == 0, 0, Native_visits/VisitNum), Other1_freq = if_else(VisitNum == 0, 0, Other1_visits/VisitNum), 
         Other2_freq = if_else(VisitNum == 0, 0, Other2_visits/VisitNum), Other1_rate = Other1_visits/totdur_min, Other2_rate = Other2_visits/totdur_min,
         APIS_visitdur = if_else(APIS_visits == 0, 0, APIS_dur/APIS_visits), APIS_visitdur2 = if_else(APIS_visits == 0, 0, APIS_dur2/APIS_visits),
         APIS_visitdur3 = if_else(APIS_visits == 0, 0, APIS_dur3/APIS_visits), APIS_visitdur4 = if_else(APIS_visits == 0, 0, APIS_dur4/APIS_visits),
         APIS_visitdur5 = if_else(APIS_visits == 0, 0, APIS_dur5/APIS_visits), BOMB_visitdur = if_else(BOMB_visits == 0, 0, BOMB_dur/BOMB_visits),
         BOMB_visitdur2 = if_else(BOMB_visits == 0, 0, BOMB_dur2/BOMB_visits), BOMB_visitdur3 = if_else(BOMB_visits == 0, 0, BOMB_dur3/BOMB_visits),
         BOMB_visitdur4 = if_else(BOMB_visits == 0, 0, BOMB_dur4/BOMB_visits), BOMB_visitdur5 = if_else(BOMB_visits == 0, 0, BOMB_dur5/BOMB_visits),
         PEPO_visitdur = if_else(PEPO_visits == 0, 0, PEPO_dur/PEPO_visits),PEPO_visitdur2 = if_else(PEPO_visits == 0, 0, PEPO_dur2/PEPO_visits), PEPO_visitdur3 = if_else(PEPO_visits == 0, 0, PEPO_dur3/PEPO_visits),
         PEPO_visitdur4 = if_else(PEPO_visits == 0, 0, PEPO_dur4/PEPO_visits), PEPO_visitdur5 = if_else(PEPO_visits == 0, 0, PEPO_dur5/PEPO_visits),
         Other1_visitdur = if_else(Other1_visits == 0, 0, Other1_dur/Other1_visits), Other1_visitdur2 = if_else(Other1_visits == 0, 0, Other1_dur2/Other1_visits),
         Other1_visitdur3 = if_else(Other1_visits == 0, 0, Other1_dur3/Other1_visits), Other1_visitdur4 = if_else(Other1_visits == 0, 0, Other1_dur4/Other1_visits), 
         Other1_visitdur5 = if_else(Other1_visits == 0, 0, Other1_dur5/Other1_visits), Other2_visitdur = if_else(Other2_visits == 0, 0, Other2_dur/Other2_visits), 
         Other2_visitdur2 = if_else(Other2_visits == 0, 0, Other2_dur2/Other2_visits), Other2_visitdur3 = if_else(Other2_visits == 0, 0, Other2_dur3/Other2_visits), 
         Other2_visitdur4 = if_else(Other2_visits == 0, 0, Other2_dur4/Other2_visits), Other2_visitdur5 = if_else(Other2_visits == 0, 0, Other2_dur5/Other2_visits))

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

# calculate the total duration per all bee visits
visitation$DurationPerVisit <- visitation$VisitDur/(visitation$VisitNum)
visitation$DurationPerVisit[is.na(visitation$DurationPerVisit)] <- 0  # convert NAs where there were zero visits to zero

View(visitation)



## Correlation between visit number and total duration (summed seconds of all visits)
cor.test(visitation$VisitNum, visitation$VisitDur)

## Correlation between visit number and total duration per visit (seconds/visits)
cor.test(visitation$VisitNum, visitation$DurationPerVisit)


# Ranges of visit number for each species
range(visitation$APIS_visits)
range(visitation$BOMB_visits)
range(visitation$PEPO_visits)
range(visitation$Other1_visits)  # with PEPO
range(visitation$Other2_visits)  # without PEPO

# Ranges of total duration per visit for each species (seconds per visit)
range(visitation$APIS_visitdur)
range(visitation$BOMB_visitdur)
range(visitation$PEPO_visitdur)
range(visitation$Other1_visitdur)  # with PEPO
range(visitation$Other2_visitdur)  # without PEPO


# Ranges of total duration per visit on petals for each species (seconds per visit)
range(visitation$APIS_visitdur2)
range(visitation$BOMB_visitdur2)
range(visitation$PEPO_visitdur2)
range(visitation$Other1_visitdur2)  # with PEPO
range(visitation$Other2_visitdur2)  # without PEPO

# Ranges of total duration per visit on pollen for each species (seconds per visit)
range(visitation$APIS_visitdur4)
range(visitation$BOMB_visitdur4)
range(visitation$PEPO_visitdur4)
range(visitation$Other1_visitdur4)  # with PEPO
range(visitation$Other2_visitdur4)  # without PEPO

# Ranges of total duration per visit on pollen+nectar for each species (seconds per visit)
range(visitation$APIS_visitdur5)
range(visitation$BOMB_visitdur5)
range(visitation$PEPO_visitdur5)
range(visitation$Other1_visitdur5)  # with PEPO
range(visitation$Other2_visitdur5)  # without PEPO


# average based on visit (1 or 2) to each site
avgs <- visitation %>%
  group_by(Site, Visit) %>%
  summarize(VisitDur = mean(VisitDur), VisitNum = mean(VisitNum), VisitRichnessPerFlower = mean(VisitRichness),
            APIS_visits = mean(APIS_visits), BOMB_visits = mean(BOMB_visits), PEPO_visits = mean(PEPO_visits), Other1_visits = mean(Other1_visits), Other2_visits = mean(Other2_visits), Native_visits = mean(Native_visits),
            APIS_freq = mean(APIS_freq), BOMB_freq = mean(BOMB_freq), PEPO_freq = mean(PEPO_freq), Other1_freq = mean(Other1_freq), Other2_freq = mean(Other2_freq), Native_freq = mean(Native_freq),
            APIS_rate = mean(APIS_rate), BOMB_rate = mean(BOMB_rate), PEPO_rate = mean(PEPO_rate), Other1_rate = mean(Other1_rate), Other2_rate = mean(Other2_rate),
            APIS_dur = mean(APIS_dur), APIS_dur2 = mean(APIS_dur2), APIS_dur3 = mean(APIS_dur3), APIS_dur4 = mean(APIS_dur4),
            APIS_dur5 = mean(APIS_dur5), APIS_visitdur = mean(APIS_visitdur), APIS_visitdur2 = mean(APIS_visitdur2), 
            APIS_visitdur3 = mean(APIS_visitdur3), APIS_visitdur4 = mean(APIS_visitdur4), APIS_visitdur5 = mean(APIS_visitdur5),
            BOMB_dur = mean(BOMB_dur), BOMB_dur2 = mean(BOMB_dur2), BOMB_dur3 = mean(BOMB_dur3), BOMB_dur4 = mean(BOMB_dur4), 
            BOMB_dur5 = mean(BOMB_dur5), BOMB_visitdur = mean(BOMB_visitdur), BOMB_visitdur2 = mean(BOMB_visitdur2), 
            BOMB_visitdur3 = mean(BOMB_visitdur3), BOMB_visitdur4 = mean(BOMB_visitdur4), BOMB_visitdur5 = mean(BOMB_visitdur5),
            PEPO_dur = mean(PEPO_dur), PEPO_dur2 = mean(PEPO_dur2), PEPO_dur3 = mean(PEPO_dur3), PEPO_dur4 = mean(PEPO_dur4), 
            PEPO_dur5 = mean(PEPO_dur5), PEPO_visitdur = mean(PEPO_visitdur), PEPO_visitdur2 = mean(PEPO_visitdur2), 
            PEPO_visitdur3 = mean(PEPO_visitdur3), PEPO_visitdur4 = mean(PEPO_visitdur4), PEPO_visitdur5 = mean(PEPO_visitdur5),
            Other1_dur = mean(Other1_dur), Other1_dur2 = mean(Other1_dur2), Other1_dur3 = mean(Other1_dur3), Other1_dur4 = mean(Other1_dur4),
            Other1_dur5 = mean(Other1_dur5), Other1_visitdur = mean(Other1_visitdur), Other1_visitdur2 = mean(Other1_visitdur2), 
            Other1_visitdur3 = mean(Other1_visitdur3), Other1_visitdur4 = mean(Other1_visitdur4), Other1_visitdur5 = mean(Other1_visitdur5),
            Other2_dur = mean(Other2_dur), Other2_dur2 = mean(Other2_dur2), Other2_dur3 = mean(Other2_dur3), Other2_dur4 = mean(Other2_dur4),
            Other2_dur5 = mean(Other2_dur5), Other2_visitdur = mean(Other2_visitdur), Other2_visitdur2 = mean(Other2_visitdur2), 
            Other2_visitdur3 = mean(Other2_visitdur3), Other2_visitdur4 = mean(Other2_visitdur4), Other2_visitdur5 = mean(Other2_visitdur5),
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
write.csv(data, "data/NosemaAnalysis.csv", quote = F)



##################################
# Appendix S1: Table S4 in manuscript
#################################

# create summary table of visitation factors for each field site
avgs_perSite <- visitation %>%
  group_by(Site) %>%
  summarize(VisitRichnessPerFlower = mean(VisitRichness), VisitRichnessPerFlower_min = min(VisitRichness), VisitRichnessPerFlower_max = max(VisitRichness),
            APIS_visits_mean = mean(APIS_visits), APIS_visits_min = min(APIS_visits), APIS_visits_max = max(APIS_visits), 
            BOMB_visits_mean = mean(BOMB_visits), BOMB_visits_min = min(BOMB_visits), BOMB_visits_max = max(BOMB_visits),
            PEPO_visits_mean = mean(PEPO_visits), PEPO_visits_min = min(PEPO_visits), PEPO_visits_max = max(PEPO_visits),
            Other_visits_mean = mean(Other2_visits), Other_visits_min = min(Other2_visits), Other_visits_max = max(Other2_visits),
            APIS_rate_mean = mean(APIS_rate), APIS_rate_min = min(APIS_rate), APIS_rate_max = max(APIS_rate), 
            BOMB_rate_mean = mean(BOMB_rate), BOMB_rate_min = min(BOMB_rate), BOMB_rate_max = max(BOMB_rate),
            PEPO_rate_mean = mean(PEPO_rate), PEPO_rate_min = min(PEPO_rate), PEPO_rate_max = max(PEPO_rate),
            Other_rate_mean = mean(Other2_rate), Other_rate_min = min(Other2_rate), Other_rate_max = max(Other2_rate),
            APIS_visitdur_mean = mean(APIS_visitdur), APIS_visitdur_min = min(APIS_visitdur), APIS_visitdur_max = max(APIS_visitdur),
            BOMB_visitdur_mean = mean(BOMB_visitdur), BOMB_visitdur_min = min(BOMB_visitdur), BOMB_visitdur_max = max(BOMB_visitdur),
            PEPO_visitdur_mean = mean(PEPO_visitdur), PEPO_visitdur_min = min(PEPO_visitdur), PEPO_visitdur_max = max(PEPO_visitdur),
            Other_visitdur_mean = mean(Other2_visitdur), Other_visitdur_min = min(Other2_visitdur), Other_visitdur_max = max(Other2_visitdur),
            APIS_visitdur2_mean = mean(APIS_visitdur2), APIS_visitdur2_min = min(APIS_visitdur2), APIS_visitdur2_max = max(APIS_visitdur2),
            APIS_visitdur3_mean = mean(APIS_visitdur3), APIS_visitdur3_min = min(APIS_visitdur3), APIS_visitdur3_max = max(APIS_visitdur3),
            APIS_visitdur4_mean = mean(APIS_visitdur4), APIS_visitdur4_min = min(APIS_visitdur4), APIS_visitdur4_max = max(APIS_visitdur4),
            APIS_visitdur5_mean = mean(APIS_visitdur5), APIS_visitdur5_min = min(APIS_visitdur5), APIS_visitdur5_max = max(APIS_visitdur5),
            BOMB_visitdur2_mean = mean(BOMB_visitdur2), BOMB_visitdur2_min = min(BOMB_visitdur2), BOMB_visitdur2_max = max(BOMB_visitdur2), 
            BOMB_visitdur3_mean = mean(BOMB_visitdur3), BOMB_visitdur3_min = min(BOMB_visitdur3), BOMB_visitdur3_max = max(BOMB_visitdur3),
            BOMB_visitdur4_mean = mean(BOMB_visitdur4), BOMB_visitdur4_min = min(BOMB_visitdur4), BOMB_visitdur4_max = max(BOMB_visitdur4),
            BOMB_visitdur5_mean = mean(BOMB_visitdur5),BOMB_visitdur5_min = min(BOMB_visitdur5), BOMB_visitdur5_max = max(BOMB_visitdur5),
            PEPO_visitdur2_mean = mean(PEPO_visitdur2), PEPO_visitdur2_min = min(PEPO_visitdur2), PEPO_visitdur2_max = max(PEPO_visitdur2), 
            PEPO_visitdur3_mean = mean(PEPO_visitdur3), PEPO_visitdur3_min = min(PEPO_visitdur3), PEPO_visitdur3_max = max(PEPO_visitdur3),
            PEPO_visitdur4_mean = mean(PEPO_visitdur4), PEPO_visitdur4_min = min(PEPO_visitdur4), PEPO_visitdur4_max = max(PEPO_visitdur4),
            PEPO_visitdur5_mean = mean(PEPO_visitdur5),PEPO_visitdur5_min = min(PEPO_visitdur5), PEPO_visitdur5_max = max(PEPO_visitdur5),
            Other_visitdur2_mean = mean(Other2_visitdur2), Other_visitdur2_min = min(Other2_visitdur2), Other_visitdur2_max = max(Other2_visitdur2),
            Other_visitdur3_mean = mean(Other2_visitdur3), Other_visitdur3_min = min(Other2_visitdur3), Other_visitdur3_max = max(Other2_visitdur3),
            Other_visitdur4_mean = mean(Other2_visitdur4), Other_visitdur4_min = min(Other2_visitdur4), Other_visitdur4_max = max(Other2_visitdur4),
            Other_visitdur5_mean = mean(Other2_visitdur5), Other_visitdur5_min = min(Other2_visitdur5), Other_visitdur5_max = max(Other2_visitdur5),)
View(avgs_perSite)



##########################
## Appendix S1, Table S4
##########################
# export file to make the table
write.csv(avgs_perSite, "tables/AppendixTableS4_MeanBeeVisitationPerSite.csv", quote = F)




##############################
# Manipulate visitation data set to analyze differences in visitation by species
#############################

# Data collected per flower
dim(visitation)
vistation_meta <- visitation %>%
  select(FlowerID:totdur_sec)

visitation_bySpp <- visitation %>%
  select(FlowerID:totdur_sec, matches("APIS"), matches("BOMB"), matches("PEPO"), matches("Other2")) %>%
  pivot_longer(APIS_visits:Other2_visitdur5, names_to = c("Genus", "visit_metric"), names_sep = "_", values_to = "value") %>%
  filter(Genus != "APISBOMB", visit_metric != "freq") %>%
  pivot_wider(names_from = "visit_metric", values_from = "value")
visitation_bySpp <- full_join(visitation_bySpp, LatLong)
dim(visitation_bySpp)
visitation_bySpp$Genus <- dplyr::recode(visitation_bySpp$Genus, Other2 = "Other")
View(visitation_bySpp)

# export file for visitation analyses
write.csv(visitation_bySpp, "data/Visitation_bySpp.csv", quote = F)

avgs_bySpp <- avgs %>%
  select(Site:Visit, matches("APIS"), matches("BOMB"), matches("PEPO"), matches("Other2")) %>%
  pivot_longer(APIS_visits:Othe2r_visitdur5, names_to = c("Genus", "visit_metric"), names_sep = "_", values_to = "value") %>%
  filter(visit_metric != "freq") %>%
  pivot_wider(names_from = "visit_metric", values_from = "value")
View(avgs_bySpp)

# export file for visitation analyses
#write.csv(avgs_bySpp, "data/VisitationAvgs_bySpp.csv", quote = F)

 

#########################
# Figures of Apis vs Bombus visitation metrics (exploration, not shown in manuscript)
########################

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


