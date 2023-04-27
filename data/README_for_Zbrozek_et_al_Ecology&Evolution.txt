Metadata for Zbrozek, Fearon, Weise, and Tibbetts. 2023. Honeybee visitation to shared floral resources increases 
Vairimorpha ceranae prevalence in bumblebees. Ecology and Evolution.

Data collected by: Michelle L. Fearon and Maryellen Zbrozek
Code written by: Michelle L. Fearon and Maryellen Zbrozek
Contact: mlfearon@umich.edu


The R code files that generates the analysis and figures for the manuscript can be used as follows: 
NosemaDataProcessing.R			Generates the calculcations and data files that are used in the analyses of V. ceranae prevalence and analyses of pollinator visitation behaviors based on pollinator group.
NosemaAnalysis_Apr.2023.R		Generates the statistical analyses and figures for V. ceranae prevalence in honeybees and bumblebees based on different visitation metrics.
VisitationAnalysis_Apr.2023.R		Generates the statistical analyses and figures of differences in visitation patterns among honeybees, bumblebees, and other pollinators.


Definitions:
Nosema or Nosema ceranae is the previous name for Vairimorpha ceranae. The data and code typically use "Nosema", while in the manuscript we updated to the current taxonomic name.
Data/Code uses "Apis" to refer to honeybees.
Data/Code uses "Bombus" or Bomb to refer to bumblebees.
Data/Code uses "Other" to refer to all other pollinators as a group.

File name:
Nosema_pos.csv

Description: 
Records of which individual bees from two focal host species (Apis mellifera and Bombus impatiens) that were positive or negative
when tested for the presence or absence of Vairimorpha (=Nosema) ceranae by PCR and additional metadata included for each sample. 
All bees were collected in or adjacent to winter squash fields. If there were more than 8 samples per host species per visit to each site
then samples were selected randomly. If less than 8 samples collected, then all samples were used in the study. Pollinator sampling overlaps with
Fearon and Tibbetts (2021) and Fearon et al. (2022), but focuses on different parasite species within the bee communities.

Headings:
Sample		Unique sample number for each individual in the study. These sample numbers match those found in PollinatorComm2015_2016_publish.csv
18S		Binary for the presence (1) or absence (0) of the 18S PCR band. NA indicates samples were not run or were contaminated.
Nosema		Binary for the presence (1) or absence (0) of Vairimorpha (=Nosema) ceranae in the sample based on PCR.
Transect.ID	Unique transect ID that indicates the site, year of collection, transect A B C or D, visit number to the site, and net or pan trap collected.  All individuals collected from the same transect share a transect ID. These unique IDs match those found in PollinatorComm2015_2016_publish.csv
Individual.ID 	Includes the transect ID with a sequential number for each sample. These unique IDs match those found in PollinatorComm2015_2016_publish.csv
Year 		Year of sample collection, all from 2016
Site		Unique site code for each of the 6 sites (see Appendix S1: Table S1 for site names, dates visited, location, etc).
Visit		Visit number to the site, each site was visited twice for this study. The first and second visit are indicated by “1” and “2”, respectively.
Site_Visit	Concatenated site code and visit number
Lat		Latitude coordinate for the Site
Long		Longitude coordinate for the Site
Date_Collected	Date the sample was collected.
Type		Either “APIS” to indicate a managed, Apis mellifera sample, or “NON” to indicate a native bee and non-Apis sample.
Family		Bee family. All were from Apidae.
Genus		Genus name of bee sample. Either Apis or Bombus.
Species		Species name of bee sample.
Code		Species identity code (code determined by the first two letters of the genus and species names)
		APME = Apis mellifera
		BOIM = Bombus impatiens
		BOFE = Bombus fervidus
		BOBI = Bombus bimaculatus
		BOPE = Bombus pensylvanicus
Sex		Sex of the bee sampled (all samples were female).
Date_ID		Date that sample was identified to species.
ID_Initials	Initials of the person that did the identification. MZ = Maryellen Zbrozek
Dissection_Date	Date that sample was dissected.
Dissect_InitialsInitials of the person that did the dissection. MZ = Maryellen Zbrozek
Notes		Additional notes and observations about the sample.


File name:
Video_Summary_2016_forNosemaAnalysis.csv

Description:
Summary of pollinator visitation from recorded approximately 30 min videos of a single flower. We visited 6 winter squash farms 
two times during the flowering period, and video recorded 8 flowers during each site visit for 30 min. We recorded the number and 
idenitity to species or morphospecies of all pollinator visitors, and the duration of their contact with each part of the inflorescence
(petals, stamen, nectar, or both stamen and nectar simultaneously). This data also includes the metadata for the time and environmental 
conditions when the videos were taken.

Headings:
FlowerID		Unique ID of each flower in the study.
Site			Unique site code for each of the 6 sites (see Appendix S1: Table S1 for site names, dates visited, location, etc).
Year			Year of sample collection, all from 2016
Site_Year		Concatenated Site and Year.
Visit			Visit number to the site, each site was visited twice for this study. The first and second visit are indicated by “1” and “2”, respectively.
Date			Date that the video recording was conducted during the visit to that site.
Lat			Latitude coordinate for the Site
Long			Longitude coordinate for the Site
StartHour		The hour that the video recording started.
TimeStart		The exact time that the video recording started (HH:MM format).
TimeStop		The exact time that the video recording ended (HH:MM format).
Duration		The total time of the video recording (H:MM:SS format)
Temperature		Temperature in degrees Fahrenheit at the start of the video recording.
WindSpeed		Wind speed in meters per second at the start of the video recording.
Clouds			Cloud cover at the start of the video recording (scale of 0 = no clouds to 7 = full cloud cover)
FlowerCollectionTime	Time that the flower was collected in a tube after the video recording.
totdur_min		Total duration of the video recording in minutes (continuous)
totdur_sec		Total duration of the video recording in seconds (continuous)
VisitDur		The duration of all pollinator visits to a flower in seconds
VisitNum		The number of all pollinator visits to a flower (number of visits per 30 min).
VisitRate		The rate (visits/min) of all pollinator visits to a flower (visit number/ total time of video recording in minutes)
VisitRichness		The number of different species (or morphospecies) that visited a single flower
Duration_2		The duration (seconds) of all pollinators that were only interacting with the petals of the flower.
Duration_3		The duration (seconds) of all pollinators that were only interacting with the nectar of the flower.
Duration_4		The duration (seconds) of all pollinators that were only interacting with the pollen of the flower.
Duration_5		The duration (seconds) of all pollinators that were interacting with the pollen and nectar simultaneously of the flower.
Duration_6		The duration (seconds) of all pollinators that were interacting directly with another pollinator on the flower.
APIS_visits		The number of Apis mellifera (honeybee) visits to a flower during the video recording duration (number of visits).
BOMB_visits		The number of Bombus spp. (bumblebee) visits to a flower during the video recording duration (number of visits).
HALI_visits		The number of small olive green halictid morphospecies (e.g. Halticus or Lasioglossum spp) visits to a flower during the video recording duration (number of visits).
HFLY_visits		The number of hoverfly visits to a flower during the video recording duration (number of visits).
PEPO_visits		The number of Eucera pruinosa (squash bee) visits to a flower during the video recording duration (number of visits).
AUGO_visits		The number of small bright green halictid morphospecies (e.g. Augochlora, Augochlorella, or Augochloropsis spp) visits to a flower during the video recording duration (number of visits).
MELI_visits		The number of Melissodes spp visits to a flower during the video recording duration (number of visits).
VESP_visits		The number of Vespula wasp spp visits to a flower during the video recording duration (number of visits).
TRIE_visits		The number of Triepeolus spp visits to a flower during the video recording duration (number of visits).
APIS_rate		The rate of Apis mellifera (honeybee) visits to a flower per minute (number of visits per min).
BOMB_rate		The rate of Bombus spp. (bumblebee) visits to a flower per minute (number of visits per min).
HALI_rate		The rate of small olive green halictid morphospecies (e.g. Halticus or Lasioglossum spp) visits to a flower per minute (number of visits per min).
HFLY_rate		The rate of hoverfly visits to a flower per minute (number of visits per min).
PEPO_rate		The rate of Eucera pruinosa (squash bee) visits to a flower per minute (number of visits per min).
AUGO_rate		The rate of small bright green halictid morphospecies (e.g. Augochlora, Augochlorella, or Augochloropsis spp) visits to a flower per minute (number of visits per min).
MELI_rate		The rate of Melissodes spp visits to a flower per minute (number of visits per min).
VESP_rate		The rate of Vespula wasp spp visits to a flower per minute (number of visits per min).
TRIE_rate		The rate of Triepeolus spp visits to a flower per minute (number of visits per min).
APISBOMB_visits		The number of Apis mellifera (honeybee) + Bombus spp. (bumblebee) visits to a flower during the video recording duration (number of visits).
Native_visits		The number of native bee (i.e. excluding honeybees) visits to a flower during the video recording duration (number of visits).
Other_visits		The number of non-honeybee and non-bumblebee visits to a flower during the video recording duration (number of visits).
APIS_dur		The duration of all Apis mellifera (honeybee) visits to a flower (any part) in seconds.
BOMB_dur		The duration of all Bombus spp. (bumblebee) visits to a flower (any part) in seconds.
Native_dur		The duration of all native bee (i.e. excluding honeybees) visits to a flower (any part) in seconds.
Other_dur		The duration of all non-honeybee and non-bumblebee visits to a flower (any part) in seconds.
APIS_dur2		The duration (seconds) of all Apis mellifera (honeybees) that were only interacting with the petals of the flower.
APIS_dur3		The duration (seconds) of all Apis mellifera (honeybees) that were only interacting with the nectar of the flower.
APIS_dur4		The duration (seconds) of all Apis mellifera (honeybees) that were only interacting with the pollen/stamen of the flower.
APIS_dur5		The duration (seconds) of all Apis mellifera (honeybees) that were interacting with the pollen/stamen and nectar simultaneously of the flower.
BOMB_dur2		The duration (seconds) of all Bombus spp. (bumblebees) that were only interacting with the petals of the flower.
BOMB_dur3		The duration (seconds) of all Bombus spp. (bumblebees) that were only interacting with the nectar of the flower.
BOMB_dur4		The duration (seconds) of all Bombus spp. (bumblebees) that were only interacting with the pollen/stamen of the flower.
BOMB_dur5		The duration (seconds) of all Bombus spp. (bumblebees) that were interacting with the pollen/stamen ad nectar simultaneously of the flower.
Other_dur2		The duration (seconds) of all non-honeybees and non-bumblebees that were only interacting with the petals of the flower.
Other_dur3		The duration (seconds) of all non-honeybees and non-bumblebees that were only interacting with the nectar of the flower.
Other_dur4		The duration (seconds) of all non-honeybees and non-bumblebees that were only interacting with the pollen/stamen of the flower.
Other_dur5		The duration (seconds) of all non-honeybees and non-bumblebees that were interacting with the pollen/stamen and nectar simultaneously of the flower.



File name:
NosemaAnalysis.csv

Description: 
Data file used for the analysis of V. ceranae prevalence in honeybees and bumblebees with each visitation metric.
Records of which individual bees from two focal host species (Apis mellifera and Bombus impatiens) that were positive or negative
when tested for the presence or absence of Vairimorpha (=Nosema) ceranae by PCR and pollinator visitation metrics averaged for each 
visit to each site. 

Headings:
Sample			Unique sample number for each individual in the study. These sample numbers match those found in PollinatorComm2015_2016_publish.csv
18S			Binary for the presence (1) or absence (0) of the 18S PCR band. NA indicates samples were not run or were contaminated.
Nosema			Binary for the presence (1) or absence (0) of Vairimorpha (=Nosema) ceranae in the sample based on PCR.
Transect.ID		Unique transect ID that indicates the site, year of collection, transect A B C or D, visit number to the site, and net or pan trap collected.  All individuals collected from the same transect share a transect ID. These unique IDs match those found in PollinatorComm2015_2016_publish.csv
Individual.ID 		Includes the transect ID with a sequential number for each sample. These unique IDs match those found in PollinatorComm2015_2016_publish.csv
Year 			Year of sample collection, all from 2016
Site			Unique site code for each of the 6 sites (see Appendix S1: Table S1 for site names, dates visited, location, etc).
Visit			Visit number to the site, each site was visited twice for this study. The first and second visit are indicated by “1” and “2”, respectively.
Site_Visit		Concatenated site code and visit number
Lat			Latitude coordinate for the Site
Long			Longitude coordinate for the Site
Date_Collected		Date the sample was collected.
Type			Either “APIS” to indicate a managed, Apis mellifera sample, or “NON” to indicate a native bee and non-Apis sample.
Family			Bee family. All were from Apidae.
Genus			Genus name of bee sample. Either Apis or Bombus.
Species			Species name of bee sample.
Code			Species identity code (code determined by the first two letters of the genus and species names)
			APME = Apis mellifera
			BOIM = Bombus impatiens
			BOFE = Bombus fervidus
			BOBI = Bombus bimaculatus
			BOPE = Bombus pensylvanicus
Sex			Sex of the bee sampled (all samples were female).
VisitDur		The average duration of all pollinator visits to a flower in seconds for each visit to each site
VisitNum		The average number of all pollinator visits to a flower for each visit to each site.
VisitRichnessPerFlower	The average number of different species (or morphospecies) that visited a single flower for each visit to each site
APIS_visits		The average number of Apis mellifera (honeybee) visits to a flower during the video recording duration for each visit to each site (number of visits).
BOMB_visits		The average number of Bombus spp. (bumblebee) visits to a flower during the video recording duration for each visit to each site (number of visits).
Other_visits		The average number of non-honeybee and non-bumblebee visits to a flower during the video recording duration for each visit to each site (number of visits).
Native_visits		The average number of native bee (i.e. excluding honeybees) visits to a flower during the video recording duration for each visit to each site (number of visits).
APIS_rate		The average rate of Apis mellifera (honeybee) visits to a flower per minute for each visit to each site (number of visits per min).
BOMB_rate		The average rate of Bombus spp. (bumblebee) visits to a flower per minute for each visit to each site (number of visits per min).
Other_rate		The average rate of non-honeybee and non-bumblebee visits to a flower per minute for each visit to each site (number of visits per min).
APIS_dur		The average duration (seconds) of all Apis mellifera (honeybees) visits that were interacting with any part of the flower for each visit to each site.
APIS_dur2		The average duration (seconds) of all Apis mellifera (honeybees) that were only interacting with the petals of the flower for each visit to each site.
APIS_dur3		The average duration (seconds) of all Apis mellifera (honeybees) that were only interacting with the nectar of the flower for each visit to each site.
APIS_dur4		The average duration (seconds) of all Apis mellifera (honeybees) that were only interacting with the pollen/stamen of the flower for each visit to each site.
APIS_dur5		The average duration (seconds) of all Apis mellifera (honeybees) that were interacting with the pollen/stamen and nectar simultaneously of the flower for each visit to each site.
APIS_visitdur		The average duration per visit (seconds per visit) of all Apis mellifera (honeybees) visits that were interacting with any part of the flower for each visit to each site.
APIS_visitdur2		The aaverage duration per visit (seconds per visit) of all Apis mellifera (honeybees) that were only interacting with the petals of the flower for each visit to each site.
APIS_visitdur3		The average duration per visit (seconds per visit) of all Apis mellifera (honeybees) that were only interacting with the nectar of the flower for each visit to each site.
APIS_visitdur4		The average duration per visit (seconds per visit) of all Apis mellifera (honeybees) that were only interacting with the pollen/stamen of the flower for each visit to each site.
APIS_visitdur5		The average duration per visit (seconds per visit) of all Apis mellifera (honeybees) that were interacting with the pollen/stamen and nectar simultaneously of the flower for each visit to each site.
BOMB_dur		The average duration (seconds) of all Bombus spp. (bumblebees) visits that were interacting with any part of the flower for each visit to each site.
BOMB_dur2		The average duration (seconds) of all Bombus spp. (bumblebees) that were only interacting with the petals of the flower for each visit to each site.
BOMB_dur3		The average duration (seconds) of all Bombus spp. (bumblebees) that were only interacting with the nectar of the flower for each visit to each site.
BOMB_dur4		The average duration (seconds) of all Bombus spp. (bumblebees) that were only interacting with the pollen/stamen of the flower for each visit to each site.
BOMB_dur5		The average duration (seconds) of all Bombus spp. (bumblebees) that were interacting with the pollen/stamen ad nectar simultaneously of the flower for each visit to each site.
BOMB_visitdur		The average duration per visit (seconds per visit) of all Bombus spp. (bumblebees) visits that were interacting with any part of the flower for each visit to each site.
BOMB_visitdur2		The average duration per visit (seconds per visit) of all Bombus spp. (bumblebees) that were only interacting with the petals of the flower for each visit to each site.
BOMB_visitdur3		The average duration per visit (seconds per visit) of all Bombus spp. (bumblebees) that were only interacting with the nectar of the flower for each visit to each site.
BOMB_visitdur4		The average duration per visit (seconds per visit) of all Bombus spp. (bumblebees) that were only interacting with the pollen/stamen of the flower for each visit to each site.
BOMB_visitdur5		The average duration per visit (seconds per visit) of all Bombus spp. (bumblebees) that were interacting with the pollen/stamen ad nectar simultaneously of the flower for each visit to each site.
Other_dur		The average duration (seconds) of all non-honeybee and non-bumblebee visits that were interacting with any part of the flower for each visit to each site.
Other_dur2		The duration (seconds) of all non-honeybees and non-bumblebees that were only interacting with the petals of the flower for each visit to each site.
Other_dur3		The duration (seconds) of all non-honeybees and non-bumblebees that were only interacting with the nectar of the flower for each visit to each site.
Other_dur4		The duration (seconds) of all non-honeybees and non-bumblebees that were only interacting with the pollen/stamen of the flower for each visit to each site.
Other_dur5		The duration (seconds) of all non-honeybees and non-bumblebees that were interacting with the pollen/stamen and nectar simultaneously of the flower for each visit to each site.
Other_visitdur		The average duration per visit (seconds per visit) of all non-honeybee and non-bumblebee visits that were interacting with any part of the flower for each visit to each site.
Other_visitdur2		The average duration per visit (seconds per visit) of all non-honeybees and non-bumblebees that were only interacting with the petals of the flower for each visit to each site.
Other_visitdur3		The average duration per visit (seconds per visit) of all non-honeybees and non-bumblebees that were only interacting with the nectar of the flower for each visit to each site.
Other_visitdur4		The average duration per visit (seconds per visit) of all non-honeybees and non-bumblebees that were only interacting with the pollen/stamen of the flower for each visit to each site.
Other_visitdur5		The average duration per visit (seconds per visit) of all non-honeybees and non-bumblebees that were interacting with the pollen/stamen and nectar simultaneously of the flower for each visit to each site.
VisitShannon		The Shannon diversity index of all pollinator visitors to flowers for each visit to a site.


File name:
Visitation_BySpp.csv

Description:
Data file used for the analysis of how each visitation metric differs based on pollinator group (i.e. "Genus"; honeybees, bumblebees, and other pollinators).
This file was modified from the VideoSummary_2016_forNosemaAnalysis.csv file so that each row has the all visitation metrics for a single pollinator group at
a single flower from a visit to each site. This allowed us to test each visitation metric 

Headings:
FlowerID	Unique ID of each flower in the study.
Site		Unique site code for each of the 6 sites (see Appendix S1: Table S1 for site names, dates visited, location, etc).
Year		Year of sample collection, all from 2016
Visit		Visit number to the site, each site was visited twice for this study. The first and second visit are indicated by “1” and “2”, respectively.
Date		Date that the video recording was conducted during the visit to that site.
StartHour	The hour that the video recording started.
totdur_min	Total duration of the video recording in minutes (continuous)
totdur_sec	Total duration of the video recording in seconds (continuous)
Genus		Code to identify the three groups that we were testing for differences in visitation behavior
		APME = Apis mellifera
		BOMB = Bombus spp
		Other = all other pollinators
visits		The number of honeybee, bumblebee, or other pollinators visits to a flower during the video recording duration (number of visits per 30 min).
rate		The rate of honeybee, bumblebee, or other pollinators visits to a flower per minute (number of visits per min).
dur		The duration (seconds) of honeybee, bumblebee, or other pollinators visits to a flower (interacting with any part).
dur2		The duration (seconds) of honeybee, bumblebee, or other pollinators that were only interacting with the petals of the flower.		
dur3		The duration (seconds) of honeybee, bumblebee, or other pollinators that were only interacting with the nectar of the flower.
dur4		The duration (seconds) of honeybee, bumblebee, or other pollinators that were only interacting with the pollen/stamen of the flower.
dur5		The duration (seconds) of honeybee, bumblebee, or other pollinators that were interacting with the pollen/stamen and nectar simultaneously of the flower.
visitdur	The duration per visit (seconds per visit) of honeybee, bumblebee, or other pollinators visits to a flower (interacting with any part).
visitdur2	The duration per visit (seconds per visit) of honeybee, bumblebee, or other pollinators that were only interacting with the petals of the flower.
visitdur3	The duration per visit (seconds per visit) of honeybee, bumblebee, or other pollinators that were only interacting with the nectar of the flower.
visitdur4	The duration per visit (seconds per visit) of honeybee, bumblebee, or other pollinators that were only interacting with the pollen/stamen of the flower.
visitdur5	The duration per visit (seconds per visit) of honeybee, bumblebee, or other pollinators that were interacting with the pollen/stamen and nectar simultaneously of the flower.
Lat			Latitude coordinate for the Site
Long			Longitude coordinate for the Site