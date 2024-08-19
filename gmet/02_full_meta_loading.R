#################################
# Full metabolomics data
#################################



### Load packages
library(tidyverse)
library(readxl)
library(readr)




# load metabolite data
full_ion_intensities <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/spreadsheets/DATA.xlsx",
  sheet="intensities1")
full_ion_intensities %>% head()

full_injection <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/spreadsheets/DATA.xlsx",
  sheet="injections")
full_injection %>% head()


# load crosswalk data
atp_id_cw <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/atp_id_crosswalk.csv",
                      col_types=cols(participantkey=col_character(), barcode=col_character()))

id_cw_dup <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/id_crosswalk_dup_info.csv",
                      col_types = cols(participantkey=col_character(), barcode=col_character()))




# # load full questionnaire data
# full_dat <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")







#======================================================================#
# add "ion" in front of each metabolite ID
# transpose full ion intensities so rows = samples, columns = metabolites
# extract sample ID and change to numeric form
full_ion_t <- full_ion_intensities %>%
  mutate(ion=paste("ion", ionIdx, sep="_")) %>% # add "ion" in front of each metabolite ID
  select(-ionIdx, -ionMz, -ionActive) %>% # remove irrelevant columns
  pivot_longer(!ion, names_to = "sampleid") %>%
  pivot_wider(names_from = ion,
              values_from = value) %>% # transpose data
  mutate(sampleid = as.numeric(str_replace_all(sampleid, "[^[:digit:]]", ""))) # change sample ID to numeric form






#======================================================================#
# separate sample IDs and cohort
# also change all ATP sample IDs to corresponding studyid
### 1.  If the cohort is ATP, remove the string "ATP_" in the front
###     and add an underscore before "a" or "b" to match BCGP strings
### 2.
full_inj2 <- full_injection %>%
  mutate(tubeid_iter = ifelse(dsUser3=="ATP",
                              gsub('^(.{4})(.*)(a|b)$', '\\2_S1\\3', dsCode),
                              dsCode)) %>%
  mutate(tubeid = str_replace_all(tubeid_iter, "(a|b)$", "")) %>%
  separate(tubeid,
           into=c("sampleid", "iter"),
           sep="_",
           convert=TRUE,
           remove=FALSE) %>%
  separate(dsUser4,
           into=c("x", "cohort", "plate"),
           sep="_",
           convert=TRUE,
           extra="merge") %>%
  mutate(plate=as.numeric(str_replace_all(plate, "[^[:digit:]]", "")),
         cohort=dsUser3,
         sampletype=dsPert)



# join injection data with crosswalk
## if sample is in the ATP cohort, then
## replace sampleid with participantkey from the crosswalk data set
full_inj3 <- full_inj2 %>%
  left_join(id_cw_dup,
            by=join_by(sampleid==barcode),
            keep=TRUE) %>%
  mutate(studyid=case_when(
    cohort=="BCGP" ~ sampleid,
    cohort=="ATP" ~ participantkey,
    TRUE ~ sampleid
  ))
dim(full_inj3)



# also label QC samples as "dup" since we are using them for ICC and CV
full_inj4 <- full_inj3 %>%
  mutate(dup=ifelse(sampletype=="QC", "YES", dup))

full_inj4 %>% filter(dup=="YES") %>%
  count(studyid, cohort) %>%
  count(cohort)





full_inj <- full_inj4



#################################################################
# join injection and ion intensities
# also remove any samples that are not in the questionnaire data
full_ion <- full_ion_t %>%
  full_join(full_inj,
            by=c("sampleid"="dsIdx"),
            keep=TRUE) %>%
  select(dsIdx, dsCode, studyid, tubeid, iter, dup, sampletype, cohort, plate, contains("ion")) %>%
  filter(sampletype=="QC" |
           sampletype=="Sample" & studyid %in% id_cw_dup$participantkey)

# double check number of study samples
n_distinct(full_ion$studyid[full_ion$sampletype=="Sample"])




#################################################################
# save data to local folder
write_csv(
  full_ion, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized.csv")

