#################################
# Full metabolomics data
#################################



### Load packages
library(tidyverse)
library(readxl)
library(readr)




# load metabolite data
ion_matrix <- readxl::read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/DATA_blood_exposeome.xlsx",
  sheet="ion_matrix")
ion_matrix %>% head()

sample_info <- readxl::read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/DATA_blood_exposeome.xlsx",
  sheet="samples")
sample_info %>% head()


# load crosswalk data
atp_id_cw <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/atp_id_crosswalk.csv",
                      col_types=cols(studyid=col_character(), barcode=col_character()))

id_cw <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/id_crosswalk_dup_info.csv",
                      col_types = cols(studyid=col_character(), barcode=col_character()))





#======================================================================#
# Transpose metabolomics ion matrix
ion_matrix %>% head()

ion_t <- ion_matrix %>%
  # rename first column to ionIdx
  rename(ionIdx = `...1`) %>%
  # remove the blank line
  filter(dsSampleCode!='ionMz') %>%
  # add the string "ion_" in front of ion ID
  mutate(ion=ifelse(ionIdx!='dsIdx', paste("ion", ionIdx, sep="_"), ionIdx)) %>%
  select(-ionIdx, -dsSampleCode) %>%
  # transpose data so rows=samples and columns=metabolites
  pivot_longer(!ion, names_to="dsSampleCode") %>%
  mutate(ion=ifelse(is.na(ion), 'dsIdx', ion)) %>%
  pivot_wider(names_from=ion,
              values_from=value) %>%
  # duplicate sample IDs for all second measurements
  mutate(dsSampleCode=ifelse(str_detect(dsSampleCode, "PB|ATP"), dsSampleCode, NA)) %>%
  fill(dsSampleCode, .direction="down")

ion_t %>% head()
ion_t %>% count(dsSampleCode)

n_distinct(ion_t$dsSampleCode)
# all good



#======================================================================#
# remove "ATP_", "_S1" and "_S2" to get barcode information and join with crosswalk data
sample2 <- sample_info %>%
  rename(cohort=dsUser3) %>%
  mutate(barcode=str_replace_all(dsSampleCode, "ATP_|_S1|_S2", "")) %>%
  right_join(id_cw)




#################################################################
# join injection and ion intensities
# also remove any samples that are not in the questionnaire data
full_ion <- ion_t %>%
  right_join(sample2) %>%
  select(dsIdx, studyid, barcode, dsPlate, dsPlateInj, dup, cohort, contains("ion"))

# double check number of study samples
n_distinct(full_ion$studyid)

# double check number of duplicates
full_ion %>% count(studyid) %>% count(n)



#################################################################
# save data to local folder
write_csv(
  full_ion,
  "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/exposome/meta_non_normalized.csv")

