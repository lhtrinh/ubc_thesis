###
library(tidyverse)
library(haven)
library(readxl)
library(lubridate)
library(psych)
library(moments)
library(readr)




#===========================================================#
# load full cohort data for imputation
bc_all <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/Multiple imputation task/Forimputation_BL_HLQ_PM_clean.csv")
colnames(bc_all) <- tolower(colnames(bc_all))
head(bc_all)

bc_all <- bc_all %>% select(-`...1`)

# load crosswalk
bc_cw <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/Multiple imputation task/XWALK_Ly_BCGP2641.csv")
colnames(bc_cw) <- tolower(colnames(bc_cw))
head(bc_cw)
bc_cw <- bc_cw %>% select(-...1)



# load original raw data
bc_raw_og <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_DATA_PT_Qx_RegistryInfo_BIO_2_0.sav")
colnames(bc_raw_og) <- tolower(colnames(bc_raw_og))

bc_raw_og <- bc_raw_og[!duplicated(bc_raw_og$studyid),]


#===========================================================#


setdiff(bc_cw$study_id, bc_all$study_id)
bc_cw %>% filter(is.na(study_id)) #one participant did not have a studyid
# full_dat %>% filter(studyid=="PB000521")
# can confirm that this participant is not in the cleaned data set
# all good



bc_proj <- bc_all %>%
  inner_join(bc_cw) %>%
  select(-study_id) %>%
  rename(studyid=bcgp2641_bhatti)



add_cols <- setdiff(colnames(bc_raw_og), colnames(bc_proj))

bc_raw_og2 <- bc_raw_og %>%
  arrange(studyid) %>%
  select(studyid, all_of(add_cols))

bc_raw_og3 <- bc_raw_og2[!duplicated(bc_raw_og2$studyid),]

bc_raw_og3 %>% filter(studyid=="PB000251")


bc_raw <- bc_proj %>%
  right_join(bc_raw_og3) %>%
  arrange(studyid)



#===========================================================#
write_csv(bc_raw,
          file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP_data.csv")
