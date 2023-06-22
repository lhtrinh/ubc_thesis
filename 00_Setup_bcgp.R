#######################################################
# In this script:
## 1. Load packages
## 2. Load BCGP data sets
## 3. Edit names and combine data
#######################################################




#========================================#
## 1. Load packages
#========================================#
library(tidyverse)
library(haven)
library(readxl)
library(lubridate)
library(psych)
library(moments)





#========================================#
## 2. Load BCGP data
#========================================#


# load data sets and change all column names to lowercase

# questionnaire data
bc_raw <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_DATA_PT_Qx_RegistryInfo_BIO_2_0.sav")
colnames(bc_raw) <- tolower(colnames(bc_raw))

# hormone receptors/HER2 status data for cases
bc_hormones <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_Y2022_BCCR_incidence_breast_hormones.sav")
colnames(bc_hormones) <- tolower(colnames(bc_hormones))

# matched pair IDs
bc_match_ind <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_MatchIndicator.csv")
colnames(bc_match_ind) <- tolower(colnames(bc_match_ind))





#========================================#
## 3. Edit variable names and combine data
#========================================#


# remove one participant from BCGP (selected in error)
bc_raw2 <- bc_raw %>%
  filter(studyid!="PB000429")



# shorten variable names for HR status
bc_hormones <- bc_hormones %>%
  mutate(er_status=tolower(bccr_er_status_final),
         pr_status=tolower(bccr_pr_status_final),
         her2_status=tolower(bccr_her2_status_final))



# add "BCGP" to match ID to distinguish from ATP match IDs
bc_match_ind$match_id <- paste("bcgp", bc_match_ind$match_id, sep="_")





# remove duplicates in bc_raw and combine BCGP data
bc_dat <- bc_raw2 %>%
  filter(dup=="NO") %>%
  inner_join(bc_match_ind) %>%
  left_join(bc_hormones)





# remove one studyid without a matched control
bc_dat %>% count(match_id) %>% count(n) # 1 studyid with no match
unmatched <- bc_dat %>% count(match_id) %>% filter(n==1)
bc_dat <- bc_dat %>% filter(!match_id %in% unmatched$match_id)


bc_dat$cohort <- "bcgp"
dim(bc_dat)



# remove pairs where at least one had no baseline questionnaire
bc_dat %>% count(is.na(adm_qx_completion)) #10 participants with no baseline qx

no_baseline <- bc_dat %>%
  filter(is.na(adm_qx_completion)) %>%
  count(match_id)


bc_dat <- bc_dat %>%
  filter(!match_id %in% no_baseline$match_id)



dim(bc_dat) #696
n_distinct(bc_dat$match_id)


###############################################################################