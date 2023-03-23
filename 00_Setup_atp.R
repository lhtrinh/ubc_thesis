library(tidyverse)
library(haven)
library(readxl)
library(lubridate)
library(psych)





#====================#
## ATP data
#====================#

# load data sets
atp_raw <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_HFRQ_BL_CPAC.csv")
colnames(atp_raw) <- tolower(colnames(atp_raw))
atp_pm_measure <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_PM_BL_CPAC.csv")
colnames(atp_pm_measure) <- tolower(colnames(atp_pm_measure))
atp_hormones <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_ACR_Breast_DX_ER_PR.csv")
colnames(atp_hormones) <- tolower(colnames(atp_hormones))



# #atp unharmonized data
atp_cases_bio <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_ACR.csv")
colnames(atp_cases_bio) <- tolower(colnames(atp_cases_bio))





############################################################

# remove duplicates
atp_hormones <- atp_hormones %>%
  filter(acr_cancer_site_aggregate=="Breast") %>%
  arrange(participantkey, acr_malignancy_number)

atp_hormones <- atp_hormones %>%
  mutate(er_status=acr_estrogen_receptor_assay,
         pr_status=acr_progesterone_receptor_assay,
         her2_status=acr_her2_assay) %>%
  mutate(across(ends_with("_status"),
                function(x) ifelse(grepl("Positive", x), "positive",
                                   ifelse(grepl("Negative", x), "negative", "unknown"))))

atp_hormones_new <- atp_hormones[!duplicated(atp_hormones$participantkey),]





atp_cases_bio <- atp_cases_bio %>%
  filter(acr_cancer_site_aggregate=="Breast") %>%
  arrange(participantkey, acr_diagnosis_year)
atp_cases_new <- atp_cases_bio[!duplicated(atp_cases_bio$participantkey),]




# combine ATP data
temp <- atp_hormones_new %>%
  inner_join(atp_cases_new)


atp_dat <- atp_raw %>%
  inner_join(atp_pm_measure) %>%
  left_join(atp_hormones_new) %>%
  left_join(atp_cases_new)




# remove 2 participants that are no longer eligible
# atp_dat <- atp_dat %>%
#   filter(!participantkey %in% c(210503833, 210524604))



# add control and case labels
atp_dat$gp <- "CNTL"
atp_dat$gp[atp_dat$participantkey %in% atp_hormones_new$participantkey] <- "CASE"



# change date form
atp_dat$adm_qx_completion <- as.Date(atp_dat$adm_qx_completion)



########################################################################
#====================#
# Replace all -7 with NAs
#====================#
atp_dat[atp_dat==-7] <- NA


#====================#
# Create new variables from data
#====================#

atp_dat$cohort <- "atp"
atp_dat$studyid <- as.character(atp_dat$participantkey)


# rename age at diagnosis to match atp
atp_dat <- rename(atp_dat, age_at_diagnosis=acr_age_at_diagnosis)
