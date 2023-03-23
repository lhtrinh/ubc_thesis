library(tidyverse)
library(haven)
library(readxl)
library(lubridate)
library(psych)


#######################################################
## LOAD DATA



#====================#
## BCGP data
#====================#

# load data sets
bc_raw <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_DATA_PT_Qx_RegistryInfo_BIO_2_0.sav")
colnames(bc_raw) <- tolower(colnames(bc_raw))
bc_hormones <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_Y2022_BCCR_incidence_breast_hormones.sav")
colnames(bc_hormones) <- tolower(colnames(bc_hormones))
bc_match_ind <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_MatchIndicator.csv")
colnames(bc_match_ind) <- tolower(colnames(bc_match_ind))


bc_hormones <- bc_hormones %>%
  mutate(er_status=tolower(bccr_er_status_final),
         pr_status=tolower(bccr_pr_status_final),
         her2_status=tolower(bccr_her2_status_final))


# remove duplicates in bc_raw and combine BCGP data
bc_dat <- bc_raw %>% filter(dup=="NO") %>%
  inner_join(bc_match_ind) %>%
  left_join(bc_hormones)