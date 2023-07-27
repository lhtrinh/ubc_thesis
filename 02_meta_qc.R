######################################################
###   QC FOR FULL METABOLOMICS SAMPLE   ###
######################################################

library(tidyverse)

# load script of pre-defined functions
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load non-normalized metabolomics data set
full_ion_plate <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized.csv")


# # QC for samples only
# full_ion_sample <- full_ion_plate %>%
#   filter(sampletype=="Sample")



#===============================================#
# CALCULATE CV% AND ICC FOR EACH RUN
# to compare iterations a and b of each tube

cv_by_tubeid <- calc_cv_med(full_ion_plate, tubeid)
head(cv_by_tubeid)
summary(cv_by_tubeid$median_cv)
# good CV% between 2 samples of each run, no CV over 20%


icc_by_tubeid <- calc_icc(full_ion_plate, tubeid)
summary(icc_by_tubeid$icc)



#===============================================#
# take mean intensity from each run
full_ion_avg <- full_ion_plate %>%
  group_by(studyid, tubeid, dup) %>%
  summarise(across(contains("ion"), mean)) %>%
  ungroup()


full_ion_avg_dup <- full_ion_avg %>%
  filter(dup=="YES")


#===============================================#
# CALCULATE CV% AND ICC FOR BLIND REPLICATES
cv_by_studyid <- calc_cv_med(full_ion_avg_dup, studyid)
head(cv_by_studyid)
summary(cv_by_studyid$median_cv)
table(cv_by_studyid$median_cv > 20)


icc_by_studyid <- calc_icc(full_ion_avg_dup, studyid)
icc_by_studyid
summary(icc_by_studyid$icc)
table(icc_by_studyid$icc>=.75)










###############################################################

## OPTION 2: CALCULATE CV AND ICC FOR ALL STUDYID

cv_2 <- calc_cv_med(full_ion_plate, studyid)
summary(cv_2)
table(cv_2$median_cv>20)


icc_2 <- calc_icc(full_ion_plate, studyid)
summary(icc_2)
