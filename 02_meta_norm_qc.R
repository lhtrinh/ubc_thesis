###############################################
###  QC FOR NORMALIZED DATA  ###
###############################################



## SETUP
library(tidyverse)
library(readr)

# load script of pre-defined functions
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load normalized metabolomics data set
full_ion_norm <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_normalized.csv")



# load questionnaire data
full_dat <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")





icc_breaks <- c(0, .4, .5, .75, .9, 1)




#===============================================#
# OPTION 1:
# 1. CALCULATE CV% AND ICC FOR EACH RUN TO COMPARE ITERATIONS FROM EACH RUN
# 2. TAKE MEAN ION LEVELS FROM EACH RUN & CALCULATE CV% AND ICC FOR BLIND REPLICATES ONLY

cv_norm_tubeid <- calc_cv_med(full_ion_norm, tubeid)
head(cv_norm_tubeid)
summary(cv_norm_tubeid$median_cv)
# good CV% between 2 samples of each run, no CV over 20%


icc_norm_tubeid <- calc_icc(full_ion_norm, tubeid)
summary(icc_norm_tubeid$icc)
proportions(table(cut(icc_norm_tubeid$icc, breaks=icc_breaks)))



#===============================================#

full_ion_norm_studyid <- full_ion_norm %>%
  group_by(studyid, tubeid, sampletype, cohort, dup) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()

full_norm_dup <- full_ion_norm_studyid %>%
  filter(dup=="YES")

cv_norm_1 <- calc_cv_med(full_norm_dup, studyid)
head(cv_norm_1)
summary(cv_norm_1$median_cv)
table(cv_norm_1$median_cv > .2)


icc_norm_1 <- calc_icc(full_norm_dup, studyid)
head(icc_norm_1)
summary(icc_norm_1$icc)
proportions(table(cut(icc_norm_1$icc, breaks=icc_breaks)))


qc_norm_1 <- cv_norm_1 %>%
  full_join(icc_norm_1)



qc_norm_1 %>%
  count(median_cv>.2, icc<.4)

qc_norm_1 %>%
  count(median_cv>.2, icc<.5)



qc_norm_1 %>%
  ggplot(aes(icc, after_stat(scaled))) +
  geom_density() +
  scale_x_continuous(breaks = seq(0,1,.1))







###############################################################
## OPTION 2: CALCULATE CV AND ICC FOR ALL STUDYID INCLUDING BOTH RUNS FOR EACH IONS

# cv_norm_2 <- calc_cv_med(full_ion_norm, studyid)
# summary(cv_norm_2)
# table(cv_norm_2$median_cv>.2)
#
#
# icc_norm_2 <- calc_icc(full_ion_norm, studyid)
# summary(icc_norm_2$icc)
# proportions(table(cut(icc_norm_2$icc, breaks=icc_breaks)))
#
#
# qc_norm_2 <- cv_norm_2 %>%
#   full_join(icc_norm_2)
#
#
# qc_norm_2 %>%
#   count(median_cv>.2, icc<.4)
#
#
# qc_norm_2 %>%
#   count(median_cv>.2, icc<.5)






###############################################################

## SELECT OPTION 1
## Remove metabolites with CV >.2 or ICC < .5
ions_to_rm_norm <- qc_norm_1$metabolite[qc_norm_1$median_cv>.2 | qc_norm_1$icc<.5]


# remove the metabolites and
# take the average metabolite levels for all duplicates
full_norm_avg <- full_ion_norm %>%
  filter(sampletype=="Sample") %>%
  group_by(studyid) %>%
  summarise(across(starts_with("ion"), mean))


full_norm_avg_qc <- full_norm_avg %>%
  select(-all_of(ions_to_rm_norm))

write_csv(full_ion_avg_cv,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_normalized_after_qc.csv")







###############################################################
###  PREPROCESSING METABOLOMICS DATA
###  NON-NORMALIZED DATA
###############################################################

## Test for normality and skewness before and after log-transformation

normality_norm_pre <- normality_skewness(full_norm_avg)
# how many metabolites have p>0.05 for each test (indicating possible normality)
normality_norm_pre %>% count(ks_pval>.05)
normality_norm_pre %>% count(shapiro_pval>.05)
normality_norm_pre %>% count(jarque_pval>.05)
describe(normality_norm_pre$skew)
proportions(table(cut(normality_norm_pre$skew, breaks=c(-2,-1,-.5,0,.5,1,2))))

# none normal


full_norm_log <- full_norm_avg %>%
  mutate(across(starts_with("ion"), log))

full_norm_log2 <- full_norm_avg %>%
  mutate(across(starts_with("ion"), log2))

full_norm_log10 <- full_norm_avg %>%
  mutate(across(starts_with("ion"), log10))

normality_norm_post <- normality_skewness(full_norm_log)
# how many metabolites have p>0.05 for each test (indicating possible)
normality_norm_post %>% count(ks_pval>.05)
normality_norm_post %>% count(shapiro_pval>.05)
normality_norm_post %>% count(jarque_pval>.05)
describe(normality_norm_post$skew)
proportions(table(cut(normality_norm_post$skew, breaks=c(-2,-1,-.5,0,.5,1,2))))

# log-transforming data only increased normality in about 21-27 metabolites
# we will keep the data not log-transformed, only scaled and mean-centered


###############################################################

# Pareto-scale ions
full_norm_sc <- full_ion_avg_cv %>%
  mutate(across(starts_with("ion"),
                function(x) (x-mean(x))/sqrt(sd(x))))



# center each pair by within-pair mean
# also merge with questionnare data
full_norm_centered <- full_norm_sc %>%
  inner_join(full_dat) %>%
  group_by(match_id) %>%
  mutate(across(starts_with("ion"),
                function(x) x-mean(x))) %>%
  ungroup()






###############################################################
write_csv(full_norm_centered,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized_preprocessed.csv")


