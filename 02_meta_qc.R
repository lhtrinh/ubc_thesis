######################################################
###   QC FOR FULL METABOLOMICS SAMPLE   ###
######################################################

library(tidyverse)

# load script of pre-defined functions
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load non-normalized metabolomics data set
full_ion_nonorm <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized.csv")


# load questionnaire data
full_dat <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")


# load ion hit frequency
ion_hit_freq_full <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/DATA_ionGapHits_NoPss.xlsx",
                          sheet="ions")


icc_breaks <- c(0, .4, .5, .75, .9, 1)




#===============================================#
# REMOVE IONS WITH LESS THAN 50% HIT FREQUENCY
dim(ion_hit_freq_full)
head(ion_hit_freq_full)


# there are 2,550 measurements from 1,275 samples
ion_freq2 <- ion_hit_freq_full %>%
  mutate(ion = paste("ion", ionIdx, sep="_"))


# calculate count of metabolites by hit frequency
for (i in seq(.5,1,.1)){
  prop_hit <- i
  ion_cnt <- length(which(ion_freq2$ionGapHits>=2550*prop_hit))
  print(paste("Detected in at least ",
            prop_hit*100,
            "% of samples: ",
            ion_cnt,
            " (",
            round(ion_cnt/1677*100, 1),
            "%)",
            sep=""))
}






thres <- .8
full_ion_trunc <- full_ion_nonorm %>%
  select(-all_of(ion_freq2$ion[ion_freq2$ionGapHits<2550*thres]))

dim(full_ion_nonorm)
dim(full_ion_trunc)


#===============================================#
# OPTION 1:
# 1. CALCULATE CV% AND ICC FOR EACH RUN TO COMPARE ITERATIONS FROM EACH RUN
# 2. TAKE MEAN ION LEVELS FROM EACH RUN & CALCULATE CV% AND ICC FOR BLIND REPLICATES ONLY

# cv_nonorm_tubeid <- calc_cv_med(full_ion_nonorm, tubeid)
# head(cv_nonorm_tubeid)
# summary(cv_nonorm_tubeid$median_cv)
# # good CV% between 2 samples of each run, no CV over 20%
#
#
# icc_nonorm_tubeid <- calc_icc(full_ion_nonorm, tubeid)
# summary(icc_nonorm_tubeid$icc)
# proportions(table(cut(icc_nonorm_tubeid$icc, breaks=icc_breaks)))



#===============================================#

full_nonorm_studyid <- full_ion_trunc %>%
  group_by(studyid, tubeid, sampletype, cohort, plate, dup) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()


full_nonorm_dup <- full_nonorm_studyid %>%
  filter(dup=="YES")

cv_nonorm_1 <- calc_cv_med(full_nonorm_dup, studyid)
head(cv_nonorm_1)
summary(cv_nonorm_1$median_cv)
table(cv_nonorm_1$median_cv > .2)


icc_nonorm_1 <- calc_icc(full_nonorm_dup, studyid)
head(icc_nonorm_1)
summary(icc_nonorm_1$icc)
proportions(table(cut(icc_nonorm_1$icc, breaks=icc_breaks)))



qc_nonorm_1 <- cv_nonorm_1 %>%
  full_join(icc_nonorm_1)

qc_nonorm_1 %>%
  count(median_cv>.2, icc<.4)

qc_nonorm_1 %>%
  count(median_cv>.2, icc<.5)


qc_nonorm_1 %>%
  ggplot(aes(icc, after_stat(scaled))) +
  geom_density() +
  scale_x_continuous(breaks = seq(0,1,.1))




###############################################################

## SELECT OPTION 1
## Remove metabolites with CV >.2 or ICC < .5
ions_to_rm_nonorm <- qc_nonorm_1$metabolite[qc_nonorm_1$median_cv>.2 | qc_nonorm_1$icc<.4]


# remove the metabolites and
# take the average metabolite levels for all duplicates
full_nonorm_avg <- full_ion_trunc %>%
  filter(sampletype=="Sample") %>%
  group_by(studyid, tubeid, sampletype, cohort, plate, dup) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()

full_nonorm_avg_qc <- full_nonorm_avg %>%
  select(-all_of(ions_to_rm_nonorm))


write_csv(full_ion_avg_qc,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized_after_qc.csv")







###############################################################
###  PREPROCESSING METABOLOMICS DATA
###  NON-NORMALIZED DATA
###############################################################

## Test for normality before and after log-transformation

normality_nonorm_pre <- normality_skewness(full_nonorm_avg)
# how many metabolites have p>0.05 for each test (indicating possible normality)
normality_nonorm_pre %>% count(ks_pval>.05)
normality_nonorm_pre %>% count(shapiro_pval>.05)
normality_nonorm_pre %>% count(jarque_pval>.05)
summary(normality_nonorm_pre$skew)

# none normal


full_nonorm_log <- full_nonorm_avg %>%
  mutate(across(starts_with("ion"), log))

full_nonorm_log2 <- full_nonorm_avg %>%
  mutate(across(starts_with("ion"), log2))

normality_nonorm_post <- normality_skewness(full_nonorm_log)
# how many metabolites have p>0.05 for each test (indicating possible)
normality_nonorm_post %>% count(ks_pval>.05)
normality_nonorm_post %>% count(shapiro_pval>.05)
normality_nonorm_post %>% count(jarque_pval>.05)
summary(normality_nonorm_post$skew)

# log-transforming data only increased normality in about 21-27 metabolites
# we will keep the data not log-transformed, only scaled and mean-centered


###############################################################

# Pareto-scale ions
full_nonorm_sc <- full_nonorm_avg_qc %>%
  mutate(across(starts_with("ion"),
                function(x) (x-mean(x))/sqrt(sd(x))))



# Within-pair mean centered
full_nonorm_centered <- full_nonorm_sc %>%
  inner_join(subset(full_dat, select=c("studyid", "match_id", "gp"))) %>%
  group_by(match_id) %>%
  mutate(across(starts_with("ion"),
                function(x) x-mean(x))) %>%
  ungroup()


###############################################################

# QC again

full_nonorm_dup2 <- full_nonorm_centered %>%
  filter(dup=="YES")

cv_nonorm_2 <- calc_cv_med(full_nonorm_dup2, studyid)
head(cv_nonorm_2)
summary(cv_nonorm_2$median_cv)
table(cv_nonorm_2$median_cv > .2)


icc_nonorm_2 <- calc_icc(full_nonorm_dup2, studyid)
head(icc_nonorm_2)
summary(icc_nonorm_2$icc)
proportions(table(cut(icc_nonorm_2$icc, breaks=icc_breaks)))



qc_nonorm_2 <- cv_nonorm_2 %>%
  full_join(icc_nonorm_2)

qc_nonorm_2 %>%
  count(median_cv>.2, icc<.4)

qc_nonorm_2 %>%
  count(median_cv>.2, icc<.5)





###############################################################
write_csv(full_nonorm_centered,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized_preprocessed.csv")
