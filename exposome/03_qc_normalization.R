# NORMALIZATION AND QC PIPELINE ####



#===============================================#
## Summary ####
#==============#
## This code details the following steps:
### 1. Filter metabolites that were missing in more than 50% of the experimental samples
### 2. Calculate CV and ICC for non-normalized samples
### 3. PCA plots for non-normalized samples (by plate and cohort)
### 4. Normalization methods
### 5. Select method
### 6. Calculate mean by study ID for each duplicated sample
### 7. Mean-centering
#===============================================#






#===============================================#
## 0. Setup ####
#==============#

library(tidyverse) # data handling
library(readxl)
library(readr)
library(factoextra) # PCA plots
library(gridExtra) # to fit multiple plots on the same page
library(vsn) # variance stabilization normalization (VSN)
library(affy) # cubic spline
library(NCStats) # cubic spline
library(preprocessCore) # quantile
library(NormalizeMets) # cyclic LOESS using QC samples (QC-RLSC)
library(pmp) # cubic spline using QC samples




# load script of pre-defined functions for ICC and CV
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/exposome/00_functions.R")


# load non-normalized metabolomics data set
full_ion_nonorm <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/exposome/meta_non_normalized.csv")


# load questionnaire data
full_dat <- import_survey()


# load ion hit frequency for missingness information
ion_hit_freq_full <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/DATA_blood_exposeome.xlsx",
  sheet="ions")





#===============================================#
## 1. Filtering ####
#### Remove metabolites with more than 50% missing
#===============================================#


dim(ion_hit_freq_full)
head(ion_hit_freq_full)


ion_freq2 <- ion_hit_freq_full %>%
  filter(!is.na(ionIdx)) %>%
  mutate(ion = paste("ion", ionIdx, sep="_"))

ion_freq2


# calculate count of metabolites by hit frequency
for (i in seq(.5,1,.1)){
  prop_hit <- i
  ion_cnt <- length(which(ion_freq2$ionGapHits>=max(ion_freq2$ionGapHits)*prop_hit))
  print(paste("Detected in at least ",
              prop_hit*100,
              "% of samples: ",
              ion_cnt,
              " (",
              round(ion_cnt/nrow(ion_freq2)*100, 1),
              "%)",
              sep=""))
}



# set cutoff threshold
thres <- .5
high_missing_ions <- ion_freq2$ion[ion_freq2$ionGapHits<2550*thres]

ion_trunc <- full_ion_nonorm %>%
  select(-all_of(high_missing_ions))

dim(full_ion_nonorm)
dim(ion_trunc)





# calculate mean of each consecutive pair of samples (sampleid)
# so that each tube id only has one measurement for each metabolite
ion_per_sample <- ion_trunc %>%
  group_by(pick(-c(dsIdx, starts_with("ion")))) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()

dim(ion_per_sample)




#===============================================#
## Average all dups to get one measurement

full_studyid <- ion_per_sample %>%
  group_by(studyid) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()


#===============================================#
## Center by within-pair mean ####
# subtract all metabolites by within-pair mean

full_centered <- full_studyid %>%
  full_join(full_dat[,colnames(full_dat) %in% c("studyid", "match_id")]) %>%
  group_by(match_id) %>%
  mutate(across(starts_with("ion"), function(x) x-mean(x))) %>%
  ungroup()


#===============================================#
## Scaling ####
full_sc <- full_centered %>%
  mutate(across(starts_with("ion"), function(x) x/sd(x)))

# save to csv
write_csv(full_sc,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data//exposome/unnormalized_metabolomics_scaled.csv")



#===============================================#
## 2. Calculate CV and ICC ####
##### for non-normalized samples
#==============#



# filter only samples with duplicates
ion_dup <- ion_per_sample %>%
  filter(dup=="YES")



## Calculate median coefficient of variation
cv_nonorm_1 <- calc_cv_med(ion_dup, studyid)

# summarize median CV
round(quantile(cv_nonorm_1$median_cv, probs=c(0,.05, .5, .95, 1)),3)

# find metabolites with high CV (>20%)
high_cv_ions <- cv_nonorm_1 %>%
  filter(median_cv>.2) %>%
  pull(metabolite)

high_cv_ions


# remove high CV ions from data
ion_cv <- ion_per_sample %>%
  select(-all_of(high_cv_ions))




#===============================================#
## 3. Non-normalized PCA plots ####
##### for non-normalized samples (by plate and cohort)
#==============#


dat_nonorm_pca <- ion_cv %>% select(starts_with("ion"))

pca_nonorm <- prcomp(dat_nonorm_pca)

# plot by plate
fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  ion_cv$dsPlate,
             alpha=.5,
             addEllipses = TRUE)+
  labs(title="PCA before normalization")

# plot by cohort
fviz_pca_ind(pca_nonorm,
             label="none",
             alpha=.5,
             habillage = ion_cv$cohort,
             palette=c("blue", "red"),
             addEllipses = TRUE) +
  labs(title="PCA before normalization")






#===============================================#
## 4. Normalization methods ####
#==============#

#===================#
### cubic splines ####
library(NCStats)
library(affy)
other_df <- ion_cv %>%
  select(!starts_with("ion")) %>%
  as.data.frame()
head(other_df)

ion_df <- ion_cv %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()
head(ion_df)

set.seed(56789)
spline_dat <- t(normalize.qspline(
  ion_df,
  samples=0.05,
  target=apply(ion_df,1,geomean),
  verbose=TRUE
))
colnames(spline_dat) <- rownames(ion_df)

ion_df[1:10,1:10]
spline_dat[1:10,1:10]


full_norm_spline <- cbind(other_df, spline_dat) %>% as_tibble()
full_norm_spline[1:10,1:10]





#===============================================#
## Select method ####

# cubic spline
full_norm <- full_norm_spline

# # median (gmet's version)
# full_norm <- full_norm_gmet

# no normalization
# full_norm <- ion_cv



#===============================================#
### EVALUATE: PCA
dat_norm_pca <- full_norm_spline %>% select(starts_with("ion"))
pca_norm <- prcomp(dat_norm_pca)

# plot by cohort
fviz_pca_ind(pca_norm,
             label="none",
             alpha=.5,
             habillage = full_norm_spline$cohort,
             palette=c("blue", "red"),
             addEllipses = TRUE,
             legend.title="Cohort") +
  labs(title="PCA after normalization")


fviz_pca_ind(pca_norm,
             label="none",
             alpha=.5,
             habillage = full_norm_spline$plate,
             addEllipses = TRUE,
             legend.title="Plate") +
  labs(title="PCA after normalization")



#===========================================#
## RLE plots ####
### EVALUATE: RLE

# randomly sample 200 samples
samp_rows <- sample(1:nrow(ion_per_sample), 200, replace=FALSE)

# before normalization
dat_nonorm_rle <- ion_per_sample[samp_rows,] %>%
  mutate(across(starts_with("ion"),
                function(x) x-median(x))) %>%
  select(studyid, cohort, starts_with("ion")) %>%
  pivot_longer(-c(studyid, cohort), names_to="ion")

dat_nonorm_rle %>%
  ggplot()+
  geom_boxplot(aes(studyid, value, fill=cohort), outlier.size=.001) +
  geom_hline(yintercept=0) +
  ylim(-20000, 20000) +
  labs(x="sample", y="relative log expression") +
  theme(axis.text.x=element_blank())



# after normalization
dat_norm_rle <- full_norm_spline[samp_rows,] %>%
  # group_by(cohort) %>%
  mutate(across(starts_with("ion"),
                function(x) x-median(x))) %>%
  select(studyid, cohort, starts_with("ion")) %>%
  pivot_longer(-c(studyid, cohort), names_to="ion")

dat_norm_rle %>%
  ggplot()+
  geom_boxplot(aes(studyid, value, fill=cohort), outlier.size=.01) +
  geom_hline(yintercept=0) +
  ylim(-20000, 20000) +
  labs(x="sample", y="relative log expression") +
  theme(axis.text.x=element_blank())







#===============================================#
## Average all dups to get one measurement

full_studyid <- full_norm %>%
  group_by(studyid) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()


#===============================================#
## Center by within-pair mean ####
# subtract all metabolites by within-pair mean

full_centered <- full_studyid %>%
  full_join(full_dat[,colnames(full_dat) %in% c("studyid", "match_id")]) %>%
  group_by(match_id) %>%
  mutate(across(starts_with("ion"), function(x) x-mean(x))) %>%
  ungroup()


#===============================================#
## Scaling ####
full_sc <- full_centered %>%
  mutate(across(starts_with("ion"), function(x) x/sd(x)))

# save to csv
write_csv(full_sc,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data//exposome/normalized_metabolomics_scaled.csv")











#
#
#
# #==========================================================================#
#
# ## Appendix: Other normalization methods ####
#
# #### other methods considered by not reported
# #===================#
#
#
#
#
# #===================#
#
# ### GMet's method ####
# ## divide by sample median and multiply by mean across all medians
# full_norm_gmet <- ion_cv %>%
#   rowwise() %>%
#   mutate(m=median(c_across(starts_with("ion")))) %>%
#   ungroup() %>%
#   group_by(cohort) %>%
#   mutate(mm=mean(m)) %>%
#   ungroup() %>%
#   mutate(across(starts_with("ion"),
#                 function(x) x*mm/m)) %>%
#   select(-m, -mm)
#
#
#
#
# #===================#
# ### linear baseline normalization ####
# other_df <- ion_cv %>%
#   select(!starts_with("ion")) %>%
#   as.data.frame()
#
# ion_df <- ion_cv %>%
#   select(starts_with("ion")) %>%
#   as.data.frame() %>%
#   t()
#
# baseline <- apply(ion_df, 1, median)
# baseline_mean <- mean(baseline)
# samp_means <- apply(ion_df, 2, mean)
# linear_fct <- baseline_mean/samp_means
#
# linear_norm_df <- t(ion_df)*linear_fct
#
# ion_df[1:10, 1:10]
# linear_norm_df[1:10, 1:10]
#
# full_norm_linear <- cbind(other_df, linear_norm_df) %>% as_tibble()
# full_norm_linear[1:10,1:10]
#
#
#
#
#
#
# #===================#
# ### non-linear baseline normalization (Li-Wong) ####
# library(affy)
#
# data <- ion_cv %>%
#   select(starts_with("ion")) %>%
#   t()
# other_df <- ion_cv %>%
#   select(!starts_with("ion"))
#
# #---First step: Find baseline sample
# average.intensity <- apply(data,2,mean)
# median.number <- round(ncol(data)/2 + 0.1)				#R has an add way of rounding.
# #the additional 0.1 ensures that it rounds properly
# ordering <- order(average.intensity)
# median.sample.number <- ordering[median.number]
# median.sample <- data[,median.sample.number]
#
# #---Apply normalization
# liwong.data=vector()
# for(i in 1:ncol(data)){
#   liwong.model <- normalize.invariantset(data=data[,i],
#                                          ref=median.sample,
#                                          prd.td=c(0.003,0.007))			#the threshold of the rank-invariant set might need to be adjusted from case to case
#   liwong.sample <- predict(liwong.model$n.curve$fit,		#chosen such that the rank-invariant set it sufficiently large
#                            data[,i])
#   liwong.data <- cbind(liwong.data,liwong.sample$y)
# }
#
# full_norm_liwong <- cbind(other_df, t(liwong.data)) %>% as_tibble()
#
#
#
#
#
#
#
# #===================#
# ### quantile normalization ####
# dat_ion <- ion_cv %>%
#   select(starts_with("ion"))
#
#
# dat_ion_quantile <- as.data.frame(t(normalize.quantiles(t(as.matrix(dat_ion)))))
# colnames(dat_ion_quantile) <- colnames(dat_ion)
#
# dat_no_ion <- ion_cv %>%
#   select(!starts_with("ion"))
#
# full_norm_qtl <- cbind(dat_no_ion, dat_ion_quantile) %>% as_tibble()
# full_norm_qtl[1:10,1:10]
#
#
#
#
#
#
#
#
#
# #===================#
# ### probabilistic quotient normalization (PQN) ####
#
# other_df <- ion_cv %>%
#   select(!starts_with("ion")) %>%
#   as.data.frame()
# ion_df <- ion_cv %>%
#   select(starts_with("ion")) %>%
#   as.data.frame() %>%
#   t()
#
# ion_df[1:10, 1:10]
#
# ref <- apply(ion_df, 1, median)
# quotient <- ion_df/ref
# quotient_med <- apply(quotient, 2, median)
# pqn_df <- t(ion_df)/quotient_med
#
# pqn_df[1:10, 1:10]
#
# full_norm_pqn <- cbind(other_df, pqn_df) %>% as_tibble()
# full_norm_pqn[1:10,1:10]
#
#
#
# #===================#
# ### probabilistic quotient normalization (PQN) with QC ####
#
# other_df <- ion_cv %>%
#   select(!starts_with("ion")) %>%
#   as.data.frame()
#
# ion_df <- ion_cv %>%
#   select(starts_with("ion")) %>%
#   as.data.frame()
# ion_df[1:10, 1:10]
#
# ion_qc_df <- nonorm_pss %>%
#   select(starts_with("ion")) %>%
#   as.data.frame()
#
# ref <- apply(ion_qc_df, 2, median)
# quotient <- ion_df/ref
# quotient_med <- apply(quotient, 1, median)
# pqn_df <- ion_df/quotient_med
#
# pqn_df[1:10, 1:10]
#
# full_norm_pqn_pss <- cbind(other_df, pqn_df) %>% as_tibble()
# full_norm_pqn_pss[1:10,1:10]
#
#
#
# # post QC
# qc_pqn_pss <- norm_qc(full_norm_pqn_pss)
#
# round(quantile(qc_pqn_pss$median_cv, probs=c(0,.05, .5, .95, 1)),3)
# qc_pqn_pss %>% count(median_cv>.2)
#
# round(quantile(qc_pqn_pss$icc, probs=c(0,.05, .5, .95, 1)), 3)
# table(cut(qc_pqn_pss$icc, breaks=icc_breaks, right=FALSE))
#
#
#
#
# #######################################################################!!!
#
#
#
#
#
#
#
#
#
# #===================#
# ### cubic splines using QC samples (version 11/01/2023) ####
# qc_spec <- ion_per_sample %>%
#   filter(cohort=="pSS") %>%
#   select(-all_of(high_cv_ions$metabolite)) %>%
#   summarise(across(starts_with("ion"), geomean)) %>%
#   as_vector()
#
#
# set.seed(30340)
# qc_spline_dat <- t(normalize.qspline(
#   ion_df,
#   samples=.05,
#   target=qc_spec,
#   verbose=TRUE
# ))
#
# colnames(qc_spline_dat) <- rownames(ion_df)
# qc_spline_dat[1:10,1:10]
#
#
# full_norm_qc_spline <- cbind(other_df, qc_spline_dat) %>% as_tibble()
# full_norm_qc_spline[1:10,1:10]
#
#
#
#
# dat_norm_pca <- full_norm_qc_spline %>% select(starts_with("ion"))
# pca_norm <- prcomp(dat_norm_pca)
#
# # plot by cohort
# fviz_pca_ind(pca_norm,
#              label="none",
#              alpha=.5,
#              habillage = full_norm_qc_spline$cohort,
#              palette=c("blue", "red"),
#              addEllipses = FALSE)
#
#
#
#
#
#
# #===================#
# ### cubic splines using QC samples (QC-RSC) ####
# # BiocManager::install("pmp")
# # library(pmp)
#
# # other non-ion features
# other_df <- ion_cv %>%
#   select(!starts_with("ion")) %>%
#   as.data.frame()
#
# # metabolite levels - transpose so that columns=samples and rows=ion intensities
# ft_df <- ion_cv %>%
#   select(starts_with("ion")) %>%
#   as.data.frame() %>%
#   t()
# ft_df[1:10,1:10]
#
#
# # add batch, class, and order information
# ## all batch=1
# ## class: QC=0, CNTL=1, CASE=2
# ## order = run order
# smp_df <- ion_cv %>%
#   left_join(subset(full_dat, select=c("studyid", "gp"))) %>%
#   mutate(batch=1,
#          class=ifelse(is.na(gp), 0, ifelse(gp=="CASE", 2, 1))) %>%
#   select(batch, class, order) %>%
#   as.data.frame()
# smp_df[1:10,]
# summary(smp_df)
#
#
# corrected_data <- QCRSC(df=ft_df,
#                         order=smp_df$order,
#                         batch=smp_df$batch,
#                         classes=smp_df$class,
#                         spar=0,
#                         log=FALSE,
#                         minQC=5,
#                         qc_label="0")
# t(corrected_data)[1:10,1:10]
#
# comb <- cbind(other_df, t(corrected_data)) %>% as_tibble()
#
# # check for NA
# ion_nmiss <- comb %>%
#   summarise(across(starts_with("ion"), function(x) sum(is.na(x)))) %>%
#   pivot_longer(cols = everything())
#
# # remove columns with missing ions
# full_norm_spline_qc <- comb %>%
#   select(-all_of(ion_nmiss$name[ion_nmiss$value>0]))
# full_norm_spline_qc[1:10,1:10]
#
#
# ###########################################################################!!!
#
#
#
# #===================#
# ### variance stabilization normalization (VSN) ####
# library(vsn)
#
# other_df <- ion_cv %>%
#   select(!starts_with("ion")) %>%
#   as.data.frame()
#
# ion_df <- ion_cv %>%
#   select(starts_with("ion")) %>%
#   as.data.frame() %>%
#   t()
#
# vsn.model <- vsn2(ion_df)
# vsn.data <- predict(vsn.model, ion_df)
#
#
# vsn.data[1:10,1:10]
#
#
# full_norm_vsn <- cbind(other_df, t(vsn.data)) %>% as_tibble()
#
# # check for missing
# full_norm_vsn %>%
#   select(starts_with("ion")) %>%
#   summarise(across(everything(), function(x) sum(is.na(x)))) %>%
#   pivot_longer(cols = everything(),
#                names_to="ion",
#                values_to="nmiss") %>%
#   count(nmiss>0)
#
#
#
#
#
#
#
# ### normalize 1: divide by cohort median and multiply by overall median
#
# # calculate median across all samples
# full_med <- ion_per_sample %>%
#   filter(sampletype=="Sample") %>%
#   summarise(across(starts_with("ion"), median)) %>%
#   as_vector()
# head(full_med)
#
# # calculate cohort-specific median
# ion_df <- ion_per_sample %>%
#   group_by(cohort) %>%
#   mutate(across(starts_with("ion"),
#                 function(x) x/median(x))) %>%
#   ungroup() %>%
#   select(starts_with("ion"))
# dim(test1)
#
# other_df <- ion_per_sample %>% select(!starts_with("ion"))
# dim(other_df)
#
#
# full_norm_batch_med <- cbind(other_df, t(apply(ion_df, 1, function(x) x*full_med))) %>% as_tibble()
#
#
#
#
# #===================#
# # normalize 1b: normalize to pSS median
# # divide by cohort median and multiply by pSS median
# # calculate median across all samples
# full_med <- ion_per_sample %>%
#   filter(cohort=="pSS") %>%
#   summarise(across(starts_with("ion"), median)) %>%
#   as_vector()
# head(full_med)
#
# # calculate cohort-specific median
# test1 <- ion_per_sample %>%
#   group_by(cohort) %>%
#   mutate(across(starts_with("ion"),
#                 function(x) x/median(x))) %>%
#   ungroup() %>%
#   select(starts_with("ion"))
# dim(test1)
#
# test2 <- ion_per_sample %>% select(!starts_with("ion"))
# dim(test2)
#
#
# full_norm_pss_med <- cbind(test2, t(apply(test1, 1, function(x) x*full_med))) %>% as_tibble()
#
#
#
#
#
#
# #===============================================#
#
# #EVALUATE:
# # choose method
# # full_norm <- full_norm_pqn
#
#
#
#
# #===============================================#
# ### EVALUATE: PCA
# dat_norm_pca <- full_norm %>% select(starts_with("ion"))
# pca_norm <- prcomp(dat_norm_pca)
#
# # plot by cohort
# fviz_pca_ind(pca_norm,
#              label="none",
#              alpha=.5,
#              habillage = full_norm_spline$cohort,
#              palette=c("blue", "red"),
#              addEllipses = TRUE,
#              legend.title="Cohort") +
#   labs(title="PCA after normalization")
#
#
# fviz_pca_ind(pca_norm,
#              label="none",
#              alpha=.5,
#              habillage = full_norm_spline$plate,
#              addEllipses = TRUE,
#              legend.title="Plate") +
#   labs(title="PCA after normalization")
#
#
#
# #===========================================#
# ## RLE plots ####
# ### EVALUATE: RLE
#
# # randomly sample 200 samples
# samp_rows <- sample(1:nrow(ion_per_sample), 200, replace=FALSE)
#
# # before normalization
# dat_nonorm_rle <- ion_per_sample[samp_rows,] %>%
#   mutate(across(starts_with("ion"),
#                 function(x) x-median(x))) %>%
#   select(studyid, cohort, starts_with("ion")) %>%
#   pivot_longer(-c(studyid, cohort), names_to="ion")
#
# dat_nonorm_rle %>%
#   ggplot()+
#   geom_boxplot(aes(studyid, value, fill=cohort), outlier.size=.001) +
#   geom_hline(yintercept=0) +
#   ylim(-20000, 20000) +
#   labs(x="sample", y="relative log expression") +
#   theme(axis.text.x=element_blank())
#
#
#
# # after normalization
# dat_norm_rle <- full_norm_spline[samp_rows,] %>%
#   # group_by(cohort) %>%
#   mutate(across(starts_with("ion"),
#                 function(x) x-median(x))) %>%
#   select(studyid, cohort, starts_with("ion")) %>%
#   pivot_longer(-c(studyid, cohort), names_to="ion")
#
# dat_norm_rle %>%
#   ggplot()+
#   geom_boxplot(aes(studyid, value, fill=cohort), outlier.size=.01) +
#   geom_hline(yintercept=0) +
#   ylim(-20000, 20000) +
#   labs(x="sample", y="relative log expression") +
#   theme(axis.text.x=element_blank())
