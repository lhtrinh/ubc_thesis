##### NORMALIZATION AND QC PIPELINE
##### KEEP ALL FEATURES DETECTED AT LEAST 20%



#===============================================#
#####=====#####
## This code details the following steps:
### 1. Filter metabolites that were missing in more than 50% of the experimental samples
### 2. Calculate CV and ICC for non-normalized samples
### 3. PCA plots for non-normalized samples (by plate and cohort)
### 4. Normalization methods. For each method:
###  a. Perform normalization
###  b. Calculate post CV and ICC
###  c. PCA plots for normalized samples (by plate and cohort)
### 5. Select method and average by studyid
### 6. Pareto scaling and within-pair mean centering
#===============================================#






#===============================================#
##########
### Setup
##########

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
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load non-normalized metabolomics data set
full_ion_nonorm <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized.csv")


# load questionnaire data
full_dat <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")


# load ion hit frequency for missingness information
ion_hit_freq_full <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/DATA_ionGapHits_NoPss.xlsx",
  sheet="ions")

# define ICC cutoff points
icc_breaks <- c(0, .4, .6, .75, 1)






#===============================================#
##########
### 1. Filter metabolites with more than 50% missing
##########

# dim(ion_hit_freq_full)
# head(ion_hit_freq_full)
#
#
# # there are 2,550 measurements from 1,275 samples
# ion_freq2 <- ion_hit_freq_full %>%
#   mutate(ion = paste("ion", ionIdx, sep="_"))
#
#
# # calculate count of metabolites by hit frequency
# for (i in seq(.5,1,.1)){
#   prop_hit <- i
#   ion_cnt <- length(which(ion_freq2$ionGapHits>=2550*prop_hit))
#   print(paste("Detected in at least ",
#               prop_hit*100,
#               "% of samples: ",
#               ion_cnt,
#               " (",
#               round(ion_cnt/1677*100, 1),
#               "%)",
#               sep=""))
# }
#
# # set cutoff threshold
# thres <- .5
# full_ion_trunc <- full_ion_nonorm %>%
#   select(-all_of(ion_freq2$ion[ion_freq2$ionGapHits<2550*thres]))
#
# dim(full_ion_nonorm)
# dim(full_ion_trunc)



# add variable to indicate run order
full_ion_nonorm <- full_ion_nonorm %>% arrange(dsIdx)
identical(sort(full_ion_nonorm$dsIdx), full_ion_nonorm$dsIdx) # data is sorted by run order

# add tubeid order
full_ion_nonorm$tubeid_lag <- lag(full_ion_nonorm$tubeid)
order <- c(1, rep(0,nrow(full_ion_nonorm)-1))
for (i in 2:nrow(full_ion_nonorm)){
  order[i] <- ifelse(full_ion_nonorm$tubeid[i]==full_ion_nonorm$tubeid_lag[i],
                     order[i-1],
                     order[i-1]+1)
}
full_ion_nonorm$order <- order


# calculate mean of each consecutive pair of samples
# so that each tube id only has one measurement for each metabolite
full_nonorm_sampleid <- full_ion_nonorm %>%
  group_by(studyid, tubeid, order, sampletype, cohort, plate, dup) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup() %>%
  filter(cohort %in% c("BCGP", "ATP", "pSS")) %>%
  arrange(order)


# # filter only samples with duplicates or QC samples
# full_nonorm_dup <- full_nonorm_sampleid %>%
#   filter(dup=="YES")





#===============================================#
##########
### 2. Calculate CV and ICC for non-normalized samples
##########
nonorm_samp <- full_nonorm_sampleid %>%
  filter(sampletype=="Sample")


# calculate CV
nonorm_samp <- full_nonorm_sampleid %>%
  filter(sampletype=="Sample")


cv_nonorm_samp <- calc_cv_med(nonorm_samp, studyid)
summary(cv_nonorm_samp)

# cv_nonorm_qc <- calc_cv_med(full_nonorm_sampleid[full_nonorm_sampleid$sampletype=="QC",], studyid)


# calculate ICC
icc_nonorm_samp <- calc_icc(nonorm_samp, studyid)

summary(icc_nonorm_samp$icc)
proportions(table(cut(icc_nonorm_samp$icc, breaks=icc_breaks)))




# combine CV and ICC results in 1 df
qc_nonorm <- cv_nonorm_samp %>%
  full_join(icc_nonorm_samp)

round(quantile(qc_nonorm$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_nonorm %>% count(median_cv>.2)

round(quantile(qc_nonorm$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_nonorm$icc, breaks=icc_breaks, right=FALSE))










#===============================================#
##########
### 3. PCA plots for non-normalized samples (by plate and cohort)
##########

dat_nonorm <- full_nonorm_sampleid %>%
  filter(sampletype=="Sample") %>%
  mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))

dat_nonorm_pca <- dat_nonorm %>% select(starts_with("ion"))

pca_nonorm <- prcomp(dat_nonorm_pca)

# plot by plate 11 vs other plates
fviz_pca_ind(pca_nonorm,
             axes = c(1,2),
             label="none",
             alpha=.5,
             habillage =  dat_nonorm$plate11_yn,
             palette=c("blue", "red"),
             addEllipses = TRUE)

# plot by plate
fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  dat_nonorm$plate,
             alpha=.5,
             addEllipses = TRUE)

# plot by cohort
fviz_pca_ind(pca_nonorm,
             label="none",
             alpha=.5,
             habillage = dat_nonorm$cohort,
             palette=c("blue", "red"),
             addEllipses = FALSE)






#===============================================#
##########
### 4. Normalization methods
##########

# normalize using only pSS and sample data
dat_before_norm <- full_nonorm_sampleid %>%
  filter(cohort %in% c("ATP", "BCGP", "pSS")) %>%
  arrange(order)



#===================
# normalize 1: divide by cohort median and multiply by overall median

# calculate median across all samples
full_med <- dat_before_norm %>%
  filter(sampletype=="Sample") %>%
  summarise(across(starts_with("ion"), median)) %>%
  as_vector()
head(full_med)

# calculate cohort-specific median
ion_df <- dat_before_norm %>%
  group_by(cohort) %>%
  mutate(across(starts_with("ion"),
                function(x) x/median(x))) %>%
  ungroup() %>%
  select(starts_with("ion"))
dim(ion_df)

other_df <- dat_before_norm %>% select(!starts_with("ion"))
dim(other_df)


full_norm_cohort_med <- cbind(other_df, t(apply(ion_df, 1, function(x) x*full_med))) %>% as_tibble()



# post QC

qc_cohort_med <- norm_qc(full_norm_cohort_med)

round(quantile(qc_cohort_med$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_cohort_med %>% count(median_cv>.2)

round(quantile(qc_cohort_med$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_cohort_med$icc, breaks=icc_breaks, right=FALSE))

#===================
# normalize 1b: normalize to pSS median
# divide by cohort median and multiply by pSS median
# calculate median across all samples
full_med <- dat_before_norm %>%
  filter(cohort=="pSS") %>%
  summarise(across(starts_with("ion"), median)) %>%
  as_vector()
head(full_med)

# calculate cohort-specific median
test1 <- dat_before_norm %>%
  group_by(cohort) %>%
  mutate(across(starts_with("ion"),
                function(x) x/median(x))) %>%
  ungroup() %>%
  select(starts_with("ion"))
dim(test1)

test2 <- dat_before_norm %>% select(!starts_with("ion"))
dim(test2)


full_norm_pss_med <- cbind(test2, t(apply(test1, 1, function(x) x*full_med))) %>% as_tibble()





# post QC
qc_pss_med <- norm_qc(full_norm_pss_med)

round(quantile(qc_pss_med$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_pss_med %>% count(median_cv>.2)

round(quantile(qc_pss_med$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_pss_med$icc, breaks=icc_breaks, right=FALSE))




#===================
# normalize 2: GMet's method
## divide by sample median and multiply by mean across all medians
full_norm_gmet <- dat_before_norm %>%
  filter(sampletype=="Sample") %>%
  rowwise() %>%
  mutate(m=median(c_across(starts_with("ion")))) %>%
  ungroup() %>%
  group_by(cohort) %>%
  mutate(mm=mean(m)) %>%
  ungroup() %>%
  mutate(across(starts_with("ion"),
                function(x) x*mm/m)) %>%
  select(-m, -mm)


###
# post QC
qc_gmet <- norm_qc(full_norm_gmet)

# summary for CV
round(quantile(qc_gmet$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_gmet %>% count(median_cv>.2)

# summary for ICC
round(quantile(qc_gmet$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_gmet$icc, breaks=icc_breaks, right=FALSE))



#===================
# normalize 2b: linear baseline normalization (similar to the GMet method?)
other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

baseline <- apply(ion_df, 1, median)
baseline_mean <- mean(baseline)
samp_means <- apply(ion_df, 2, mean)
linear_fct <- baseline_mean/samp_means

linear_norm_df <- t(ion_df)*linear_fct

ion_df[1:10, 1:10]
linear_norm_df[1:10, 1:10]

full_norm_linear <- cbind(other_df, linear_norm_df) %>% as_tibble()
full_norm_linear[1:10,1:10]





qc_linear <- norm_qc(full_norm_linear)

round(quantile(qc_linear$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_linear %>% count(median_cv>.2)

round(quantile(qc_linear$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_linear$icc, breaks=icc_breaks, right=FALSE))





# #===================
# # normalize 2c: non-linear baseline normalization (Li-Wong)
# library(affy)
#
#
# data <- dat_before_norm %>%
#   select(starts_with("ion")) %>%
#   t()
#
# other_df <- dat_before_norm %>%
#   select(!starts_with("ion"))
#
#
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
#
# full_norm_liwong <- cbind(other_df, t(liwong.data)) %>% as_tibble()
#
#
#
# qc_linear <- norm_qc(full_norm_linear)
#
# round(quantile(qc_linear$median_cv, probs=c(0,.05, .5, .95, 1)),3)
# qc_linear %>% count(median_cv>.2)
#
# round(quantile(qc_linear$icc, probs=c(0,.05, .5, .95, 1)), 3)
# table(cut(qc_linear$icc, breaks=icc_breaks, right=FALSE))




#===================
# normalize 3: divide by median of all metabolites across each sample
full_norm_samp_med <- dat_before_norm %>%
  rowwise() %>%
  mutate(m=median(c_across(starts_with("ion")))) %>%
  ungroup() %>%
  mutate(across(starts_with("ion"), function(x) x/m)) %>%
  select(-m)




#===================
## normalize 4: quantile normalization for each metabolite
library(preprocessCore)

dat_ion <- dat_before_norm %>%
  select(starts_with("ion"))


dat_ion_quantile <- as.data.frame(normalize.quantiles.robust(as.matrix(test)))
colnames(testt) <- colnames(test)

dat_no_ion <- dat_before_norm %>%
  select(!starts_with("ion"))

full_norm_qtl_meta <- cbind(dat_no_ion, dat_ion_quantile) %>% as_tibble()
full_norm_qtl_meta[1:10,1:10]





#===================
## normalize 4b: quantile normalization for each sample (pick this)
dat_ion <- dat_before_norm %>%
  select(starts_with("ion"))


dat_ion_quantile <- as.data.frame(t(normalize.quantiles(t(as.matrix(dat_ion)))))
colnames(dat_ion_quantile) <- colnames(dat_ion)

dat_no_ion <- dat_before_norm %>%
  select(!starts_with("ion"))

full_norm_qtl <- cbind(dat_no_ion, dat_ion_quantile) %>% as_tibble()
full_norm_qtl[1:10,1:10]


# post QC
qc_qtl <- norm_qc(full_norm_qtl)

round(quantile(qc_qtl$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_qtl %>% count(median_cv>.2)

round(quantile(qc_qtl$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_qtl$icc, breaks=icc_breaks, right=FALSE))




#===================
# normalize 5: cyclic loess using QC samples (QC-RLSC) from Dunn et al

# install.packages("devtools")
# devtools::install_github("metabolomicstats/NormalizeMets")
library(NormalizeMets)

other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ft_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame()

rownames(ft_df) <- paste("tube", dat_before_norm$order, sep="_")
ft_df[1:10,1:10]

smp_df <- dat_before_norm %>%
  left_join(subset(full_dat, select=c("studyid", "gp"))) %>%
  mutate(batch=ifelse(cohort %in% c("BCGP", "ATP"), cohort, "QC"),
         class=ifelse(is.na(gp), 0, ifelse(gp=="CASE", 2, 1)),
         order=1:nrow(dat_before_norm)) %>%
  select(batch, class, order) %>%
  as.data.frame()
smp_df[1:10,]

rownames(smp_df) <- paste("tube", dat_before_norm$order, sep="_")



test_out <- NormQcsamples(featuredata=ft_df, sampledata=smp_df, span=.75, lg=FALSE)

full_norm_loess <- cbind(other_df, test_out$featuredata) %>%
  as_tibble()

full_norm_loess[1:10,1:15]



# post QC
qc_loess <- norm_qc(full_norm_loess)

round(quantile(qc_loess$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_loess %>% count(median_cv>.2)

round(quantile(qc_loess$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_loess$icc, breaks=icc_breaks, right=FALSE))




#===================
# normalize 6: probabilistic quotient normalization (PQN)

other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>% t()

ion_df[1:10, 1:10]


# ref <- apply(ion_df, 1, median)
# quotient <- ion_df/ref
# quotient_med <- apply(quotient, 2, median)
# pqn_df <- t(ion_df)/quotient_med
#
# pqn_df[1:10, 1:10]
#
# full_norm_pqn <- cbind(other_df, pqn_df) %>% as_tibble()
# full_norm_pqn[1:10,1:10]



######
# use QC samples as reference
ion_qc_df <- dat_before_norm %>%
  filter(sampletype=="Sample") %>%
  select(starts_with("ion")) %>% t()


ref <- apply(ion_qc_df, 1, median)
quotient <- ion_df/ref
quotient_med <- apply(quotient, 2, median)

pqn_df_from_qc <- t(ion_df/quotient_med)
pqn_df_from_qc[1:10,1:10]
dim(pqn_df_from_qc)


full_norm_pqn_qc <- cbind(other_df, pqn_df_from_qc) %>% as_tibble()



# post QC
qc_pqn_from_qc <- norm_qc(full_norm_pqn_qc)

round(quantile(qc_pqn_from_qc$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_pqn_from_qc %>% count(median_cv>.2)

round(quantile(qc_pqn_from_qc$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_pqn_from_qc$icc, breaks=icc_breaks, right=FALSE))





#===================
# normalize 7: cubic spline
library(NCStats)
library(affy)
other_df <- dat_before_norm %>%
  filter(sampletype=="Sample") %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- dat_before_norm %>%
  filter(sampletype=="Sample") %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

set.seed(14101)
spline_dat <- t(normalize.qspline(
  ion_df,
  samples=.05,
  target=apply(ion_df,1,geomean),
  verbose=TRUE
))
colnames(spline_dat) <- rownames(ion_df)

ion_df[1:10,1:10]
spline_dat[1:10,1:10]


full_norm_spline <- cbind(other_df, spline_dat) %>% as_tibble()
full_norm_spline[1:10,1:10]
dim(full_norm_spline)


# post QC
qc_spline <- norm_qc(full_norm_spline)

round(quantile(qc_spline$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_spline %>% count(median_cv>.2)

round(quantile(qc_spline$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_spline$icc, breaks=icc_breaks, right=FALSE))






#===================
# normalize 7b: cubic spline using QC samples (QC-RSC)
# BiocManager::install("pmp")
# library(pmp)

# other non-ion features
other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

# metabolite levels - transpose so that columns=samples and rows=ion intensities
ft_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()
ft_df[1:10,1:10]


# add batch, class, and order information
## all batch=1
## class: QC=0, CNTL=1, CASE=2
## order = run order
smp_df <- dat_before_norm %>%
  left_join(subset(full_dat, select=c("studyid", "gp"))) %>%
  mutate(batch=1,
         class=ifelse(is.na(gp), 0, ifelse(gp=="CASE", 2, 1))) %>%
  select(batch, class, order) %>%
  as.data.frame()
smp_df[1:10,]
summary(smp_df)


corrected_data <- QCRSC(df=ft_df,
                        order=smp_df$order,
                        batch=smp_df$batch,
                        classes=smp_df$class,
                        spar=0,
                        log=FALSE,
                        minQC=5,
                        qc_label="0")
t(corrected_data)[1:10,1:10]

comb <- cbind(other_df, t(corrected_data)) %>% as_tibble()

# check for NA
ion_nmiss <- comb %>%
  summarise(across(starts_with("ion"), function(x) sum(is.na(x)))) %>%
  pivot_longer(cols = everything())

# remove columns with missing ions
full_norm_spline_qc <- comb %>%
  select(-all_of(ion_nmiss$name[ion_nmiss$value>0]))
full_norm_spline_qc[1:10,1:10]



# post QC
qc_spline_qc <- norm_qc(full_norm_spline_qc)

round(quantile(qc_spline_qc$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_spline_qc %>% count(median_cv>.2)

round(quantile(qc_spline_qc$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_spline_qc$icc, breaks=icc_breaks, right=FALSE))



#===================
# normalize 8: variance stabilization normalization (VSN)
library(vsn)

other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

vsn.model <- vsn2(ion_df)
vsn.data <- predict(vsn.model, ion_df)


vsn.data[1:10,1:10]


full_norm_vsn <- cbind(other_df, t(vsn.data)) %>% as_tibble()

# check for missing
full_norm_vsn %>%
  select(starts_with("ion")) %>%
  summarise(across(everything(), function(x) sum(is.na(x)))) %>%
  pivot_longer(cols = everything(),
               names_to="ion",
               values_to="nmiss") %>%
  count(nmiss>0)



# post QC
qc_vsn <- norm_qc(full_norm_vsn)

round(quantile(qc_vsn$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_vsn %>% count(median_cv>.2)

round(quantile(qc_vsn$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_vsn$icc, breaks=icc_breaks, right=FALSE))






#===============================================#
##########
### 5. Select method and average by studyid
###    to get only 1 measurement per sample
##########

# cubic spline
full_norm <- full_norm_spline %>%
  filter(sampletype=="Sample") %>%
  group_by(studyid, dup) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()

# # median (gmet's version)
# full_norm <- full_norm_gmet %>%
#   filter(sampletype=="Sample") %>%
#   group_by(studyid, dup) %>%
#   summarise(across(starts_with("ion"), mean)) %>%
#   ungroup()

# # no normalization
# full_norm <- full_nonorm_sampleid %>%
#   filter(sampletype=="Sample") %>%
#   group_by(studyid, dup) %>%
#   summarise(across(starts_with("ion"), mean)) %>%
#   ungroup()
# dim(full_norm)



# vsn
full_norm <- full_norm_vsn %>%
  filter(sampletype=="Sample") %>%
  group_by(studyid, dup) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup()


#===============================================#
##########
### 6. Finish preprocessing
###    - Calculate mean of duplicate samples, Pareto-scale, and within-pair mean centering
###    - Merge with questionnaire data
##########

full_sc <- full_norm %>%
  mutate(across(starts_with("ion"),
                function(x) x-mean(x)/sqrt(sd(x))))


# Calculate mean of duplicate samples
full_centered <- full_sc %>%
  full_join(full_dat[,colnames(full_dat) %in% c("studyid", "match_id")]) %>%
  group_by(match_id) %>%
  mutate(across(starts_with("ion"), mean)) %>%
  ungroup()






