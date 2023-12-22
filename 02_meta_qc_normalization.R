# NORMALIZATION AND QC PIPELINE ####



#===============================================#
## Summary ####
#==============#
## This code details the following steps:
### 1. Filter metabolites that were missing in more than 50% of the experimental samples
### 2. Calculate CV and ICC for non-normalized samples
### 3. PCA plots for non-normalized samples (by plate and cohort)
### 4. Normalization methods. For each method:
  ###  a. Perform normalization
  ###  b. Calculate post CV and ICC
  ###  c. PCA plots for normalized samples (by plate and cohort)
### 5. Select method
### 6. Calculate mean by study ID for each duplicated sample
### 7. Mean-centering
### 8. Merge with questionnaire data
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
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load non-normalized metabolomics data set
full_ion_nonorm <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized.csv")


# load questionnaire data
full_dat <- import_survey()


# load ion hit frequency for missingness information
ion_hit_freq_full <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/DATA_ionGapHits_NoPss.xlsx",
  sheet="ions")

# define ICC cutoff points
icc_breaks <- c(0, .4, .6, .75, 1)






#===============================================#
## 1. Filtering ####
#### Remove metabolites with more than 50% missing
#==============#



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

# set cutoff threshold
thres <- .5
full_ion_trunc <- full_ion_nonorm %>%
  select(-all_of(ion_freq2$ion[ion_freq2$ionGapHits<2550*thres]))

dim(full_ion_nonorm)
dim(full_ion_trunc)



# add variable to indicate run order
full_ion_trunc <- full_ion_trunc %>% arrange(dsIdx)
identical(sort(full_ion_trunc$dsIdx), full_ion_trunc$dsIdx) # data is sorted by run order

# add tubeid order
full_ion_trunc$tubeid_lag <- lag(full_ion_trunc$tubeid)
order <- c(1, rep(0,nrow(full_ion_trunc)-1))
for (i in 2:nrow(full_ion_trunc)){
  order[i] <- ifelse(full_ion_trunc$tubeid[i]==full_ion_trunc$tubeid_lag[i],
                     order[i-1],
                     order[i-1]+1)
}
full_ion_trunc$order <- order


# calculate mean of each consecutive pair of samples (sampleid)
# so that each tube id only has one measurement for each metabolite
full_nonorm_sampleid <- full_ion_trunc %>%
  group_by(studyid, tubeid, order, sampletype, cohort, plate, dup) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup() %>%
  arrange(order)








#===============================================#
## 2. Calculate CV and ICC ####
##### for non-normalized samples
#==============#



# filter only samples with duplicates
full_nonorm_dup <- full_nonorm_sampleid %>%
  filter(sampletype=="Sample" & dup=="YES")



## CV
cv_nonorm_1 <- calc_cv_med(full_nonorm_dup, studyid)

round(quantile(cv_nonorm_1$median_cv, probs=c(0,.05, .5, .95, 1)),3)
high_cv_ions <- cv_nonorm_1 %>%
  filter(median_cv>.2)

high_cv_ions

# qc_nonorm_1 <- cv_nonorm_1 %>%
#   full_join(icc_nonorm_1)


## ICC
# icc_nonorm_1 <- calc_icc(full_nonorm_dup, studyid)

# round(quantile(qc_nonorm_1$icc, probs=c(0,.05, .5, .95, 1)), 3)
# table(cut(icc_nonorm_1$icc, breaks=icc_breaks, right=FALSE))





# # calculate count of metabolites with low ICC or high CV
# qc_nonorm_1 %>%
#   count(median_cv>.2, icc<.5)
#
# # plot ICC distribution
# qc_nonorm_1 %>%
#   ggplot(aes(icc, after_stat(scaled))) +
#   geom_density() +
#   scale_x_continuous(breaks = seq(0,1,.1))




# filter out ions with CV > 20%
# then average metabolite values of samples
full_nonorm_cv <- full_nonorm_sampleid %>%
  select(-all_of(high_cv_ions$metabolite))




#===============================================#
## 3. Non-normalized PCA plots ####
##### for non-normalized samples (by plate and cohort)
#==============#


dat_nonorm <- full_nonorm_cv %>%
  filter(sampletype=="Sample") %>%
  mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))

dat_nonorm_pca <- dat_nonorm %>% select(starts_with("ion"))

pca_nonorm <- prcomp(dat_nonorm_pca)

# plot by plate 11 vs other plates
# fviz_pca_ind(pca_nonorm,
#              axes = c(1,2),
#              label="none",
#              alpha=.5,
#              habillage =  dat_nonorm$plate11_yn,
#              palette=c("blue", "red"),
#              addEllipses = TRUE)

# plot by plate
fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  dat_nonorm$plate,
             alpha=.5,
             addEllipses = TRUE)+
  labs(title="PCA before normalization")

# plot by cohort
fviz_pca_ind(pca_nonorm,
             label="none",
             alpha=.5,
             habillage = dat_nonorm$cohort,
             palette=c("blue", "red"),
             addEllipses = TRUE) +
  labs(title="PCA before normalization")






#===============================================#
## 4. Normalization methods ####
#==============#

# normalize using only pSS and sample data


nonorm_samp <- full_nonorm_cv %>%
  filter(sampletype=="Sample")
nonorm_pss <- full_nonorm_cv %>%
  filter(cohort=="pSS")

nonorm_samp_pss <- full_nonorm_cv %>%
  filter(cohort %in% c("BCGP", "ATP", "pSS"))



#===================#
### cubic splines ####
library(NCStats)
library(affy)
other_df <- nonorm_samp %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- nonorm_samp %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

set.seed(30340)
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


# post QC
# qc_spline <- norm_qc(full_norm_spline)
#
# round(quantile(qc_spline$median_cv, probs=c(0,.05, .5, .95, 1)),3)
# qc_spline %>% count(median_cv>.2)
#
# round(quantile(qc_spline$icc, probs=c(0,.05, .5, .95, 1)), 3)
# table(cut(qc_spline$icc, breaks=icc_breaks, right=FALSE))



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







#===============================================#
## Select method ####

# cubic spline
full_norm <- full_norm_spline

# # median (gmet's version)
# full_norm <- full_norm_gmet

# no normalization
# full_norm <- nonorm_samp




#===========================================#
## RLE plots ####
### EVALUATE: RLE


# before normalization
ids <- nonorm_samp$tubeid

set.seed(46867)
samp_id <- sample(ids, 200, replace=FALSE)


dat_nonorm_rle <- nonorm_samp %>%
  mutate(across(starts_with("ion"),
                function(x) x-median(x))) %>%
  select(tubeid, cohort, starts_with("ion")) %>%
  pivot_longer(-c(tubeid, cohort), names_to="ion") %>%
  filter(tubeid %in% samp_id)


dat_nonorm_rle %>%
  ggplot()+
  geom_boxplot(aes(tubeid, value, fill=cohort), outlier.size=.001) +
  geom_hline(yintercept=0) +
  ylim(-20000, 20000) +
  labs(x="sample", y="relative log expression") +
  theme(axis.text.x=element_blank())




# after normalization
# ions <- names(full_sc)[grepl("^ion_", names(full_sc))]
# samp_ion <- sample(ions, 30, replace=FALSE)

dat_norm_rle <- full_norm_spline %>%
  # group_by(cohort) %>%
  mutate(across(starts_with("ion"),
                function(x) x-median(x))) %>%
  select(tubeid, cohort, starts_with("ion")) %>%
  pivot_longer(-c(tubeid, cohort), names_to="ion") %>%
  filter(tubeid %in% samp_id)


dat_norm_rle %>%
  ggplot()+
  geom_boxplot(aes(tubeid, value, fill=cohort), outlier.size=.01) +
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



# save to csv
write_csv(full_centered,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/normalized_metabolomics_pair_centered.csv")





#===============================================#
## Scaling ####
full_sc <- full_centered %>%
  mutate(across(starts_with("ion"), function(x) x/sd(x)))

# save to csv
write_csv(full_sc,
          "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/normalized_metabolomics_scaled.csv")














#==========================================================================#

## Appendix: Other normalization methods ####

#### other methods considered by not reported
#===================#




#===================#

### GMet's method ####
## divide by sample median and multiply by mean across all medians
full_norm_gmet <- nonorm_samp %>%
  rowwise() %>%
  mutate(m=median(c_across(starts_with("ion")))) %>%
  ungroup() %>%
  group_by(cohort) %>%
  mutate(mm=mean(m)) %>%
  ungroup() %>%
  mutate(across(starts_with("ion"),
                function(x) x*mm/m)) %>%
  select(-m, -mm)



# post QC
qc_gmet <- norm_qc(full_norm_gmet)

round(quantile(qc_gmet$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_gmet %>% count(median_cv>.2)

round(quantile(qc_gmet$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_gmet$icc, breaks=icc_breaks, right=FALSE))



#===================#
### linear baseline normalization ####
other_df <- nonorm_samp %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- nonorm_samp %>%
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





#===================#
### non-linear baseline normalization (Li-Wong) ####
library(affy)


data <- nonorm_samp %>%
  select(starts_with("ion")) %>%
  t()

other_df <- nonorm_samp %>%
  select(!starts_with("ion"))



#---First step: Find baseline sample
average.intensity <- apply(data,2,mean)
median.number <- round(ncol(data)/2 + 0.1)				#R has an add way of rounding.
#the additional 0.1 ensures that it rounds properly
ordering <- order(average.intensity)
median.sample.number <- ordering[median.number]
median.sample <- data[,median.sample.number]

#---Apply normalization
liwong.data=vector()
for(i in 1:ncol(data)){
  liwong.model <- normalize.invariantset(data=data[,i],
                                         ref=median.sample,
                                         prd.td=c(0.003,0.007))			#the threshold of the rank-invariant set might need to be adjusted from case to case
  liwong.sample <- predict(liwong.model$n.curve$fit,		#chosen such that the rank-invariant set it sufficiently large
                           data[,i])
  liwong.data <- cbind(liwong.data,liwong.sample$y)
}


full_norm_liwong <- cbind(other_df, t(liwong.data)) %>% as_tibble()


# post QC
qc_liwong <- norm_qc(full_norm_liwong)

round(quantile(qc_liwong$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_liwong %>% count(median_cv>.2)

round(quantile(qc_liwong$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_liwong$icc, breaks=icc_breaks, right=FALSE))






#===================#
### quantile normalization ####
dat_ion <- nonorm_samp %>%
  select(starts_with("ion"))


dat_ion_quantile <- as.data.frame(t(normalize.quantiles(t(as.matrix(dat_ion)))))
colnames(dat_ion_quantile) <- colnames(dat_ion)

dat_no_ion <- nonorm_samp %>%
  select(!starts_with("ion"))

full_norm_qtl <- cbind(dat_no_ion, dat_ion_quantile) %>% as_tibble()
full_norm_qtl[1:10,1:10]


# post QC
qc_qtl <- norm_qc(full_norm_qtl)

round(quantile(qc_qtl$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_qtl %>% count(median_cv>.2)

round(quantile(qc_qtl$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_qtl$icc, breaks=icc_breaks, right=FALSE))






#===================#
### cyclic loess using QC samples (QC-RLSC) ####

# install.packages("devtools")
# devtools::install_github("metabolomicstats/NormalizeMets")
library(NormalizeMets)

other_df <- nonorm_samp %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ft_df <- nonorm_samp %>%
  select(starts_with("ion")) %>%
  as.data.frame()

rownames(ft_df) <- paste("tube", nonorm_samp$order, sep="_")
ft_df[1:10,1:10]

smp_df <- nonorm_samp %>%
  left_join(subset(full_dat, select=c("studyid", "gp"))) %>%
  mutate(batch=ifelse(cohort %in% c("BCGP", "ATP"), cohort, "QC"),
         class=ifelse(is.na(gp), 0, ifelse(gp=="CASE", 2, 1)),
         order=1:nrow(nonorm_samp)) %>%
  select(batch, class, order) %>%
  as.data.frame()
smp_df[1:10,]

rownames(smp_df) <- paste("tube", nonorm_samp$order, sep="_")



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








#===================#
### probabilistic quotient normalization (PQN) ####

other_df <- nonorm_samp %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- nonorm_samp %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

ion_df[1:10, 1:10]


ref <- apply(ion_df, 1, median)
quotient <- ion_df/ref
quotient_med <- apply(quotient, 2, median)
pqn_df <- t(ion_df)/quotient_med

pqn_df[1:10, 1:10]

full_norm_pqn <- cbind(other_df, pqn_df) %>% as_tibble()
full_norm_pqn[1:10,1:10]



# post QC
qc_pqn <- norm_qc(full_norm_pqn)

round(quantile(qc_pqn$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_pqn %>% count(median_cv>.2)

round(quantile(qc_pqn$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_pqn$icc, breaks=icc_breaks, right=FALSE))




#===================#
### probabilistic quotient normalization (PQN) with QC ####

other_df <- nonorm_samp %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- nonorm_samp %>%
  select(starts_with("ion")) %>%
  as.data.frame()
ion_df[1:10, 1:10]

ion_qc_df <- nonorm_pss %>%
  select(starts_with("ion")) %>%
  as.data.frame()

ref <- apply(ion_qc_df, 2, median)
quotient <- ion_df/ref
quotient_med <- apply(quotient, 1, median)
pqn_df <- ion_df/quotient_med

pqn_df[1:10, 1:10]

full_norm_pqn_pss <- cbind(other_df, pqn_df) %>% as_tibble()
full_norm_pqn_pss[1:10,1:10]



# post QC
qc_pqn_pss <- norm_qc(full_norm_pqn_pss)

round(quantile(qc_pqn_pss$median_cv, probs=c(0,.05, .5, .95, 1)),3)
qc_pqn_pss %>% count(median_cv>.2)

round(quantile(qc_pqn_pss$icc, probs=c(0,.05, .5, .95, 1)), 3)
table(cut(qc_pqn_pss$icc, breaks=icc_breaks, right=FALSE))




#######################################################################!!!









#===================#
### cubic splines using QC samples (version 11/01/2023) ####
qc_spec <- full_nonorm_sampleid %>%
  filter(cohort=="pSS") %>%
  select(-all_of(high_cv_ions$metabolite)) %>%
  summarise(across(starts_with("ion"), geomean)) %>%
  as_vector()


set.seed(30340)
qc_spline_dat <- t(normalize.qspline(
  ion_df,
  samples=.05,
  target=qc_spec,
  verbose=TRUE
))

colnames(qc_spline_dat) <- rownames(ion_df)
qc_spline_dat[1:10,1:10]


full_norm_qc_spline <- cbind(other_df, qc_spline_dat) %>% as_tibble()
full_norm_qc_spline[1:10,1:10]




dat_norm_pca <- full_norm_qc_spline %>% select(starts_with("ion"))
pca_norm <- prcomp(dat_norm_pca)

# plot by cohort
fviz_pca_ind(pca_norm,
             label="none",
             alpha=.5,
             habillage = full_norm_qc_spline$cohort,
             palette=c("blue", "red"),
             addEllipses = FALSE)






#===================#
### cubic splines using QC samples (QC-RSC) ####
# BiocManager::install("pmp")
# library(pmp)

# other non-ion features
other_df <- nonorm_samp %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

# metabolite levels - transpose so that columns=samples and rows=ion intensities
ft_df <- nonorm_samp %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()
ft_df[1:10,1:10]


# add batch, class, and order information
## all batch=1
## class: QC=0, CNTL=1, CASE=2
## order = run order
smp_df <- nonorm_samp %>%
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

###########################################################################!!!



#===================#
### variance stabilization normalization (VSN) ####
library(vsn)

other_df <- nonorm_samp %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- nonorm_samp %>%
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














### normalize 1: divide by cohort median and multiply by overall median

# calculate median across all samples
full_med <- full_nonorm_sampleid %>%
  filter(sampletype=="Sample") %>%
  summarise(across(starts_with("ion"), median)) %>%
  as_vector()
head(full_med)

# calculate cohort-specific median
ion_df <- full_nonorm_sampleid %>%
  group_by(cohort) %>%
  mutate(across(starts_with("ion"),
                function(x) x/median(x))) %>%
  ungroup() %>%
  select(starts_with("ion"))
dim(test1)

other_df <- full_nonorm_sampleid %>% select(!starts_with("ion"))
dim(other_df)


full_norm_batch_med <- cbind(other_df, t(apply(ion_df, 1, function(x) x*full_med))) %>% as_tibble()




#===================#
# normalize 1b: normalize to pSS median
# divide by cohort median and multiply by pSS median
# calculate median across all samples
full_med <- full_nonorm_sampleid %>%
  filter(cohort=="pSS") %>%
  summarise(across(starts_with("ion"), median)) %>%
  as_vector()
head(full_med)

# calculate cohort-specific median
test1 <- full_nonorm_sampleid %>%
  group_by(cohort) %>%
  mutate(across(starts_with("ion"),
                function(x) x/median(x))) %>%
  ungroup() %>%
  select(starts_with("ion"))
dim(test1)

test2 <- full_nonorm_sampleid %>% select(!starts_with("ion"))
dim(test2)


full_norm_pss_med <- cbind(test2, t(apply(test1, 1, function(x) x*full_med))) %>% as_tibble()


qc_pss_med <- norm_qc(full_norm_pss_med)