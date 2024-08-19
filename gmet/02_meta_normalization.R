# normalize using only pSS and sample data
dat_before_norm <- full_nonorm_sampleid %>%
  filter(cohort %in% c("ATP", "BCGP", "pSS")) %>%
  arrange(order)



####################################################

# normalize 1: divide by cohort median and multiply by overall median

# calculate median across all samples
full_med <- dat_before_norm %>%
  filter(sampletype=="Sample") %>%
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


full_ion_norm <- cbind(test2, t(apply(test1, 1, function(x) x*full_med))) %>% as_tibble()



###
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


full_ion_norm <- cbind(test2, t(apply(test1, 1, function(x) x*full_med))) %>% as_tibble()


###############################################

# normalize 2: GMet's method
## divide by sample median and multiply by mean across all medians
full_ion_norm <- dat_before_norm %>%
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
# normalize 2b: linear baseline normalization (similarto the GMet method?)
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

full_ion_norm <- cbind(other_df, linear_norm_df) %>% as_tibble()






#############################################################

# normalize 3: divide by median of all metabolites across each sample
full_ion_norm <- dat_before_norm %>%
  rowwise() %>%
  mutate(m=median(c_across(starts_with("ion")))) %>%
  ungroup() %>%
  mutate(across(starts_with("ion"), function(x) x/m)) %>%
  select(-m)


#############################################################
# normalize 4: divide by median NIST, multiply by mean of median NISTs
nist_med <- dat_before_norm %>%
  filter(cohort=="NIST") %>%
  summarise(across(starts_with("ion"), median)) %>%
  as_vector()


test1 <- dat_before_norm %>% select(starts_with("ion"))
test2 <- dat_before_norm %>% select(!starts_with("ion"))
dim(test2)



full_ion_norm <- cbind(test2, t(apply(test1, 1, function(x) x/nist_med*mean(nist_med)))) %>% as_tibble()


#############################################################
# normalize 5: divide by cohort median only
full_ion_norm <- dat_before_norm %>%
  group_by(cohort) %>%
  mutate(across(starts_with("ion"),
                function(x) x/median(x))) %>%
  ungroup()




#############################################################
# normalize 6: quantile normalization for each metabolite
library(preprocessCore)

test <- dat_before_norm %>%
  select(starts_with("ion"))


testt <- as.data.frame(normalize.quantiles.robust(as.matrix(test)))
colnames(testt) <- colnames(test)

test2 <- dat_before_norm %>%
  select(!starts_with("ion"))

full_ion_norm <- cbind(test2, testt) %>% as_tibble()
full_ion_norm[1:10,1:10]


#############################################################
# normalize 7: quantile normalization for each sample
test <- dat_before_norm %>%
  select(starts_with("ion"))


testt <- as.data.frame(normalize.quantiles(as.matrix(test)))

testt <- as.data.frame(t(normalize.quantiles(t(as.matrix(test)))))
colnames(testt) <- colnames(test)

test2 <- dat_before_norm %>%
  select(!starts_with("ion"))

full_ion_norm <- cbind(test2, testt) %>% as_tibble()



#############################################################
# normalize 8: cyclic loess
install.packages("devtools")
devtools::install_github("metabolomicstats/NormalizeMets")
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

full_ion_norm <- cbind(other_df, test_out$featuredata) %>%
  as_tibble()

full_ion_norm[1:10,1:15]
# test <- cbind(other_df, ft_df - test_out$featuredata) %>%
#   as_tibble()




#############################################################
# normalize 9: probabilistic quotient normalization (PQN)

other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

ion_df[1:10, 1:10]


ref <- apply(ion_df, 1, median)
quotient <- ion_df/ref
quotient_med <- apply(quotient, 2, median)
pqn_df <- t(ion_df)/quotient_med

pqn_df[1:10, 1:10]

full_ion_norm <- cbind(other_df, pqn_df) %>% as_tibble()
full_ion_norm[1:10,1:10]



#############################################################
# normalize 10a: cubic spline
library(NCStats)
library(affy)
other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

ion_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

spline_dat <- t(normalize.qspline(
  ion_df,
  samples=.05,
  target=apply(ion_df,1,geomean),
  verbose=TRUE
))

colnames(spline_dat) <- rownames(ion_df)

ion_df[1:10,1:10]
spline_dat[1:10,1:10]


full_ion_norm <- cbind(other_df, spline_dat) %>% as_tibble()



###
# normalize 10b: cubic spline - quality control (QC-RSC)
BiocManager::install("pmp")
library(pmp)


df <- dat_before_norm %>%
  filter(cohort %in% c("ATP" , "BCGP", "pSS"))

# other non-ion features
other_df <- dat_before_norm %>%
  select(!starts_with("ion")) %>%
  as.data.frame()

# metabolite levels - transpose so that columns=samples and rows=ion intensities
ft_df <- dat_before_norm %>%
  select(starts_with("ion")) %>%
  as.data.frame() %>%
  t()

# rownames(ft_df) <- paste("tube", dat_before_norm$order, sep="_")
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


# rownames(smp_df) <- paste("tube", dat_before_norm$order, sep="_")



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
full_ion_norm <- comb %>%
  select(-all_of(ion_nmiss$name[ion_nmiss$value>0]))






#############################################################
# normalize 11: variance stabilization normalization (VSN)
# library(vsn)
#
# vsn_df <- justvsn(ft_df)
# vsn_df[1:10,1:10]
#
#
# full_ion_norm <- cbind(other_df, t(vsn_df)) %>% as_tibble()
#
#
# # combine with median normalization??




##############################################################
##############################################################
# pareto-scaling

full_nonorm_sc <- full_ion_norm %>%
  mutate(across(starts_with("ion"),
                function(x) (x-mean(x))/sqrt(sd(x))))



# mean-centered
full_centered <- full_nonorm_sc %>%
  mutate(cohort=tolower(cohort)) %>%
  full_join(full_dat[c("studyid", "match_id")]) %>%
  group_by(match_id) %>%
  mutate(across(starts_with("ion"),
                function(x) x-mean(x))) %>%
  ungroup()


# %>%
#   mutate(across(c(gp, cohort,
#                   menopause_stt,
#                   ethnicity,
#                   edu_level,
#                   income_level,
#                   wh_contraceptives_ever,
#                   wh_hrt_ever,
#                   fam_hist_breast, fam_hist_breast_cat,
#                   alc_ever, alc_cur_freq_cat, alc_binge_cat,
#                   smk_cig_stt,bmi_cat,
#                   er_status, pr_status, her2_status, hr_subtype,
#                   hist_subtype),
#                 as_factor))

full_final <- full_centered %>%
  filter(sampletype=="Sample") %>%
  group_by(studyid) %>%
  summarise(across(starts_with("ion"), mean)) %>%
  ungroup() %>%
  full_join(full_dat)

