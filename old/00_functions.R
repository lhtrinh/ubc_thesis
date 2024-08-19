#======================================================#
#
###    MANUALLY DEFINED FUNCTIONS FOR METABOLOMCIS   ###
#
#======================================================#
library(tidyverse)
library(lme4)
library(performance)
library(readr)


#=====================#
# GROUP VARIABLES ####

y_col <- "gp"

id_cols <- c("studyid", "match_id")

match_cols <-  c("cohort", "collect_yr_cat", "collect_age", "menopause_stt")

context_cols <- c("sdc_age_calc", "ethnicity",
                  "edu_level", "income_level",
                  "wh_menstruation_age_cat", "wh_menopause_age_cat",
                  "wh_contraceptives_ever",
                  "wh_gravidity", "wh_live_births", "wh_preg_first_age_cat",
                  "wh_breastfeeding_cat",
                  "wh_hrt_ever", "wh_hrt_duration_yr",
                  "fam_hist_breast",
                  "bmi", "bmi_cat",
                  "alc_ever", "alc_cur_freq_cat", "alc_binge_cat",
                  "smk_cig_status", "smk_cig_dur", "smk_first_age_cat",
                  "pse_childhood_duration", "pse_adult_home_duration", "pse_adult_wrk_duration")

imp_cols <- c(y_col, match_cols, context_cols)






fct_cols <- c("gp", "cohort",
              "menopause_stt", "collect_yr_cat",
              "ethnicity",
              "edu_level",
              "income_level",
              "wh_menstruation_age_cat", "wh_menopause_age_cat",
              "wh_contraceptives_ever",
              "wh_breastfeeding_cat",
              "wh_preg_first_age_cat",
              "wh_hrt_ever",
              "fam_hist_breast",
              "alc_ever", "alc_cur_freq_cat", "alc_binge_cat",
              "smk_cig_status", "smk_first_age_cat",
              "bmi_cat",
              "er_status", "pr_status", "her2_status", "hr_subtype",
              "site_code", "hist_code", "hist_subtype")

int_cols <- c("sdc_age_calc", "wh_gravidity", "wh_live_births")


num_cols <- c("sdc_age_calc", "collect_age", "wh_gravidity",
              "wh_live_births", "wh_hrt_duration_yr", "smk_cig_dur",
              "pse_childhood_duration", "pse_adult_home_duration", "pse_adult_wrk_duration",
              "bmi", "age_at_diagnosis")



n_imp <- 10



#=====================#




# IMPORT DATA ####
import_survey <- function(){
  full_dat <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv") %>%
    mutate(across(all_of(fct_cols), as.factor)) %>%
    mutate(across(all_of(int_cols), as.integer))

  full_dat
}



# Import normalized metabolomics data
import_meta_paired <- function(){
  full_centered <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/normalized_metabolomics.csv")

  full_centered
}

import_meta_sc <- function(){
  full_sc <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/normalized_metabolomics_scaled.csv")
}



# Combine full questionnaire data with metabolomics data
# combine_data <- function() full_join(import_survey(), import_meta())





## Import imputed data

import_imp <- function() {
  ## Impute data
  path <- Sys.glob("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/*.csv")
  boot_list_import <- import_list(path)


  for (i in 1:length(boot_list_import)){
    boot_list_import[[i]][boot_list_import[[i]]==""] <- NA
    boot_list_import[[i]][fct_cols] <- lapply(boot_list_import[[i]][fct_cols], as.factor)
    boot_list_import[[i]][int_cols] <- lapply(boot_list_import[[i]][int_cols], as.integer)
  }

  boot_list_import
}



#==================================================================#

# CALCULATE CV #########

# cv_calc: calculate coefficient of variation (CV%)
# per metabolite per participant by dividing each
# metabolite's SD by its mean and multiply by 100%
calc_cv <- function(dat, group_var){

  group_var <- enquo(group_var)

  # calculate CV per metabolite per participant
  replicate_cv <- dat %>%
    filter(dup=="YES") %>%
    group_by(!!group_var) %>%
    summarise(across(starts_with("ion"),
                     function(x) sd(x)/mean(x))) %>%
    ungroup()

  return(replicate_cv)
}





#=====================#
# cv_calc_med
# calculate median coefficient of variation
# across each metabolite
calc_cv_med <- function(dat, group_var){
  group_var <- enquo(group_var)

  cv_dat <- calc_cv(dat, !!group_var)


  replicate_cv_med <- cv_dat %>%
    summarise(across(starts_with("ion"), function(x) median(x, na.rm=TRUE))) %>%
    t()

  # record the CV into a dataset
  ion_cv_dat <- data.frame(metabolite=rownames(replicate_cv_med),
                           median_cv=replicate_cv_med[,1])
  rownames(ion_cv_dat) <- NULL
  return(ion_cv_dat)
}






#=====================#
# CALCULATE ICC ########


calc_icc <- function(dat, group_var){

  # select data with duplicates
  dat_dup <- dat %>% filter(dup=="YES")

  # define list of ion names
  ion_cols <- names(dat_dup)[str_detect(names(dat_dup), "ion")]

  # create placeholder data frame to record ICC values
  icc_dat <- data.frame(
    metabolite = ion_cols,
    icc = 0)

  # calculate ICC for each metabolite
  for (i in 1:length(ion_cols)){
    ion_col <- ion_cols[i]
    ion_col <- enquo(ion_col)
    group_var <- enquo(group_var)

    ion_dat <- dat_dup %>%
      select(!!group_var, !!ion_col)
    colnames(ion_dat) <- c("id", "y")
    my_mod <- lme4::lmer(y ~ 1 + (1|id), data=ion_dat)
    my_summ <- summary(my_mod)

    my_icc <- my_summ$varcor$id[1] / (my_summ$varcor$id[1] + my_summ$sigma^2)

    # icc_dat[i,1] <- ion_col
    icc_dat[i,2] <- my_icc
  }
  icc_dat
}






#=====================#
# TEST FOR NORMALITY OF METABOLITE LEVELS #######



# test for normality
normality_skewness <- function(df) {
  norm_df <- data.frame(metabolite=rep(0, ncol(df)),
                        shapiro_pval=rep(0,ncol(df)),
                        jarque_pval=rep(0,ncol(df)),
                        ks_pval=rep(0, ncol(df)),
                        skew=rep(0, ncol(df)))
  ion_df <- df %>% select(starts_with("ion"))
  for (i in 1:ncol(ion_df)){
    x <- ion_df[,i] %>% as_vector()

    # for all three normality tests, p>.05 indicates normality
    norm_df[i,1] <- colnames(ion_df)[i]
    norm_df[i,2] <- shapiro.test(x)$p.value
    norm_df[i,3] <- jarque.test(x)$p.value
    norm_df[i,4] <- ks.test(x, "pnorm")$p.value
    norm_df[i,5] <- 3*(mean(x) - median(x))/sd(x)
  }
  norm_df
}




#=====================#
# QC FOR NORMALIZATION ####



# norm_qc function:
# - calculates median CV and ICC and saves in a data set
# - does PCA and saves biplot by cohort and study plate and saves image
norm_qc <- function(df){
  df_dup <- df %>%
    filter(dup=="YES")

  cv_norm <- calc_cv_med(df_dup, studyid)
  icc_norm <- calc_icc(df_dup, studyid)

  # save median CV and ICC for each metabolite and add tertile info
  qc_norm <- cv_norm %>%
    full_join(icc_norm)

  # PCA
  df_for_pca <- df %>%
    filter(sampletype=="Sample") %>%
    mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))

  ion_df_pca <- df_for_pca %>% select(starts_with("ion"))
  pca_mod <- prcomp(ion_df_pca)

  fig_cohort <- fviz_pca_ind(pca_mod,
                             label="none",
                             alpha=.5,
                             habillage = df_for_pca$cohort,
                             palette=c("blue", "red"),
                             addEllipses = TRUE)

  fig_plate <- fviz_pca_ind(pca_mod,
                            label="none",
                            habillage =  df_for_pca$plate,
                            alpha=.5,
                            addEllipses = TRUE)

  # save PCA biplots
  figs <- grid.arrange(fig_cohort, fig_plate, nrow=1)
  ggsave(filename=paste("pca_", deparse(substitute(df)), ".jpg", sep=""), plot=figs)

  return(qc_norm)
}





#=====================#
# CALCULATE IMPUTED CI (RUBIN'S RULE) ####
imputed_ci <- function(x){
  # x is a vector of estimates


  # calculate pooled estimate
  pooled_est <- mean(x)

  # calculate pooled SE
  m <- length(x)
  se <- sd(x)/sqrt(m)

  within_var <- mean(se^2)
  btwn_var <- sum((x-pooled_est)^2)/(m-1)
  total_var <- within_var + btwn_var + btwn_var/m

  pooled_se <- sqrt(total_var)

  # calculate t-statistic for 95% CI
  alpha <- .05
  deg_free <- m-1
  t.score <- qt(p=alpha/2, df=deg_free, lower.tail = FALSE)

  # calculate 95% CI lower and upper bounds
  ci_lb <- pooled_est - t.score*pooled_se
  ci_ub <- pooled_est + t.score*pooled_se

  c(pooled_est, pooled_se, ci_lb, ci_ub)
}