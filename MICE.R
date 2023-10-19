#======================================================#
library(mice)
library(tidyverse)


full_all <- full_all %>%
  mutate(
    across(c(
      gp, studyid, match_id, cohort,
      collect_yr_cat,
      menopause_stt,
      ethnicity,
      edu_level,
      income_level,
      wh_menstruation_age_cat, wh_menopause_age_cat,
      wh_contraceptives_ever,
      wh_preg_first_age_cat,
      wh_breastfeeding_cat,
      wh_hrt_ever,
      fam_hist_breast,
      alc_ever, alc_cur_freq_cat, alc_binge_cat,
      smk_cig_status, smk_first_age_cat,
      bmi_cat,
      er_status, pr_status, her2_status, hr_subtype,
      site_code, hist_code, hist_subtype),
    as.factor)) %>%
  mutate(
    across(c(
      sdc_age_calc, collect_age,
      wh_gravidity, wh_live_births, wh_hrt_duration_yr,
      smk_cig_dur,
      bmi,
      pse_childhood_duration, pse_adult_home_duration, pse_adult_wrk_duration,
      age_at_diagnosis,
      starts_with("ion")),
    as.numeric))

dat_pre_imp_id <- dat_pre_imp %>% select(all_of(id_cols))


# summary(dat_pre_imp)

# remove information not for imputing
dat_for_imp <- dat_pre_imp %>%
  select(gp, all_of(match_cols), all_of(context_cols))
summary(dat_for_imp)

#======================================================#

# imputation

init <- mice(dat_for_imp,
             m=5,
             maxit=5,
             seed=292920)


# record all imputed data sets
dat_imp <- list()
for (i in 1:5){
  full_dat <- cbind(dat_pre_imp_id, complete(init, action=i))
  dat_imp[[i]] <- full_join(full_centered, full_dat)
}