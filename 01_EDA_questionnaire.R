####################################################
###  CALCULATE SUMMARY STATISTICS FOR DATA
####################################################



source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load full questionnaire data
full_dat <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")



full_dat <- full_dat %>%
  mutate(across(c(gp, cohort, baseline_info_avail,
                  menopause_stt,
                  ethnicity,
                  edu_level, sdc_edu_level,
                  sdc_income, income_level,
                  wh_contraceptives_ever,
                  wh_hft_ever, wh_hrt_ever,
                  fam_hist_breast,
                  alc_ever, alc_cur_freq, alc_cur_freq_cat,
                  alc_binge_freq_female, alc_binge_cat,
                  smk_cig_status,bmi_cat,
                  er_status, pr_status, her2_status, hr_subtype,
                  hist_subtype),
  as_factor))


full_dat %>% count(gp)

full_dat %>% group_by(gp) %>% count(cohort) %>% mutate(p=n/sum(n))




#=================================================================#
# summary for numerical variables
# mean, SD, and count of missing
full_num_summ <- full_dat %>%
  group_by(gp) %>%
  summarise(across(c(sdc_age_calc, bmi,
         wh_menstruation_age,
         wh_gravidity,
         wh_preg_first_age,
         wh_live_births,
         wh_contraceptives_age,
         wh_contraceptives_duration_yr,
         wh_hrt_age,
         alc_drink_per_week,
         nut_veg_qty,
         nut_fruits_qty),
         c(avg=function(x) round(mean(x, na.rm=TRUE), 1),
           sd=function(x) round(sd(x, na.rm=TRUE), 1),
           nmiss=function(x) round(sum(is.na(x)), 1)))) %>%
  pivot_longer(-gp, names_to="var") %>%
  arrange(var, gp)
View(full_num_summ)


#=================================================================#
# summary for categorical variables
# non-missing count and percentage, plus count of missing with no pct
full_cat1 <- full_dat %>%
  select(gp,
         cohort,
         ethnicity, income_level, edu_level,
         bmi_cat,
         fam_hist_breast,
         menopause_stt,
         wh_contraceptives_ever,
         wh_hrt_ever,
         alc_cur_freq_cat,
         alc_binge_cat,
         smk_cig_status) %>%
  pivot_longer(-gp, names_to="var") %>%
  group_by(gp, var) %>%
  count(value)

full_cat_summ <- full_cat1 %>%
  filter(!is.na(value)) %>%
  mutate(pct=round(100*n/sum(n),1)) %>%
  ungroup() %>%
  full_join(full_cat1) %>%
  arrange(var, value, gp)

View(full_cat_summ)




#=================================================================#
# cancer characteristics
library(psych)

full_case <- full_cat %>%
  filter(gp=="CASE")

describe(full_case$age_at_diagnosis)

full_case %>%
  count(hist_subtype) %>%
  mutate(freq=round(100*n/sum(n), 1))


full_case %>%
  count(hr_subtype) %>%
  mutate(freq=round(100*n/sum(n), 1))
