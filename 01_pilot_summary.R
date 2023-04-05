#===========================================#
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/01_pilot_pretreatment.R")





table(pm_dat$gp)


# age at baseline
with(pm_dat, by(sdc_age_calc, gp, describe))


#ethnicity
pm_dat %>% count(gp, ethnicity)

# Menopause status
with(pm_dat, by(wh_menopause_ever, gp, table))


# age at menopause
with(pm_dat, by(wh_menopause_age, gp, describe))


# Age at cancer diagnosis
with(pm_dat, by(age_at_diagnosis, gp, describe))


# hormone receptor status subgroups
pm_dat %>% count(gp, hr_subgroup_alt)


# age at menarche
with(pm_dat, by(wh_menstruation_age, gp, describe))


# contraceptives
pm_dat %>% count(gp, wh_contraceptives_ever)

# age at first contraceptive use
with(pm_dat, by(wh_contraceptives_age, gp, describe))

# months of contraceptive use
with(pm_dat, by(wh_contraceptives_duration, gp, describe))


# number or pregnancies
pm_dat %>% count(gp, wh_gravidity>0)
with(pm_dat, by(wh_gravidity, gp, describe))
with(pm_dat[pm_dat$wh_gravidity>0,], by(wh_preg_first_age, gp, describe))
with(pm_dat[pm_dat$wh_gravidity>0,], by(wh_live_births, gp, describe))


# hrt use
pm_dat %>% count(gp, wh_hrt_ever)


# age at hrt
with(pm_dat[pm_dat$wh_hrt_ever==1,], by(wh_hrt_age, gp, describe))


# family history of breast cancer
pm_dat %>%
  group_by(gp, dis_cancer_fam_ever) %>%
  count(dis_cancer_m_breast)

pm_dat %>%
  group_by(gp, dis_cancer_fam_ever) %>%
  count(dis_cancer_sib_breast)

pm_dat %>%
  count(gp, fam_hist_breast)


# BMI
with(pm_dat, by(bmi, gp, describe))

pm_dat %>% count(gp, bmi_cat)


# alcohol use
pm_dat %>% count(gp, alc_cur_cat)


# binge drinking
pm_dat %>% count(gp, alc_binge_cat)


# smoking status
pm_dat %>% count(gp, smk_cig_status)

# duration of heavy smoking
with(pm_dat, by(smk_cig_dur, gp, describe))

# cigarettes per day during heaviest smoking duration
pm_dat %>% count(gp, smk_cig_qty)

pm_dat %>% filter(smk_cig_status!=0) %>%
  count(gp, smk_cig_qty)


# diet
with(pm_dat, by(nut_fruits_qty, gp, describe))
with(pm_dat, by(nut_veg_qty, gp, describe))


# education levels
pm_dat %>% count(gp, sdc_edu_level)

# income levels
pm_dat %>% count(gp, sdc_income)

#######################################
































