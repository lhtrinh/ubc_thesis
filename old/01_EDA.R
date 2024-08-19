source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup_both.R")



##############################################################
table(dat$gp)

# BMI
with(dat, by(bmi, gp, mean, na.rm=TRUE))
with(dat, by(bmi, gp, describe))

with(dat, by(bmi_cat, gp, summary))



# Menopause status
with(dat, by(wh_menopause_ever, gp, table))

# age at menopause
with(dat, by(wh_menopause_age, gp, describe))


# Age at cancer diagnosis
with(dat, by(age_at_diagnosis, gp, describe))



# hormone receptor status subgroups
with(dat, by(sdc_age_calc, gp, describe))

# age at menarche
with(dat, by(wh_menstruation_age, gp, describe))

# number or pregnancies
with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_gravidity, gp, summary))
with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_preg_first_age, gp, summary))
with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_live_births, gp, summary))


# hrt use
with(dat, by(wh_hrt_ever, gp, table))

# age at hrt
with(dat[dat$wh_hrt_ever==1,], by(wh_hrt_age, gp, describe))

# family history of breast cancer
dat %>%
  group_by(gp, dis_cancer_fam_ever) %>%
  count(dis_cancer_m_breast)

dat %>%
  group_by(gp, dis_cancer_fam_ever) %>%
  count(dis_cancer_sib_breast)

dat %>%
  count(gp, fam_hist_breast)


# alcohol use
dat %>% count(gp, alc_cur_cat)

dat %>% count(gp, alc_binge_cat)

# smoking
dat %>% count(gp, smk_cig_status)
dat %>% count(gp, smk_cig_daily_cur_qty)
dat %>% count(gp, smk_cig_daily_avg_qty)

with(dat, by(smk_cig_dur, gp, describe))
dat %>% count(gp, smk_cig_qty)

dat %>% filter(smk_cig_status %in% c(1,2,3)) %>%
  count(gp, smk_cig_status)

with(dat[dat$smk_cig_status %in% c(1,2,3),], by(smk_cig_dur, gp, summary))


dat %>% filter(smk_cig_status!=0) %>%
  count(gp, smk_cig_qty)


#ethnicity
dat %>% count(gp, ethnicity)

# contraceptives
dat %>% count(gp, wh_contraceptives_ever)
with(dat, by(wh_contraceptives_duration, gp, describe))

# diet
with(dat, by(nut_fruits_qty, gp, describe))
with(dat, by(nut_veg_qty, gp, describe))