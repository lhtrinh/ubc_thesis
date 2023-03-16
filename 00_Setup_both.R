source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup_v2.R")
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup_atp.R")


dat_raw <- bc_dat %>%
  full_join(atp_dat)




#====================#
# Replace all -7 with NAs
#====================#
dat_raw[dat_raw==-7] <- NA


#====================#
# Create new variables from data
#====================#



### BMI and BMI categories
temp1 <- dat_raw %>%
  mutate(bmi=as.numeric(coalesce(pm_tanitabmi, pm_bioimped_bmi, pm_bmi_sr))) %>%
  mutate(bmi_cat=factor(findInterval(bmi, c(18.5, 25, 30)),
                        labels=c("underweight", "normal", "overweight", "obese")))


### Molecular subtype
temp2 <- temp1 %>%
  mutate(hr_subgroup = case_when(
    (bccr_er_status_final=="Positive" | bccr_pr_status_final=="Positive") & bccr_her2_status_final=="Negative" ~ "luminal_a",
    (bccr_er_status_final=="Positive" | bccr_pr_status_final=="Positive") & bccr_her2_status_final=="Positive" ~ "luminal_b",
    bccr_er_status_final=="Negative" & bccr_pr_status_final=="Negative" & bccr_her2_status_final=="Positive" ~ "her2_positive",
    bccr_er_status_final=="Negative" & bccr_pr_status_final=="Negative" & bccr_her2_status_final=="Negative" ~ "triple_negative",
    TRUE ~ "unknown"
  ))


### Family history of breast cancer
temp1 <- temp2 %>%
  mutate(
    fam_hist_breast1=ifelse(is.na(dis_cancer_m_breast) & is.na(dis_cancer_sib_breast), 0, 1)) %>%
  mutate(
    fam_hist_breast=ifelse(fam_hist_breast1==1, 1,
                           ifelse(!is.na(dis_cancer_fam_ever), 0, NA)))

### alcohol consumption
temp2 <- temp1 %>%
  mutate(alc_cur_cat = case_when(
    alc_cur_freq == 0 ~ "never",
    alc_cur_freq %in% c(1,2,3,4) ~ "1 or less a week",
    alc_cur_freq == 5 ~ "2-3 times a week",
    alc_cur_freq == 6 ~ "4-5 times a week",
    alc_cur_freq == 7 ~ "6-7 times a week"
  ))

temp1 <- temp2 %>%
  mutate(alc_binge_cat = case_when(
    alc_binge_freq_female == 0 ~ "never",
    alc_binge_freq_female %in% c(1,2,3,4,5) ~ "1 or less a week",
    alc_binge_freq_female == 6 ~ "2-3 times a week",
    alc_binge_freq_female == 7 ~ "4-5 times a week",
    alc_binge_freq_female == 8 ~ "6-7 times a week"
  ))


# ethnicity
temp2 <- temp1 %>%
  mutate(ethnicity = ifelse(sdc_eb_white==1, "white", "other"))

# smoking
temp1 <- temp2 %>%
  mutate(
    smk_cig_dur = as.numeric(coalesce(
      smk_cig_daily_cur_dur,
      smk_cig_former_daily_dur,
      smk_cig_heaviest_dur)),
    smk_cig_qty = coalesce(
      smk_cig_daily_avg_qty,
      smk_cig_former_daily_qty,
      smk_cig_heaviest_qty
    ))




dat <- temp1



##############################################################
# temp1 %>% count(smk_cig_status)
#
# with(bc_dat, by(wh_menopause_age, gp, describe))
# with(dat, by(sdc_age_calc, gp, describe))
#
# with(dat, by(wh_menopause_age, gp, describe))
# with(dat, by(age_at_diagnosis, gp, describe))
#
#
#
# with(dat, by(wh_menstruation_age, gp, describe))
#
# # number or pregnancies
# with(dat, by(wh_gravidity, gp, describe))
# with(dat, by(wh_preg_first_age, gp, describe))
# with(dat, by(wh_live_births, gp, describe))
#
#
# with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_gravidity, gp, summary))
# with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_preg_first_age, gp, summary))
# with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_live_births, gp, summary))
#
# # hrt use
# with(temp1, by(wh_hrt_ever, gp, table))
#
# # age at hrt
# with(dat[dat$wh_hrt_ever==1,], by(wh_hrt_age, gp, describe))
#
# # contraceptives
# temp1 %>% count(gp, wh_contraceptives_ever)
# with(dat[dat$wh_contraceptives_ever==1,], by(wh_contraceptives_duration, gp, summary))
#
# # diet
# with(dat, by(nut_fruits_qty, gp, describe))
# with(dat, by(nut_veg_qty, gp, describe))
#
#
# with(dat, by(bmi, gp, describe))







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