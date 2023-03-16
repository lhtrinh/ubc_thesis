source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup_v2.R")
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup_atp.R")


dat <- bc_dat %>%
  full_join(atp_dat)

with(bc_dat, by(wh_menopause_age, gp, describe))
with(dat, by(sdc_age_calc, gp, describe))

with(dat, by(wh_menopause_age, gp, describe))
with(dat, by(age_at_diagnosis, gp, describe))



with(dat, by(wh_menstruation_age, gp, describe))

# number or pregnancies
with(dat, by(wh_gravidity, gp, describe))
with(dat, by(wh_preg_first_age, gp, describe))
with(dat, by(wh_live_births, gp, describe))


with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_gravidity, gp, summary))
with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_preg_first_age, gp, summary))
with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_live_births, gp, summary))

# hrt use
with(temp1, by(wh_hrt_ever, gp, table))

# age at hrt
with(dat[dat$wh_hrt_ever==1,], by(wh_hrt_age, gp, describe))

# contraceptives
temp1 %>% count(gp, wh_contraceptives_ever)
with(dat[dat$wh_contraceptives_ever==1,], by(wh_contraceptives_duration, gp, summary))

# diet
with(dat, by(nut_fruits_qty, gp, describe))
with(dat, by(nut_veg_qty, gp, describe))

# bmi
dat <- dat %>%
  mutate(bmi=as.numeric(coalesce(pm_tanitabmi, pm_bioimped_bmi, pm_bmi_sr))) %>%
  mutate(bmi_cat=factor(findInterval(bmi, c(18.5, 25, 30)),
                        labels=c("underweight", "normal", "overweight", "obese")))
with(dat, by(bmi, gp, describe))
