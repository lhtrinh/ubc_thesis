### Descriptive statistics ###


num_summ <- function(x){
  x <- enquo(x)
  full_dat %>%
    filter(!is.na(!!x)) %>%
    group_by(gp) %>%
    summarise(meanx = mean(!!x),
              sdx = sd(!!x))
}

cat_summ <- function(x) {
  x <- enquo(x)
  full_dat %>%
    filter(!is.na(!!x)) %>%
    group_by(gp) %>%
    count(!!x) %>%
    mutate(p=proportions(n))
}


### MATCHED VARIABLES
# This section checks matching variables between cases and controls
# to make sure they're matched with no significant differences

### Cohort count
full_dat %>% count(cohort, gp)
#227 pairs from ATP, 348 pairs from BCGP

### menopause status
full_dat %>% count(menopause_stt) %>% mutate(p=proportions(n))
# all matched




### DEMOGRAPHICS

# Age
with(full_dat, by(sdc_age_calc, gp, describe))
num_summ(sdc_age_calc)

# ethnicity
full_dat %>% count(ethnicity) %>% mutate(p=proportions(n))
cat_summ(ethnicity)

# income level
cat_summ(income_level)

# education level
cat_summ(edu_level)




### REPRODUCTIVE AND HEALTH

# BMI
num_summ(bmi)
cat_summ(bmi_cat)

# family history
cat_summ(fam_hist_breast)

# menopause status
cat_summ(menopause_stt)

# age at menarche
num_summ(wh_menstruation_age)

# gravidity
num_summ(wh_gravidity)

# live births
num_summ(wh_live_births)

# contraceptive use
cat_summ(wh_contraceptives_ever)

# use of hrt
cat_summ(wh_hrt_ever)



### LIFESTYLE BEHAVIORS

# alcohol
cat_summ(alc_ever)

cat_summ(alc_cur_freq)
cat_summ(alc_cur_freq_cat)

num_summ(alc_drink_per_week)

cat_summ(alc_binge_cat)


# smoking
cat_summ(smk_cig_status)

# nutrition
num_summ(nut_veg_qty)

num_summ(nut_fruits_qty)
