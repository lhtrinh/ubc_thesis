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


# ### Molecular subtype
# temp1 <- temp1 %>%
#   mutate(hr_subgroup = case_when(
#     (er_status=="positive" | pr_status=="positive") & her2_status=="negative" ~ "luminal_a",
#     (er_status=="positive" | pr_status=="positive") & her2_status=="positive" ~ "luminal_b",
#     er_status=="negative" & pr_status=="negative" & her2_status=="positive" ~ "her2_positive",
#     er_status=="negative" & pr_status=="negative" & her2_status=="negative" ~ "triple_negative",
#     TRUE ~ "unknown"
#   ))


# another breast cancer subtype
temp1 <- temp1 %>%
  mutate(hr_subgroup_alt=case_when(
    er_status=="negative" & pr_status=="negative" & her2_status=="negative" ~ "triple_negative",
    her2_status=="positive" ~ "her2_positive",
    (er_status=="positive" | pr_status=="positive") ~ "er/pr_positive",
    gp=='CNTL' ~ "control",
    TRUE ~ "unknown"
  ))


### Participant 210503316's tumour grade is false
temp1$hr_subgroup_alt[temp1$studyid=="210503316"]
### Currently marked as unknown, will keep that way


### Family history of breast cancer
temp1 <- temp1 %>%
  mutate(
    fam_hist_breast1=ifelse(is.na(dis_cancer_m_breast) & is.na(dis_cancer_sib_breast), 0, 1)) %>%
  mutate(
    fam_hist_breast=ifelse(fam_hist_breast1==1, 1,
                           ifelse(!is.na(dis_cancer_fam_ever), 0, NA)))

### alcohol consumption
temp1 <- temp1 %>%
  mutate(alc_cur_cat = case_when(
    alc_cur_freq == 0 ~ "never",
    alc_cur_freq %in% c(1,2,3,4) ~ "1 or less a week",
    alc_cur_freq == 5 ~ "2-3 times a week",
    alc_cur_freq == 6 ~ "4-5 times a week",
    alc_cur_freq == 7 ~ "6-7 times a week"
  ))

temp1 <- temp1 %>%
  mutate(alc_binge_cat = case_when(
    alc_binge_freq_female == 0 ~ "never",
    alc_binge_freq_female %in% c(1,2,3,4,5) ~ "1 or less a week",
    alc_binge_freq_female == 6 ~ "2-3 times a week",
    alc_binge_freq_female == 7 ~ "4-5 times a week",
    alc_binge_freq_female == 8 ~ "6-7 times a week"
  ))


# ethnicity
temp1 <- temp1 %>%
  mutate(ethnicity = ifelse(sdc_eb_white==1, "white", "other"))

# smoking
temp1 <- temp1 %>%
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

