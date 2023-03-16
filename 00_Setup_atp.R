#====================#
## ATP data
#====================#

# load data sets
atp_raw <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_HFRQ_BL_CPAC.csv")
colnames(atp_raw) <- tolower(colnames(atp_raw))
atp_pm_measure <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_PM_BL_CPAC.csv")
colnames(atp_pm_measure) <- tolower(colnames(atp_pm_measure))
atp_hormones <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_ACR_Breast_DX_ER_PR.csv")
colnames(atp_hormones) <- tolower(colnames(atp_hormones))





# remove duplicates

er <- atp_hormones %>%
  select(participantkey, acr_estrogen_receptor_assay) %>%
  group_by(participantkey) %>%
  count(acr_estrogen_receptor_assay) %>%
  mutate(
    er_status = ifelse(
      str_detect(acr_estrogen_receptor_assay, "Positive"), "positive",
      ifelse(str_detect(acr_estrogen_receptor_assay, "Negative"), "negative",
             NA))) %>%
  select(participantkey, er_status) %>%
  arrange(participantkey, er_status) %>%
  fill(er_status, .direction = "down") %>%
  mutate(er_status=replace_na(er_status, "unknown")) %>%
  unique() %>%
  ungroup()


pr <- atp_hormones %>%
  select(participantkey, acr_progesterone_receptor_assay) %>%
  group_by(participantkey) %>%
  count(acr_progesterone_receptor_assay) %>%
  mutate(
    pr_status = ifelse(
      str_detect(acr_progesterone_receptor_assay, "Positive"), "positive",
      ifelse(str_detect(acr_progesterone_receptor_assay, "Negative"), "negative",
             NA))) %>%
  select(participantkey, pr_status) %>%
  arrange(participantkey, pr_status) %>%
  fill(pr_status, .direction = "down") %>%
  mutate(pr_status=replace_na(pr_status, "unknown")) %>%
  unique() %>%
  ungroup()


her2 <- atp_hormones %>%
  select(participantkey, acr_her2_assay) %>%
  group_by(participantkey) %>%
  count(acr_her2_assay) %>%
  mutate(
    her2_status = ifelse(
      str_detect(acr_her2_assay, "Positive"), "positive",
      ifelse(str_detect(acr_her2_assay, "Negative"), "negative",
             NA))) %>%
  select(participantkey, her2_status) %>%
  arrange(participantkey, her2_status) %>%
  fill(her2_status, .direction = "down") %>%
  mutate(her2_status=replace_na(her2_status, "unknown")) %>%
  unique()


atp_hormones_new <- er %>%
  inner_join(pr) %>%
  inner_join(her2)







# combine ATP data
atp_dat <- atp_raw %>%
  inner_join(atp_pm_measure) %>%
  left_join(atp_hormones_new)




# remove 2 participants that are no longer eligible
# atp_dat <- atp_dat %>%
#   filter(!participantkey %in% c(210503833, 210524604))



# add control and case labels
atp_dat$gp <- "CNTL"
atp_dat$gp[atp_dat$participantkey %in% atp_hormones_new$participantkey] <- "CASE"



# change date form
atp_dat$adm_qx_completion <- as.Date(atp_dat$adm_qx_completion)



########################################################################
#====================#
# Replace all -7 with NAs
#====================#
atp_dat[atp_dat==-7] <- NA


#====================#
# Create new variables from data
#====================#

atp_dat$cohort <- "atp"


### BMI and BMI categories
temp1 <- atp_dat %>%
  mutate(bmi=as.numeric(coalesce(pm_tanitabmi, pm_bioimped_bmi, pm_bmi_sr))) %>%
  mutate(bmi_cat=factor(findInterval(bmi, c(18.5, 25, 30)),
                        labels=c("underweight", "normal", "overweight", "obese")))

### Molecular subtype
atp_dat <- temp1 %>%
  mutate(hr_subgroup = case_when(
    (er_status=="positive" | pr_status=="positive") & her2_status=="negative" ~ "luminal_a",
    (er_status=="positive" | pr_status=="positive") & her2_status=="positive" ~ "luminal_b",
    er_status=="negative" & pr_status=="negative" & her2_status=="positive" ~ "her2_positive",
    er_status=="negative" & pr_status=="negative" & her2_status=="negative" ~ "triple_negative",
    TRUE ~ "unknown"
  ))

### Participant 210503316's tumour grade is false
atp_dat$hr_subgroup[atp_dat$participantkey==210503316]
### Currently marked as unknown, will keep that way



# ### Family history of breast cancer
# temp1 <- temp2 %>%
#   mutate(
#     fam_hist_breast1=ifelse(is.na(dis_cancer_m_breast) & is.na(dis_cancer_sib_breast), 0, 1)) %>%
#   mutate(
#     fam_hist_breast=ifelse(fam_hist_breast1==1, 1,
#                            ifelse(!is.na(dis_cancer_fam_ever), 0, NA)))
#
# ### alcohol consumption
# temp2 <- temp1 %>%
#   mutate(alc_cur_cat = case_when(
#     alc_cur_freq == 0 ~ "never",
#     alc_cur_freq %in% c(1,2,3,4) ~ "1 or less a week",
#     alc_cur_freq == 5 ~ "2-3 times a week",
#     alc_cur_freq == 6 ~ "4-5 times a week",
#     alc_cur_freq == 7 ~ "6-7 times a week"
#   ))
#
# temp1 <- temp2 %>%
#   mutate(alc_binge_cat = case_when(
#     alc_binge_freq_female == 0 ~ "never",
#     alc_binge_freq_female %in% c(1,2,3,4,5) ~ "1 or less a week",
#     alc_binge_freq_female == 6 ~ "2-3 times a week",
#     alc_binge_freq_female == 7 ~ "4-5 times a week",
#     alc_binge_freq_female == 8 ~ "6-7 times a week"
#   ))
#
#
# # ethnicity
# temp2 <- temp1 %>%
#   mutate(ethnicity = ifelse(sdc_eb_white==1, "white", "other"))
#
# # smoking
# temp1 <- temp2 %>%
#   mutate(
#     smk_cig_dur = as.numeric(coalesce(
#       smk_cig_daily_cur_dur,
#       smk_cig_former_daily_dur,
#       smk_cig_heaviest_dur)),
#     smk_cig_qty = coalesce(
#       smk_cig_daily_avg_qty,
#       smk_cig_former_daily_qty,
#       smk_cig_heaviest_qty
#     ))
#
#
#
#
#
#
#
#
#
#
# ###############################################################
#
# table(atp_dat$gp)
#
# # BMI
# with(temp1, by(bmi, gp, mean, na.rm=TRUE))
# with(temp1, by(bmi, gp, describe))
#
# with(temp1, by(bmi_cat, gp, summary))
#
#
#
# # Menopause status
# with(temp1, by(wh_menopause_ever, gp, table))
#
# # age at menopause
# with(temp1, by(wh_menopause_age, gp, describe))
#
#
# # Age at cancer diagnosis
# with(temp1, by(age_at_diagnosis, gp, describe))
#
#
#
# # hormone receptor status subgroups
# with(temp2, by(sdc_age_calc, gp, describe))
#
# # age at menarche
# with(temp1, by(wh_menstruation_age, gp, describe))
#
# # number or pregnancies
# with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_gravidity, gp, summary))
# with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_preg_first_age, gp, summary))
# with(atp_dat[atp_dat$wh_gravidity>0,], by(wh_live_births, gp, summary))
#
#
# # hrt use
# with(temp1, by(wh_hrt_ever, gp, table))
#
# # age at hrt
# with(temp1[temp1$wh_hrt_ever==1,], by(wh_hrt_age, gp, describe))
#
# # family history of breast cancer
# temp1 %>%
#   group_by(gp, dis_cancer_fam_ever) %>%
#   count(dis_cancer_m_breast)
#
# temp1 %>%
#   group_by(gp, dis_cancer_fam_ever) %>%
#   count(dis_cancer_sib_breast)
#
# temp1 %>%
#   count(gp, fam_hist_breast)
#
#
# # alcohol use
# temp2 %>% count(gp, alc_cur_cat)
#
# temp1 %>% count(gp, alc_binge_cat)
#
# # smoking
# temp1 %>% count(gp, smk_cig_status)
# temp1 %>% count(gp, smk_cig_daily_cur_qty)
# temp1 %>% count(gp, smk_cig_daily_avg_qty)
#
# with(temp1, by(smk_cig_dur, gp, describe))
# temp1 %>% count(gp, smk_cig_qty)
#
# temp1 %>% filter(smk_cig_status %in% c(1,2,3)) %>%
#   count(gp, smk_cig_status)
#
# with(temp1[temp1$smk_cig_status %in% c(1,2,3),], by(smk_cig_dur, gp, summary))
#
#
# temp1 %>% filter(smk_cig_status!=0) %>%
#   count(gp, smk_cig_qty)
#
#
# #ethnicity
# temp2 %>% count(gp, ethnicity)
#
# # contraceptives
# temp1 %>% count(gp, wh_contraceptives_ever)
# with(temp1, by(wh_contraceptives_duration, gp, describe))
#
# # diet
# with(temp1, by(nut_fruits_qty, gp, describe))
# with(temp1, by(nut_veg_qty, gp, describe))

