library(tidyverse)
library(haven)
library(readxl)
library(lubridate)
library(psych)


#######################################################
## LOAD DATA



#====================#
## BCGP data
#====================#

# load data sets
bc_raw <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_DATA_PT_Qx_RegistryInfo_BIO_2_0.sav")
colnames(bc_raw) <- tolower(colnames(bc_raw))
bc_hormones <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_Y2022_BCCR_incidence_breast_hormones.sav")
colnames(bc_hormones) <- tolower(colnames(bc_hormones))
bc_match_ind <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_MatchIndicator.csv")
colnames(bc_match_ind) <- tolower(colnames(bc_match_ind))


# remove duplicates in bc_raw and combine BCGP da=ta
bc_dat <- bc_raw %>% filter(dup=="NO") %>%
  inner_join(bc_match_ind) %>%
  left_join(bc_hormones)




#######################################################
###### CLEANING ###########




# #====================#
# # Format SPSS date formats
# #====================#
# temp1 <- bc_dat %>%
#   mutate(across(.cols=c(adm_qx_completion, f1_adm_qx_completion),
#                 .fns = function(x) as.Date(x/86400, origin = "1582-10-14")))
#
# head(temp1)

#====================#
# Replace all -7 with NAs
#====================#
bc_dat[bc_dat==-7] <- NA


#====================#
# Create new variables from data
#====================#



### BMI and BMI categories
temp1 <- bc_dat %>%
  mutate(bmi=as.numeric(coalesce(pm_tanitabmi, pm_bioimped_bmi, pm_bmi_sr))) %>%
  mutate(bmi_cat=factor(findInterval(bmi, c(18.5, 25, 30)),
                        labels=c("underweight", "normal", "overweight", "obese")))

### Molecular subtype
bc_dat <- temp1 %>%
  mutate(hr_subgroup = case_when(
    (bccr_er_status_final=="Positive" | bccr_pr_status_final=="Positive") & bccr_her2_status_final=="Negative" ~ "luminal_a",
    (bccr_er_status_final=="Positive" | bccr_pr_status_final=="Positive") & bccr_her2_status_final=="Positive" ~ "luminal_b",
    bccr_er_status_final=="Negative" & bccr_pr_status_final=="Negative" & bccr_her2_status_final=="Positive" ~ "her2_positive",
    bccr_er_status_final=="Negative" & bccr_pr_status_final=="Negative" & bccr_her2_status_final=="Negative" ~ "triple_negative",
    TRUE ~ "unknown"
  ))


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
# table(bc_dat$gp)
#
# # BMI
# with(temp1, by(bmi, gp, mean, na.rm=TRUE))
# with(temp1, by(bmi, gp, describe))
#
# with(temp1, by(bmi_cat, gp, summary))
# with(temp1, proportions(table(gp, bmi_cat), 1))
#
#
# # Menopause status
# with(temp1, by(menopause_status, gp, table))
#
#
# # Age at cancer diagnosis
# with(temp1, by(age_at_diagnosis, gp, describe))
#
#
#
# # hormone receptor status subgroups
# with(temp2, by(sdc_age_calc, gp, summary))
#
# # age at menarche
# with(temp1, by(wh_menstruation_age, gp, describe))
#
# # number or pregnancies
# with(temp1, by(wh_gravidity, gp, describe))
# with(temp1, by(wh_preg_first_age, gp, describe))
# with(temp1, by(wh_live_births, gp, describe))
#
# # hrt use
# with(temp1, by(wh_hrt_ever, gp, summary))
#
# # age at hrt
# with(temp1, by(wh_hrt_age, gp, describe))
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
# #ethnicity
# temp2 %>% count(gp, ethnicity)
#
# # contraceptives
# temp1 %>% count(gp, wh_contraceptives_ever)
# with(temp1, by(wh_contraceptives_duration, gp, summary))
#
# # diet
# with(temp1, by(nut_veg_qty, gp, describe))
# with(temp1, by(nut_fruits_qty, gp, describe))
