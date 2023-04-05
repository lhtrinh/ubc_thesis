library(tidyverse)
library(haven)
library(readxl)
library(lubridate)
library(psych)
library(moments)

#######################################################



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



bc_hormones <- bc_hormones %>%
  mutate(er_status=tolower(bccr_er_status_final),
         pr_status=tolower(bccr_pr_status_final),
         her2_status=tolower(bccr_her2_status_final))




# add "BCGP" to match ID to distinguish from ATP match IDs
bc_match_ind$match_id <- paste("bcgp", bc_match_ind$match_id, sep="_")



##############################################################
##############################################################
##############################################################
# remove duplicates in bc_raw and combine BCGP data
# bc_dat <- bc_raw %>% filter(dup=="NO") %>%
#   inner_join(bc_match_ind) %>%
#   left_join(bc_hormones)
##############################################################
##############################################################
##############################################################


bc_dat$cohort <- "bcgp"






###############################################################################






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



# atp unharmonized data
atp_cases_bio <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_ACR.csv")
colnames(atp_cases_bio) <- tolower(colnames(atp_cases_bio))




# match ID
atp_matched_pairs <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Pilot/2105_Study2105_03_29_2023_CaseControl_ParticipantKeyMatchedList.xlsx",
                               sheet="Sheet1")



#====================#
# add match IDs to matched pairs
atp_matched_pairs %>% head()
atp_match_ids <- atp_matched_pairs %>%
  rename(CNTL=CONTROL) %>%
  mutate(match_id=1:nrow(atp_matched_pairs)) %>%
  pivot_longer(!match_id,
               names_to="gp",
               values_to="participantkey") %>%
  mutate(match_id=paste("atp", match_id, sep="_"))
atp_match_ids %>% head()




#====================#
# remove duplicates
#====================#

atp_hormones <- atp_hormones %>%
  filter(acr_cancer_site_aggregate=="Breast") %>%
  arrange(participantkey, acr_malignancy_number)

# change variable names for HR/HER2 status and change all values to lowercase
atp_hormones <- atp_hormones %>%
  mutate(er_status=acr_estrogen_receptor_assay,
         pr_status=acr_progesterone_receptor_assay,
         her2_status=acr_her2_assay) %>%
  mutate(across(ends_with("_status"),
                function(x) ifelse(grepl("Positive", x), "positive",
                                   ifelse(grepl("Negative", x), "negative", NA))))

atp_hormones_new <- atp_hormones[!duplicated(atp_hormones$participantkey),]




# for cases, only select the first breast cancer diagnosis
atp_cases_bio <- atp_cases_bio %>%
  filter(acr_cancer_site_aggregate=="Breast") %>%
  arrange(participantkey, acr_diagnosis_year)
atp_cases_new <- atp_cases_bio[!duplicated(atp_cases_bio$participantkey),]




# combine ATP data
temp <- atp_hormones_new %>%
  inner_join(atp_cases_new)


atp_dat <- atp_raw %>%
  inner_join(atp_pm_measure) %>%
  left_join(atp_hormones_new) %>%
  left_join(atp_cases_new) %>%
  left_join(atp_match_ids)






#====================#
# remove 2 participants that are no longer eligible
# atp_dat <- atp_dat %>%
#   filter(!participantkey %in% c(210503833, 210524604))



# add control and case labels
atp_dat$gp <- "CNTL"
atp_dat$gp[atp_dat$participantkey %in% atp_hormones_new$participantkey] <- "CASE"



# change date form for questionnaire completion date
atp_dat$adm_qx_completion <- as.Date(atp_dat$adm_qx_completion)



#====================#
# Replace all -7 with NAs
#====================#
atp_dat[atp_dat==-7] <- NA


#====================#
# Create new variables from data
#====================#

atp_dat$cohort <- "atp"
atp_dat$studyid <- as.character(atp_dat$participantkey)
atp_dat$participantkey <- NULL


# rename age at diagnosis to match bcgp
atp_dat <- rename(atp_dat, age_at_diagnosis=acr_age_at_diagnosis)











###############################################################################



# COMBINE THE TWO DATA SETS

dat_raw <- bc_dat %>%
  full_join(atp_dat)




# change all match IDs








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
    TRUE ~ NA_character_
  ))


### Participant 210503316's tumour grade is false
temp1$hr_subgroup_alt[temp1$studyid=="210503316"]
### Currently marked as unknown, will keep that way


### Family history of breast cancer
temp1 <- temp1 %>%
  mutate(
    fam_hist_breast = ifelse(
      dis_cancer_fam_ever==0, 0,
      ifelse(dis_cancer_m_breast==1 | dis_cancer_sib_breast==1 | dis_cancer_child_breast, 1,
             NA)))






### alcohol consumption
temp1 <- temp1 %>%
  mutate(alc_cur_cat = case_when(
    alc_cur_freq == 0 ~ "never",
    alc_cur_freq %in% c(1,2,3,4) ~ "1 or less a week",
    alc_cur_freq == 5 ~ "2-3 times a week",
    alc_cur_freq == 6 ~ "4-5 times a week",
    alc_cur_freq == 7 ~ "6-7 times a week",
    TRUE ~ NA_character_
  ))

temp1 <- temp1 %>%
  mutate(alc_binge_cat = case_when(
    alc_binge_freq_female == 0 ~ "never",
    alc_binge_freq_female %in% c(1,2,3,4,5) ~ "1 or less a week",
    alc_binge_freq_female == 6 ~ "2-3 times a week",
    alc_binge_freq_female == 7 ~ "4-5 times a week",
    alc_binge_freq_female == 8 ~ "6-7 times a week",
    TRUE ~ NA_character_
  ))


# ethnicity
# temp1 <- temp1 %>%
#   mutate(ethnicity = do.call(coalesce,
#                              across(c(starts_with("sdc_eb_white"),
#                                       starts_with("sdc_eb_")
#                                       ))
#                              ))

temp1 <- temp1 %>%
  mutate(ethnicity=ifelse(sdc_eb_white==1, "white", "non_white"))




# education levels
temp1 <- temp1 %>%
  mutate(edu_level = case_when(
    sdc_edu_level %in% 1:2 ~ "high_school_or_less",
    sdc_edu_level %in% 3:5 ~ "trade_community_certificate",
    sdc_edu_level %in% 6:7 ~ "uni_graduate",
    TRUE ~ NA_character_
  ))


# income level
temp1 <- temp1 %>%
  mutate(income_level = case_when(
    sdc_income %in% 1:3 ~ "50k_or_less",
    sdc_income %in% 4:5 ~ "50k_to_100k",
    sdc_income %in% 6:8 ~ "100k_or_more",
    TRUE ~ NA_character_
  ))




# alcohol

temp1 <- temp1 %>%
  mutate(across(
    c(alc_cur_freq, alc_binge_freq_female, ends_with("week")),
    function(x) ifelse(alc_ever==0, 0, x)))




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



# set live births to 0 if gravidity=0
temp1 <- temp1 %>%
  mutate(wh_live_births = ifelse(wh_gravidity==0, 0, wh_live_births))


# set contraceptive duration to 0 if never used contraceptives
temp1 <- temp1 %>%
  mutate(wh_contraceptives_duration=ifelse(wh_contraceptives_ever==0, 0, wh_contraceptives_duration))



# impute menopause status for atp and combine with bc
temp1 <- temp1 %>%
  mutate(ms_bc=ifelse(str_detect(menopause_status, "Post"), 1, 0),
         ms_atp=ifelse(!is.na(wh_menopause_ever), wh_menopause_ever,
                       ifelse(sdc_age_calc>=51, 1, 0))
  ) %>%
  mutate(menopause_stt = coalesce(ms_bc, ms_atp))



###############################################

dat <- temp1
