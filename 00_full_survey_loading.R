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


# remove one participant from BCGP (selected in error)
bc_raw2 <- bc_raw %>%
  filter(studyid!="PB000429")



# shorten variable names for HR status
bc_hormones <- bc_hormones %>%
  mutate(er_status=tolower(bccr_er_status_final),
         pr_status=tolower(bccr_pr_status_final),
         her2_status=tolower(bccr_her2_status_final))



# add "BCGP" to match ID to distinguish from ATP match IDs
bc_match_ind$match_id <- paste("bcgp", bc_match_ind$match_id, sep="_")
# remove duplicates
bc_match_ind <- bc_match_ind[!duplicated(bc_match_ind),]




# remove one set of duplicates in bc_raw and combine BCGP data
bc_raw2 <- bc_raw2 %>% arrange(studyid)
dup_ind <- which(duplicated(bc_raw2$studyid))
dup_ind


bc_raw3 <- bc_raw2[-dup_ind,]



bc_raw4 <- bc_raw3 %>%
  inner_join(bc_match_ind) %>%
  left_join(bc_hormones)


bc_raw4$cohort <- "bcgp"
dim(bc_raw4)



# remove pairs where at least one had no baseline questionnaire
bc_raw4 %>% count(is.na(adm_qx_completion)) #10 participants with no baseline qx

no_baseline <- bc_raw4 %>%
  filter(is.na(adm_qx_completion)) %>%
  count(match_id)


bc_raw5 <- bc_raw4 %>%
  mutate(baseline_info_avail=ifelse(match_id %in% no_baseline$match_id, "NO", "YES"))


bc_dat <- bc_raw5


dim(bc_dat)
n_distinct(bc_dat$match_id)


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
atp_cases <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_ACR.csv")
colnames(atp_cases) <- tolower(colnames(atp_cases))


# match ID
atp_matched_pairs <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/spreadsheets/2105_Study2105_03_29_2023_CaseControl_ParticipantKeyMatchedList.xlsx",
                               sheet="Sheet1")


# biospecimens data
atp_bio <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_BIO.csv")


#====================#
# add match IDs to matched pairs
atp_matched_pairs %>% head()
atp_match_id <- atp_matched_pairs %>%
  rename(CNTL=CONTROL) %>%
  mutate(match_id=1:nrow(atp_matched_pairs)) %>%
  pivot_longer(!match_id,
               names_to="gp",
               values_to="participantkey") %>%
  mutate(match_id=paste("atp", match_id, sep="_"))
atp_match_id %>% head()




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
atp_cases <- atp_cases %>%
  filter(acr_cancer_site_aggregate=="Breast") %>%
  arrange(participantkey, acr_diagnosis_year)
atp_cases_new <- atp_cases[!duplicated(atp_cases$participantkey),]




# combine ATP data
atp_raw2 <- atp_raw %>%
  inner_join(atp_pm_measure) %>%
  left_join(atp_hormones_new) %>%
  left_join(atp_cases_new) %>%
  left_join(atp_match_id)






#====================#
# remove 2 participants that are no longer eligible
# atp_dat %>%
#   filter(!participantkey %in% c(210503833, 210524604))






# change date form for questionnaire completion date
atp_raw2$adm_qx_completion <- as.Date(atp_dat$adm_qx_completion)



#====================#
# Replace all -7 with NAs
#====================#
atp_raw2[atp_raw2==-7] <- NA


#====================#
# Create new variables from data
#====================#

atp_raw2$cohort <- "atp"
atp_raw2$studyid <- as.character(atp_raw2$participantkey)
atp_raw2$participantkey <- NULL


# rename age at diagnosis to match bcgp
atp_raw2 <- rename(atp_raw2, age_at_diagnosis=acr_age_at_diagnosis)





atp_dat <- atp_raw2





###############################################################################



# COMBINE THE TWO DATA SETS

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
temp1 %>% count(bmi_cat) #84 NAs
summary(temp1$bmi)



# another breast cancer subtype
temp1 <- temp1 %>%
  mutate(hr_subtype=case_when(
    er_status=="negative" & pr_status=="negative" & her2_status=="negative" ~ "triple_negative",
    her2_status=="positive" ~ "her2_positive",
    (er_status=="positive" | pr_status=="positive") ~ "hr_positive",
    gp=='CNTL' ~ "control",
    TRUE ~ NA_character_
  ))
temp1 %>% count(hr_subtype)





### Family history of breast cancer
temp1 <- temp1 %>%
  mutate(
    cancer_fam_breast = ifelse(dis_cancer_m_breast==1 | dis_cancer_sib_breast==1 | dis_cancer_child_breast==1, 1, 0)) %>%
  mutate(
    fam_hist_breast = case_when(
      cancer_fam_breast==1 ~ 1,
      is.na(cancer_fam_breast) & !is.na(dis_cancer_fam_ever) ~ 0,
      TRUE ~ NA
    )
  )

temp1 %>% count(fam_hist_breast)
#############





### alcohol consumption
temp1 <- temp1 %>%
  mutate(alc_cur_freq_cat = case_when(
    alc_cur_freq == 0 ~ "never",
    alc_cur_freq %in% c(1,2,3,4) ~ "1 time a week or less",
    alc_cur_freq == 5 ~ "2-3 times a week",
    alc_cur_freq %in% c(6,7) ~ "4 times a week or more",
    TRUE ~ NA_character_
  ))
temp1 %>% count(alc_cur_freq_cat)

summary(temp1$alc_cur_freq)

# binge drinking
temp1 <- temp1 %>%
  mutate(alc_binge_cat = case_when(
    alc_binge_freq_female == 0 ~ "never",
    alc_binge_freq_female %in% c(1,2,3,4,5) ~ "1 or less a week",
    alc_binge_freq_female == 6 ~ "2-3 times a week",
    alc_binge_freq_female == 7 ~ "4-5 times a week",
    alc_binge_freq_female == 8 ~ "6-7 times a week",
    TRUE ~ NA_character_
  ))


# if alc_ever=0, set alcohol frequency and binging frequency to 0
# temp1 <- temp1 %>%
#   mutate(across(
#     c(alc_cur_freq, alc_binge_freq_female, ends_with("week")),
#     function(x) ifelse(alc_ever==0, 0, x)))


temp1 <- temp1 %>%
  mutate(alc_cur_freq = ifelse(alc_ever==0, 0, alc_cur_freq)) %>%
  mutate(across(c(alc_binge_freq_female, ends_with("week")),
                function(x) ifelse(alc_cur_freq==0, 0, x)))

temp1 %>% count(alc_cur_freq, alc_binge_freq_female)

# calculate number of drinks per week
temp1 <- temp1 %>%
  mutate(alc_drink_per_week = rowSums(across(ends_with("week")), na.rm=TRUE))





# ethnicity
temp1 <- temp1 %>%
  mutate(ethnicity = coalesce(sdc_eb_white, !!!select(., starts_with("sdc_eb_"))))
temp1 %>% count(ethnicity)


# education levels
temp1 <- temp1 %>%
  mutate(edu_level = case_when(
    sdc_edu_level %in% 1:2 ~ "high_school_or_less",
    sdc_edu_level %in% 3:5 ~ "trade_community_certificate",
    sdc_edu_level %in% 6:7 ~ "uni_or_more",
    TRUE ~ NA_character_
  ))
temp1 %>% count(edu_level)
temp1 %>% count(sdc_edu_level)


# income level
temp1 <- temp1 %>%
  mutate(income_level = case_when(
    sdc_income %in% 1:3 ~ "50k_or_less",
    sdc_income %in% 4:5 ~ "50k_to_100k",
    sdc_income %in% 6:8 ~ "100k_or_more",
    TRUE ~ NA_character_
  ))






### NEED LITERATURE ###
# smoking
# temp1 <- temp1 %>%
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



# set live births to 0 if gravidity=0
temp1 %>% count(wh_live_births)
temp1 <- temp1 %>%
  mutate(wh_live_births = ifelse(wh_gravidity==0, 0, wh_live_births))


# set contraceptive duration to 0 if never used contraceptives
temp1 %>% count(wh_contraceptives_ever)
summary(temp1$wh_contraceptives_duration)
temp1 <- temp1 %>%
  mutate(wh_contraceptives_duration=ifelse(wh_contraceptives_ever==0, 0, wh_contraceptives_duration))

# set hrt age to 0 if never used HRT
temp1 %>% count(wh_hrt_ever)
with(temp1, by(wh_hrt_age, wh_hrt_ever, summary))



# impute menopause status for ATP and combine with BCGP
temp1 <- temp1 %>%
  mutate(ms_bc=ifelse(str_detect(menopause_status, "Post"), 1, 0),
         ms_atp=ifelse(!is.na(wh_menopause_ever), wh_menopause_ever,
                       ifelse(sdc_age_calc>=51, 1, 0))
  ) %>%
  mutate(menopause_stt = coalesce(ms_bc, ms_atp))






# Cancer descriptions

# Combine topology variables
# BCGP = site, site_desc; ATP = acr_icd_o_topography
temp1 %>% count(site)
temp1 %>% count(site_desc)
temp1 %>% count(acr_icd_o_topography)

# Combine site codes for BCGP and ATP
# Extract the code part of acr_icd_o_topography
# remove "." and convert "C" to uppercase
temp1 <- temp1 %>%
  mutate(
    site_code=ifelse(
      !is.na(site), site,
      str_to_upper(
        str_replace(
          str_sub(acr_icd_o_topography, start=1L, end=5L),
          "\\.", "")
      )
    ))


# Combine morphology
# BCGP = hist1, hist1_desc; ATP = acr_icd_o_morphology
temp1 %>% filter(cohort=="bcgp") %>% count(gp, hist1, hist1_desc)
temp1 %>% filter(cohort=="atp") %>% count(gp, acr_icd_o_morphology)


#
temp1 <- temp1 %>%
  mutate(
    hist_subtype = ifelse(
      !is.na(hist1), hist1,
      str_replace(
        str_sub(acr_icd_o_morphology, start=1L, end=6L),
        "/", ""
      )
    )
  )

table(temp1$hist_subtype)


## remove all in situ carcinomas (morphology ends with 2) and their controls
temp1 %>% count(str_detect(hist_subtype, "2$")) #3 cases
temp1 %>% filter(str_detect(hist_subtype, "2$")) %>% select(hist_subtype)

pair_in_situ <- temp1$match_id[str_which(temp1$hist_subtype, "2$")]

# dim(temp1)
# temp1 <- temp1 %>% filter(!match_id %in% pair_in_situ)
# dim(temp1)


# assign histological subtype
temp1 <- temp1 %>%
  mutate(
    hist_subtype=case_when(
      str_detect(hist_subtype, "^850") ~ "ductal",
      str_detect(hist_subtype, "85203") ~ "lobular",
      str_detect(hist_subtype, "8522|8523|8524|8255") ~ "mixed",
      is.na(hist_subtype) ~ NA_character_,
      TRUE ~ "other"
    )
  )


temp1 %>% filter(gp=="CASE") %>% count(hist_subtype) %>% mutate(p=proportions(n))




###############################################

full_dat <- temp1 %>%
  select(gp, studyid, match_id, cohort, dup, baseline_info_avail,
         sdc_age_calc, #collect_year, followuptime,
         menopause_stt,
         ethnicity,
         edu_level, sdc_edu_level,
         sdc_income, income_level,
         wh_menstruation_age,
         wh_contraceptives_ever, # wh_contraceptives_age, wh_contraceptives_duration,
         wh_gravidity, wh_live_births, # wh_preg_first_age,
         # wh_breastfeeding_duration,
         wh_hft_ever, wh_hrt_ever, # wh_hrt_age, wh_hrt_duration,
         fam_hist_breast,
         nut_veg_qty, nut_fruits_qty,
         alc_ever, alc_cur_freq, alc_cur_freq_cat,
         alc_binge_freq_female, alc_binge_cat,
         alc_drink_per_week,
         smk_cig_status, # smk_cig_dur, smk_cig_qty, # need to find more smoking variables
         bmi, bmi_cat,
         age_at_diagnosis,
         er_status, pr_status, her2_status, hr_subtype,
         site_code, hist_subtype)


full_dat <- full_dat %>%
  mutate(across(c(gp, cohort, menopause_stt, dup,
                  ethnicity, edu_level, sdc_edu_level,
                  sdc_income, income_level,
                  wh_contraceptives_ever, wh_hft_ever, wh_hrt_ever,
                  fam_hist_breast,
                  alc_ever, alc_cur_freq, alc_cur_freq_cat,
                  alc_binge_freq_female, alc_binge_cat,
                  smk_cig_status, bmi_cat
  ),
  as_factor))




write.csv(full_dat, file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")

