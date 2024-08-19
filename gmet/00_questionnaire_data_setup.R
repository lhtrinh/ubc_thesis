library(tidyverse)
library(haven)
library(readxl)
library(lubridate)
library(psych)
library(moments)
library(readr)

#=========================================================#




#====================#
# BCGP data ####
#====================#

# load data sets

# this code loads BCGP survey data
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/gmet/00_setup_add_var_bc_dat.R")
colnames(bc_raw) <- tolower(colnames(bc_raw))
dim(bc_raw)

# hormone receptor status data
bc_hormones <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_Y2022_BCCR_incidence_breast_hormones.sav")
colnames(bc_hormones) <- tolower(colnames(bc_hormones))

# matching information data
bc_match_ind <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP/BCGP2641_MatchIndicator.csv")
colnames(bc_match_ind) <- tolower(colnames(bc_match_ind))





# shorten variable names for HR status
bc_hormones <- bc_hormones %>%
  mutate(er_status=tolower(bccr_er_status_final),
         pr_status=tolower(bccr_pr_status_final),
         her2_status=tolower(bccr_her2_status_final))



# add "BCGP" to match ID to distinguish from ATP match IDs
bc_match_ind$match_id <- paste("BCGP", bc_match_ind$match_id, sep="_")
# remove duplicates
bc_match_ind <- bc_match_ind[!duplicated(bc_match_ind),]




# remove one set of duplicates in bc_raw and combine BCGP data
bc_raw2 <- bc_raw %>%
  inner_join(bc_match_ind) %>%
  left_join(bc_hormones)


bc_raw2$cohort <- "BCGP"
dim(bc_raw2)




# if a participant has no baseline questionnaire info,
# mark them and their matched pairs as "no baseline avail"
no_baseline <- bc_raw2 %>%
  filter(is.na(adm_qx_completion)) %>%
  select(studyid, match_id)


bc_raw3 <- bc_raw2 %>%
  filter(!match_id %in% no_baseline$match_id)



# Calculate age at blood collection
## by calculating difference between year of survey and year of blood collection
bc_raw4 <- bc_raw3 %>%
  mutate(collect_age = (sdc_age_calc+(collect_year-year(adm_qx_completion))))



# !!! one patient had the wrong age
# participant PB000321's age is 52, but written as 62
bc_raw5 <- bc_raw4 %>%
  mutate(sdc_age_calc=ifelse(studyid=="PB000321", 52, sdc_age_calc))
bc_raw5$sdc_age_calc[bc_raw5$studyid=="PB000321"]


bc_dat <- bc_raw5


dim(bc_dat)
n_distinct(bc_dat$match_id)





#=========================================================#






#====================#
# ATP data ####
#====================#

# survey data
atp_raw <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_HFRQ_BL_CPAC.csv")
colnames(atp_raw) <- tolower(colnames(atp_raw))

# physical measurements
atp_pm_measure <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_PM_BL_CPAC.csv")
colnames(atp_pm_measure) <- tolower(colnames(atp_pm_measure))

# hormone receptor status
atp_hormones <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_ACR_Breast_DX_ER_PR.csv")
colnames(atp_hormones) <- tolower(colnames(atp_hormones))



# unharmonized data, contains case information
atp_cases <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_ACR.csv")
colnames(atp_cases) <- tolower(colnames(atp_cases))


# matching ID data
atp_matched_pairs <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/2105_Study2105_03_29_2023_CaseControl_ParticipantKeyMatchedList.xlsx",
                               sheet="Sheet1")


# biospecimens data
atp_bio <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/vw2105_BIO.csv")
colnames(atp_bio) <- tolower(colnames(atp_bio))




#====================#
# add match IDs to matched pairs
atp_matched_pairs %>% head()
atp_match_id <- atp_matched_pairs %>%
  rename(CNTL=CONTROL) %>%
  mutate(match_id=1:nrow(atp_matched_pairs)) %>%
  pivot_longer(!match_id,
               names_to="gp",
               values_to="participantkey") %>%
  mutate(match_id=paste("ATP", match_id, sep="_"))
atp_match_id %>% head()





# atp_hormones: change variable names for HR/HER2 status and change all values to lowercase
atp_hormones <- atp_hormones %>%
  filter(acr_cancer_site_aggregate=="Breast") %>%
  arrange(participantkey, acr_malignancy_number)


atp_hormones <- atp_hormones %>%
  mutate(er_status=acr_estrogen_receptor_assay,
         pr_status=acr_progesterone_receptor_assay,
         her2_status=acr_her2_assay) %>%
  mutate(across(ends_with("_status"),
                function(x) ifelse(grepl("Positive", x), "positive",
                                   ifelse(grepl("Negative", x), "negative", NA))))

# remove duplicates
atp_hormones_new <- atp_hormones[!duplicated(atp_hormones$participantkey),]




# atp_cases: only select the first breast cancer diagnosis
atp_cases <- atp_cases %>%
  filter(acr_cancer_site_aggregate=="Breast") %>%
  arrange(participantkey, acr_diagnosis_year) %>%
  rename(age_at_diagnosis = acr_age_at_diagnosis)

atp_cases_new <- atp_cases[!duplicated(atp_cases$participantkey),]


# atp_bio: remove duplicates with NA in age at blood collection
# and change name to collect_age
atp_bio %>% count(participantkey) %>% filter(n>1) %>% left_join(atp_bio)
## 15 cases had duplicates, seems like age at blood collection is NA for each dup
atp_bio %>% count(is.na(bio_age_at_cptp_blood))
# only the duplicates had NAs for age, remove

atp_bio_new <- atp_bio %>%
  filter(!is.na(bio_age_at_cptp_blood)) %>%
  mutate(collect_age = bio_age_at_cptp_blood)


# combine ATP data
atp_raw2 <- atp_raw %>%
  inner_join(atp_pm_measure) %>%
  left_join(atp_hormones_new) %>%
  left_join(atp_cases_new) %>%
  left_join(atp_match_id) %>%
  left_join(atp_bio_new)




# change date form for questionnaire completion date
# also add year of blood collection
atp_raw2$adm_qx_completion <- as.Date(atp_raw2$adm_qx_completion)


# calculate year at blood collection
atp_raw3 <- atp_raw2 %>%
  mutate(adm_qx_completion = as.Date(atp_raw2$adm_qx_completion),
         collect_year = round(year(adm_qx_completion) + (collect_age - sdc_age_calc), 0))
summary(atp_raw3$collect_year)


#====================#
# Replace all -7 with NAs
#====================#
# atp_raw3[atp_raw3==-7] <- NA



# Rename ATP identifier variable (participantkey) to studyid

atp_raw4 <- atp_raw3 %>%
  mutate(cohort="ATP",
         studyid=as.character(participantkey))


# remove work variables
atp_raw5 <- atp_raw4 %>%
  select(-starts_with("wrk"))




atp_dat <- atp_raw5





#=========================================================#
# COMBINE DATA SETS ####

dat_raw <- bc_dat %>%
  full_join(atp_dat)





#====================#
# Replace all -7 with NAs
#====================#
dat_raw[dat_raw==-7] <- NA










#=========================================================#
# CREATE NEW VARIABLES ####



#====================#
## Matching factors ####
### Categorize collect_year
temp1 <- dat_raw %>%
  mutate(collect_yr_cat = case_when(
    between(collect_year,2009,2011) ~ "yr_09_11",
    between(collect_year,2012,2013) ~ "yr_12_13",
    between(collect_year,2014,2016) ~ "yr_14_16",
    TRUE ~ NA_character_
  ))






# impute menopause status for ATP, then combine with BCGP
temp1 <- temp1 %>%
  mutate(ms_bc=ifelse(str_detect(menopause_status, "Post"), 1, 0),
         ms_atp=ifelse(!is.na(wh_menopause_ever), wh_menopause_ever,
                       ifelse(sdc_age_calc>=51, 1, 0))
  ) %>%
  mutate(menopause_stt = coalesce(ms_bc, ms_atp))








#====================#
## Demographics ####


## ethnicity
# create white/non-white categories for ethnicity
temp1 <- temp1 %>%
  mutate(ethnicity = coalesce(sdc_eb_white, !!!select(., starts_with("sdc_eb_"))))
temp1 %>% count(ethnicity)


# education levels
temp1 <- temp1 %>%
  mutate(edu_level = case_when(
    sdc_edu_level %in% 1:2 ~ "1_high_school_or_less",
    sdc_edu_level %in% 3:5 ~ "2_some_uni",
    sdc_edu_level %in% 6:7 ~ "3_bachelor_or_more",
    TRUE ~ NA_character_
  ))
temp1 %>% count(edu_level)


# income level
temp1 <- temp1 %>%
  mutate(income_level = case_when(
    sdc_income %in% 1:3 ~ "1_50k_or_less",
    sdc_income %in% 4:5 ~ "2_50k_to_100k",
    sdc_income %in% 6:8 ~ "3_100k_or_more",
    TRUE ~ NA_character_
  ))
temp1 %>% count(income_level)










#====================#
## Family/Reproductive factors ####


### Family history of breast cancer
temp1 <- temp1 %>%
  mutate(
    cancer_fam_breast = ifelse(dis_cancer_m_breast==1 |
                                 dis_cancer_sib_breast==1 |
                                 dis_cancer_child_breast==1, 1, 0)) %>%
  mutate(
    fam_hist_breast = case_when(
      cancer_fam_breast==1 ~ 1,
      is.na(cancer_fam_breast) & !is.na(dis_cancer_fam_ever) ~ 0,
      TRUE ~ NA
    )
  )

temp1 %>% count(fam_hist_breast)


# add number of relatives with breast cancer (numeric and categorical count)
temp1 <- temp1 %>%
  mutate(
    fam_hist_breast_num=case_when(
      fam_hist_breast==0 ~ 0,
      fam_hist_breast==1 ~ rowSums(across(c(dis_cancer_m_breast, dis_cancer_sib_breast_num, dis_cancer_child_breast_num)), na.rm=TRUE),
      TRUE ~ NA_integer_
    )) %>%
  mutate(fam_hist_breast_cat=case_when(
    fam_hist_breast_num==0 ~ "0",
    fam_hist_breast_num==1 ~ "1",
    fam_hist_breast_num>=2 ~ "2_or_more",
    TRUE ~ NA_character_
  ))

temp1 %>% count(fam_hist_breast_num)
temp1 %>% count(gp, fam_hist_breast_cat) %>% mutate(p=n/sum(n))







# menstruation categories
temp1 <- temp1 %>%
  mutate(wh_menstruation_age_cat=case_when(
    wh_menstruation_age<12 ~ "cat1_under_12",
    between(wh_menstruation_age, 12, 13) ~ "cat2_12_13",
    wh_menstruation_age>=14 ~ "cat3_14_or_older"
  ))



# set live births to 0 if gravidity=0
temp1 %>% count(wh_live_births)
temp1 <- temp1 %>%
  mutate(wh_live_births = ifelse(wh_gravidity==0, 0, wh_live_births))


# make categories for breastfeeding
## summarize breastfeeding by number of live births
with(temp1, by(wh_breastfeeding_duration, wh_live_births, summary))


# calculate breastfeeding duration per live birth
temp1 <- temp1 %>%
  mutate(wh_breastfeeding_per_birth = ifelse(
    wh_live_births==0, NA,
    wh_breastfeeding_duration/wh_live_births)) %>%
  mutate(wh_breastfeeding_cat = case_when(
    wh_live_births==0 ~ "cat1_no_live_births",
    between(wh_breastfeeding_per_birth, 0, 5) ~ "cat2_0_5_months",
    between(wh_breastfeeding_per_birth, 6, 11) ~ "cat3_6_11_months",
    wh_breastfeeding_per_birth > 12 ~ "cat4_12_months_or_more",
    TRUE ~ NA_character_
  ))
temp1 %>% group_by(gp) %>%
  count(wh_breastfeeding_cat) %>%
  filter(!is.na(wh_breastfeeding_cat)) %>%
  mutate(p=proportions(n))


# categorize age at first pregnancy
summary(temp1$wh_preg_first_age)

temp1 <- temp1 %>%
  mutate(wh_preg_first_age_cat = case_when(
    wh_gravidity==0 ~ "cat1_no_preg",
    wh_preg_first_age < 20 ~ "cat2_under_20",
    between(wh_preg_first_age, 20, 24) ~ "cat3_20_24",
    between(wh_preg_first_age, 25, 29) ~ "cat4_25_29",
    wh_preg_first_age >= 30 ~ "cat5_over_30",
    TRUE ~ NA_character_
  ))
temp1 %>% group_by(gp) %>% count(wh_preg_first_age_cat) %>%
  filter(!is.na(wh_preg_first_age_cat)) %>%
  mutate(p=proportions(n))




# contraceptives
temp1 %>% count(wh_contraceptives_ever)
# summary(temp1$wh_contraceptives_duration)




# set hrt duration if never used HRT
temp1 %>% count(gp, wh_hrt_ever)
with(temp1, by(wh_hrt_age, wh_hrt_ever, summary))

temp1 <- temp1 %>%
  mutate(wh_hrt_duration=ifelse(wh_hrt_ever==0,0, wh_hrt_duration)) %>%
  mutate(wh_hrt_duration_yr=wh_hrt_duration/12)
with(temp1, by(wh_hrt_duration, wh_hrt_ever, summary))



# categorize HRT duration
temp1 <- temp1 %>%
  mutate(wh_hrt_dur_cat=case_when(
    wh_hrt_ever == 0 ~ "cat1_no_hrt",
    wh_hrt_duration < 12 ~ "cat2_less_than_1_year",
    between(wh_hrt_duration, 12, 59) ~ "cat3_1_5_years",
    between(wh_hrt_duration, 60, 119) ~ "cat4_5_10_years",
    wh_hrt_duration >= 120 ~ "cat5_10_years_or_more",
    TRUE ~ NA_character_
  ))
temp1 %>% group_by(gp) %>%
  count(wh_hrt_dur_cat) %>%
  filter(!is.na(wh_hrt_dur_cat)) %>%
  mutate(p=proportions(n))


# categorize age at menopause
temp1 <- temp1 %>%
  mutate(wh_menopause_age_cat = case_when(
    menopause_stt==0 ~ "cat0_pre",
    wh_menopause_age<45 ~ "cat1_under_45",
    between(wh_menopause_age, 45, 54) ~ "cat2_45_54",
    wh_menopause_age >= 55 ~ "cat3_55_or_over",
    TRUE ~ NA_character_
  ))










#====================#
## Lifestyle factors ####

### BMI and BMI categories
temp1 <- temp1 %>%
  mutate(bmi=as.numeric(coalesce(pm_tanitabmi, pm_bioimped_bmi, pm_bmi_sr))) %>%
  mutate(bmi_cat=factor(findInterval(bmi, c(18.5, 25, 30)),
                        labels=c("underweight", "normal", "overweight", "obese")))
temp1 %>% count(bmi_cat)
summary(temp1$bmi)





### alcohol consumption

# if alc_ever=0, set alcohol frequency and binge frequency to 0
temp1 <- temp1 %>%
  mutate(alc_cur_freq = ifelse(alc_ever==0, 0, alc_cur_freq)) %>%
  mutate(across(c(alc_binge_freq_female, ends_with("week")),
                function(x) ifelse(alc_cur_freq==0, 0, x)))


# create new categories for alcohol frequencies
temp1 <- temp1 %>%
  mutate(alc_cur_freq_cat = case_when(
    alc_cur_freq == 0 ~ "cat0_never",
    between(alc_cur_freq, 1, 3) ~ "cat1_1_3_per_month",
    between(alc_cur_freq, 4, 5) ~ "cat2_1_3_per_week",
    between(alc_cur_freq, 6, 7) ~ "cat3_4_or_more_per_week",
    TRUE ~ NA_character_
  ))
temp1 %>% count(gp, alc_cur_freq_cat) %>% arrange(alc_cur_freq_cat, gp) %>% mutate(p=n/sum(n))




# create new categories for binge drinking
temp1 <- temp1 %>%
  mutate(alc_binge_cat = case_when(
    alc_binge_freq_female == 0 | alc_ever==0 ~ "cat0_never",
    between(alc_binge_freq_female, 1, 2) ~ "cat1_1_11_per_year",
    between(alc_binge_freq_female, 3, 4) ~ "cat2_1_3_per_month",
    between(alc_binge_freq_female, 5, 8) ~ "cat3_1_or_more_per_week",
    TRUE ~ NA_character_
  ))



# # calculate number of drinks per week
# temp1 <- temp1 %>%
#   mutate(alc_drink_per_week = rowSums(across(ends_with("week")), na.rm=TRUE))





# smoking

## add smoking duration
temp1 <- temp1 %>%
  mutate(
    smk_cig_dur = as.numeric(coalesce(
        smk_cig_daily_cur_dur,
        smk_cig_former_daily_dur,
        smk_cig_heaviest_dur))) %>%
  mutate(smk_cig_dur = ifelse(smk_cig_status==0, 0, smk_cig_dur)) %>%
  mutate(smk_cig_dur_cat = case_when(
    smk_cig_status==0 ~ "cat0_never_smk",
    smk_cig_dur < 10 ~ "cat1_less_than_10_years",
    between(smk_cig_dur, 10, 19) ~ "cat2_10_19_years",
    smk_cig_dur>=20 ~ "cat5_20_years_or_more",
    TRUE ~ NA_character_
  ))


# age at smoking initiation
temp1 <- temp1 %>%
  mutate(smk_first_age = coalesce(smk_cig_daily_cur_onset, smk_cig_former_daily_onset)) %>%
  mutate(smk_first_age_cat = case_when(
    smk_cig_status==0 ~ "cat0_never",
    smk_first_age < 15 ~ "cat1_younger_than_15",
    between(smk_first_age, 15, 19) ~ "cat2_15_to_19",
    between(smk_first_age, 20, 24) ~ "cat3_20_to_24",
    smk_first_age>=25 ~ "cat4_25_or_older",
    TRUE ~ NA_character_
  )
  )

# # age at smoking relative to menarche and pregnancy
# temp1 <- temp1 %>%
#   mutate(smk_before_first_preg = case_when(
#     smk_cig_status==0 ~ "cat0_never",
#     wh_gravidity==0 ~ "never_preg",
#     smk_first_age < wh_preg_first_age ~ "before",
#     smk_first_age >= wh_preg_first_age ~ "after",
#     TRUE ~ NA_character_
#   ))













#====================#
## Cancer-specific variables ####


# time from enrolment to diagnosis
summary(temp1$age_at_diagnosis - temp1$sdc_age_calc)
temp1 %>% filter(age_at_diagnosis<sdc_age_calc) %>% select(match_id)
#2 cases were diagnosed before enrolment?!


# molecular subtype
temp1 <- temp1 %>%
  mutate(hr_subtype=case_when(
    er_status=="negative" & pr_status=="negative" & her2_status=="negative" ~ "triple_negative",
    her2_status=="positive" ~ "her2_positive",
    (er_status=="positive" | pr_status=="positive") ~ "hr_positive",
    gp=='CNTL' ~ "control",
    TRUE ~ NA_character_
  ))
temp1 %>% count(hr_subtype)



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
temp1 %>% filter(cohort=="BCGP") %>% count(gp, hist1, hist1_desc)
temp1 %>% filter(cohort=="ATP") %>% count(gp, acr_icd_o_morphology)



# remove "/" from ATP morphology codes
temp1 <- temp1 %>%
  mutate(
    hist_code = ifelse(
      !is.na(hist1), hist1,
      str_replace(
        str_sub(acr_icd_o_morphology, start=1L, end=6L),
        "/", ""
      )
    )
  )

table(temp1$hist_code)


## find all morphologies that need to be removed
# in situ carcinomas (morphology ends with 2) and their controls
temp1 %>% count(str_detect(hist_code, "2$")) #3 cases
temp1 %>% filter(str_detect(hist_code, "2$")) %>% select(hist_code)



# assign histological subtype
temp1 <- temp1 %>%
  mutate(
    hist_subtype=case_when(
      str_detect(hist_code, "^850") ~ "ductal",
      str_detect(hist_code, "85203") ~ "lobular",
      str_detect(hist_code, "8522|8523|8524|8255") ~ "mixed",
      is.na(hist_code) ~ NA_character_,
      TRUE ~ "other"
    )
  )

temp1 %>% count(hist_subtype)








#=========================================================#

# REMOVE INELIGIBLE PARTICIPANTS ####

# Participant PB000429 was selected in error
# 8 pairs had at least one participant with no baseline information
# Pair of 210503833 and 210524604 bc one was selected as a control but is now a case
# 3 pairs have cases with in-situ tumors instead of invasive
# Pair of participant 210503316 because of mismatched information
# 1 participant had incorrect enrolment age


# find matching IDs of pairs to remove
rm_pair <- temp1 %>%
  filter(studyid %in% c("210503833", "210524604", "210503316") |
           str_detect(hist_code, "2$")) %>%
  select(match_id) %>%
  inner_join(subset(temp1, select=c("match_id", "studyid", "hist_code"))) %>%
  pull(studyid)

rm_id <- c(rm_pair, "PB000429")



# filter out pairs

temp2 <- temp1 %>%
  filter(!(studyid %in% rm_id))






dim(temp2)



#=========================================================#
# FINAL DATA SET ####

full_dat <- temp2 %>%
  select(gp, studyid, match_id, cohort, dup,
         sdc_age_calc, collect_year, collect_yr_cat, collect_age,
         menopause_stt,
         ethnicity,
         edu_level,
         income_level,
         wh_menstruation_age_cat, wh_menopause_age_cat,
         wh_contraceptives_ever,
         wh_gravidity, wh_live_births, wh_preg_first_age_cat,
         wh_breastfeeding_cat,
         wh_hrt_ever, wh_hrt_duration_yr,
         fam_hist_breast,
         alc_ever, alc_cur_freq_cat, alc_binge_cat,
         smk_cig_status, smk_cig_dur, smk_first_age_cat,
         pse_childhood_duration, pse_adult_home_duration, pse_adult_wrk_duration,
         bmi, bmi_cat,
         age_at_diagnosis,
         er_status, pr_status, her2_status, hr_subtype,
         site_code, hist_code, hist_subtype)




full_dat %>% count(gp)
full_dat %>% count(cohort, gp)



write_csv(full_dat, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/gmet/data_with_missing.csv")








#======================================================================#
#======================================================================#
# CREATE SURVEY-METABOLOMICS CROSSWALK DATA ####


## load crosswalk between participant and sample IDs for ATP
atp_id_crosswalk <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/Bhatti samples pulled and duplicate IDs 2023-06-02.xlsx",
  sheet="231Cs, 231Cntr w duplc"
)
names(atp_id_crosswalk)


# 1. rename column "tube barcode of duplicate sent" to dup_barcode
# 2. change all column names to lowercase
# 3. keep only columns participantkey, aliquotbarcode, and dup_barcode
atp_id_cw2 <- atp_id_crosswalk %>%
  rename(dup_barcode = `tube barcode of duplicate sent`) %>%
  rename_with(tolower) %>%
  select(participantkey, aliquotbarcode, dup_barcode)
head(atp_id_cw2)

# one participant has a matching aliquot, but it's not in the crosswalk
# so we'll add crosswalk record to data set
atp_id_cw3 <- rbind(atp_id_cw2,
                    data.frame(participantkey="210509417",
                               aliquotbarcode="0080196574",
                               dup_barcode=NA_character_))

# combine aliquot and duplicate barcodes into one column named "barcode"
atp_id_cw4 <- atp_id_cw3 %>%
  mutate(across(c(aliquotbarcode, dup_barcode),
                function(x) as.character(as.numeric(x))),
         dup=ifelse(!is.na(dup_barcode), "YES", "NO")) %>%
  pivot_longer(-c(participantkey, dup),
               values_to="barcode",
               values_drop_na=TRUE)

head(atp_id_cw4)
summary(atp_id_cw4)


# double check for duplicates
atp_id_cw4 %>%
  count(barcode) %>%
  arrange(desc(n)) #2 duplicates

atp_id_cw4 %>%
  filter(barcode==170965728 | barcode==81869160) #same participant

# there is one duplicate barcode that appeared twice
# looks like the same ID throughout so we can remove the twice duplicated barcode
atp_id_cw5 <- atp_id_cw4[!duplicated(atp_id_cw4),]

# remove two participants with no record in questionnaire data
atp_id_cw6 <- atp_id_cw5 %>%
  filter(!participantkey %in% c("210540416", "210533315"))

# final crosswalk table
atp_id_cw <- atp_id_cw6 %>%
  mutate(studyid=participantkey, cohort="ATP") %>%
  select(studyid, cohort, dup, barcode) %>%
  filter(!(studyid %in% rm_id))


head(atp_id_cw)
dim(atp_id_cw)

#======================================================================#
# ID crosswalk for both cohorts
bcgp_id_cw <- temp2 %>%
  filter(cohort=="BCGP") %>%
  mutate(barcode=studyid) %>%
  select(studyid, cohort, dup, barcode)

bcgp_id_cw %>% count(dup)

id_cw <- rbind(atp_id_cw, bcgp_id_cw) %>%
  filter(!(studyid %in% rm_id))
head(id_cw)
dim(id_cw)
n_distinct(id_cw$studyid)




#======================================================================#


write_csv(atp_id_cw,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/atp_id_crosswalk.csv")

write_csv(id_cw,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/id_crosswalk_dup_info.csv")
