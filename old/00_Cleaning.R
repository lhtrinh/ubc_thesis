source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup.R")


###### CLEANING ###########


#=========================================================#
# Remove extra IDs
#=========================================================#

# Remove one ID that was selected in error
temp1 <- bc_dat %>%
  filter(studyid!="PB000429")

# Remove duplicate IDs
temp2 <- temp1[!duplicated(temp1$studyid),]

#=========================================================#
# Format SPSS date formats
#=========================================================#
temp1 <- temp2 %>%
  mutate(across(.cols=c(adm_qx_completion, f1_adm_qx_completion),
                .fns = function(x) as.Date(x/86400, origin = "1582-10-14")))

head(temp1)

#=========================================================#
# Replace all -7 with NAs
#=========================================================#
temp1[temp1==-7] <- NA


#=========================================================#
# Create new variables from data
#=========================================================#
# follow-up time

# age
summary(temp1$sdc_age_calc)
# no edit

# missing: fam hist of cancer

# ethnicity: white/non-white
# BMI: calculated first, then fill in the blank with self-reported
dat_noimp <- temp1 %>%
  mutate(
    ethnicity=coalesce(sdc_eb_white, !!! select(., matches("eb"))),
    bmi=coalesce(pm_tanitabmi, pm_bioimped_bmi, pm_bmi_sr)
  )
# smoking status: current daily, current occasional, past, never
# use smk_cig_status
temp1 %>%
  select(contains("smk")) %>%
  head()

# age at diagnosis
















###############################################################################




#=========================================================#
# REMOVE FOLLOW-UP COLUMNS
# Some variables that were missing from the original df
# were added using follow-up ("f1_").
# We need to combine these variables
#=========================================================#

f1_cols <- temp1 %>% select(starts_with("f1_")) %>% names()
length(f1_cols)

# 42 columns from follow-up questionnaire

og_cols <- str_replace_all(f1_cols, "f1_", "")
dat_cols <- temp1[,names(temp1) %in% og_cols] %>% names()
length(dat_cols)
# 27 column names are duplicated
# 15 columns from follow-up that do not have matching names in original

# F1_ variables that have duplicates by the same name: 27 variables
# F1_ variables that do not have the same name in the original questionnaire:
### wh_contraceptives_duration <--> f1_wh_contraceptive_dur_mo
### wh_hrt_duration <--> f1_wh_hrt_dur_mo
### pm_weight_sr_avg <--> f1_pm_weight_kg

# F1_ variables that do not have duplicates: 12 variables
### f1_sdc_marital_status
### f1_wh_contraceptive_dur_yr
### f1_wh_hrt_dur_yr
### f1_wh_hrt_type
### f1_dis_cancer_sib_ever
### f1_dis_cancer_sib_num
### f1_alc_red_wine
### f1_alc_white_wine
### f1_alc_beer
### f1_alc_liquor
### f1_alc_other_alc
### f1_pm_weight_lb
### f1_adm_qx_completion and adm_qx_completion
# Keep as is

# Final list of duplicates
dups1 <- data.frame(
  og = og_cols[which(og_cols %in% dat_cols)],
  dup = f1_cols[which(og_cols %in% dat_cols)])

dups2 <- data.frame(matrix(
  c(
    c("wh_contraceptives_duration", "f1_wh_contraceptive_dur_mo"),
    c("wh_hrt_duration", "f1_wh_hrt_dur_mo"),
    c("pm_weight_sr_avg", "f1_pm_weight_kg")
  ), ncol=2, byrow=TRUE)
)

dupnames_df <- rbind(dups1, setNames(dups2, names(dups1)))
dupnames_df

# Remove variables with duplicate names but are actually different measures

# Leave alone:
### adm_qx_completion && f1_adm_qx_completion
### sdc_age_calc && f1_sdc_age_calc
### wh_preg_cur and f1_wh_preg_cur
### wh_preg_cur_wk and f1_wh_cur_wk
###

# Coalesce:
### sdc_sex and f1_sdc_sex
### pm_weight_sr_avg and f1_pm_weight_kg OR maybe just need BMI


# Others:
### wh_contraceptives_ever and f1_wh_contraceptives_ever:
### if og=1 then 1
### if og=0 then 0
### if og is NA and f1=1 then 1
### if og is NA and f1 is NA then ??????!!!!
### sdc_income and f1_sdc_income
### wh_contraceptives_age and f1_wh_contraceptives_age: why???
### wh_gravidity and f1_wh_gravidity
### wh_live_births and f1_wh_live_births
### if og=-7 then 0
### if og!=dup then og
### wh_preg_last_age --> could maybe do a binary of late (>30) vs early pregnancies?
### wh_menopause_ever and f1_menopause_ever
### delete both, already synthesized with menopause_status
### wh_menopause_age and f1_wh_menopause_age
### do we need menopause age?
### maybe delete, but double check to make sure menopause_status is accurate
### wh_hrt_ever and f1_wh_hrt_ever
### if og=0 then 0
### if og=1 then 1
### if og is NA and f1=0 then 0
### if og is NA and f1=1 then ???
### all variables with "dis_" (family hist of cancer)
### combine to make Y/N family hist of cancer
### but use original or follow-up?
### alc_ever and f1_alc_ever, smk_cig_ever and f1_smk_cig_ever
### if og=0 then 0
### if og=1 then 1
### if og=NA and f1=0 then 0
### if og=NA and f1=1 then ???
### if follow-up in the same year as sample collection then yes???
### Contraceptives variables: just binary Y/N?
### BMI variables: pm_bmi_sr, pm_bioimped_bmi, pm_tanitabmi
### coalesce --> still have 95 NAs???????

temp1 %>%
  mutate(wt = coalesce(pm_weight_sr_avg,
                       f1_pm_weight_kg),
         ht = coalesce(pm_stand_height_sr_avg,
                       pm_standing_height_avg)) %>%
  mutate(bmi = wt/((ht/100)^2)) %>%
  mutate(bmi_calc = coalesce(pm_bioimped_bmi,
                             pm_tanitabmi,
                             pm_bmi_sr)) %>%
  count(is.na(bmi_calc))

temp1 %>% select(contains("bmi")) %>% head()


#=========================================================#
# Find index date #
#=========================================================#

# Some potential index dates:
# Date at which questionnaire is completed
# Date at which sample is collected
summary(temp1)
temp1 %>% count(gp)
temp1$sdc_age_calc

cases <- temp1 %>% filter(gp=='CASE')
cntls <- temp1 %>% filter(gp=='CNTL')
head(cases)

summary(cases$sdc_age_calc)
summary(cases$f1_sdc_age_calc)
temp1 %>%
mutate(age=coalesce(sdc_age_calc, f1_sdc_age_calc)) %>%
count(is.na(age))


temp1 %>%
  filter(is.na(sdc_age_calc)) %>%
  select(studyid, collect_year, followuptime, f1_adm_qx_completion, f1_sdc_age_calc)


temp1 %>% select(adm_qx_completion, f1_adm_qx_completion, followuptime) %>%
  mutate(og_yr = year(adm_qx_completion),
         f1_yr = year(f1_adm_qx_completion)) %>%
  count((f1_yr-og_yr)-followuptime==0)


# Index year = year of baseline questionnaire completion
# else use age at follow-up minus follow-up time
# else use



# Age at dx (for cases)
summary(cases$age_at_diagnosis)







