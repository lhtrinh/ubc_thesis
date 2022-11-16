source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup.R")
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Cleaning.R")



#=========================================================#
# Summary baseline statistics BEFORE IMPUTATION
#=========================================================#

dat_noimp %>%
  select(gp,
         sdc_age_calc,
         bmi,
         followuptime,
         age_at_diagnosis,
         wh_menopause_age,
         wh_menstruation_age,
         wh_hrt_age,
         wh_contraceptives_age,
         ) %>%
  group_by(gp) %>%
  summarize(across(
    everything(),
    .fns=c(avg=function(x) mean(x,na.rm=T), stdev=function(x) sd(x,na.rm=T))
  )) %>%
  pivot_longer(
    -gp,
    names_to="var_stat",
    values_to="n") %>%
  mutate(stat=gsub("^.*_", "", var_stat),
         var=gsub("_[^_]*$", "", var_stat)) %>%
  select(gp, var, stat, n)

write.csv(num_summ, file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/numeric_summary.csv")


dat_noimp %>%
  select(gp, ethnicity, menopause_status, smk_cig_status,
         wh_hrt_ever,
         wh_contraceptives_ever) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(
    -gp,
    names_to="var",
    values_to="level"
  ) %>%
  count(gp, var, level) %>% View()
