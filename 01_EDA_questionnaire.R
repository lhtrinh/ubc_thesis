####################################################
###  CALCULATE SUMMARY STATISTICS FOR DATA
####################################################
library(tidyverse)


source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load full questionnaire data
full_dat <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")

full_dat <- full_dat %>%
  mutate(across(all_of(fct_cols),
                as.factor)) %>%
  mutate(across(all_of(int_cols),
                as.integer))







full_dat %>% count(gp)

full_dat %>% group_by(gp) %>% count(cohort) %>% mutate(p=n/sum(n))




#=================================================================#
# summary for numerical variables
# mean, SD, and count of missing
full_num_summ <- full_dat %>%
  group_by(gp) %>%
  summarise(across(
    all_of(num_cols),
    c(avg=function(x) round(mean(x, na.rm=TRUE), 1),
           sd=function(x) round(sd(x, na.rm=TRUE), 1),
           nmiss=function(x) round(sum(is.na(x)), 1)))) %>%
  pivot_longer(-gp, names_to="var") %>%
  arrange(var, gp)
View(full_num_summ)

# # for parous women only
# full_dat %>% filter(wh_gravidity>0) %>%
#   group_by(gp) %>%
#   summarise(across(c(wh_preg_first_age, wh_live_births, wh_breastfeeding_duration),
#                    c(avg=function(x) round(mean(x, na.rm=TRUE), 1),
#                      sd=function(x) round(sd(x, na.rm=TRUE), 1),
#                      nmiss=function(x) round(sum(is.na(x)), 1)))) %>%
#   pivot_longer(-gp, names_to="var") %>%
#   arrange(var, gp)


#=================================================================#
# summary for categorical variables
# non-missing count and percentage, plus count of missing with no pct
full_dat_cat <- full_dat %>%
  select(all_of(fct_cols)) %>%
  pivot_longer(-gp, names_to="var") %>%
  group_by(gp, var) %>%
  count(value)

full_cat_summ <- full_dat_cat %>%
  filter(!is.na(value)) %>%
  mutate(pct=round(100*n/sum(n),1)) %>%
  ungroup() %>%
  full_join(full_dat_cat) %>%
  arrange(var, value, gp)

View(full_cat_summ)






#=================================================================#
# number of missing

full_dat %>%
  summarize(across(all_of(context_cols),
                   function(x) round(100*sum(is.na(x))/nrow(full_dat),1))) %>%
  pivot_longer(everything(), names_to="var", values_to="prop_missing") %>%
  arrange(prop_missing) %>%
  View()
