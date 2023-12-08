####################################################
###  CALCULATE SUMMARY STATISTICS FOR DATA
####################################################
library(tidyverse)


source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


# load full questionnaire data
full_dat <- import_survey()







full_dat %>% count(gp)

full_dat %>% group_by(gp) %>% count(cohort) %>% mutate(p=n/sum(n))




#=================================================================#
# summary for numerical variables
# mean, SD, and count of missing
full_num_summ <- full_dat %>%
  select(gp, all_of(num_cols)) %>%
  pivot_longer(-gp, names_to="var") %>%
  group_by(gp, var) %>%
  summarize(mean=round(mean(value, na.rm=TRUE), 1),
            sd=round(sd(value, na.rm=TRUE), 1),
            nmiss=round(sum(is.na(value)), 1)) %>%
  arrange(var, gp)
View(full_num_summ)




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
