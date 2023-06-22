### Descriptive statistics ###


num_summ <- function(x){
  x <- enquo(x)
  full_dat %>%
    filter(!is.na(!!x)) %>%
    group_by(gp) %>%
    summarise(meanx = mean(!!x),
              sdx = sd(!!x))
}

cat_summ <- function(x) {
  x <- enquo(x)
  full_dat %>%
    filter(!is.na(!!x)) %>%
    group_by(gp) %>%
    count(!!x) %>%
    mutate(p=proportions(n))
}


### MATCHED VARIABLES
# This section checks matching variables between cases and controls
# to make sure they're matched with no significant differences

### Cohort count
full_dat %>% count(cohort, gp)
#227 pairs from ATP, 348 pairs from BCGP

### menopause status
full_dat %>% count(menopause_stt) %>% mutate(p=proportions(n))
# all matched




### DEMOGRAPHICS

# Age
with(full_dat, by(sdc_age_calc, gp, describe))
num_summ(sdc_age_calc)

# ethnicity
full_dat %>% count(ethnicity) %>% mutate(p=proportions(n))
cat_summ(ethnicity)

# income level
cat_summ(income_level)

# education level
cat_summ(edu_level)




### REPRODUCTIVE
