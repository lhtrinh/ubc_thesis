########################################################
#
###    MANUALLY DEFINED FUNCTIONS FOR METABOLOMCIS   ###
#
########################################################
library(tidyverse)
library(lme4)
library(performance)


#=====================#


##################
# CALCULATE CV
##################

# cv_calc: calculate coefficient of variation (CV%)
# per metabolite per participant by dividing each
# metabolite's SD by its mean and multiply by 100%
calc_cv <- function(dat, group_var){

  group_var <- enquo(group_var)

  # calculate CV per metabolite per participant
  replicate_cv <- dat %>%
    group_by(!!group_var) %>%
    summarise(across(contains("ion"),
                     function(x) sd(x)/mean(x)*100)) %>%
    ungroup()

  return(replicate_cv)
}





#=====================#
# cv_calc_med
# calculate median coefficient of variation
# across each metabolite
calc_cv_med <- function(dat, group_var){
  group_var <- enquo(group_var)

  cv_dat <- calc_cv(dat, !!group_var)


  replicate_cv_med <- cv_dat %>%
    summarise(across(contains("ion"), median)) %>%
    t()

  # record the CV into a dataset
  ion_cv_dat <- data.frame(metabolite=rownames(replicate_cv_med),
                           median_cv=replicate_cv_med[,1])
  rownames(ion_cv_dat) <- NULL
  return(ion_cv_dat)
}






#=====================#
##################
# CALCULATE ICC
##################



# Use ICC function from psych package
# Calculate ICC3 to measure test-retest reliability

calc_icc <- function(dat, group_var){
  ion_cols <- names(dat)[str_detect(names(dat), "ion")]
  icc_dat <- data.frame(
    metabolite = ion_cols,
    icc = 0)

  for (i in 1:length(ion_cols)){
    ion_col <- ion_cols[i]
    ion_col <- enquo(ion_col)
    group_var <- enquo(group_var)

    ion_dat <- dat %>% select(!!group_var, !!ion_col)
    colnames(ion_dat) <- c("id", "y")
    my_mod <- lme4::lmer(y ~ 1 + (1|id), data=ion_dat)
    my_icc <- performance::icc(my_mod)$ICC_unadjusted

    icc_dat[i,2] <- my_icc
  }
  icc_dat
}
