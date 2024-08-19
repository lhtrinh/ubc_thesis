#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION ####
#=============================================================#
library(tidyverse)
library(foreach)
library(doParallel)
library(rio)
library(mice) # for pooling results
# library(survival) # for conditional logistic regression




source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/exposome/00_functions.R")


full_dat <- import_survey()
full_meta <- import_meta_sc()
full_imp <- import_imp()

# combine questionnaire and metabolomics data
full_all <- full_dat %>% full_join(full_meta)


# change outcome variable to numeric 0/1 with 0=control, 1=case
full_all %>% count(gp)


#=======================================================================#
## ADJUST FOR CONTEXTUAL VARIABLES ####

ulr_results <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/exposome/overall_lr_results.csv")

sig_ions_df <- ulr_results %>% filter(qval<=0.1)

ion_cols <- sig_ions_df$ion
ion_cols

names(ion_cols) <- ion_cols
ion_cols <- as.list(ion_cols)

n_imp <- 10




## Function for adjustment ####
confounder_ulr <- function(ion) {
  dfs <- foreach(var=context_cols, .combine="rbind", .verbose=TRUE) %:%
    foreach(i=1:n_imp, .combine="rbind", .verbose=TRUE) %do% {
      dat <- full_imp[[i]]

      adj_lr_formula <- as.formula(paste0(
        "gp ~ ",
        ion,
        " + menopause_stt + collect_year + collect_age + cohort + ",
        var))

      lr_fit <- glm(adj_lr_formula, data=dat, family=binomial(link="logit"))
      lr_coef <- summary(lr_fit)$coefficients

      coef <- lr_coef[rownames(lr_coef)==ion, 1]
      pval <- lr_coef[rownames(lr_coef)==ion, 4]

      data.frame(imp=i, ion=ion, context_var=var, beta_adj=coef, pval_adj=pval)
    }
}


confounder_ulr <- lapply(ion_cols, confounder_ulr) %>% bind_rows
confounder_ulr %>% count(imp, ion, context_var) %>% count(n)
head(confounder_ulr)

confounder_df <- confounder_ulr %>%
  mutate(or_adj=exp(beta_adj)) %>%
  select(ion, context_var, beta_adj, or_adj) %>%
  full_join(sig_ions_df) %>%
  mutate(beta_unadj = beta,
         or_unadj = or) %>%
  mutate(beta_change=1-beta_adj/beta_unadj,
         or_change = 1-or_adj/or_unadj)

head(confounder_df)


summary(confounder_df$beta_change)
summary(confounder_df$or_change)


# how many are significant confounders?
confounder_df %>% count(abs(beta_change)>=0.1)
confounder_df %>% count(abs(or_change)>=0.1)

# no confounding detected

write_csv(confounder_df,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/exposome/overall_confounder_results.csv")
