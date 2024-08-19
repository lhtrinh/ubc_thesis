#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION (DUCTAL) ####
#=============================================================#
library(tidyverse)
library(foreach)
library(doParallel)
library(rio)





source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


full_dat <- import_survey()
full_meta <- import_meta_sc()


full_all <- full_dat %>% full_join(full_meta)



# import imputed data (with metabolites)
full_imp <- import_imp()


#=======================================================================#
## Adjust only for matching factors ####

### Define logreg function ####
unadj_logreg <- function(dat){
  ion_cols <- colnames(dat)[grepl("^ion", colnames(dat))]
  n <- length(ion_cols)

  lr_out <-
    foreach (i=1:n, .combine="rbind", .verbose=TRUE) %do% {
      ion <- as.character(ion_cols[i])

      # create data set with only required column names
      # (vectors of column names defined in "00_functions.R")
      cols_ <- c(y_col, match_cols, ion)
      lr_dat <- dat[,colnames(dat) %in% cols_]

      # fit logistic regression model
      lr_mod <- glm(gp~., data=lr_dat, family=binomial(link="logit"))
      lr_confint <- confint(lr_mod)

      # extract p-value for metabolite
      lr_coef <- summary(lr_mod)$coefficients

      coef <- lr_coef[rownames(lr_coef)==ion, 1]
      pval <- lr_coef[rownames(lr_coef)==ion, 4]
      ci_lb <- lr_confint[rownames(lr_confint)==ion, 1]
      ci_ub <- lr_confint[rownames(lr_confint)==ion, 2]

      data.frame(metabolite=ion,
                 beta_unadj=coef,
                 ci_lb=ci_lb,
                 ci_ub=ci_ub,
                 pval=pval)
    }
}




#===================================#
### Apply on full data set ####
match_id_5yr <- full_all$match_id[full_all$age_at_diagnosis - full_all$sdc_age_calc<=5]
full_5yr <- full_all %>% filter(match_id %in% match_id_5yr)


pvals_5yr <- unadj_logreg(full_5yr)
head(pvals_5yr)


#===================================#
### Correction for multiple testing ####
pvals_5yr$q_fdr <- p.adjust(pvals_5yr$pval, method="fdr")
head(pvals_5yr)








write_csv(pvals_5yr,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/pvals_5yr.csv")



# x_ductal <- sig_ions_ductal %>%
#   mutate(across(c(beta_unadj, ci_lb, ci_ub), function(x) round(exp(x), 2))) %>%
#   mutate(or_ductal = paste(beta_unadj, " (", ci_lb, "-", ci_ub, ")", sep="")) %>%
#   select(metabolite,or_ductal)
#
# x_ductal
#
# write_csv(x_ductal, file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Output/x_ductal.csv")
#


#===================================#
# pvals_5yr <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/pvals_5yr.csv")

summary(pvals_5yr)
summary(pvals_5yr$q_fdr)

table(pvals_5yr$q_fdr<=0.2)
table(pvals_5yr$q_fdr<=0.1)



# extract all signficant metabolites (q<=0.1)
sig_ions_5yr <- pvals_5yr %>%
  filter(q_fdr<=0.1) %>%
  arrange(beta_unadj)

sig_ions_5yr


write_csv(sig_ions_5yr, file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions_5yr.csv")



# create OR tables for reporting
sig_ions_5yr %>%
  mutate(or_unadj=round(exp(beta_unadj),2),
         across(c(ci_lb, ci_ub), function(x) round(exp(x), 2)),
         across(c(pval, q_fdr), function(x) round(x,3))) %>%
  mutate(txt=paste(or_unadj, "(", ci_lb, "-", ci_ub, ")", sep="")) %>%
  select(metabolite, txt, pval, q_fdr)





#=======================================================================#
## Adjust for each contextual variables ####

# sig_ions_5yr <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions_5yr.csv")

ion_cols_5yr <- sig_ions_5yr$metabolite
names(ion_cols_5yr) <- ion_cols_5yr
ion_cols_5yr <- as.list(ion_cols_5yr)
ion_cols_5yr

n_imp <- 10




imp_glm_sep_5yr <- function(ion) {
  dfs <- foreach(var=context_cols, .combine="rbind", .verbose=TRUE) %:%
    foreach(i=1:n_imp, .combine="rbind", .verbose=TRUE) %do% {
      cols_ <- c(y_col, match_cols, var, ion)
      dat <- full_imp[[i]]

      dat_5yr <- dat[dat$match_id %in% dat$match_id[dat$age_at_diagnosis-dat$sdc_age_calc<=5],]

      lr_dat <- dat_5yr[,colnames(dat_5yr) %in% cols_]
      lr_fit <- glm(gp~., data=lr_dat, family=binomial(link="logit"))
      lr_coef <- summary(lr_fit)$coefficients

      coef <- lr_coef[rownames(lr_coef)==ion, 1]
      pval <- lr_coef[rownames(lr_coef)==ion, 4]

      data.frame(metabolite=ion, context_var=var, beta_adj=coef, pval_adj=pval)
    }
}




bio_5yr <- lapply(ion_cols_5yr, imp_glm_sep_5yr) %>% bind_rows

bio_5yr_comb <- bio_5yr %>%
  full_join(sig_ions_5yr) %>%
  mutate(or_adj=exp(beta_adj),
         or_unadj=exp(beta_unadj)) %>%
  mutate(beta_change=1-beta_adj/beta_unadj,
         or_change = 1-or_adj/or_unadj)

head(bio_5yr_comb)


summary(bio_5yr_comb$beta_change)
summary(bio_5yr_comb$or_change)



# calculate average change in OR
bio_5yr_comb_avg <- bio_5yr_comb %>%
  group_by(metabolite, context_var) %>%
  summarise(or_change_avg=mean(or_change)) %>%
  ungroup()

# filter variables that changed ORs by more than 10%
sig_ion_change_10pct <- bio_5yr_comb_avg %>%
  filter(abs(or_change_avg)>=.1)

sig_ion_change_10pct


sig_ion_change_10pct %>% count(metabolite)

# write_csv(bio_df_comb,
#           file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")


# bio_df_comb <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")







#==================================================#
## Plot distribution of beta change for each metabolite ####


### point plots ####

# try on one metabolite
bio_5yr_comb %>%
  filter(metabolite=="ion_63") %>%
  ggplot()+
  geom_point(aes(context_var, beta_change)) +
  geom_hline(yintercept=0.1) +
  geom_hline(yintercept=-0.1) +
  theme_minimal() +
  theme(axis.text.x=element_blank())




# point plots for all metabolites
bio_5yr_comb %>%
  mutate(metabolite=fct_reorder(metabolite, beta_unadj)) %>%
  ggplot()+
  geom_point(aes(context_var, or_change), size=.8) +
  geom_hline(yintercept=.1) +
  geom_hline(yintercept=-.1) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  scale_y_continuous(breaks=c(-.1,0,.1),
                     labels = scales::percent)+
  facet_wrap(~metabolite, ncol=3) +
  labs(x="Risk factors",
       y="Change in odds ratio estimate (%)")




### boxplots ####


bio_5yr_comb %>%
  mutate(metabolite=fct_reorder(metabolite, beta_unadj)) %>%
  ggplot()+
  geom_boxplot(aes(metabolite, or_change)) +
  geom_hline(yintercept=.1) +
  geom_hline(yintercept=-.1) +
  theme_minimal() +
  scale_x_discrete(guide = guide_axis(angle=90))+
  scale_y_continuous(breaks=c(-.1,0,.1),
                     labels = scales::percent) +
  labs(x="Metabolite",
       y="Change in odds ratio estimate (%)")












###################====================##############################===#

#==================================================#
## Adjust for ALL contextual variables ####



# select list of significant metabolites
ion_cols_5yr <- unique(sig_ions_5yr$metabolite)
names(ion_cols_5yr) <- ion_cols_5yr
ion_cols_5yr <- as.list(ion_cols_5yr)


# function to run through all significant confounders
imp_glm_all <- function(ion) {
  # select columns for logistic fit
  cols_ <- c(y_col, match_cols, context_cols, ion)

  # for each ion, run separate glm models on imputed data
  dfs <- lapply(full_imp, function(dat){

    dat_5yr
    lr_dat <- dat[,colnames(dat) %in% cols_]
    lr_fit <- glm(gp~., data=lr_dat, family=binomial(link="logit"))
    lr_fit
  }
  )
}




# run all glm models on list of metabolites
lr_adj_all <- lapply(ion_cols_5yr, imp_glm_all)

# pool all results
pool_adj_all <- lapply(lr_adj_all, function(x) summary(pool(x), conf.int=TRUE))

# select pooled results for metabolites only
pool_adj_ion <- lapply(pool_adj_all, function(x) x[grepl("^ion", x$term),]) %>%
  bind_rows()

head(pool_adj_ion)



# select pooled results for metabolites only
# also exponentiate coefficient and 95% estimates
pool_adj_ion <- lapply(pool_adj_all, function(x) x[grepl("^ion", x$term),]) %>%
  bind_rows() %>%
  mutate(or_adj=exp(estimate),
         ci_lb=exp(`2.5 %`),
         ci_ub=exp(`97.5 %`)) %>%
  arrange(estimate)


pool_adj_ion




pool_adj_ion %>%
  mutate(across(c(or_adj, ci_lb, ci_ub), function(x) round(x,2))) %>%
  rename(metabolite=term) %>%
  mutate(str_ = paste(or_adj, " (", ci_lb, "-", ci_ub, ")", sep="")) %>%
  select(metabolite, str_)





write_csv(bio_comb_all, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_all_vars.csv")





