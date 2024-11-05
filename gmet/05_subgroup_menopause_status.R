#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION (Postmenopausal) ####
#=============================================================#
library(tidyverse)
library(foreach)
library(doParallel)
library(rio)
# library(mice)





source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


full_dat <- import_survey()
full_meta <- import_meta_sc()


full_all <- full_dat %>% full_join(full_meta)


# import imputed data (with metabolites)
full_imp <- import_imp()


#=======================================================================#
## Adjust only for matching factors ####

### Define logreg function ####
unadj_logreg_mnp <- function(dat){
  ion_cols <- colnames(dat)[grepl("^ion", colnames(dat))]
  n <- length(ion_cols)

  lr_out <-
    foreach (i=1:n, .combine="rbind", .verbose=TRUE) %do% {
      ion <- as.character(ion_cols[i])

      # create data set with only required column names
      # (vectors of column names defined in "00_functions.R")
      cols_ <- setdiff(c(y_col, match_cols, ion), "menopause_stt")
      lr_dat <- dat[,colnames(dat) %in% cols_]

      # fit logistic regression model
      lr_mod <- glm(formula=gp~., data=lr_dat, family=binomial(link="logit"))
      lr_confint <- confint.default(lr_mod)

      # extract p-value for metabolite
      lr_coef <- summary(lr_mod)$coefficients

      coef <- lr_coef[rownames(lr_coef)==ion, 1]
      pval <- lr_coef[rownames(lr_coef)==ion, 4]
      ci_lb <- lr_confint[rownames(lr_confint)==ion, 1]
      ci_ub <- lr_confint[rownames(lr_confint)==ion, 2]
      or <- exp(coef)
      or_lb <- exp(ci_lb)
      or_ub <- exp(ci_ub)

      data.frame(metabolite=ion,
                 beta_unadj=coef,
                 ci_lb=ci_lb,
                 ci_ub=ci_ub,
                 or=or,
                 or_lb=or_lb,
                 or_ub=or_ub,
                 pval=pval)
    }
}




#===================================#
### Apply on menopausal set ####
full_mnp <- full_all %>% filter(menopause_stt==1)
all_pvals_mnp <- unadj_logreg_mnp(full_mnp)
head(all_pvals_mnp)





#===================================#
### Correction for multiple testing ####
all_pvals_mnp$q_fdr <- p.adjust(all_pvals_mnp$pval, method="fdr")
head(all_pvals_mnp)




write_csv(all_pvals_mnp,
          file = "Data/gmet/all_pvals_mnp_wald.csv")





#===================================#
# all_pvals_mnp <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/all_pvals_mnp.csv")

summary(all_pvals_mnp)
summary(all_pvals_mnp$q_fdr)


all_pvals_mnp %>% count(q_fdr<=0.1)




sig_ions_mnp <- all_pvals_mnp %>%
  filter(q_fdr<=0.1) %>%
  arrange(beta_unadj)


sig_ions_mnp





# save data for significant metabolites at fdr 0.1
write_csv(sig_ions_mnp,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions_mnp.csv")



sig_ions_mnp %>%
  mutate(or_unadj=exp(beta_unadj),
         lb_unadj=exp(ci_lb),
         ub_unadj=exp(ci_ub)) %>%
  mutate(across(ends_with("unadj"), function(x) round(x,2))) %>%
  mutate(unadj=paste(or_unadj, "(", lb_unadj, "-", ub_unadj, ")", sep="")) %>%
  select(metabolite, unadj)





#=======================================================================#
## Adjust for each contextual variables ####



# select list of significant metabolites
ion_cols_mnp <- sig_ions_mnp$metabolite
names(ion_cols_mnp) <- ion_cols_mnp
ion_cols_mnp <- as.list(ion_cols_mnp)


imp_glm_sep_mnp <- function(ion) {
  dfs <- foreach(var=context_cols, .combine="rbind", .verbose=TRUE) %:%
    foreach(i=1:n_imp, .combine="rbind", .verbose=TRUE) %do% {
      # remove menopause_stt from list of matching vars
      cols_ <- setdiff(c(y_col, match_cols, var, ion), c("menopause_stt"))
      dat <- full_imp[[i]]

      lr_dat <- dat[dat$menopause_stt==1,colnames(dat) %in% cols_]
      lr_fit <- glm(formula=gp~., data=lr_dat, family=binomial(link="logit"))
      lr_coef <- summary(lr_fit)$coefficients

      coef <- lr_coef[rownames(lr_coef)==ion, 1]
      pval <- lr_coef[rownames(lr_coef)==ion, 4]

      data.frame(metabolite=ion, context_var=var, beta_adj=coef, pval_adj=pval)
    }
}




# dfs <- lapply(ion_cols_mnp, imp_glm_sep_mnp)

bio_mnp <- lapply(ion_cols_mnp, imp_glm_sep_mnp) %>% bind_rows

# combine with unadjusted coefficients
# and calculate beta change ratio
bio_mnp_comb <- bio_mnp %>%
  full_join(sig_ions_mnp) %>%
  mutate(or_adj=exp(beta_adj),
         or_unadj=exp(beta_unadj)) %>%
  mutate(beta_change=1-beta_adj/beta_unadj,
         or_change = 1-or_adj/or_unadj)

head(bio_mnp_comb)


summary(bio_mnp_comb$beta_change)
summary(bio_mnp_comb$or_change)


# how many are significant confounders?
bio_mnp_comb %>% filter(abs(or_change)>=0.1) %>% count(metabolite)


# calculate average change in OR
bio_mnp_comb_avg <- bio_mnp_comb %>%
  group_by(metabolite, context_var) %>%
  summarise(or_change_avg=mean(or_change)) %>%
  ungroup()

# filter variables that changed ORs by more than 10%
sig_ion_change_10pct <- bio_mnp_comb_avg %>%
  filter(abs(or_change_avg)>=.1)

sig_ion_change_10pct %>% count(metabolite)


# write_csv(bio_df_comb,
#           file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")


# bio_df_comb <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")


















#======================================================#
## Adjust for ALL contextual variables ####
sig_ions_mnp <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions_mnp.csv")

sig_ions_mnp


# select list of significant metabolites
ion_cols_mnp <- sig_ions_mnp$metabolite
names(ion_cols_mnp) <- ion_cols_mnp
ion_cols_mnp <- as.list(ion_cols_mnp)
ion_cols_mnp

# function to run through all significant confounders
imp_glm_mnp <- function(ion) {
  # select columns for logistic fit
  cols_ <- setdiff(c(y_col, match_cols, context_cols, ion), c("menopause_stt"))


  # for each ion, run separate glm models on imputed data
  dfs <- lapply(full_imp, function(dat){
    dat_mnp <- dat[dat$menopause_stt==1,]
    lr_dat <- dat_mnp[,colnames(dat_mnp) %in% cols_]
    lr_fit <- glm(gp~., data=lr_dat, family=binomial(link="logit"))
    lr_fit
  }
  )
}



# run all glm models on list of metabolites
lr_adj_all <- lapply(ion_cols_mnp, imp_glm_mnp)

# pool all results
pool_adj_all <- lapply(lr_adj_all, function(x) summary(pool(x), conf.int=TRUE))
head(pool_adj_all)


# select pooled results for metabolites only
# also exponentiate coefficient and 95% estimates
pool_adj_ion <- lapply(pool_adj_all, function(x) x[grepl("^ion", x$term),]) %>%
  bind_rows() %>%
  full_join(sig_ions_mnp,
            by=join_by(term==metabolite),
            keep=TRUE)

pool_adj_ion




pool_adj_ion %>%
  mutate(or_unadj=exp(beta_unadj),
         lb_unadj=exp(ci_lb),
         ub_unadj=exp(ci_ub),
         or_adj=exp(estimate),
         lb_adj=exp(`2.5 %`),
         ub_adj=exp(`97.5 %`)) %>%
  mutate(across(ends_with("adj"), function(x) round(x,2))) %>%
  mutate(
    unadj = paste(or_unadj, "(", lb_unadj, "-", ub_unadj, ")", sep=""),
    adj = paste(or_adj, "(", lb_adj, "-", ub_adj, ")", sep="")) %>%
  select(term, unadj, adj)




write_csv(bio_comb_all, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_all_vars.csv")






#==================================================#
## Plot distribution of beta change for each metabolite ####


### point plots ####

# try on one metabolite
bio_df_comb %>%
  filter(metabolite=="ion_63") %>%
  ggplot()+
  geom_point(aes(context_var, beta_change)) +
  geom_hline(yintercept=0.1) +
  geom_hline(yintercept=-0.1) +
  theme_minimal() +
  theme(axis.text.x=element_blank())




# point plots for all metabolites
bio_df_comb %>%
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


bio_df_comb %>%
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
