#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION ####
#=============================================================#
library(tidyverse)
library(foreach)
library(doParallel)
library(rio)
library(mice)





source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/gmet/00_functions.R")


full_dat <- import_survey()
full_meta <- import_meta_sc()
full_meta$match_id <- toupper(full_meta$match_id)

full_all <- full_dat %>% full_join(full_meta)



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
      lr_formula <- as.formula(paste0(
        "gp~menopause_stt+cohort+collect_age+collect_year+",ion
      ))

      # fit logistic regression model
      lr_mod <- glm(lr_formula, data=dat, family=binomial(link="logit"))
      lr_confint <- confint(lr_mod)

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
                 beta=coef,
                 beta_lb=ci_lb,
                 beta_ub=ci_ub,
                 or=or,
                 or_lb=or_lb,
                 or_ub=or_ub,
                 pval=pval)
    }
  }




#===================================#
### Apply on full data set ####

unadjusted_results_all <- unadj_logreg(full_all)
head(unadjusted_results_all)


#===================================#
### Correction for multiple testing ####
unadjusted_results_all$q_fdr <- p.adjust(unadjusted_results_all$pval, method="fdr")
head(unadjusted_results_all)










write_csv(unadjusted_results_all,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/unadjusted_results_all.csv")





#===================================#
# unadjusted_results_all <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/unadjusted_results_all.csv")

summary(unadjusted_results_all)
summary(unadjusted_results_all$q_fdr)

table(unadjusted_results_all$q_fdr<=0.2)
table(unadjusted_results_all$q_fdr<=0.1)




sig_ions_df <- unadjusted_results_all %>%
  filter(q_fdr<=0.1) %>%
  arrange(beta)

sig_ions_0.1 <- sig_ions_df$metabolite
nonsig_ions_0.1 <- unadjusted_results_all$metabolite[unadjusted_results_all$q_fdr>0.1]

# sig_ions_df





# record all vars for FDR level 0.2
sig_ions_0.2 <- unadjusted_results_all$metabolite[unadjusted_results_all$q_fdr<0.2]



# save data for significant metabolites at fdr 0.1
write_csv(sig_ions_df,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/gmet/sig_ions.csv")






#=======================================================================#
## Adjust for each contextual variables ####


# full_imp <- import_list("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imp_list_%s.csv")
ion_cols <- unadjusted_results_all$metabolite[unadjusted_results_all$q_fdr<0.1]
ion_cols

names(ion_cols) <- ion_cols
ion_cols <- as.list(ion_cols)



n_imp <- 10

adjust_each_var <- function(ion) {
  dfs <- foreach(var=context_cols, .combine="rbind", .verbose=TRUE) %:%
    foreach(i=1:n_imp, .combine="rbind", .verbose=TRUE) %do% {
      cols_ <- c(y_col, match_cols, var, ion)
      dat <- full_imp[[i]]

      lr_formula <- as.formula(paste0(
        "gp ~ cohort + menopause_stt + collect_year + collect_age + ",
        ion, " + ", var
      )
      )

      lr_fit <- glm(lr_formula, data=dat, family=binomial(link="logit"))
      lr_coef <- summary(lr_fit)$coefficients

      coef <- lr_coef[rownames(lr_coef)==ion, 1]
      pval <- lr_coef[rownames(lr_coef)==ion, 4]

      data.frame(metabolite=ion,
                 context_var=var,
                 beta_adj=coef,
                 pval_adj=pval)
    }
}



bio_all <- lapply(ion_cols, adjust_each_var) %>% bind_rows


bio_all_comb <- bio_all %>%
  full_join(sig_ions_df) %>%
  mutate(or_adj=exp(beta_adj),
         or_unadj=exp(beta_unadj)) %>%
  mutate(beta_change=1-beta_adj/beta_unadj,
         or_change = 1-or_adj/or_unadj)

head(bio_all_comb)


summary(bio_all_comb$beta_change)
summary(bio_all_comb$or_change)


# how many are significant confounders?
bio_all_comb %>% count(abs(beta_change)>=0.1 | abs(or_change)>=0.1)



write_csv(bio_all_comb,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")


# bio_all_comb <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")







#==================================================#
## Plot distribution of beta change for each metabolite ####


### point plots ####

# try on one metabolite
bio_all_comb %>%
  filter(metabolite=="ion_63") %>%
  ggplot()+
  geom_point(aes(context_var, beta_change)) +
  geom_hline(yintercept=0.1) +
  geom_hline(yintercept=-0.1) +
  theme_minimal() +
  theme(axis.text.x=element_blank())




# point plots for all metabolites
bio_all_comb %>%
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


bio_all_comb %>%
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
ion_cols <- sig_ions_df$metabolite
names(ion_cols) <- ion_cols
ion_cols <- as.list(ion_cols)


imp_glm_all <- function(ion) {
  # select columns for logistic fit
  cols_ <- c(y_col, match_cols, context_cols, ion)

  # for each ion, run separate glm models on imputed data
  dfs <- lapply(full_imp, function(dat){
    lr_dat <- dat[,colnames(dat) %in% cols_]
    lr_fit <- glm(gp~., data=lr_dat, family=binomial(link="logit"))
    lr_fit
  }
  )
}

# run all glm models on list of metabolites
lr_adj_all <- lapply(ion_cols, imp_glm_all)

# pool all results
pool_adj_all <- lapply(lr_adj_all, function(x) summary(pool(x), conf.int=TRUE))

# select pooled results for metabolites only
pool_adj_ion <- lapply(pool_adj_all, function(x) x[grepl("^ion", x$term),]) %>%
  bind_rows

head(pool_adj_ion)



bio_comb_all <- pool_adj_ion %>%
  full_join(sig_ions_df,
            by=join_by(term==metabolite),
            keep=TRUE) %>%
  mutate(or_unadj=exp(beta_unadj),
         or_adj=exp(estimate)) %>%
  mutate(or_change = 1-or_adj/or_unadj)

bio_comb_all %>% filter(abs(or_change)>=.1)


bio_comb_all %>%
  mutate(
    or_unadj=round(exp(beta_unadj),2),
    lb_unadj=round(exp(ci_lb),2),
    ub_unadj=round(exp(ci_ub),2),
    or_adj=round(exp(estimate),2),
    lb_adj=round(exp(`2.5 %`),2),
    ub_adj=round(exp(`97.5 %`),2)) %>%
  mutate(
    unadj=paste(or_unadj, "(", lb_unadj, "-", ub_unadj, ")", sep=""),
    adj = paste(or_adj, "(", lb_adj, "-", ub_adj, ")", sep="")) %>%
  select(metabolite, unadj, adj)



bio_comb_all %>% filter(abs(beta_change)>=0.1)



write_csv(bio_comb_all, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_all_vars.csv")




bio_comb_all %>%
  mutate(lb_unadj=exp(ci_lb),
         ub_unadj=exp(ci_ub),
         lb_adj=exp(`2.5 %`),
         ub_adj=exp(`97.5 %`)) %>%
  mutate(across(c(ends_with("unadj"), ends_with("adj")), function(x) round(x,2))) %>%
  mutate(ion=term,
         term_unadj=paste(or_unadj, "(", lb_unadj, "-", ub_unadj, ")", sep=""),
         term_adj=paste(or_adj, "(", lb_adj, "-", ub_adj, ")", sep="")) %>%
  select(ion, term_unadj, term_adj)

