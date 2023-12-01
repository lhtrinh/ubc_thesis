#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION (Postmenopausal) ####
#=============================================================#
library(tidyverse)
library(foreach)
library(doParallel)
library(rio)
library(mice)





source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


full_dat <- import_survey()
full_meta <- import_meta_sc()


full_all <- full_dat %>% full_join(full_meta)



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
      cols_ <- c(y_col, "collect_age", "collect_year_cat", "cohort", ion)
      lr_dat <- dat[,colnames(dat) %in% cols_]

      # fit logistic regression model
      lr_mod <- glm(formula=gp~., data=lr_dat, family=binomial(link="logit"))
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
full_postmnp <- full_all %>% filter(menopause_stt==1)

all_pvals_mnp <- unadj_logreg_mnp(full_postmnp)
head(all_pvals_mnp)


#===================================#
### Correction for multiple testing ####
all_pvals_mnp$q_fdr <- p.adjust(all_pvals$pval, method="fdr")
head(all_pvals_mnp)



#===================================#
### Manual calculation for BH correction ####


pvals <- all_pvals_mnp %>% select(metabolite, pval)

# pvalues <- pvals$pval


# assign variables to calculate critical value
ranks <- rank(pvals$pval, ties.method = "last")
m <- length(pvals$pval)
Q <- 0.1

# calculate critical value
critical_val <- (ranks/m)*Q
all_pvals_mnp$critical_val <- critical_val



# check: how many p-values are smaller than their critical values?
which(pvals$pval<critical_val)

# find the largest pvalue that is smaller than its critical value
largest_pval <- max(pvals$pval[pvals$pval<critical_val])

# how many are there?
length(pvals$pval[pvals$pval<=largest_pval])

pvals[pvals$pval<=largest_pval,]
critical_val[pvals$pval<=largest_pval]

#===================================#




write_csv(all_pvals,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/all_pvals.csv")





#===================================#
# all_pvals <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/all_pvals.csv")

summary(all_pvals_mnp)
summary(all_pvals_mnp$q_fdr)


all_pvals_mnp %>% count(q_fdr<=0.1 | q_fdr<=0.2)




sig_ions_mnp <- all_pvals_mnp %>%
  filter(q_fdr<=0.1) %>%
  arrange(beta_unadj)


sig_ions_mnp





# save data for significant metabolites at fdr 0.1
write_csv(sig_ions_mnp,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions_df.csv")





# ## ONLY FOR WRITING
# x_mnp <- sig_ions_mnp %>%
#   mutate(across(c(beta_unadj, ci_lb, ci_ub), function(x) round(exp(x), 2))) %>%
#   mutate(ci = paste(ci_lb, "-", ci_ub, sep="")) %>%
#   select(metabolite, beta_unadj, ci, pval, q_fdr, critical_val)
#
# write_csv(x_mnp, file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Output/x_mnp.csv")



#=======================================================================#
## Adjust for each contextual variables ####


# full_imp <- import_list("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imp_list_%s.csv")



# select list of significant metabolites
ion_cols <- sig_ions_mnp$metabolite
names(ion_cols) <- ion_cols
ion_cols <- as.list(ion_cols)


imp_glm_sep_mnp <- function(ion) {
  dfs <- foreach(var=context_cols, .combine="rbind", .verbose=TRUE) %:%
    foreach(i=1:n_imp, .combine="rbind", .verbose=TRUE) %do% {
      cols_ <- setdiff(c(y_col, match_cols, var, ion), c("menopause_stt"))
      dat <- full_imp[[i]]

      lr_dat <- dat[,colnames(dat) %in% cols_]
      lr_fit <- glm(formula=gp~., data=lr_dat, family=binomial(link="logit"))
      lr_coef <- summary(lr_fit)$coefficients

      coef <- lr_coef[rownames(lr_coef)==ion, 1]
      pval <- lr_coef[rownames(lr_coef)==ion, 4]

      data.frame(metabolite=ion, context_var=var, beta_adj=coef, pval_adj=pval)
    }
}




dfs <- lapply(ion_cols, imp_glm_sep_mnp)

bio_df <- dfs %>% bind_rows %>% rename(pval_adj=pval)

bio_df_comb <- bio_df %>%
  full_join(sig_ions_mnp) %>%
  mutate(or_adj=exp(beta_adj),
         or_unadj=exp(beta_unadj)) %>%
  mutate(beta_change=1-beta_adj/beta_unadj,
         or_change = 1-or_adj/or_unadj)

head(bio_df_comb)


summary(bio_df_comb$beta_change)
summary(bio_df_comb$or_change)


# how many are significant confounders?
bio_df_comb %>% count(abs(beta_change)>=0.1 | abs(or_change)>=0.1)



write_csv(bio_df_comb,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")


# bio_df_comb <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")







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












###################====================##############################===#

#==================================================#
## Adjust for ALL contextual variables ####

# select list of significant metabolites
ion_cols <- as.list(sig_ions_0.1)
names(ion_cols) <- sig_ions_0.1


imp_glm_all <- function(ion) {
  # select columns for logistic fit
  cols_ <- setdiff(c(y_col, match_cols, context_cols, ion), c("menopause_stt"))

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
  bind_rows()

head(pool_adj_ion)



bio_comb_all <- pool_adj_ion %>%
  rename(metabolite=term,
         beta_adj=estimate) %>%
  full_join(sig_ions_df) %>%
  mutate(or_unadj=exp(beta_unadj),
         or_adj=exp(beta_adj)) %>%
  mutate(or_change = 1-or_adj/or_unadj)

bio_comb_all %>% filter(abs(or_change)>=.1)


bio_comb_all %>%
  mutate(or_adj=round(or_adj,2),
         lb_adj=round(exp(`2.5 %`),2),
         ub_adj=round(exp(`97.5 %`),2)) %>%
  mutate(str_ = paste(or_adj, " (", lb_adj, "-", ub_adj, ")", sep="")) %>%
  select(metabolite, str_)



bio_comb_all %>% filter(abs(beta_change)>=0.1)



write_csv(bio_comb_all, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_all_vars.csv")





