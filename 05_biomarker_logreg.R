#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION ####
#=============================================================#
library(tidyverse)
library(foreach)



## Variable categories ####

id_cols <- c("studyid", "match_id")

y_col <- "gp"

match_cols <-  c("cohort", "collect_yr_cat", "collect_age", "menopause_stt")

context_cols <- c("sdc_age_calc", "ethnicity",
                  "edu_level", "income_level",
                  "wh_menstruation_age_cat", "wh_menopause_age_cat",
                  "wh_contraceptives_ever",
                  "wh_gravidity", "wh_live_births", "wh_preg_first_age_cat",
                  "wh_breastfeeding_cat",
                  "wh_hrt_ever", "wh_hrt_duration_yr",
                  "fam_hist_breast",
                  "bmi", "bmi_cat",
                  "alc_ever", "alc_cur_freq_cat", "alc_binge_cat",
                  "smk_cig_status", "smk_cig_dur", "smk_first_age_cat",
                  "pse_childhood_duration", "pse_adult_home_duration", "pse_adult_wrk_duration")




## Transform variables to factor ####
full_dat <- full_dat %>%
  mutate(across(c(gp, cohort,
                  menopause_stt, collect_yr_cat,
                  ethnicity,
                  edu_level,
                  income_level,
                  wh_menstruation_age_cat, wh_menopause_age_cat,
                  wh_contraceptives_ever,
                  wh_breastfeeding_cat,
                  wh_preg_first_age_cat,
                  wh_hrt_ever,
                  fam_hist_breast,
                  alc_ever, alc_cur_freq_cat, alc_binge_cat,
                  smk_cig_status, smk_cig_dur, smk_first_age_cat,
                  pse_childhood_duration, pse_adult_home_duration, pse_adult_wrk_duration,
                  bmi_cat,
                  er_status, pr_status, her2_status, hr_subtype,
                  site_code, hist_code, hist_subtype),
                as_factor))


## Merge questionnaire and metabolomics data ####
full_all <- full_centered %>%
  full_join(full_dat)

write_csv(full_all, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/full_data_with_normalized_metabolomics.csv")



full_all <- read_csv(
  "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/full_data_with_normalized_metabolomics.csv"
) %>%
  mutate(across(c(gp, cohort,
                  menopause_stt, collect_yr_cat,
                  ethnicity,
                  edu_level,
                  income_level,
                  wh_menstruation_age_cat, wh_menopause_age_cat,
                  wh_contraceptives_ever,
                  wh_breastfeeding_cat,
                  wh_preg_first_age_cat,
                  wh_hrt_ever,
                  fam_hist_breast,
                  alc_ever, alc_cur_freq_cat, alc_binge_cat,
                  smk_cig_status, smk_first_age_cat,
                  bmi_cat,
                  er_status, pr_status, her2_status, hr_subtype,
                  site_code, hist_code, hist_subtype),
                as_factor))



#=======================================================================#
## Define unadjusted logistic regression function ####
unadj_logistic <- function(dat){
  ion_cols <- colnames(dat)[grepl("^ion", colnames(dat))]
  n <- length(ion_cols)

  logreg_pvals <- data.frame(
    metabolite=ion_cols,
    beta_unadj = rep(0, n),
    pval=rep(0, n)
  )

  foreach (i=1:n) %do% {
    ion <- ion_cols[i]

    cols_ <- c(y_col, match_cols, ion)
    lr_dat <- dat[,colnames(dat) %in% cols_]

    # fit logistic regression model
    lr_mod <- glm(gp~., data=lr_dat, family=binomial(link="logit"))

    # extract p-value for metabolite
    lr_coef <- summary(lr_mod)$coefficients
    coef <- lr_coef[rownames(lr_coef)==ion, 1]
    pval <- lr_coef[rownames(lr_coef)==ion, 4]

    logreg_pvals[i,1] <- ion
    logreg_pvals[i,2] <- coef
    logreg_pvals[i,3] <- pval
  }


  # adjust for p-values with FDR
  logreg_pvals$p_adjust_fdr <- p.adjust(logreg_pvals$pval, method="BH")

  logreg_pvals
}





#=======================================================================#
## Full data set ####
#=======================================================================#

### Unadjusted ####

all_pvals <- unadj_logistic(full_all)


summary(all_pvals)
summary(all_pvals$p_adjust_fdr)

table(all_pvals$p_adjust_fdr<=0.2)



sig_ions <- all_pvals %>%
  filter(p_adjust_fdr<=0.2)


write_csv(sig_ions,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions.csv")



#=======================================================================#
### Adjust for each contextual variables ####

###
# list of imputed data
dat_imp <- list()
for (i in 1:5){
  full_dat <- cbind(dat_pre_imp_id, complete(init, action=i))
  dat_imp[[i]] <- full_join(full_centered, full_dat)
}
###



# select list of significant metabolites
ion_cols <- sig_ions$metabolite


bio_df <- data.frame(
  df_ver = rep(1:5, each=length(ion_cols)*length(context_cols)),
  metabolite = rep(ion_cols, each=length(context_cols)),
  context_var = context_cols
)
head(bio_df)
dim(bio_df)





### Option a: loop without function ####
start_time=Sys.time()
for (i in 1:nrow(bio_df)){

  # extract imputed data set
  dat <- dat_imp[[bio_df$df_ver[i]]]

  # select metabolite variable
  ion <- bio_df$metabolite[i]
  context_col <- bio_df$context_var[i]

  # select column names for logistic regression model
  cols_adj <- c(y_col, match_cols, context_col, ion)

  # subset data to select only relevant columns
  lr_dat <- dat[,colnames(full_all) %in% cols_adj]

  # logreg fit
  adj_fit <- glm(gp~., data=lr_dat, family=binomial(link="logit"))
  lr_coef <- summary(adj_fit)$coefficients

  # report adjusted coefficient and pval for metabolite
  coef <- lr_coef[rownames(lr_coef)==ion, 1]
  pval <- lr_coef[rownames(lr_coef)==ion, 4]


  # save coefficient and pval in df
  bio_df$beta_adj[i] <- coef
  bio_df$pval_adj[i] <- pval
}
end_time=Sys.time()
end_time-start_time






# calculate adjusted to unadjusted ratio
bio_df_comb <- bio_df %>%
  full_join(sig_ions) %>%
  mutate(beta_ratio = abs(1-beta_adj/beta_unadj))

bio_df_comb[bio_df_comb$beta_ratio>.10,]



write_csv(bio_df_comb,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results.csv")


#==================================================#
### Adjust for ALL contextual variables ####

# select list of significant metabolites
ion_cols <- sig_ions$metabolite


bio_df <- data.frame(
  df_ver = rep(1:5, each=length(ion_cols)),
  metabolite = rep(ion_cols, 5)
)


start_time=Sys.time()
for (i in 1:nrow(bio_df)){

  # extract imputed data set
  dat <- dat_imp[[bio_df$df_ver[i]]]

  # select metabolite variable
  ion <- bio_df$metabolite[i]

  # select column names for logistic regression model
  cols_adj <- c(y_col, match_cols, context_cols, ion)

  # subset data to select only relevant columns
  lr_dat <- dat[,colnames(full_all) %in% cols_adj]

  # logreg fit
  adj_fit <- glm(gp~., data=lr_dat, family=binomial(link="logit"))
  lr_coef <- summary(adj_fit)$coefficients

  # report adjusted coefficient and pval for metabolite
  coef <- lr_coef[rownames(lr_coef)==ion, 1]
  pval <- lr_coef[rownames(lr_coef)==ion, 4]


  # save coefficient and pval in df
  bio_df$beta_adj[i] <- coef
  bio_df$pval_adj[i] <- pval
}
end_time=Sys.time()
end_time-start_time



bio_df_comb <- bio_df %>%
  full_join(sig_ions) %>%
  mutate(beta_ratio = abs(1-beta_adj/beta_unadj))

bio_df_comb %>% filter(beta_ratio<.05) %>%
  count(df_ver)

