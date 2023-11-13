#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION ####
#=============================================================#
library(tidyverse)
library(foreach)

source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")


full_dat <- import_survey()
full_meta <- import_meta_sc()


full_all <- full_dat %>% full_join(full_meta)



#=======================================================================#
## Define unadjusted logistic regression function ####
unadj_logistic <- function(dat){
  ion_cols <- colnames(dat)[grepl("^ion", colnames(dat))]
  n <- length(ion_cols)


  logreg_pvals <- data.frame(
    metabolite=rep("", n),
    or_unadj = rep(0, n),
    ci_lb = rep(0, n),
    ci_ub = rep(0, n),
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
    lr_ci <- confint(lr_mod)

    coef <- lr_coef[rownames(lr_coef)==ion, 1]
    pval <- lr_coef[rownames(lr_coef)==ion, 4]
    ci_lb <- lr_ci[rownames(lr_ci)==ion, 1]
    ci_ub <- lr_ci[rownames(lr_ci)==ion, 2]

    logreg_pvals$metabolite[i] <- ion
    logreg_pvals$or_unadj[i] <- exp(coef)
    logreg_pvals$ci_lb[i] <- exp(ci_lb)
    logreg_pvals$ci_ub[i] <- exp(ci_ub)
    logreg_pvals$pval[i] <- pval
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
table(all_pvals$p_adjust_fdr<=0.1)


all_pvals %>% filter(p_adjust_fdr<=0.2)


sig_ions <- all_pvals %>%
  filter(p_adjust_fdr<=0.1) %>%
  arrange(or_unadj)

sig_ions


write_csv(sig_ions,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions.csv")



#=======================================================================#
### Adjust for each contextual variables ####


# select list of significant metabolites
ion_cols <- sig_ions$metabolite


bio_df <- data.frame(
  df_ver = rep(1:5, each=length(ion_cols)*length(context_cols)),
  metabolite = rep(ion_cols, each=length(context_cols)),
  context_var = context_cols
)
head(bio_df)
dim(bio_df)



# reiteratively fit logistic regression functions for each risk factor
start_time=Sys.time()
for (i in 1:nrow(bio_df)){

  # extract imputed data set
  dat <- full_imp[[bio_df$df_ver[i]]]

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
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")


bio_df_comb1 <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_single_var.csv")


head(bio_df_comb1)

# those altered by more than 10%
bio_df_comb1 %>% filter(beta_ratio>.10)
bio_df_comb1 %>% filter(beta_ratio<.05)

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

  # extract imputed data version
  dat <- full_imp[[bio_df$df_ver[i]]]

  # select metabolite variable
  ion <- bio_df$metabolite[i]

  # select column names for logistic regression model
  cols_adj <- c(y_col, match_cols, context_cols, ion)

  # subset data to select only relevant columns
  lr_dat <- dat[,colnames(dat) %in% cols_adj]

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



bio_df_comb2 <- bio_df %>%
  full_join(sig_ions) %>%
  mutate(beta_ratio = abs(1-beta_adj/beta_unadj))

bio_df_comb2 %>% filter(beta_ratio>.10) %>% count(metabolite) %>% arrange(desc(n))
bio_df_comb2 %>% filter(beta_ratio<.05) %>% count(metabolite) %>% arrange(desc(n))



write_csv(bio_df_comb2, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/adjusted_logreg_results_all_vars.csv")
