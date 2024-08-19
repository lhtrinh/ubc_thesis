#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION ####
#=============================================================#
library(tidyverse)
library(foreach)
library(doParallel)
library(rio)
library(mice) # for pooling results
library(survival) # for conditional logistic regression


source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/exposome/00_functions.R")

full_dat <- import_survey()
full_meta <- import_meta_sc()
full_imp <- import_imp()


# combine questionnaire and metabolomics data
full_all <- full_dat %>% full_join(full_meta)


full_all %>% count(gp)


#=======================================================================#
## UNCONDITIONAL LOGISTIC REGRESSION ####


### Define unconditional logistic function ####
unadj_ulr <- function(dat){

  # define ion columns
  ion_cols <- colnames(dat)[grepl("^ion", colnames(dat))]
  n <- length(ion_cols)

  # loop run conditional logistic regression for each metabolite
  lr_out <-
    foreach (i=1:n, .combine="rbind", .verbose=TRUE) %dopar% {
      ion <- as.character(ion_cols[i])

      # create model formula
      lr_formula <- as.formula(paste0(
        "gp ~ ",
        ion,
        " + menopause_stt + collect_year + collect_age + cohort"))

      # fit CLR model
      lr_mod <- glm(lr_formula, data=dat, family = binomial(link="logit"))
      lr_summ <- summary(lr_mod)
      lr_confint <- confint(lr_mod)

      beta=coef(lr_mod)[names(coef(lr_mod))==ion]
      beta_lb=lr_confint[rownames(lr_confint)==ion][1]
      beta_ub=lr_confint[rownames(lr_confint)==ion][2]
      or=exp(beta)
      or_lb=exp(beta_lb)
      or_ub=exp(beta_ub)
      pval=lr_summ$coefficients[rownames(lr_summ$coefficients)==ion][4]


      data.frame(ion=ion,
                 beta=beta,
                 beta_lb=beta_lb,
                 beta_ub=beta_ub,
                 or=or,
                 or_lb=or_lb,
                 or_ub=or_ub,
                 pval=pval)
    }
  }



#===================================#
### Apply on full data set ####

clust <- makeCluster(4)
setDefaultCluster(clust)
registerDoParallel(clust)

clusterEvalQ(clust, {
  library(doParallel)
}
)

clusterExport(clust, c("full_all"), envir = .GlobalEnv)


start_t <- Sys.time()
ulr_results <- unadj_ulr(full_all)
end_t <- Sys.time()


end_t - start_t
stopCluster(clust)



#===================================#
### View results ####
head(ulr_results)

#===================================#
### Correction for multiple testing ####
ulr_results$qval <- p.adjust(ulr_results$pval, method="fdr")
head(ulr_results)

# significant ions at 0.1
ulr_results %>%
  filter(qval<=0.1) %>%
  select(ion, or, or_lb, or_ub, pval, qval)

length(which(ulr_results$qval<0.1))
length(which(ulr_results$qval<0.2))




write_csv(ulr_results,
          file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/exposome/overall_lr_results.csv")

