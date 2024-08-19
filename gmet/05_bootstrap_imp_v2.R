library(parallel)
library(mice)
library(doParallel)
library(tidyverse)
library(readr)
library(rio)



#=============================================================#
# Setup ####




sig_all <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions_full_data.csv") %>%
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
                as.factor)) %>%
  mutate(across(c(sdc_age_calc, wh_gravidity, wh_live_births),
                as.integer))



#====================================#
## Assign variable groups ####
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

imp_cols <- c(y_col, match_cols, context_cols)



fct_cols <- c("gp", "cohort",
              "menopause_stt", "collect_yr_cat",
              "ethnicity",
              "edu_level",
              "income_level",
              "wh_menstruation_age_cat", "wh_menopause_age_cat",
              "wh_contraceptives_ever",
              "wh_breastfeeding_cat",
              "wh_preg_first_age_cat",
              "wh_hrt_ever",
              "fam_hist_breast",
              "alc_ever", "alc_cur_freq_cat", "alc_binge_cat",
              "smk_cig_status", "smk_first_age_cat",
              "bmi_cat",
              "er_status", "pr_status", "her2_status", "hr_subtype",
              "site_code", "hist_code", "hist_subtype")

int_cols <- c("sdc_age_calc", "wh_gravidity", "wh_live_births")


#====================================#
## Import bootstrapped data ####
path <- Sys.glob("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/*.csv")
boot_list_import <- import_list(path)




for (i in 1:length(boot_list_import)){
  boot_list_import[[i]][boot_list_import[[i]]==""] <- NA
  boot_list_import[[i]][fct_cols] <- lapply(boot_list_import[[i]][fct_cols], as_factor)
  boot_list_import[[i]][int_cols] <- lapply(boot_list_import[[i]][int_cols], as.integer)
}




#=============================================================#
# Impute bootstrapped data ####
clust <- makeCluster(5)
setDefaultCluster(clust)

registerDoParallel(clust)

# clusterSetRNGStream(clust, 9967)

clusterEvalQ(clust, {
  library(mice);
  library(doParallel)
}
)

# test_list <- boot_list_import[1:2]

clusterExport(clust, c("test_list", "imp_cols"), envir = .GlobalEnv)




#====================================#
### method 2: use foreach loop ####

#====================================#
#### data 1-25 ####
clust <- makeCluster(3)
setDefaultCluster(clust)
registerDoParallel(clust)


boot_chunk_1_25 <- boot_list_import[1:25]


start_t <- Sys.time()
boot_imp_1_25 <-
  foreach(i=1:length(boot_chunk_1_25), .packages="mice", .combine="c", .errorhandling = "remove", .verbose=TRUE) %dopar% {
    df <- boot_chunk_1_25[[i]]

    pre_imp <- df[,colnames(df) %in% imp_cols]
    pre_imp_id <- df[,!colnames(df) %in% imp_cols]

    imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 3)
    post_imp <- complete(imp_pars, action = "all")

    full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))
    names(full_imp) <- paste("imp", i, "_dat", names(full_imp), sep="")
    full_imp
  }
end_t <- Sys.time()


end_t - start_t
stopCluster(clust)


export_list(boot_imp_1_25,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imp_1_25/%s.csv")



#====================================#
#### data 26-50 ####
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)


boot_chunk <- boot_list_import[26:50]

start_t <- Sys.time()
boot_imp <-
  foreach(i=1:length(boot_chunk),
          .packages="mice",
          .combine="c",
          .errorhandling = "remove",
          .verbose=TRUE) %dopar% {
            df <- boot_chunk[[i]]

            pre_imp <- df[,colnames(df) %in% imp_cols]
            pre_imp_id <- df[,!colnames(df) %in% imp_cols]

            imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 3)
            post_imp <- complete(imp_pars, action = "all")

            full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))
            names(full_imp) <- paste("imp", i+25, "_dat", names(full_imp), sep="")
            full_imp
            }
end_t <- Sys.time()
end_t - start_t


stopCluster(clust)


export_list(boot_imp,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imp_26_50/%s.csv")



#====================================#
#### data 51-75 ####
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)


boot_chunk <- boot_list_import[51:75]


start_t <- Sys.time()

boot_imp <-
  foreach(i=1:length(boot_chunk),
          .packages="mice",
          .combine="c",
          .errorhandling = "remove",
          .verbose=TRUE) %dopar% {
            df <- boot_chunk[[i]]

            pre_imp <- df[,colnames(df) %in% imp_cols]
            pre_imp_id <- df[,!colnames(df) %in% imp_cols]

            # imputation
            imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 3)
            post_imp <- complete(imp_pars, action = "all")

            # combine with other data
            full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))
            names(full_imp) <- paste("imp", i+50, "_dat", names(full_imp), sep="")
            full_imp
            }

end_t <- Sys.time()
end_t - start_t


stopCluster(clust)


export_list(boot_imp,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imputed_boot_data/%s.csv")




#====================================#
#### data 76-100 ####
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)


boot_chunk <- boot_list_import[76:100]

start_t <- Sys.time()
boot_imp <-
  foreach(i=1:length(boot_chunk),
          .packages="mice",
          .combine="c",
          .errorhandling = "remove",
          .verbose=TRUE) %dopar% {
            df <- boot_chunk[[i]]

            pre_imp <- df[,colnames(df) %in% imp_cols]
            pre_imp_id <- df[,!colnames(df) %in% imp_cols]

            imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 5)
            post_imp <- complete(imp_pars, action = "all")

            full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))

            # change i accordingly
            names(full_imp) <- paste("imp", i+75, "_dat", names(full_imp), sep="")
            full_imp
          }


end_t <- Sys.time()
end_t - start_t


stopCluster(clust)


export_list(boot_imp,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imputed_boot_data/%s.csv")


#====================================#
#### data 101-125 ####
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)


boot_chunk <- boot_list_import[101:125]

start_t <- Sys.time()
boot_imp <-
  foreach(i=1:length(boot_chunk),
          .packages="mice",
          .combine="c",
          .errorhandling = "remove",
          .verbose=TRUE) %dopar% {
            df <- boot_chunk[[i]]

            pre_imp <- df[,colnames(df) %in% imp_cols]
            pre_imp_id <- df[,!colnames(df) %in% imp_cols]

            imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 5)
            post_imp <- complete(imp_pars, action = "all")

            full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))
            names(full_imp) <- paste("imp", i+100, "_dat", names(full_imp), sep="")
            full_imp
          }

end_t <- Sys.time()
end_t - start_t


stopCluster(clust)


export_list(boot_imp,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imputed_boot_data/%s.csv")




#====================================#
#### data 126-150 ####
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)


boot_chunk <- boot_list_import[126:150]

start_t <- Sys.time()
boot_imp <-
  foreach(i=1:length(boot_chunk),
          .packages="mice",
          .combine="c",
          .errorhandling = "remove",
          .verbose=TRUE) %dopar% {
            df <- boot_chunk[[i]]

            pre_imp <- df[,colnames(df) %in% imp_cols]
            pre_imp_id <- df[,!colnames(df) %in% imp_cols]

            imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 5)
            post_imp <- complete(imp_pars, action = "all")

            full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))
            names(full_imp) <- paste("imp", i+125, "_dat", names(full_imp), sep="")
            full_imp
          }

end_t <- Sys.time()
end_t - start_t


stopCluster(clust)


export_list(boot_imp,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imputed_boot_data/%s.csv")



#====================================#
#### data 151-175 ####
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)


boot_chunk <- boot_list_import[151:175]

start_t <- Sys.time()
boot_imp <-
  foreach(i=1:length(boot_chunk),
          .packages="mice",
          .combine="c",
          .errorhandling = "remove",
          .verbose=TRUE) %dopar% {
            df <- boot_chunk[[i]]

            pre_imp <- df[,colnames(df) %in% imp_cols]
            pre_imp_id <- df[,!colnames(df) %in% imp_cols]

            imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 5)
            post_imp <- complete(imp_pars, action = "all")

            full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))
            names(full_imp) <- paste("imp", i+150, "_dat", names(full_imp), sep="")
            full_imp
          }

end_t <- Sys.time()
end_t - start_t


stopCluster(clust)


export_list(boot_imp,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imputed_boot_data/%s.csv")




#====================================#
#### data 176-200 ####
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)


boot_chunk <- boot_list_import[176:200]

start_t <- Sys.time()
boot_imp <-
  foreach(i=1:length(boot_chunk),
          .packages="mice",
          .combine="c",
          .errorhandling = "remove",
          .verbose=TRUE) %dopar% {
            df <- boot_chunk[[i]]

            pre_imp <- df[,colnames(df) %in% imp_cols]
            pre_imp_id <- df[,!colnames(df) %in% imp_cols]

            imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 5)
            post_imp <- complete(imp_pars, action = "all")

            full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))
            names(full_imp) <- paste("imp", i+175, "_dat", names(full_imp), sep="")
            full_imp
          }

end_t <- Sys.time()
end_t - start_t


stopCluster(clust)


export_list(boot_imp,
            file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imputed_boot_data/%s.csv")



#=============================================================#
### DEBUGGING ####


# ISSUES: these do not impute: 5, 11, 22, 25, 26, 31, 40, 41, 45, 46
# ISSUES: these do not impute:







### check on individual data sets
df <- boot_list_import[[26]]

pre_imp <- df[,colnames(df) %in% imp_cols]
pre_imp_id <- df[,!colnames(df) %in% imp_cols]

# imp_pars <- futuremice(pre_imp, m=5, parallelseed=393857, n.core = 5)


test <- pre_imp %>% select(-wh_hrt_ever, -wh_contraceptives_ever)

# Rm each variable and rerun mice - does it solve the issue?
# alc_ever: maybe
# wh_hrt_ever: no
# wh_contraceptives_ever: worked for df 5, 11; no for 22
# imp_pars <- futuremice(test, m=5, parallelseed=393857, n.core = 5)


# 5,11: wh_contraceptives_ever
# 22: wh_contraceptives_ever
# 25: wh_contraceptives_ever, alc_ever
# 26: wh_hrt_ever, wh_contraceptives_ever
imp_pars <- mice(test, m=5, seed=10395)



post_imp <- complete(imp_pars, action = "all")

full_imp <- lapply(post_imp, function(imp_df) cbind(pre_imp_id, imp_df))




stopCluster(clust)




for (i in 1:50){
  print(i)
  print(summary(boot_list_import[[i]]$wh_contraceptives_ever))
}

