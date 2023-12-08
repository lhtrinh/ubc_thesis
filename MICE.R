#======================================================#
library(mice)
library(tidyverse)
# library(foreach)
# library(doParallel)
library(rio)


source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")




full_all <- import_survey() %>%
  full_join(import_meta_sc())




#======================================================#
# MICE imputation ####

# separate columns into data for imp and other columns
pre_imp <- full_all %>%
  select(all_of(imp_cols))

pre_imp_id <- full_all %>% select(-all_of(imp_cols))


# imputation (use futuremice function for parallel processing, speed up the process)
n_imp <- 10
imp <- futuremice(pre_imp,
             m=n_imp,
             maxit=30,
             parallelseed = 292920)

plot(imp, layout=c(5,5))


# record all imputed data sets
post_imp <- complete(imp, action="all")

# combine imputed data with other variables (studyid, match_id, cancer characteristics, and ions)
full_imp <- lapply(post_imp, function(dat) cbind(pre_imp_id, dat))


export_list(full_imp, file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/imp_list_%s.csv")



# #======================================================#
# # Bootstrap ####
# n_fold <- 10
#
# ## boot function ####
# boot_fn <- function(dat){
#   # Bootstrap data
#   match_ids <- unique(dat$match_id)
#   n_boot <- 200
#
#   boot_list <- foreach(i=1:n_boot, .packages = "dplyr", .verbose=TRUE) %dopar% {
#     # sample match IDs with replacement from list
#     boot_match <- sample(match_ids, length(match_ids), replace=T)
#     boot_fold <- sample(rep(sample(1:n_fold), length=length(boot_match)), replace=F)
#
#     # create bootstrapped data by selecting columns from full data
#     boot_df <- data.frame(match_id=boot_match, fold=boot_fold) %>% left_join(dat)
#     boot_df
#   }
# }
#
#
#
# # create cluster for parallel processing
# clust <- makeCluster(3)
# setDefaultCluster(clust)
#
# clusterEvalQ(clust, {
#   library(mice);
#   library(doParallel)
# }
# )
#
# clusterExport(clust, c("full_all", "boot_fn"), envir = .GlobalEnv)
#
#
# boot_list <- list_flatten(lapply(full_imp, boot_fn))
#
#
# stopCluster(clust)
#
#
#
# length(boot_list)
#
#
#
#
#
#
#
#
# library(rio)
# export_list(boot_list, file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/boot_list_%s.csv")
