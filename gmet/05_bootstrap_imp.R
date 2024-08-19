# Bootstrap data ####

# create cluster for parallel processing
# clust <- makeCluster(4)
# setDefaultCluster(clust)
#
#
#
# clusterEvalQ(clust, {
#   library(mice);
#   library(doParallel)
# }
# )
#
# clusterExport(clust, c("sig_all", "imp_fn"), envir = .GlobalEnv)
#
#
#
# start_time <- Sys.time()
# match_ids <- unique(sig_all$match_id)
# n_boot <- 200
#
# boot_list <- foreach(i=1:n_boot, .packages = "dplyr", .verbose=TRUE) %dopar% {
#   boot_match <- sample(match_ids, length(match_ids), replace=T)
#   boot_df <- data.frame(match_id=boot_match) %>% left_join(sig_all)
# }
#
# end_time <- Sys.time()
# end_time - start_time
#
# stopCluster(clust)
#
#
#
# print("BOOTSTRAPPING DONE.")


export_list(boot_list, file="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Bootstrapped_data/boot_list_%s.csv")