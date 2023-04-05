#############################################
# preprocessing
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/01_pilot_qc.R")



#===========================================#
# remove metabolites with CV > 20%
ions_to_discard <- names(which(blind_dup_cvs > 20))
pm_cv <- pilot_meta %>% select(!all_of(ions_to_discard))





#===========================================#
# Use average concentrations for each run
pm_avg <- pm_cv %>%
  group_by(type, sample_id, id, participantkey) %>%
  summarise(across(contains("ion"),
                   function(x) round(mean(x)))) %>%
  ungroup() %>%
  left_join(atp_match_ids)




#===========================================#
# Randomly select samples from duplicate pairs to be included in analysis
set.seed(123)
dups_to_discard <- dups_meta %>%
  group_by(participantkey) %>%
  filter(tubeid==tubeid[sample(length(tubeid), 1)]) %>%
  select(participantkey, tubeid)



# remove duplicates
pm_dedup <- pm_avg %>%
  filter(!sample_id %in% dups_to_discard$tubeid)





#===========================================#
# log-transform
pm_log <- pm_dedup %>%
  mutate(across(contains("ion"),
                log))



# test for normality before and after log-transforming
metabolite_norm <- function(df) {
  norm_df <- data.frame(metabolite=rep(0, ncol(df)),
                        shapiro_pval=rep(0,ncol(df)),
                        jarque_pval=rep(0,ncol(df)))
  ion_df <- df %>% select(starts_with("ion"))
  for (i in 1:ncol(ion_df)){
    x <- ion_df[,i] %>% as_vector()
    norm_df[i,1] <- colnames(ion_df)[i]
    norm_df[i,2] <- shapiro.test(x)$p.value
    norm_df[i,3] <- jarque.test(x)$p.value
  }
  norm_df
}

# before log-transforming
norm <- metabolite_norm(pm_dedup %>% filter(type=="Sample"))
proportions(table(norm$shapiro_pval>0.05))
proportions(table(norm$jarque_pval>0.05))


# after log-transforming
norm_log <- metabolite_norm(pm_log %>% filter(type=="Sample"))
proportions(table(norm_log$shapiro_pval>0.05))
proportions(table(norm_log$jarque_pval>0.05))





#===========================================#
# center each pair by within-pair mean
pm_pair_centered <- pm_log %>%
  filter(type=="Sample") %>%
  arrange(match_id) %>%
  group_by(match_id) %>%
  mutate(across(contains("ion"),
                function(x) x-mean(x))) %>%
  ungroup()





#===========================================#
# pareto scaling
# (center by mean and divide by square root of standard deviation)
pm_mean_centered <- pm_pair_centered %>%
  mutate(across(contains("ion"),
                function(x) x-mean(x)))

pm_sc <- pm_mean_centered %>%
  mutate(across(contains("ion"),
                function(x) x/sqrt(sd(x))))






#===========================================#
# Assign pre-treated data name
pm_treated <- pm_sc




#===========================================#
# Merge with contextual data
pm_treated2 <- pm_treated %>%
  mutate(studyid=as.character(participantkey)) %>%
  select(-type, -sample_id, -id, -participantkey)


dat2 <- dat %>% select(-c(participantkey, match_id))

pm_dat <- pm_treated2 %>%
  left_join(dat2,
            by=c("studyid", "gp"))



