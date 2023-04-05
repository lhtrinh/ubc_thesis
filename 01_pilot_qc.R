#=======================================#
## Check assay reproducibility
#=======================================#
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/01_Pilot.R")







#=======================================#
# Calculate coefficient of variation (CV)
#=======================================#

# Each column has a technical replicate
# Six samples also have blind replicates

# Calculate CV for each technical replicate
calc_cv <- function(dat, group_var){
  group_var <- enquo(group_var)
  replicate_cv <- dat %>%
    group_by(!!group_var) %>%
    summarise(across(contains("ion"),
                     function(x) sd(x)/mean(x)*100)) %>%
    ungroup() %>%
    summarise(across(contains("ion"),
                     median)) %>%
    as_vector()
  replicate_cv
}

summary(calc_cv(pilot_meta, id))
# median CVs range from 0.32% to 10.70%
# all good


# Calculate CV for the blind duplicates
head(pilot_ids)

dups_meta <- pilot_ids %>%
  count(participantkey) %>%
  filter(n>=2) %>%
  left_join(pilot_ids) %>%
  left_join(pilot_meta,
            by=c("tubeid"="sample_id", "participantkey"))

dups_meta %>% count(participantkey)


blind_dup_cvs <- calc_cv(dups_meta, participantkey)
summary(blind_dup_cvs)
table(blind_dup_cvs>=20)
proportions(table(blind_dup_cvs>=20))

# 13 ions have median CVs > 20%

ions_to_discard <- names(which(blind_dup_cvs >= 20))










##########################################
# # Compare metabolites between cases and controls
#
#
# # histogram of distribution for a few metabolites
# ions_df <- pm_dedup %>% select(starts_with("ion"))
#
# nrow <- 3
# ncol <- 8
# dim <- nrow*ncol
# par(mfrow=c(nrow, ncol))
# for (i in 1:dim){
#   x <- ions_df[,i]
#   xval <- as_vector(as.matrix(x))
#   xname <- colnames(x)
#   boxplot(xval, xlab="", main=xname)
# }





#=======================================#
# Calculate ICC
#=======================================#