#======================================================#
library(mice)
library(tidyverse)


source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")




full_all



# summary(dat_pre_imp)

# remove information not for imputing
pre_imp <- full_all %>%
  select(all_of(imp_cols))

pre_imp_id <- full_all %>% select(-all_of(imp_cols))


#======================================================#

# imputation

init <- futuremice(pre_imp,
             m=5,
             maxit=5,
             parallelseed = 292920)


# record all imputed data sets
post_imp <- complete(init, action="all")

full_imp <- lapply(post_imp, function(df_n) cbind(pre_imp_id, df_n))




