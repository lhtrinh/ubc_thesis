############################################################
# Multiple Imputation (MICE) of BCGP data
############################################################
library(mice)
library(VIM)

names(full_dat)




# remove diagnosis-specific variables
# (molecular and histological subtypes)
pre_imp <- full_dat %>%
  select(-c(age_at_diagnosis,
            er_status, pr_status, her2_status, hr_subtype,
            site_code, hist_subtype))

pre_imp_cln <- zap_labels(pre_imp)



# pre_imp_bcgp <- pre_imp_cln %>%
#   filter(cohort=="bcgp")
# dim(pre_imp_bcgp)
#
# pre_imp_atp <- pre_imp_cln %>%
#   filter(cohort=="atp")
# dim(pre_imp_atp)




#====================================#
# MICE FOR BCGP #



# Count % of NAs for each column
summary(pre_imp_cln)
apply(pre_imp_cln, 2, function(x) round(sum(is.na(x))/length(x), 3))

# plot of missing data
aggr(pre_imp_cln, col=c("navyblue", "red"), numbers=TRUE)


# impute bcgp data with 5 complete datasets and 5 iterations
mice_mod <- mice(data=pre_imp_cln, seed=123)

# check distributions of numeric variables
pre_imp_cln %>%
  group_by(gp) %>%
  summarize(across(where(is.numeric),
                   c(mean=function(x) mean(x, na.rm=TRUE),
                     sd=function(x) sd(x, na.rm=TRUE)))) %>%
  as.data.frame()

complete(mice_mod,1) %>%
  group_by(gp) %>%
  summarize(across(is.numeric,
                   c(mean=function(x) mean(x, na.rm=TRUE),
                     sd=function(x) sd(x, na.rm=TRUE)))) %>%
  as.data.frame()
complete(mice_bcgp,2) %>%
  group_bymice_mod %>%
  summarize(across(is.numeric,
                   c(mean=function(x) mean(x, na.rm=TRUE),
                     sd=function(x) sd(x, na.rm=TRUE)))) %>%
  as.data.frame()
complete(mice_mod,3)  %>%
  group_by(gp) %>%
  summarize(across(is.numeric,
                   c(mean=function(x) mean(x, na.rm=TRUE),
                     sd=function(x) sd(x, na.rm=TRUE)))) %>%
  as.data.frame()
complete(mice_mod,4) %>%
  group_by(gp) %>%
  summarize(across(is.numeric,
                   c(mean=function(x) mean(x, na.rm=TRUE),
                     sd=function(x) sd(x, na.rm=TRUE)))) %>%
  as.data.frame()
complete(mice_mod,5) %>%
  group_by(gp) %>%
  summarize(across(is.numeric,
                   c(mean=function(x) mean(x, na.rm=TRUE),
                     sd=function(x) sd(x, na.rm=TRUE)))) %>%
  as.data.frame()


# check distributions of categorical variables
pre_imp_cln %>%
  select(where(is.factor), -cohort) %>%
  pivot_longer(!gp) %>%
  group_by(gp, name) %>%
  count(value) %>%
  filter(!is.na(value)) %>%
  mutate(p=proportions(n)) %>%
  View()



imp_dat1 <- complete(mice_mod, 1)
imp_dat2 <- complete(mice_mod, 2)
imp_dat3 <- complete(mice_mod, 3)
imp_dat4 <- complete(mice_mod, 4)
imp_dat5 <- complete(mice_mod, 5)





