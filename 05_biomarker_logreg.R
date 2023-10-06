#=============================================================#
# BIOMARKER DISCOVERY WITH LOGISTIC REGRESSION ####
#=============================================================#


lr_bio_fn <- function(dat){
  ion_cols <- colnames(dat)[grepl("^ion", colnames(dat))]
  length(ion_cols)

  logreg_pvals <- data.frame(
    metabolite=ion_cols,
    pval=rep(0, length(ion_cols)))


  for (i in 1:length(ion_cols)){
    ion <- ion_cols[i]

    cols_ <- c("gp", "cohort", "collect_age", "collect_yr_cat", "menopause_stt", ion)
    lr_dat <- dat[,colnames(dat) %in% cols_]

    # fit logistic regression model
    lr_mod <- glm(gp~., data=lr_dat, family=binomial(link="logit"))

    # extract p-value for metabolite
    lr_coef <- summary(lr_mod)$coefficients
    pval <- lr_coef[rownames(lr_coef)==ion, 4]

    logreg_pvals[i,1] <- ion
    logreg_pvals[i,2] <- pval
  }


  # adjust for p-values with FDR
  logreg_pvals$p_adjust_fdr <- p.adjust(logreg_pvals$pval, method="BH")

  logreg_pvals
}





#=======================================================================#
## Full data set ####
full_all <- full_centered %>%
  full_join(full_dat) %>%
  mutate(across(c(gp, collect_yr_cat, menopause_stt, cohort), as.factor))


all_pvals <- lr_bio_fn(full_all)


summary(all_pvals$pval)
summary(all_pvals$p_adjust_fdr)

table(all_pvals$p_adjust_fdr<=0.2)





#=======================================================================#
## Restrict to non-HRT users ####
full_all %>% count(wh_hrt_ever)


full_nohrt <- full_all %>%
  filter(wh_hrt_ever==0)


nohrt_pvals <- lr_bio_fn(full_nohrt)

summary(nohrt_pvals$pval)
summary(nohrt_pvals$p_adjust_fdr)

table(nohrt_pvals$pval<=.05)
table(nohrt_pvals$p_adjust_fdr<=0.2)





#=======================================================================#
## Restrict to ER+ ####
hrpos_match <- full_all$match_id[which(full_all$hr_subtype=="hr_positive")]

full_hrpos <- full_all %>%
  filter(match_id %in% hrpos_match)


hrpos_pvals <- lr_bio_fn(full_hrpos)

summary(hrpos_pvals$pval)
summary(hrpos_pvals$p_adjust_fdr)

table(hrpos_pvals$pval<=.05)
table(hrpos_pvals$p_adjust_fdr<=0.2)





#=======================================================================#
## Restrict to postmenopausal ####

full_postmnp <- full_all %>%
  filter(menopause_stt==1)


postmnp_pvals <- lr_bio_fn(full_postmnp)

summary(postmnp_pvals$pval)
summary(postmnp_pvals$p_adjust_fdr)

table(postmnp_pvals$pval<=.05)
table(postmnp_pvals$p_adjust_fdr<=0.2)





#=======================================================================#
## Restrict to ductal ####

full_postmnp <- full_all %>%
  filter(menopause_stt==1)


postmnp_pvals <- lr_bio_fn(full_postmnp)

summary(postmnp_pvals$pval)
summary(postmnp_pvals$p_adjust_fdr)

table(postmnp_pvals$pval<=.05)
table(postmnp_pvals$p_adjust_fdr<=0.2)
