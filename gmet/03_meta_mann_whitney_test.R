###################################################
###  MANN-WHITNEY U TEST
###  WITH FDR CORRECTION
###################################################

source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")




###########################################
# run Mann-Whitney U-test for each metabolite
wilcox_pval_func <- function(df) {
  ion_list <- df %>% select(contains("ion")) %>% names()
  wilcox_pvals <- data.frame(ion=ion_list,
                             pval=rep(0,length(ion_list)))

  # loop through list of ions
  for (i in 1:length(ion_list)){
    ion <- ion_list[i]
    case <- df[df$gp=="CASE",ion] %>% as_vector()
    cntl <- df[df$gp=="CNTL",ion] %>% as_vector()

    pval <- wilcox.test(case, cntl)$p.value

    wilcox_pvals[i,2] <- pval
    wilcox_pvals$p_adjust_fdr <- p.adjust(wilcox_pvals$pval, method="BH")
  }
  wilcox_pvals
}




# adjust for p-values with FDR
wilcox_pvals <- wilcox_pval_func(full_all)
# wilcox_pvals$p_adjust_fdr <- p.adjust(wilcox_pvals$pval, method="BH")

summary(wilcox_pvals$pval)
summary(wilcox_pvals$p_adjust_fdr)

wilcox_pvals[wilcox_pvals$p_adjust_fdr<=0.2,]

# 14 metabolites were determined to be significant in univariate analysis




###########################################################
## Test for normality before and after log-transformation

normality_pre <- normality_skewness(full_all)
# how many metabolites have p>0.05 for each test (indicating possible normality)
normality_nonorm_pre %>% count(ks_pval>.05)
normality_nonorm_pre %>% count(shapiro_pval>.05)
normality_nonorm_pre %>% count(jarque_pval>.05)
summary(normality_nonorm_pre$skew)

# none normal


full_nonorm_log <- full_all %>%
  mutate(across(starts_with("ion"), log))

full_nonorm_log2 <- full_nonorm_avg %>%
  mutate(across(starts_with("ion"), log2))

normality_nonorm_post <- normality_skewness(full_nonorm_log)
# how many metabolites have p>0.05 for each test (indicating possible)
normality_nonorm_post %>% count(ks_pval>.05)
normality_nonorm_post %>% count(shapiro_pval>.05)
normality_nonorm_post %>% count(jarque_pval>.05)
summary(normality_nonorm_post$skew)






################################################
## TRY T-TESTS

t_pval_func <- function(df) {
  ion_list <- df %>% select(contains("ion")) %>% names()
  t_pvals <- data.frame(ion=ion_list,
                        pval=rep(0,length(ion_list)))

  # loop through list of ions
  for (i in 1:length(ion_list)){
    ion <- ion_list[i]
    case <- df[df$gp=="CASE",ion] %>% as_vector()
    cntl <- df[df$gp=="CNTL",ion] %>% as_vector()

    pval <- t.test(case, cntl)$p.value

    t_pvals[i,2] <- pval
  }
  t_pvals$p_adjust_fdr <- p.adjust(t_pvals$pval, method="BH")
  t_pvals
}


t_df <- t_pval_func(full_nonorm_log)

summary(t_df)
