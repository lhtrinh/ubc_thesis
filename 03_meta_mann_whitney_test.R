###################################################
###  MANN-WHITNEY U TEST
###  WITH FDR CORRECTION
###################################################

source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")





# run t-tests for each metabolite
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
  }
  wilcox_pvals
}




# adjust for p-values with FDR
pvals <- wilcox_pval_func(full_nonorm_centered)
pvals$p_adjust_fdr <- p.adjust(pvals$pval, method="BH")

summary(pvals$p_adjust_fdr)
table(pvals$p_adjust_fdr<=0.05)
table(pvals$pval<=0.05)

pvals[pvals$p_adjust_fdr<=0.05,]

# 14 metabolites were determined to be significant in univariate analysis