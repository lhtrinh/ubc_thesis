source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/01_pilot_pretreatment.R")


# Paired t-tests with adjustment using the false discovery rate



#==================================#
# Compare metabolites between cases and controls
# Append group status to metabolite dataset
pm_sc_groups <- pm_treated %>%
  filter(type=="Sample")


# run t-tests for each metabolite
t_test_func <- function(df) {
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
  t_pvals
}




# adjust for p-values with FDR
pvals <- t_test_func(pm_sc_groups)
pvals$p_adjust_fdr <- p.adjust(pvals$pval, method="BH")

summary(pvals$p_adjust_fdr)
table(pvals$p_adjust_fdr<=0.05)
table(pvals$pval<=0.05)

pvals[pvals$p_adjust_fdr<=0.05,]



