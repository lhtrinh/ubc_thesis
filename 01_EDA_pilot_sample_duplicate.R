##########################################
## Check assay reproducibility
##########################################
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/01_Pilot.R")





# select metabolomics data of samples with blind replicates
dups_meta <- pilot_ids %>%
  count(participantkey) %>%
  filter(n>=2) %>%
  left_join(pilot_ids) %>%
  left_join(pilot_meta) %>%
  select(-n)

dups_meta[,1:10]


# calculate CV per metabolite per pair
# calculate min, median, and max CV for each metabolite
dups_cv <- dups_meta %>%
  select(-c(gp, tubeid)) %>%
  group_by(participantkey) %>%
  summarise(across(.cols=everything(),
                   .fns=function(x) sd(x)/mean(x)*100)) %>%
  ungroup() %>%
  summarise(across(-participantkey,
                   median)) %>%
  as_vector()

summary(dups_cv)
table(dups_cv>20)
proportions(table(dups_cv>20))
proportions(table(dups_cv>10))







##########################################
# Randomly select samples from duplicate pairs to be included in analysis
dups_meta$tubeid
dups_meta$participantkey

set.seed(123)
dups_to_discard <- dups_meta %>%
  group_by(participantkey) %>%
  filter(tubeid==tubeid[sample(length(tubeid), 1)]) %>%
  select(participantkey, tubeid)

# remove duplicates
# scale to mean 0 and sd 1
pilot_postqc <- pilot_meta %>%
  filter(!tubeid %in% dups_to_discard$tubeid)






#====================#
# !!! TO DO: center by within-pair mean before scaling
#====================#











##########################################
# Compare metabolites between cases and controls


# histogram of distribution for a few metabolites
ions_df <- pilot_postqc %>% select(starts_with("ion"))

nrow <- 6
ncol <- 6
dim <- nrow*ncol
par(mfrow=c(nrow, ncol))
for (i in 1:dim){
  x <- ions_df[,i]
  xval <- as_vector(as.matrix(x))
  xname <- colnames(x)
  hist(xval, xlab="", main=xname)
}
par(mfrow=c(1,1))



# test for normality
metabolite_norm <- function(df) {
  norm_df <- data.frame(metabolite=rep(0, ncol(df)),
                        norm_pval=rep(0,ncol(df)))
  for (i in 1:ncol(df)){
    norm_df[i,1] <- colnames(df)[i]
    norm_df[i,2] <- shapiro.test(df[,i])$p.value
  }

  out <- prop.table(table(norm_df$norm_pval>0.05))
  names(out) <- c("not_normal", "normal")
  out
}


ions_df <- pilot_postqc %>% select(starts_with("ion"))
metabolite_norm(ions_df)



# 3/4 of the metabolites' distributions are significantly different from a normal distribution;
# 1/4 are not

pilot_postqc[1:6, 1:10]






# check for missing values in pilot data
test <- pilot_meta %>%
  summarise(across(contains("ion"),
                   function(x) sum(is.na(x)))) %>%
  as_vector()
summary(test) # no missing values


#############################################
# log-transform and t-test OR leave alone and wilcoxon test

# log-transform
pilot_logmeta <- pilot_postqc %>%
  mutate(across(contains("ion"),
                log))

pilot_logmeta[,1:10] %>% head()


# check for normality
df <- pilot_logmeta %>% select(starts_with("ion"))
metabolite_norm(df)
# 57% log-transformed metabolites are normally distributed
# 43% log-transformed metabolites are not
# t-tests not a good idea



## wilcoxon test


# run wilcoxon signed rank test for each metabolite
ion_list <- pilot_postqc %>% select(starts_with("ion")) %>% names()
wilcox_pvals <- data.frame(ion=ion_list,
                             pval=rep(0,length(ion_list)))

for (i in 1:length(ion_list)){
  ion <- ion_list[i]
  cases <- pilot_postqc[,ion][pilot_postqc$gp=="CASE"]
  controls <- pilot_postqc[,ion][pilot_postqc$gp=="CNTL"]

  pval <- t.test(cases, controls)$p.value
  t_pvals[i,2] <- pval
}


# adjust for p-values with FDR
wilcox_adjust_fdr <- p.adjust(wilcox_pvals$pval, method="BH")
summary(wilcox_adjust_fdr)
table(wilcox_adjust_fdr<=0.05)
table(wilcox_pvals$pval<=0.05)
