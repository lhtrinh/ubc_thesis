##########################################
## Check assay reproducibility
##########################################




# select metabolomics data of samples with blind replicates
dups_meta <- pilot_ids %>%
  count(participantkey) %>%
  filter(n>=2) %>%
  left_join(pilot_ids) %>%
  left_join(pilot_meta,
            by=c("tubeid"="sampleid"))

dups_meta[,1:10]


# calculate CV per metabolite per pair
# calculate min, median, and max CV for each metabolite
dups_cv <- dups_meta %>%
  select(-c(n, pilot, arm, tubeid)) %>%
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

pilot_postqc <- pilot_meta %>%
  filter(!sampleid %in% dups_to_discard$tubeid)
dim(pilot_postqc)






##########################################
# Compare metabolites between cases and controls



nrow <- 7
ncol <- 5
dim <- nrow*ncol
par(mfrow=c(nrow, ncol))
for (i in 1:dim){
  x <- pilot_postqc[,i+1]
  xval <- as_vector(as.matrix(x))
  xname <- colnames(x)
  hist(xval, xlab="", main=xname)
}
par(mfrow=c(1,1))


# test for normality
norm_df <- data.frame(metabolite=rep(0, 1120), norm_pval=rep(0,1120))
for (i in 1:(ncol(pilot_meta)-2)){
  norm_df[i,1] <- colnames(pilot_meta)[i+2]
  norm_df[i,2] <- shapiro.test(pilot_meta[,i+2])$p.value
}

table(norm_df$norm_pval>0.05)


# log-transform and t-test OR leave alone and wilcoxon test