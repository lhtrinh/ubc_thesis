###################################################
###  PCA ON NON-NORMALIZED DATA
###################################################
library(factoextra)

#non-normalized
dat <- full_nonorm_sampleid %>%
  filter(sampletype=="Sample") %>%
  mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))


####
# # QC samples and dat
# dat <- full_nonorm_sampleid %>%
#   filter(cohort %in% c("ATP", "BCGP", "pSS")) %>%
#   mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))
####

#normalized
dat <- full_ion_norm %>%
  # inner_join(full_dat[,c("studyid", "gp")]) %>%
  filter(sampletype=="Sample") %>%
  mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))


# normalized and pareto scaled
# dat <- full_nonorm_sc %>%
#   # inner_join(full_dat[,c("studyid", "gp")]) %>%
#   filter(sampletype=="Sample") %>%
#   mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))


# normalized, pareto scaled, and mean-centered
dat <- full_centered %>%
  # inner_join(full_dat[,c("studyid", "gp")]) %>%
  filter(sampletype=="Sample") %>%
  mutate(plate11_yn=ifelse(plate==11, "Plate 11", "Other plates"))





####################################################
# isolate data for PCA
dat_for_pca <- dat %>% select(starts_with("ion"))

pca_nonorm <- prcomp(dat_for_pca)



# # plot percentage of variance explained
# fviz_eig(pca_nonorm, addlabels = TRUE)

# plot by plate 11 vs other plates
fviz_pca_ind(pca_nonorm,
             axes = c(1,2),
             label="none",
             alpha=.5,
             habillage =  dat$plate11_yn,
             # palette=c("blue", "red"),
             addEllipses = TRUE) #+ xlim(-1*10^4, 5*10^3)


# plot by plate
fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  dat$plate,
             alpha=.5,
             addEllipses = TRUE) #+ xlim(-2*10^4, 5*10^3)


# plot by cohort
fviz_pca_ind(pca_nonorm,
             label="none",
             alpha=.5,
             habillage = dat$cohort,
             # palette=c("blue", "red"),
             addEllipses = FALSE) #+ xlim(-10^7, 5*10^6)



# plot components against one another

# fviz_pca_ind(pca_nonorm,
#              axes=c(1,2),
#              label="none",
#              habillage = dat$gp,
#              palette=c("blue", "red"),
#              addEllipses = TRUE)
#
#
# fviz_pca_ind(pca_nonorm,
#              axes=c(1,3),
#              label="none",
#              habillage =  dat$gp,
#              palette=c("blue", "red"),
#              addEllipses = TRUE)
#
#
# fviz_pca_ind(pca_nonorm,
#              label="none",
#              habillage =  dat$gp,
#              palette=c("blue", "red"),
#              addEllipses = TRUE)











# 10 components explain all variance in the data
# there is batch effect when showing plate 11 vs other plates
# maybe look into batch correction methods