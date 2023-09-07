###################################################
###  PCA ON NON-NORMALIZED DATA
###################################################
library(factoextra)


dat <- full_nonorm_centered %>%
  inner_join(full_dat[,c("studyid", "gp")]) %>%
  mutate(plate11_yn=ifelse(plate==11, 1, 0))

dat_for_pca <- dat %>% select(starts_with("ion"))

pca_nonorm <- prcomp(subset(dat, select=-gp), scale=FALSE)


fviz_pca_ind(pca_nonorm,
             axes=c(1,2),
             label="none",
             habillage = dat$gp,
             palette=c("blue", "red"),
             addEllipses = TRUE)


fviz_pca_ind(pca_nonorm,
             axes=c(1,3),
             label="none",
             habillage =  dat$gp,
             palette=c("blue", "red"),
             addEllipses = TRUE)


fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  dat$gp,
             palette=c("blue", "red"),
             addEllipses = TRUE)


# plot by plate
fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  dat$plate,
             addEllipses = TRUE)

fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  dat$plate11_yn,
             addEllipses = TRUE)


# plot by cohort
fviz_pca_ind(pca_nonorm,
             label="none",
             habillage =  dat$cohort,
             addEllipses = TRUE)


# plot percentage of variance explained
fviz_eig(pca_nonorm, addlabels = TRUE)


# 10 components explain all variance in the data
# there is batch effect when showing plate 11 vs other plates
# maybe look into batch correction methods