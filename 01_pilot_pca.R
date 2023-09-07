#==================================#
# PCA
pm_pca <- pm_treated %>% select(contains("ion"))
meta_pca <- prcomp(pm_pca, scale=FALSE)
summary(meta_pca)


# Visualize the percentage of variance explained
library(factoextra)
fviz_eig(meta_pca, ncp=20)
get_eig(meta_pca)



# meta_pca_summary <- summary(meta_pca)
# ggplot() +
#   geom_point(aes(1:46, meta_pca_summary$importance[3,]), color="blue") +
#   geom_line(aes(1:46, meta_pca_summary$importance[3,r]), color="blue") +
#   geom_point(aes(1:46, meta_pca_summary$importance[2,]), color="red") +
#   geom_line(aes(1:46, meta_pca_summary$importance[2,]), color="red") +
#   labs(x="principal component",
#        y="proportion of variance explained") +
#   theme_classic()


# biplot of the attributes
fviz_pca_var(meta_pca, col.var="cos2",
             gradient.cols=c("black", "orange", "green"),
             repel=TRUE)
fviz_pca(meta_pca)




# plot of data points by group for PC1 vs PC2
par(mfrow=c(2,2))
fviz_pca_ind(meta_pca,
             axes=c(1,2),
             label="none",
             habillage=pm_treated$gp,
             palette=c("blue", "red"),
             addEllipses = TRUE
             )

fviz_pca_ind(meta_pca,
             axes=c(2,3),
             label="none",
             habillage=pm_treated$gp,
             palette=c("blue", "red"),
             addEllipses = TRUE
             )


fviz_pca_ind(meta_pca,
             axes=c(1,3),
             label="none",
             habillage=pm_treated$gp,
             palette=c("blue", "red"),
             addEllipses = TRUE
)
par(mfrow=c(1,1))



#==================================#

# Correlate metabolites with risk factors
dups <- atp_pilot$participantkey[duplicated(atp_pilot$participantkey)]
dups_tube <- atp_pilot$tubeid[atp_pilot$participantkey %in% dups]






