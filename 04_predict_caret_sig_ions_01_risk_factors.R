library(doParallel)
library(rio)

# # If error in parallel processing, use this function to clean up
# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }





source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")




#=============================#
## Assign modeling and holdout set ####



full_dat <- import_survey()
full_sc <- import_meta_sc()

# full_all <- full_dat %>% full_join(full_sc)


full_imp <- import_imp()



all_pvals <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/all_pvals.csv")


ion_cols <- all_pvals$metabolite
# ion_cols <- all_pvals$metabolite[all_pvals$q_fdr<=0.1]
# ion_cols <- all_pvals$metabolite[all_pvals$q_fdr<=0.2]





#=============================#
## Set up parallel processing ##
clust <- makeCluster(3)
setDefaultCluster(clust)

registerDoParallel(clust)

clusterEvalQ(clust, {
  library(pROC);
  library(caret);
  library(randomForest);
  library(kernlab);
  library(tidyverse);
  library(doParallel)
}
)

clusterExport(clust, c("full_imp"), envir = .GlobalEnv)




mods <- foreach(i=1:10, .errorhandling = "stop", .combine="rbind", .verbose=TRUE)%:%
  foreach(j=1:10, .errorhandling="stop", .combine="rbind", .verbose=TRUE) %dopar% {

    my_df <- full_imp[[i]] %>%
      select(gp, studyid, match_id, all_of(c(match_cols, context_cols, ion_cols)))


    # Randomly sample by match ID
    matchIds <- unique(my_df$match_id)

    # set.seed(2931)
    trainMatchId <- sample(matchIds, length(matchIds)*.80, replace = FALSE)


    # Split data training (90%) and validation (10%) sets
    model_dat <- my_df %>% filter(match_id %in% trainMatchId) %>% select(-studyid, -match_id)
    holdout_dat <- my_df %>% filter(!match_id %in% trainMatchId) %>% select(-studyid, -match_id)



    #=============================#
    # Create fold indices for 10-fold CV
    ## randomly split match IDs into 10 equal-sized folds
    n_fold <- 10
    # set.seed(39671)
    foldInd <- sample(rep(sample(1:n_fold), length=length(trainMatchId)), replace=FALSE)
    # table(foldInd)

    ## divide indices in training data by folds
    ctrlIndexOut <- split(1:nrow(model_dat), foldInd)

    ## Specify train control to be 10-fold CV
    ctrl <- trainControl(
      method="cv",
      number=n_fold,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      indexOut = ctrlIndexOut,
      savePredictions = TRUE
    )
    #=============================#


#
#     #=============================================================#
#     ## Lasso ####
#     # Create grid of hyper parameters
#     lasso_grid <- expand.grid(alpha=1,
#                               lambda=seq(.001,1,length=100))
#
#     lasso_train <- train(gp~.,
#                          data=model_dat,
#                          method="glmnet",
#                          trControl=ctrl,
#                          tuneGrid=lasso_grid,
#                          metric="ROC")
#
#     # lasso_train
#
#     # Predict on holdout test set
#     lasso_pred <- predict(lasso_train,
#                           ncomp=lasso_train$bestTune,
#                           newdata=holdout_dat,
#                           type="prob")
#     lasso_pred_prob <- lasso_pred$CASE
#     lasso_pred_label <- as.factor(ifelse(lasso_pred_prob>=0.5,"CASE", "CNTL"))
#     true_label <- holdout_dat$gp
#
#
#     lasso_roc <- roc(true_label, lasso_pred_prob)
#
#     lasso_mods <- data.frame(
#       dfver=i,
#       iter=j,
#       method="lasso",
#       auc=lasso_roc$auc,
#       sen=sensitivity(lasso_pred_label, true_label),
#       spe=specificity(lasso_pred_label, true_label),
#       ppv=posPredValue(factor(lasso_pred_label), factor(true_label)))
#
#
#
#
#
#     #=============================================================#
#     ## Partial least squares ####
#     plsda_train <- train(gp~.,
#                          data=model_dat,
#                          method="pls",
#                          scale=TRUE,
#                          trControl=ctrl,
#                          tuneLength=10,
#                          # tuneGrid=rf_grid,
#                          metric="ROC")
#
#     # plsda_train
#
#     # Predict on holdout test set
#     pls_pred <- predict(plsda_train, ncomp=plsda_train$bestTune, newdata=holdout_dat, type="prob")
#     pls_pred_prob <- pls_pred$CASE
#     pls_pred_label <- as.factor(ifelse(pls_pred_prob>=0.5,"CASE", "CNTL"))
#     true_label <- holdout_dat$gp
#
#
#     pls_roc <- roc(true_label, pls_pred_prob)
#
#     pls_mods <- data.frame(
#       dfver=i,
#       iter=j,
#       method="pls-da",
#       auc=pls_roc$auc,
#       sen=sensitivity(pls_pred_label, true_label),
#       spe=specificity(pls_pred_label, true_label),
#       ppv=posPredValue(factor(pls_pred_label), factor(true_label)))
#
#


 #=============================================================#
    # # Random forest ####
    # rf_grid <- expand.grid(.mtry=1:15)
    #
    # rf_train <- train(gp~.,
    #                   data=model_dat,
    #                   method="rf",
    #                   trControl=ctrl,
    #                   tuneGrid=rf_grid,
    #                   metric="ROC")
    #
    # # Predict on holdout test set
    # rf_test <- predict(rf_train, mtry=rf_train$bestTune, newdata=holdout_dat, type="prob")
    # rf_pred_prob <- rf_test$CASE
    # rf_pred_label <- as.factor(ifelse(rf_pred_prob>=0.5,"CASE", "CNTL"))
    # true_label <- holdout_dat$gp
    #
    #
    # rf_roc <- roc(true_label, rf_pred_prob)
    #
    # rf_mods <- data.frame(
    #   dfver=i,
    #   iter=j,
    #   method="rf",
    #   auc=rf_roc$auc,
    #   sen=sensitivity(rf_pred_label, true_label),
    #   spe=specificity(rf_pred_label, true_label),
    #   ppv=posPredValue(factor(rf_pred_label), factor(true_label))
    # )


#
#
#
#     #===============================================================#
#     ## SVM linear kernel ####
#     svm_train <- train(gp~.,
#                        data=model_dat,
#                        method="svmLinear",
#                        trControl=ctrl,
#                        tuneLength=10)
#
#     svm_train$bestTune
#
#
#     # Predict on holdout test set
#     svm_test <- predict(svm_train, newdata=holdout_dat, type="prob")
#     svm_pred_prob <- svm_test$CASE
#     svm_pred_label <- as.factor(ifelse(svm_pred_prob>=0.5,"CASE", "CNTL"))
#     true_label <- holdout_dat$gp
#
#     table(svm_pred_label, true_label)
#
#
#     svm_roc <- roc(true_label, svm_pred_prob)
#
#     svm_mods <- data.frame(
#       dfver=i,
#       iter=j,
#       method="svm_linear",
#       auc=svm_roc$auc,
#       sen=sensitivity(svm_pred_label, true_label),
#       spe=specificity(svm_pred_label, true_label),
#       ppv=posPredValue(factor(svm_pred_label), factor(true_label))
#     )
#




    #===============================================================#
    ## SVM radial kernel ####
    svmr_train <- train(gp~.,
                        data=model_dat,
                        method="svmRadial",
                        trControl=ctrl,
                        tuneLength=10)

    svmr_train$bestTune


    # Predict on holdout test set
    svmr_test <- predict(svmr_train, newdata=holdout_dat, type="prob")
    svmr_pred_prob <- svmr_test$CASE
    svmr_pred_label <- as.factor(ifelse(svmr_pred_prob>=0.5,"CASE", "CNTL"))
    true_label <- holdout_dat$gp

    table(svmr_pred_label, true_label)


    svmr_roc <- roc(true_label, svmr_pred_prob)

    svmr_mods <- data.frame(
      dfver=i,
      iter=j,
      method="svm_radial",
      auc=svmr_roc$auc,
      sen=sensitivity(svmr_pred_label, true_label),
      spe=specificity(svmr_pred_label, true_label),
      ppv=posPredValue(factor(svmr_pred_label), factor(true_label))
    )



    # all_mods <- rbind(
      # lasso_mods,
      # pls_mods)
      # ,
                      # rf_mods
    # ,
                      # svm_mods
    # ,
                      svmr_mods
    # )
}



# Clear up all parallel clusters
stopCluster(clust)
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()



write_csv(mods, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/holdout_svmr_all_risk_100times.csv")

print("Done with running all metabolites with risk factors.")
