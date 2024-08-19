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
full_all <- full_dat %>% full_join(full_sc)



all_pvals <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/all_pvals.csv")


# ion_cols <- all_pvals$metabolite
ion_cols <- all_pvals %>% filter(q_fdr<=0.1) %>% pull(metabolite)
# ion_cols <- all_pvals$metabolite[all_pvals$q_fdr<=0.2]

n_ion <- length(ion_cols)
n_ion



# postmenopausal
# my_df <- full_all %>%
#   filter(menopause_stt==1) %>%
#   select(gp, studyid, match_id, all_of(ion_cols))

# ductal
my_df <- full_all %>%
  filter(match_id %in% full_all$match_id[full_all$hist_subtype=="ductal"]) %>%
  select(gp, studyid, match_id, all_of(ion_cols))
#
# # er/pr
# my_df <- full_all %>%
#   filter(match_id %in% full_all$match_id[full_all$er_status=="positive" | full_all$pr_status=="positive"])




names(my_df)
dim(my_df)


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
  library(rminer)
  library(doParallel)
}
)
clusterExport(clust, c("my_df"), envir = .GlobalEnv)




mods <- foreach(i=1:100, .errorhandling = "stop", .combine="rbind", .verbose=TRUE) %dopar% {
  # Randomly sample by match ID
  matchIds <- unique(my_df$match_id)

  # set.seed(2931)
  trainMatchId <- sample(matchIds, length(matchIds)*.80, replace = FALSE)


  # Split data training (90%) and validation (10%) sets
  model_dat <- my_df %>% filter(match_id %in% trainMatchId) %>% select(-studyid, -match_id) %>% as.data.frame
  holdout_dat <- my_df %>% filter(!match_id %in% trainMatchId) %>% select(-studyid, -match_id) %>% as.data.frame



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



  #=============================================================#
  ## Lasso ####
  # Create grid of hyper parameters
  lasso_grid <- expand.grid(alpha=1,
                            lambda=seq(.001,1,length=100))

  # train
  lasso_train <- train(gp~.,
                       data=model_dat,
                       method="glmnet",
                       trControl=ctrl,
                       tuneGrid=lasso_grid,
                       metric="ROC")



  # Predict on holdout test set
  lasso_pred <- predict(lasso_train, ncomp=lasso_train$bestTune, newdata=holdout_dat, type="prob")
  lasso_pred_prob <- lasso_pred$CASE
  lasso_pred_label <- as.factor(ifelse(lasso_pred_prob>=0.5,"CASE", "CNTL"))
  true_label <- holdout_dat$gp

  # get prediction metrics
  lasso_auc <- roc(true_label, lasso_pred_prob)$auc
  lasso_sen <- sensitivity(lasso_pred_label, true_label)
  lasso_spe <- specificity(lasso_pred_label, true_label)
  lasso_ppv <- posPredValue(factor(lasso_pred_label), factor(true_label))
  lasso_fscore <- 2*lasso_ppv*lasso_sen/(lasso_ppv+lasso_sen)

  # get variable importance
  lasso_imp <- varImp(lasso_train, scale=FALSE)$importance

  # record in data table
  lasso_mods <- data.frame(iter=rep(i, n_ion),
                           method=rep("lasso", n_ion),
                           auc=rep(lasso_auc, n_ion),
                           sen=rep(lasso_sen, n_ion),
                           spe=rep(lasso_spe, n_ion),
                           ppv=rep(lasso_ppv, n_ion),
                           fscore=rep(lasso_fscore, n_ion),
                           var=rownames(lasso_imp),
                           varImp=lasso_imp[[1]])






  #=============================================================#
  ## Partial least squares ####
  plsda_train <- train(gp~.,
                       data=model_dat,
                       method="pls",
                       scale=TRUE,
                       trControl=ctrl,
                       tuneLength=10,
                       metric="ROC")


  # Predict on holdout test set
  plsda_pred <- predict(plsda_train, ncomp=plsda_train$bestTune, newdata=holdout_dat, type="prob")
  plsda_pred_prob <- plsda_pred$CASE
  plsda_pred_label <- as.factor(ifelse(plsda_pred_prob>=0.5,"CASE", "CNTL"))
  true_label <- holdout_dat$gp


  # get prediction metrics
  plsda_auc <- roc(true_label, plsda_pred_prob)$auc
  plsda_sen <- sensitivity(plsda_pred_label, true_label)
  plsda_spe <- specificity(plsda_pred_label, true_label)
  plsda_ppv <- posPredValue(factor(plsda_pred_label), factor(true_label))
  plsda_fscore <- 2*plsda_ppv*plsda_sen/(plsda_ppv+plsda_sen)

  # get variable importance
  plsda_imp <- varImp(plsda_train, scale=FALSE)$importance

  # record in data table
  plsda_mods <- data.frame(iter=rep(i, n_ion),
                           method=rep("plsda", n_ion),
                           auc=rep(plsda_auc, n_ion),
                           sen=rep(plsda_sen, n_ion),
                           spe=rep(plsda_spe, n_ion),
                           ppv=rep(plsda_ppv, n_ion),
                           fscore=rep(plsda_fscore, n_ion),
                           var=rownames(plsda_imp),
                           varImp=plsda_imp[[1]])




  #=============================================================#
  ## Random forest ####
  rf_grid <- expand.grid(.mtry=1:15)

  rf_train <- train(gp~.,
                    data=model_dat,
                    method="rf",
                    trControl=ctrl,
                    tuneGrid=rf_grid,
                    metric="ROC")


  ## Predict on holdout test set
  rf_test <- predict(rf_train, mtry=rf_train$bestTune, newdata=holdout_dat, type="prob")
  rf_pred_prob <- rf_test$CASE
  rf_pred_label <- as.factor(ifelse(rf_pred_prob>=0.5,"CASE", "CNTL"))
  true_label <- holdout_dat$gp


  # get prediction metrics
  rf_auc <- roc(true_label, rf_pred_prob)$auc
  rf_sen <- sensitivity(rf_pred_label, true_label)
  rf_spe <- specificity(rf_pred_label, true_label)
  rf_ppv <- posPredValue(factor(rf_pred_label), factor(true_label))
  rf_fscore <- 2*rf_ppv*rf_sen/(rf_ppv+rf_sen)

  # get variable importance
  rf_imp <- varImp(rf_train, scale=FALSE)$importance

  # record in data table
  rf_mods <- data.frame(iter=rep(i, n_ion),
                        method=rep("rf", n_ion),
                        auc=rep(rf_auc, n_ion),
                        sen=rep(rf_sen, n_ion),
                        spe=rep(rf_spe, n_ion),
                        ppv=rep(rf_ppv, n_ion),
                        fscore=rep(rf_fscore, n_ion),
                        var=rownames(rf_imp),
                        varImp=rf_imp[[1]])



    #===============================================================#
    ## SVM linear kernel ####

    # 10-fold cross-validation
    svmlin_train <-   train(gp~.,
                            data=model_dat,
                            method="svmLinear",
                            trControl=ctrl,
                            tuneLength=10,
                            metric="ROC")

    # fit final model using rminer
    svmlin_final_fit <- rminer::fit(
      gp~.,
      data=model_dat,
      task="prob",
      model="ksvm",
      kernel="vanilladot",
      C=svmlin_train$bestTune$C
    )


    # Predict on holdout test set
    svmlin_test <- predict(svmlin_final_fit, newdata=holdout_dat, type="prob")
    svmlin_pred_prob <- svmlin_test[,colnames(svmlin_test)=="CASE"]
    svmlin_pred_label <- as.factor(ifelse(svmlin_pred_prob>=0.5,"CASE", "CNTL"))
    true_label <- holdout_dat$gp


    # get prediction metrics
    svmlin_auc <- roc(true_label, svmlin_pred_prob)$auc
    svmlin_sen <- sensitivity(svmlin_pred_label, true_label)
    svmlin_spe <- specificity(svmlin_pred_label, true_label)
    svmlin_ppv <- posPredValue(factor(svmlin_pred_label), factor(true_label))
    svmlin_fscore <- 2*svmlin_ppv*svmlin_sen/(svmlin_ppv+svmlin_sen)

    # get variable importance and corresponding var names
    svmlin_imp <- Importance(svmlin_final_fit, model_dat)$imp[which(names(model_dat)!="gp")]
    svmlin_var <- names(model_dat)[which(names(model_dat)!="gp")]

    # record in data table
    svmlin_mods <- data.frame(iter=rep(i, n_ion),
                            method=rep("svmlin", n_ion),
                            auc=rep(svmlin_auc, n_ion),
                            sen=rep(svmlin_sen, n_ion),
                            spe=rep(svmlin_spe, n_ion),
                            ppv=rep(svmlin_ppv, n_ion),
                            fscore=rep(svmlin_fscore, n_ion),
                            var=svmlin_var,
                            varImp=svmlin_imp)





    #===============================================================#
    ## SVM radial kernel ####

    # 10-fold cross-validation
    svmrad_train <-   train(gp~.,
                            data=model_dat,
                            method="svmRadial",
                            trControl=ctrl,
                            tuneLength=10,
                            metric="ROC")

    # fit final model using rminer
    svmrad_final_fit <- rminer::fit(
      gp~.,
      data=model_dat,
      task="prob",
      model="ksvm",
      kernel="rbfdot",
      C=svmrad_train$bestTune$C,
      kpar=list(sigma=svmrad_train$bestTune$sigma)
    )


    # Predict on holdout test set
    svmrad_test <- predict(svmrad_final_fit, newdata=holdout_dat, type="prob")
    svmrad_pred_prob <- svmrad_test[,colnames(svmrad_test)=="CASE"]
    svmrad_pred_label <- as.factor(ifelse(svmrad_pred_prob>=0.5,"CASE", "CNTL"))
    true_label <- holdout_dat$gp


    # get prediction metrics
    svmrad_auc <- roc(true_label, svmrad_pred_prob)$auc
    svmrad_sen <- sensitivity(svmrad_pred_label, true_label)
    svmrad_spe <- specificity(svmrad_pred_label, true_label)
    svmrad_ppv <- posPredValue(factor(svmrad_pred_label), factor(true_label))
    svmrad_fscore <- 2*svmrad_ppv*svmrad_sen/(svmrad_ppv+svmrad_sen)

    # get variable importance and corresponding var names
    svmrad_imp <- Importance(svmrad_final_fit, model_dat)$imp[which(names(model_dat)!="gp")]
    svmrad_var <- names(model_dat)[which(names(model_dat)!="gp")]

    # record in data table
    svmrad_mods <- data.frame(iter=rep(i, n_ion),
                              method=rep("svmrad", n_ion),
                              auc=rep(svmrad_auc, n_ion),
                              sen=rep(svmrad_sen, n_ion),
                              spe=rep(svmrad_spe, n_ion),
                              ppv=rep(svmrad_ppv, n_ion),
                              fscore=rep(svmrad_fscore, n_ion),
                              var=svmrad_var,
                              varImp=svmrad_imp)



  #===============================================================#
  # Combine all outputs from the five sets of models

  all_mods <- rbind(
    lasso_mods,
    plsda_mods,
    rf_mods,
    svmlin_mods,
    svmrad_mods
  )

  all_mods
}




# Clear up all parallel clusters
stopCluster(clust)
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()



write_csv(mods, "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/holdout_duct_meta_100times_varimp.csv")

print("Done with running all metabolites without risk factors.")


# remove(mods)
