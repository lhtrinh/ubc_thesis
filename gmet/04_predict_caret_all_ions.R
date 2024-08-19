#=================================================#
#=================================================#


#  LASSO FOR META AND SURVEY DATA ####


#=================================================#
#=================================================#



source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_functions.R")




#=============================#
## Setup ####
library(caret)
library(randomForest)
library(glmnet)
library(pROC)
library(ROCR)
library(pls)
library(rio)



full_dat <- import_survey()
full_sc <- import_meta_sc()

full_all <- full_dat %>% full_join(full_sc)


all_pvals <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/all_pvals.csv")
ion_cols <- all_pvals$metabolite

sig_ions_0.1 <- all_pvals$metabolite[all_pvals$q_fdr<=0.1]
sig_ions_0.2 <- all_pvals$metabolite[all_pvals$q_fdr<=0.2]
#=============================#




#=============================#
## Assign modeling and holdout set ####



my_df <- full_all %>%
  select(gp, studyid, match_id, starts_with("ion"))
names(my_df)


# Randomly sample by match ID
matchIds <- unique(my_df$match_id)
length(matchIds)

set.seed(2931)
trainMatchId <- sample(matchIds, length(matchIds)*.80, replace = FALSE)


# Split data training (90%) and validation (10%) sets
model_dat <- my_df %>% filter(match_id %in% trainMatchId) %>% select(-studyid, -match_id)
holdout_dat <- my_df %>% filter(!match_id %in% trainMatchId) %>% select(-studyid, -match_id)

#=============================#





#=============================#
# Create fold indices for 10-fold CV
## randomly split match IDs into 10 equal-sized folds
n_fold <- 10
# set.seed(39671)
foldInd <- sample(rep(sample(1:n_fold), length=length(trainMatchId)), replace=FALSE)
table(foldInd)

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





#=============================#
## Lasso ####

# Create grid of hyper parameters
lasso_grid <- expand.grid(alpha=1,
                          lambda=seq(.001,1,length=100))


# Fit model
lasso_train <- train(gp~.,
                     data=model_dat,
                     method="glmnet",
                     trControl=ctrl,
                     tuneGrid=lasso_grid,
                     metric="ROC")

lasso_train

# extract fold results
lasso_train$results %>%
  filter(lambda==lasso_train$bestTune$lambda)


lasso_train$bestTune




#===============================================================#


# Predict on holdout test set
lasso_pred <- predict(lasso_train, lambda=lasso_train$bestTune, newdata=holdout_dat, type="prob")
lasso_pred_prob <- lasso_pred$CASE
lasso_pred_label <- as.factor(ifelse(lasso_pred_prob>=0.5,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(lasso_pred_label, true_label)


lasso_roc <- roc(true_label, lasso_pred_prob)
lasso_roc$auc

plot(lasso_roc)

sensitivity(lasso_pred_label, true_label)
specificity(lasso_pred_label, true_label)
posPredValue(factor(lasso_pred_label), factor(true_label))


pred <- prediction(lasso_pred_label, true_label)




#=============================================================#
## Random forest ####
rf_grid <- expand.grid(.mtry=1:15)

rf_train <- train(gp~.,
                  data=model_dat,
                  method="rf",
                  trControl=ctrl,
                  tuneGrid=rf_grid,
                  metric="ROC")


# extract fold results
rf_train$results %>%
  filter(mtry==rf_train$bestTune$mtry)


rf_train$bestTune




#=============================#
# Predict on holdout test set
rf_pred <- predict(rf_train, newdata=holdout_dat, type="prob")
rf_pred_prob <- rf_pred$CASE
rf_pred_label <- as.factor(ifelse(rf_pred_prob>=0.6,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(rf_pred_label, true_label)


rf_roc <- roc(true_label, rf_pred_prob)
rf_roc
plot(rf_roc)


sensitivity(rf_pred_label, true_label)
specificity(rf_pred_label, true_label)
posPredValue(factor(rf_pred_label), factor(true_label))






#=============================================================#
## Partial least squares ####
plsda_train <- train(gp~.,
                     data=model_dat,
                     method="pls",
                     scale=TRUE,
                     trControl=ctrl,
                     tuneLength=10,
                     # tuneGrid=rf_grid,
                     metric="ROC")

plsda_train

# Predict on holdout test set
pls_pred <- predict(plsda_train, ncomp=plsda_train$bestTune, newdata=holdout_dat, type="prob")
pls_pred_prob <- pls_pred$CASE
pls_pred_label <- as.factor(ifelse(pls_pred_prob>=0.5,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(pls_pred_label, true_label)


pls_roc <- roc(true_label, pls_pred_prob)
pls_roc
plot(pls_roc)


sensitivity(pls_pred_label, true_label)
specificity(pls_pred_label, true_label)
posPredValue(factor(pls_pred_label), factor(pls_pred_label))





#=============================================================#
## Support vector machine (SVM) ####
svm_grid <- expand.grid(.mtry=1:15)

svm_train <- train(gp~.,
                  data=model_dat,
                  method="svm",
                  trControl=ctrl,
                  tuneGrid=svm_grid,
                  metric="ROC")


# extract fold results
svm_train$results %>%
  filter(mtry==svm_train$bestTune$mtry)


svm_train$bestTune




#=============================#
# Predict on holdout test set
svm_pred <- predict(svm_train, newdata=holdout_dat, type="prob")
svm_pred_prob <- svm_pred$CASE
svm_pred_label <- as.factor(ifelse(svm_pred_prob>=0.6,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(svm_pred_label, true_label)


svm_roc <- roc(true_label, svm_pred_prob)
svm_roc
plot(svm_roc)


sensitivity(svm_pred_label, true_label)
specificity(svm_pred_label, true_label)
posPredValue(factor(svm_pred_label), factor(true_label))

