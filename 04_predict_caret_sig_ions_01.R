#=================================================#
#=================================================#


#  LASSO FOR METABOLITES, ONLY SIGNIFICANT ####


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
  select(gp, studyid, match_id, all_of(match_cols), all_of(sig_ions_0.1))
names(my_df)


# Randomly sample by match ID
matchIds <- unique(my_df$match_id)
length(matchIds)

set.seed(2931)
trainMatchId <- sample(matchIds, length(matchIds)*.8, replace = FALSE)


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

# extract fold results
lasso_train$results %>%
  filter(lambda==lasso_train$bestTune$lambda)


lasso_train$bestTune




#===============================================================#


# Predict on holdout test set
lasso_test <- predict(lasso_train, lambda=lambda_train$bestTune, newdata=holdout_dat, type="prob")
lasso_pred_prob <- lasso_test$CASE
lasso_pred_label <- as.factor(ifelse(lasso_pred_prob>=0.5,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(lasso_pred_label, true_label)


lasso_roc <- roc(true_label, lasso_pred_prob)


lasso_roc$auc

sensitivity(lasso_pred_label, true_label)
specificity(lasso_pred_label, true_label)
posPredValue(factor(lasso_pred_label), factor(true_label))





#===============================================================#

## Random forest ####
rf_grid <- expand.grid(.mtry=1:15)

rf_train <- train(gp~.,
                  data=model_dat,
                  method="rf",
                  trControl=ctrl,
                  tuneGrid=rf_grid,
                  metric="ROC")
rf_train


#=============================#
# Predict on holdout test set
test <- predict(rf_train, mtry=rf_train$bestTune, newdata=holdout_dat, type="prob")
pred_prob <- test$CASE
pred_label <- as.factor(ifelse(pred_prob>=0.5,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(pred_label, true_label)


roc_ <- roc(true_label, pred_prob)
roc_
plot(roc_)


sensitivity(pred_label, true_label)
specificity(pred_label, true_label)
posPredValue(factor(pred_label), factor(true_label))






#===============================================================#
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
test <- predict(plsda_train, newdata=holdout_dat, type="prob")
pred_prob <- test$CASE
pred_label <- as.factor(ifelse(pred_prob>=0.5,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(pred_label, true_label)


roc_ <- roc(true_label, pred_prob)
roc_
plot(roc_)


sensitivity(pred_label, true_label)
specificity(pred_label, true_label)
posPredValue(factor(pred_label), factor(true_label))






#===============================================================#
## SVM ####
library(kernlab)
# svm_grid  <- expand.grid(sigma=seq(0,2,length=20))

svm_train <- train(gp~.,
                   data=model_dat,
                   method="svmPoly",
                   trControl=ctrl,
                   tuneLength=10)

svm_train$bestTune


# Predict on holdout test set
svm_test <- predict(svm_train, newdata=holdout_dat, type="prob")
svm_pred_prob <- svm_test$CASE
svm_pred_label <- as.factor(ifelse(svm_pred_prob>=0.5,"CASE", "CNTL"))
true_label <- holdout_dat$gp

table(svm_pred_label, true_label)


roc_ <- roc(true_label, svm_pred_prob)
roc_
plot(roc_)


sensitivity(svm_pred_label, true_label)
specificity(svm_pred_label, true_label)
posPredValue(factor(svm_pred_label), factor(true_label))

