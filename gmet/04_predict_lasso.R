###################################################
###  LASSO FOR META AND SURVEY DATA
###################################################


library(glmnet)
library(pROC)
library(ROCR)
library(caret)




#
#
# #===========================================#
# # select variables for modeling


# metabolite data: full_nonorm_centered
# survey data with missingness: full_dat




#===========================================#
# Train-test split
ids <- unique(full_nonorm_centered$match_id)
set.seed(2345)
samp <- sample(ids, size=length(ids)*.7, replace=FALSE)

# split data
train_df <- full_nonorm_centered %>%
  filter(match_id %in% samp) %>%
  select(gp, starts_with("ion"))
train_x <- model.matrix(gp~., train_df)[,-1]
train_y <- train_df$gp

test_df <- full_nonorm_centered %>%
  filter(!match_id %in% samp) %>%
  select(gp, starts_with("ion"))
test_x <- model.matrix(gp~., test_df)[,-1]
test_y <- test_df$gp







#===========================================#
# lasso with metabolites only

set.seed(12345)
mod2_cv <- cv.glmnet(train_x, train_y,
                     family = "binomial",
                     alpha=1)
plot(mod2_cv)

mod2_cv$lambda.min
mod2_cv$lambda.1se

bestlam <- mod2_cv$lambda.min

mod2_lasso <- glmnet(train_x, train_y,
                     alpha=1,
                     family="binomial",
                     lambda=bestlam)

# predicted probabilities
pred2 <- predict(mod2_lasso, s=bestlam, newx=test_x, type="response")
pred2_gp <- ifelse(pred2>=0.5, "CASE", "CNTL")
table(test_y, pred2_gp)




#####
# calculate predicted sensitivity, specificity, and AUC
mod2_roc <- roc(test_y, pred2)
mod2_roc # auc 0.593
plot(mod2_roc)

# accuracy
sum(test_y==pred2_gp)



# select important features
pred2_coef <- predict(mod2_lasso, s=bestlam, newx=test_x, type="coefficients")
rownames(pred2_coef)[pred2_coef[,"s1"]!=0]









#===========================================#
# lasso with metabolites and contextual factors


# split data
train_both <- pm_dat_for_logistic %>%
  filter(match_id %in% samp) %>%
  select(-match_id, -studyid)
train_x_both <- model.matrix(gp~., train_both)
train_y_both <- train_both$gp

test_both <- pm_dat_for_logistic %>%
  filter(!match_id %in% samp) %>%
  select(-c(match_id, studyid))
test_x_both <- model.matrix(gp~., test_both)[,-1]
test_y_both <- test_both$gp



set.seed(2456)
mod3_cv <- cv.glmnet(train_x_both, train_y_both,
                     family = "binomial",
                     alpha=1)
plot(mod3_cv)

lam_min <- mod3_cv$lambda.min

mod3_lasso <- glmnet(train_x_both, train_y_both,
                     alpha=1,
                     family="binomial",
                     lambda=lam_min)

# predicted probabilities
pred3 <- predict(mod3_lasso, s=bestlam, newx=test_x_both, type="response")
pred3_gp <- ifelse(pred3>=0.5, 1, 0)
table(test_y_both, pred3_gp)

# select important features
pred3_coef <- predict(mod3_lasso, s=bestlam, newx=test_x, type="coefficients")
rownames(pred3_coef)[pred3_coef[,"s1"]!=0]

round(pred3_coef[pred3_coef[,"s1"]!=0], 3)










#======================================================#
# lasso regression with each df


imp_lasso_df <- data.frame(
  imp_id <- rep(0, 10)

)






# i=10
lasso_coef <- tibble(iter=rep(0,1),
                     var=rep("char", 1))

for (i in 1:10){
  # extract imputed data
  full_contextual_imp <- cbind(dat_pre_impute[,"studyid"], complete(init, action=i))

  # join meta with imputed data
  full_dat_lasso <- full_nonorm_centered %>%
    left_join(full_contextual_imp) %>%
    mutate(gp=as.factor(gp))


  # split data
  train_both <- full_dat_lasso %>%
    filter(match_id %in% samp) %>%
    select(-match_id, -studyid, -dup, -sampletype, -cohort)
  train_x_both <- model.matrix(gp~., train_both)[,-1]
  train_y_both <- train_both$gp

  test_both <- full_dat_lasso %>%
    filter(!match_id %in% samp) %>%
    select(-match_id, -studyid, -dup, -sampletype, -cohort)
  test_x_both <- model.matrix(gp~., test_both)[,-1]
  test_y_both <- test_both$gp


  # cv for lambda
  set.seed(39393)
  mod_cv <- cv.glmnet(train_x_both, train_y_both,
                      family = "binomial",
                      alpha=1)
  lam <- mod_cv$lambda.min


  # fit model with lambda min
  mod <- glmnet(train_x_both, train_y_both,
                alpha=1,
                family="binomial",
                lambda=lam)


  # predicted probabilities
  preds <- predict(mod, s=lam, newx=test_x_both, type="response")
  preds_gp <- ifelse(preds>=0.5, "CASE", "CNTL")
  table(test_y_both, preds_gp)

  # evaluate model performance
  mod2_roc <- roc(test_y, pred2)
  mod2_roc # auc 0.593
  plot(mod2_roc)

  # accuracy
  mean(test_y==pred2_gp)

  # select important features
  preds_coef <- predict(mod, s=lam, newx=test_x_both, type="coefficients")
  preds_var <- rownames(preds_coef)[preds_coef[,"s1"]!=0]

  round(preds_coef[preds_coef[,"s1"]!=0], 3)



  lasso_coef <- rbind(lasso_coef,
                      data.frame(iter=rep(i,length(preds_var)),
                                 var=preds_var))
}

lasso_coef <- lasso_coef %>%
  filter(iter!=0 & var!="(Intercept)")

lasso_coef %>%
  count(var) %>%
  filter(n==10) %>%
  View()