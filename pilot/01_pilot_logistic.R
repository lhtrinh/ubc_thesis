source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/01_pilot_pretreatment.R")
library(glmnet)




#===========================================#
# select variables for modeling

pm_dat_for_logistic <- pm_dat %>%
  select(
    studyid, match_id, gp,
    sdc_age_calc,
    # ethnicity,
    edu_level, income_level,
    bmi, wh_menopause_ever,
    fam_hist_breast,
    wh_menstruation_age,
    wh_contraceptives_ever, wh_contraceptives_age, wh_contraceptives_duration,
    wh_gravidity, wh_preg_first_age, wh_live_births,
    wh_hrt_ever, wh_hrt_age,
    fam_hist_breast, bmi,
    alc_cur_cat, alc_binge_cat,
    smk_cig_status, smk_cig_dur, smk_cig_qty,
    nut_fruits_qty, nut_veg_qty,
    starts_with("ion")
    ) %>%
  mutate(across(c(gp,
                  income_level,
                  wh_menopause_ever,
                  wh_contraceptives_ever,
                  wh_hrt_ever,
                  fam_hist_breast,
                  smk_cig_status,
                  alc_cur_cat,
                  alc_binge_cat),
                as.factor)) %>%
  ungroup()









#===========================================#
# Train-test split
ids <- unique(pm_dat_for_logistic$match_id)
set.seed(345)
samp <- sample(ids, size=length(ids)*.7, replace=TRUE)

# split data
train_df <- pm_dat_for_logistic %>%
  filter(match_id %in% samp) %>%
  select(gp, starts_with("ion"))
train_x <- model.matrix(gp~., train_df)[,-1]
train_y <- train_df$gp

test_df <- pm_dat_for_logistic %>%
  filter(!match_id %in% samp) %>%
  select(gp, starts_with("ion"))
test_x <- model.matrix(gp~., test_df)[,-1]
test_y <- test_df$gp




#===========================================#
# logistic regression

train1 <- train_df %>%
  select(gp, starts_with("ion"))
test1 <- test_df %>% select(gp, starts_with("ion"))

mod1 <- glm(gp~., data=test1, family=binomial(link="logit"))
pred1 <- predict(mod1, newdata=test1, type="response")
pred_gp1 <- ifelse(pred1>=0.5, 1, 0)
table(test1$gp, pred_gp1)

# perfect separation, definitely overfitting








#===========================================#
# lasso with metabolites only

set.seed(1245)
mod2_cv <- cv.glmnet(train_x, train_y,
                  family = "binomial",
                  alpha=1)
plot(mod2_cv)

bestlam <- mod2_cv$lambda.min

mod2_lasso <- glmnet(train_x, train_y,
                     alpha=1,
                     family="binomial",
                     lambda=bestlam)

# predicted probabilities
pred2 <- predict(mod2_lasso, s=bestlam, newx=test_x, type="response")
pred2_gp <- ifelse(pred2>=0.5, 1, 0)
table(test_y, pred2_gp)

# select important features
pred2_coef <- predict(mod2_lasso, s=bestlam, newx=test_x, type="coefficients")
rownames(pred2_coef)[pred2_coef[,"s1"]!=0]

round(pred2_coef[pred2_coef[,"s1"]!=0], 3)









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

