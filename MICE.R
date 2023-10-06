#======================================================#
library(mice)


dat_pre_impute <- full_final %>%
  select(
    gp, cohort, studyid, match_id, #dup,
    sdc_age_calc, collect_year, collect_age,
    menopause_stt,
    ethnicity,
    edu_level,
    income_level,
    wh_menstruation_age, wh_menopause_age,
    wh_contraceptives_ever,
    wh_gravidity, wh_live_births, wh_preg_first_age,
    wh_breastfeeding_duration,
    wh_hrt_ever, wh_hrt_duration_yr,
    fam_hist_breast, fam_hist_breast_cat,
    #nut_veg_qty, nut_fruits_qty,
    alc_cur_freq_cat,
    alc_binge_cat,
    smk_cig_stt, # need to find more smoking variables
    bmi, bmi_cat #,
    #starts_with("ion")
  ) %>%
  mutate(across(c(gp, cohort,
                  ethnicity,
                  menopause_stt,
                  edu_level, income_level,
                  wh_contraceptives_ever,
                  wh_hrt_ever,
                  fam_hist_breast,
                  bmi_cat,
                  alc_cur_freq_cat,
                  alc_binge_cat,
                  smk_cig_stt),
                as.factor)) %>%
  mutate(across(c(sdc_age_calc,
                wh_menstruation_age, wh_menopause_age,
                wh_gravidity, wh_live_births, wh_preg_first_age,
                bmi),
                as.numeric)) %>%
  ungroup()




summary(dat_pre_impute)


dat_for_impute <- dat_pre_impute %>% select(-studyid)

#======================================================#

# imputation

set.seed(292920)
init <- mice(subset(dat_pre_impute, select=-studyid),
             m=5,
             maxit=10,
             seed=292920)

comp_impute <- complete(init)

comp_impute





#======================================================#
# logistic regression with each df


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
