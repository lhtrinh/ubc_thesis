#======================================================#
# try mice with atp data
library(mice)


dat_pre_impute <- dat %>%
  filter(cohort=="atp") %>%
  select(
    studyid, # match_id,
    gp,
    sdc_age_calc,
    ethnicity,
    sdc_edu_level, sdc_income,
    bmi,
    menopause_stt,
    wh_menstruation_age,
    starts_with("wh_contraceptives"),
    wh_gravidity, wh_preg_first_age, wh_live_births,
    wh_hrt_ever, # wh_hrt_age,
    fam_hist_breast,
    bmi,
    alc_cur_freq, alc_binge_freq_female,
    # alc_cur_cat, alc_binge_cat,
    smk_cig_status, # smk_cig_dur, smk_cig_qty,
    nut_fruits_qty, nut_veg_qty,
    starts_with("ion")
  ) %>%
  mutate(across(c(gp,
                  ethnicity,
                  sdc_edu_level, sdc_income,
                  menopause_stt,
                  wh_contraceptives_ever,
                  wh_hrt_ever,
                  fam_hist_breast,
                  alc_cur_freq,
                  alc_binge_freq_female,
                  smk_cig_status),
                as.factor)) %>%
  mutate(across(c(sdc_age_calc,
                wh_menstruation_age,
                wh_contraceptives_age,
                wh_gravidity,
                wh_preg_first_age,
                nut_fruits_qty,
                nut_veg_qty),
                as.numeric)) %>%
  ungroup()




summary(dat_pre_impute)


dat_for_impute <- dat_pre_impute %>% select(-studyid)

#======================================================#

# imputation

set.seed(1034)
init <- mice(subset(dat_pre_impute, select=-studyid),
             m=10,
             maxit=10)

comp_impute <- complete(ini)

comp_impute





#======================================================#
# logistic regression with each df


imp_lasso_df <- data.frame(
  imp_id <- rep(0, 10)

)






i=10
# for (i in 1:10){
  # extract imputed data
  imp_dat <- complete(init, action = i)
  pm_contextual <- cbind(dat_pre_impute[,"studyid"], imp_dat)

  pm_dat_lasso <- pm_treated2 %>% left_join(pm_contextual) %>%
    mutate(gp=as.factor(gp))


  # split data
  train_both <- pm_dat_lasso %>%
    filter(match_id %in% samp) %>%
    select(-match_id, -studyid)
  train_x_both <- model.matrix(gp~., train_both)[,-1]
  train_y_both <- train_both$gp

  test_both <- pm_dat_lasso %>%
    filter(!match_id %in% samp) %>%
    select(-c(match_id, studyid))
  test_x_both <- model.matrix(gp~., test_both)[,-1]
  test_y_both <- test_both$gp


  # cv for lambda
  set.seed(3939)
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
  preds_gp <- ifelse(preds>=0.5, 1, 0)
  table(test_y_both, preds_gp)

  # select important features
  preds_coef <- predict(mod, s=lam, newx=test_x_both, type="coefficients")
  rownames(preds_coef)[preds_coef[,"s1"]!=0]

  round(preds_coef[preds_coef[,"s1"]!=0], 3)
  # }

