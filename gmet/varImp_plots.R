##### Varimp plots #########

varimp1 <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/holdout_results_FDR01_100times_no_svm_unscale_varimp.csv") %>%
  select(-varImp) %>%
  rename(varImp=varimp2)
varimp2 <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/holdout_results_FDR01_100times_only_svm.csv") %>%
  select(-iter)

varimp <- rbind(varimp1, varimp2)
varimp %>% head()



ion_name=as.data.frame(rbind(
  c("ion_63", "Monomethyl sulfate"),
  c("ion_123", "Oxalurate"),
  c("ion_147", "2,3,5-Trihydroxytoluene"),
  c("ion_184", "Cinnamate"),
  c("ion_200", "2-Isopropyl-1,4-benzenediol"),
  c("ion_227", "Isopropylmaleate"),
  c("ion_630", "10-Hydroxymyristic acid methyl ester"),
  c("ion_634", "Celerin"),
  c("ion_689", "N-[4''-hydroxy-(E)-cinnamoyl]-L-aspartic acid"),
  c("ion_724", "16a-Hydroxyestrone"),
  c("ion_729", "Lauryl diethanolamide"),
  c("ion_809", "Arachidonate"),
  c("ion_1110", "Corchoionoside A"),
  c("ion_1139", "N-Benzoyl-D-arginine-4-nitroanilide"),
  c("ion_1158", "Acetyl tributyl citrate"),
  c("ion_1194", "Atractyloside G"),
  c("ion_1260", "Fuscaxanthone C"),
  c("ion_1282", "1-a,24R,25-Trihydroxyvitamin D2"),
  c("ion_1370", "11-Oxo-androsterone glucuronide"),
  c("ion_1377", "Fusicoccin H"),
  c("ion_1600", "PE(18:1/18:1)"),
  c("ion_1612", "PE(16:1/22:6)"),
  c("ion_1624", "PE(22:6/18:1)"),
  c("ion_1625", "PE(22:5/18:1)")
))
colnames(ion_name) <- c("metabolite", "name")



imp_m <- varimp %>%
  mutate(method=case_when(
    method=="svmlin" ~ "svm_linear",
    method=="svmrad" ~ "svm_radial",
    method=="plsda" ~ "pls-da",
    .default=method)) %>%
  group_by(method, var) %>%
  summarize(
    meanImp=mean(varImp),
    lb=quantile(varImp, probs = 0.025),
    ub=quantile(varImp, probs=0.975),
    .groups = "drop") %>%
  inner_join(ion_name,
             by=join_by(var==metabolite), keep=TRUE)


#================================================================#
library(scales)
show_col(hue_pal()(5))


## Plot for lasso
imp_m %>%
  filter(method=="lasso") %>%
  mutate(name=fct_reorder(name, meanImp)) %>%
  ggplot(aes(name, meanImp)) +
  geom_col(fill="#F8766D") +
  geom_errorbar(aes(x=name, ymin=lb, ymax=ub), width=0.5)+
  labs(x="", y="Absolute penalized coefficients")+
  ylim(0,0.6)+
  coord_flip()+
  theme_classic()



## Plot for PLS-DA
imp_m %>%
  filter(method=="pls-da") %>%
  mutate(name=fct_reorder(name, meanImp)) %>%
  ggplot(aes(name, meanImp)) +
  geom_col(fill="#A3A500") +
  geom_errorbar(aes(x=name, ymin=lb, ymax=ub), width=0.5)+
  labs(x="", y="Weighted sum of absolute coefficients")+
  coord_flip()+
  theme_classic()



## Plot for RF
imp_m %>%
  filter(method=="rf") %>%
  mutate(name=fct_reorder(name, meanImp)) %>%
  ggplot(aes(name, meanImp)) +
  geom_col(fill="#00BF7D") +
  geom_errorbar(aes(x=name, ymin=lb, ymax=ub), width=0.5)+
  labs(x="", y="Loss in OOB accuracy")+
  coord_flip()+
  theme_classic()



## Plot for SVM_linear
imp_m %>%
  filter(method=="svm_linear") %>%
  mutate(name=fct_reorder(name, meanImp)) %>%
  ggplot(aes(name, meanImp)) +
  geom_col(fill="#00B0F6") +
  geom_errorbar(aes(x=name, ymin=lb, ymax=ub), width=0.5)+
  labs(x="", y="Average AAD from median")+
  coord_flip()+
  theme_classic()




## Plot for SVM_radial
imp_m %>%
  filter(method=="svm_radial") %>%
  mutate(name=fct_reorder(name, meanImp)) %>%
  ggplot(aes(name, meanImp)) +
  geom_col(fill="#E76BF3") +
  geom_errorbar(aes(x=name, ymin=lb, ymax=ub), width=0.5)+
  labs(x="", y="Average AAD from median")+
  coord_flip()+
  theme_classic()





# One plot with all five

imp_m %>%
  mutate(name=fct_reorder(name, meanImp)) %>%
  # group_by(method) %>%
  ggplot(aes(name, meanImp)) +
  geom_col(aes(fill=method)) +
  geom_errorbar(aes(x=name, ymin=lb, ymax=ub), width=0.5)+
  coord_flip()+
  theme_classic()+
  facet_wrap(~method, ncol=3)
