require(tidyverse)
require(readr)
require(forcats)

pvals <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/sig_ions.csv") %>%
  mutate(across(c(beta_unadj, ci_lb, ci_ub), exp))
head(pvals)

mnp <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/all_pvals_mnp.csv") %>%
  mutate(across(c(beta_unadj, ci_lb, ci_ub), exp))
ductal <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ductal_pvals.csv") %>%
  mutate(across(c(beta_unadj, ci_lb, ci_ub), exp))
erpr <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/erpr_pvals.csv") %>%
  mutate(across(c(beta_unadj, ci_lb, ci_ub), exp))




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



## Plot
pvals %>%
  inner_join(ion_name) %>%
  mutate(name = fct_reorder(name, desc(beta_unadj))) %>%
  ggplot(aes(name, beta_unadj))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_lb, ymax=ci_ub), width=0.5)+
  geom_hline(yintercept=1, linetype="dashed")+
  labs(x="", y="Odds ratio")+
  coord_flip()+
  theme_bw()


or_all <- pvals %>% transmute(metabolite=metabolite, or_all=beta_unadj)

## Subgroup plots
mnp %>%
  inner_join(ion_name) %>%
  inner_join(or_all) %>%
  mutate(name = fct_reorder(name, desc(or_all)),
         sig=as.factor(ifelse(q_fdr<0.1,1,0))) %>%
  ggplot(aes(name, beta_unadj, color=sig))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_lb, ymax=ci_ub), width=0.5)+
  geom_hline(yintercept=1, linetype="dashed")+
  labs(x="", y="Odds ratio")+
  scale_color_manual(values=c("0"="grey", "1"="black"))+
  coord_flip()+
  theme_bw()



ductal %>%
  inner_join(ion_name) %>%
  inner_join(or_all) %>%
  mutate(name = fct_reorder(name, desc(or_all)),
         sig=as.factor(ifelse(q_fdr<0.1,1,0))) %>%
  ggplot(aes(name, beta_unadj, color=sig))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_lb, ymax=ci_ub), width=0.5)+
  geom_hline(yintercept=1, linetype="dashed")+
  labs(x="", y="Odds ratio")+
  scale_color_manual(values=c("0"="grey", "1"="black"))+
  coord_flip()+
  theme_bw()



erpr %>%
  inner_join(ion_name) %>%
  inner_join(or_all) %>%
  mutate(name = fct_reorder(name, desc(or_all)),
         sig=as.factor(ifelse(q_fdr<0.1,1,0))) %>%
  ggplot(aes(name, beta_unadj, color=sig))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_lb, ymax=ci_sub), width=0.5)+
  geom_hline(yintercept=1, linetype="dashed")+
  labs(x="", y="Odds ratio")+
  scale_color_manual(
    name="significant",
    values=c("0"="grey", "1"="black"),
    labels=c("no", "yes"))+
  coord_flip()+
  theme_bw()
