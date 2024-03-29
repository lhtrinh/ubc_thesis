df1 <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Preds/fdr01_norisk.csv") %>%
  mutate(f1=2*sen*ppv/(sen+ppv))

df2 <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Preds/fdr01_risk.csv") %>%
  mutate(f1=2*sen*ppv/(sen+ppv))


labs = c("AUC", "F1 score", "sensitivity", "specificity")
names(labs)=c("auc", "f1", "sen", "spe")


df1 %>%
  select(-iter, -ppv) %>%
  pivot_longer(-method,
               names_to="metric") %>%
  ggplot()+
  geom_boxplot(aes(y=value, fill=method)) +
  ylim(0.4, 0.75)+
  facet_wrap(
    ~metric, nrow=4,
    labeller=labeller(metric=labs))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


df2 %>%
  select(-iter, -dfver, -ppv) %>%
  pivot_longer(-method,
               names_to="metric") %>%
  ggplot()+
  geom_boxplot(aes(y=value, fill=method)) +
  ylim(0.4, 0.75)+
  facet_wrap(
    ~metric, nrow=4,
    labeller=labeller(metric=labs))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())




#===================================================================#
## Subgroups ####
df <- read_csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Preds/subgroup_5yr.csv") %>%
  mutate(f1=2*sen*ppv/(sen+ppv))


df %>%
  select(-iter, -ppv) %>%
  pivot_longer(-method,
               names_to="metric") %>%
  ggplot()+
  geom_boxplot(aes(y=value, fill=method)) +
  ylim(0.4, 0.75)+
  facet_wrap(
    ~metric, nrow=4,
    labeller=labeller(metric=labs))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
