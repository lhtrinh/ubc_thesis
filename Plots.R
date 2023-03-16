my_pink <- rgb(red=255, green=51, blue=153, maxColorValue = 255)
my_grey <- rgb(red=118, green=113, blue=113, maxColorValue = 255)
my_color <- c(my_pink, my_grey)

my_theme <- theme(
  line=element_line(
    color=my_color,
    linewidth=2
  ),
  panel.background = element_rect(
    fill="white",
    linetype=NULL
  )
)

# plot age by group
dat %>%
  ggplot(aes(sdc_age_calc, y=..count.., group=gp, color=gp)) +
  geom_density(size=.9) +
  labs(x="age at baseline", color="group") +
  theme_classic() +
  ylim(0, 40) +
  scale_color_manual(values=my_color)


# plot hr presence
dat %>% filter(gp=='CASE') %>%
  pivot_longer(
  cols=c('er_status', 'pr_status', 'her2_status'),
  names_to='hr',
  values_to='status'
) %>%
  mutate(hr=factor(hr, levels=c("er_status", "pr_status", "her2_status"))) %>%
  select(hr, status) %>%
  ggplot(aes(hr, fill=status)) +
  labs(x="", fill="") +
  geom_bar() +
  scale_fill_brewer(palette = "Blues") +
  theme_classic()


# plot hr subgroups
dat %>% filter(gp=='CASE') %>%
  # mutate(
  #   hr_subgroup_alt=factor(hr_subgroup_alt,
  #                      levels=c("er/pr_positive", "her2_positive", "triple_negative", "unknown"))) %>%
  ggplot(aes(hr_subgroup_alt, fill=hr_subgroup_alt)) +
  labs(x="") +
  geom_bar() +
  scale_fill_brewer(palette = "Purples") +
  theme_classic()




# plot family hist of breast cancer
dat %>%
  ggplot(aes(gp, fill=as.factor(fam_hist_breast)))+
  geom_bar() +
  labs(x="", fill="")+
  theme_classic() +
  scale_fill_manual(
    labels=c("No family history of breast cancer", "Mother and/or sibling with breast cancer"),
    values=c("darkgrey", "pink"),
    na.value = "lightgrey") +
  coord_flip()



# plot age at diag
dat %>%
  mutate(
    dx_age_cat=findInterval(age_at_diagnosis,
                   c(40, 50, 60, 70))) %>% count(dx_age_cat)
  ggplot(aes(dx_age_cat))+
  geom_bar() +
  theme_classic()

