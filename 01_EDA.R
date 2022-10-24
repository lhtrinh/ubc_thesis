source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/00_Setup.R")

## Explore data


## How is data coded?

dat_raw %>% slice_sample(n=20) %>% View()

## Empty rows have "." connotations
## Rows with -7 values are not applicable

## Replace "." with empty cells and -7 with ?


clean_func <- function(x) {
  x1 <- str_replace_all(x, "\\.", "")
  x2 <- ifelse(x1=="-7", NA, x1)
  y <- str_trim(x2)
  return(y)
}

temp1 <- data.frame(lapply(dat_raw, clean_func))

View(temp1)
str(temp1)


temp2 <- lapply(temp1,
              function(x) ifelse(str_detect(x, "^[:digit:]+$") |
                                   x=="",
                                 as.numeric(x),
                                 x)) %>%
  as_tibble()

summary(temp2[,81:110])
table(temp1$first_diag_cancer)
temp1 %>% count(gp, followuptime) #all case no
temp1 %>% filter(followuptime=='31') %>% count(collect_year)

dat_raw %>% count(dis_cancer_sib_breast)

# adm_qx_completion is supposed to be date type. Not sure why it shows as numbers.
table(temp1$dis_cancer_m_breast)

# some indicators are only marked for 1. The empty cells are either 0 or NA.
# double check


dim(temp2)

