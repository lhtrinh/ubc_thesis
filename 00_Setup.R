library(tidyverse)
library(haven)
library(readxl)
library(lubridate)


# load data
dat_raw <- read_sav("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP2641_DATA_PT_Qx_RegistryInfo_BIO.sav")
# change all variable names to lowercase
colnames(dat_raw) <- tolower(colnames(dat_raw))
colnames(dat_raw)


# load data dictionary files


path_dict <- file.path("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/BCGP2641_Data_Dictionary.xlsx")
dat_dict <- read_excel(path_dict, sheet="variables")
head(dat_dict)
colnames(dat_dict)
dat_dict$variable_name <- tolower(dat_dict$variable_name)

dat_cat <- read_excel(path_dict, sheet = "categories")
head(dat_cat)
colnames(dat_cat)
dat_cat$variable_name <- tolower(dat_cat$variable_name)



# add underscore to varnames and change all to lowercase
colnames(dat_dict) <- tolower(str_replace_all(colnames(dat_dict), " ", "_"))
colnames(dat_cat) <- tolower(str_replace_all(colnames(dat_cat), " ", "_"))


names(dat_dict)
names(dat_cat)






