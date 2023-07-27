##################################
# Set up full metabolomics data
# Normalized version
##################################



# Load packages
library(tidyverse)
library(readxl)



# load metabolite data
full_ion_norm_raw_raw <- read_excel(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/spreadsheets/DATA_NORM.xlsx",
  sheet="ion_matrix",
  .name_repair="universal")
full_ion_norm_raw[1:10, 1:10]
names(full_ion_norm_raw)

full_inj_norm <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/spreadsheets/DATA_NORM.xlsx",
  sheet="samples")
full_inj_norm %>% head()




# load full questionnaire data
full_dat <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")



#=====================================================#
# format names of ion columns and transpose data
# so that columns are ions and rows are participants

temp1 <- full_ion_norm_raw[-2,-2]
temp1[1:10, 1:10]

names(temp1)[1] <- "ion"
temp1$ion[is.na(temp1$ion)] <- "dsIdx"


temp2 <- temp1 %>%
  pivot_longer(-ion)
temp2_idx <- temp2 %>% filter(ion=="dsIdx") %>%
  select(-ion) %>%
  rename(ionIdx=value)
temp2b <- temp2 %>% filter(ion!="dsIdx")
dim(temp2b)
dim(temp2_idx)


setdiff(temp2_idx$name, temp2b$name)
setdiff(temp2b$name, temp2_idx$name)
# names are identical


temp3 <- temp2 %>%
  mutate(ion_id = ifelse(ion=="dsIdx", ion,
                         paste("ion", ion, sep="_"))) %>%
  select(-ion) %>%
  pivot_wider(names_from="ion_id", values_from="value") %>%
  select(-name)


full_ion <- temp3
full_ion[1:10, 1:10]


#=====================================================#
# merge with injection information
full_ion_merge <- full_inj_norm %>%
  full_join(full_ion)

full_1 <- full_ion_merge %>%
  rename(sampletype=dsSampleType,
         cohort=dsUser3,
         plate=dsPlate,
         iter=dsPlateInj) %>%
  mutate(
    sampleid = ifelse(cohort=="ATP",
                         gsub('^(.{4})(.*)$', '\\2', dsSampleCode),
                         gsub('(.*)_(.*)$', '\\1', dsSampleCode))) %>%
  mutate(tubeid=paste(sampleid, iter, sep="_"))



full_2 <- full_1 %>%
  left_join(atp_id_cw,
            by=c("sampleid"="barcode"),
            keep=TRUE) %>%
  mutate(studyid=case_when(
    cohort=="BCGP" ~ sampleid,
    cohort=="ATP" ~ participantkey,
    TRUE ~ sampleid
  ))





### add duplicates indicator for BCGP and QC samples
full_3 <- full_2 %>%
  mutate(
    dup=ifelse(sampletype=="QC" |
                 studyid %in% full_dat$studyid[full_dat$dup=="YES"],
                                                    "YES", "NO"))


full_3 %>% count(dup, cohort)






#=====================================================#

full_ion_norm <- full_2 %>%
  select(studyid, tubeid, iter, plate, cohort, sampletype, dup, contains("ion"))


full_ion_norm[1:10, 1:10]
