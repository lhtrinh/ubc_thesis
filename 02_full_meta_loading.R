#################################
# Full metabolomics data
#################################



### Load packages
library(tidyverse)
library(readxl)

# load full questionnaire data
full_dat <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/data_with_missing.csv")


# load metabolite data
full_ion_intensities <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/spreadsheets/DATA.xlsx",
  sheet="intensities1")
full_ion_intensities %>% head()

full_injection <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/spreadsheets/DATA.xlsx",
  sheet="injections")
full_injection %>% head()







#======================================================================#
# load crosswalk between participant and sample IDs for ATP
atp_id_crosswalk <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Metabolomics/Full/Bhatti samples pulled and duplicate IDs 2023-06-02.xlsx",
  sheet="231Cs, 231Cntr w duplc"
)
names(atp_id_crosswalk)


# 1. fix column names and remove the last column (blank column with one note only)
# 2. remove leading 00s from barcode columns
# 3. combine original and duplicated barcodes into one single column named "sampleid"
# 4. mark columns with replicate barcodes as dup
atp_id_cw2 <- atp_id_crosswalk %>%
  rename(dup_barcode = `tube barcode of duplicate sent`) %>%
  rename_with(tolower) %>%
  select(participantkey, aliquotbarcode, dup_barcode) %>%
  mutate(across(c(aliquotbarcode, dup_barcode),
                function(x) as.character(as.numeric(x))),
         dup=ifelse(!is.na(dup_barcode), "YES", "NO")) %>%
  pivot_longer(-c(participantkey, dup),
               values_to="barcode",
               values_drop_na=TRUE)

head(atp_id_cw2)
summary(atp_id_cw2)


# double check for duplicates
atp_id_cw2 %>%
  count(barcode) %>%
  arrange(desc(n)) #2 duplicates

atp_id_cw2 %>%
  filter(barcode==170965728 | barcode==81869160) #same participant

# there is one duplicate barcode that appeared twice
# looks like the same ID throughout
# we can remove the "double-duplicated" barcode
atp_id_cw3 <- atp_id_cw2[!duplicated(atp_id_cw2),]



# final crosswalk IDs
atp_id_cw <- atp_id_cw3
head(atp_id_cw)


#======================================================================#
# add "ion" in front of each metabolite ID
# transpose full ion intensities so rows = samples, columns = metabolites
# extract sample ID and change to numeric form
full_ion_t <- full_ion_intensities %>%
  mutate(ion=paste("ion", ionIdx, sep="_")) %>% # add "ion" in front of each metabolite ID
  select(-ionIdx, -ionMz, -ionActive) %>% # remove irrelevant columns
  pivot_longer(!ion, names_to = "sampleid") %>%
  pivot_wider(names_from = ion,
              values_from = value) %>% # transpose data
  mutate(sampleid = as.numeric(str_replace_all(sampleid, "[^[:digit:]]", ""))) # change sample ID to numeric form




# separate sample IDs and cohort
# also change all ATP sample IDs to corresponding studyid
### 1.  If the cohort is ATP, remove the string "ATP_" in the front
###     and add an underscore before "a" or "b" to match BCGP strings
### 2.
full_inj2 <- full_injection %>%
  mutate(tubeid_iter = ifelse(dsUser3=="ATP",
                              gsub('^(.{4})(.*)(a|b)$', '\\2_S1\\3', dsCode),
                              dsCode)) %>%
  mutate(tubeid = str_replace_all(tubeid_iter, "(a|b)$", "")) %>%
  separate(tubeid,
           into=c("sampleid", "iter"),
           sep="_",
           convert=TRUE,
           remove=FALSE) %>%
  separate(dsUser4,
           into=c("x", "cohort", "plate"),
           sep="_",
           convert=TRUE,
           extra="merge") %>%
  mutate(plate=as.numeric(str_replace_all(plate, "[^[:digit:]]", "")),
         cohort=dsUser3,
         sampletype=dsPert)



# join injection data with crosswalk
# if sample is in the ATP cohort, replace sampleid with participantkey
full_inj3 <- full_inj2 %>%
  left_join(atp_id_cw,
            by=join_by(sampleid==barcode),
            keep=TRUE) %>%
  mutate(studyid=case_when(
    cohort=="BCGP" ~ sampleid,
    cohort=="ATP" ~ participantkey,
    TRUE ~ sampleid
  ))


# add duplicate info from BCGP
# We also label QC samples as "dup" since we are using them for QC

dup_id_bcgp <- unique(full_dat$studyid[full_dat$dup=="YES"])
dup_id_bcgp <- dup_id_bcgp[!is.na(dup_id_bcgp)]
dup_id_bcgp

full_inj4 <- full_inj3 %>%
  mutate(dup=ifelse(studyid %in% dup_id_bcgp |
                      sampletype=="QC", "YES", dup))



full_inj4 %>% filter(dup=="YES") %>% count(studyid, cohort) %>% count(n)


full_inj <- full_inj4



#################################################################
full_ion_plate <- full_ion_t %>%
  full_join(full_inj,
            by=c("sampleid"="dsIdx")) %>%
  select(studyid, tubeid, iter, dup, sampletype, plate, contains("ion"))

# number of study samples
length(unique(full_ion_plate$studyid[full_ion_plate$sampletype=="Sample"]))



write.csv(full_ion_plate, file = "C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/meta_non_normalized.csv")

