####################################################
#  PILOT SAMPLE ANALYSIS  #
####################################################
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup.R")



# load pilot IDs
pilot_ids <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/2105_Study2105_02_23_2023_Bhatti_Pilot_Box_Pulls_Plating_2105_Researcher.csv") %>%
  transmute(
    tubeid=as.character(Tube.ID),
    gp=ifelse(ARM=="Case", "CASE", "CNTL"),
    participantkey=participant.Key)

# load metabolomics data
# pilot_ions <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Pilot/DATA.xlsx",
#                         sheet="ions")
pilot_ion_matrix <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Pilot/DATA.xlsx",
  sheet="ion_matrix")
pilot_annotation <- read_xlsx(
  path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Pilot/DATA.xlsx",
  sheet="annotation")








####################################################
# transpose metabolites
head(pilot_ion_matrix)

# remove blank rows
temp <- pilot_ion_matrix[-(1:2),-2]
head(temp)

names(temp)[1] <- "ion"

# each sample is measured twice
# Excel reads the second column of a merged cell as blank
# therefore, all odd-numbered columns are blanked and automatically given a name by read_xlsx
# add names back
ids <- 1:ncol(temp) %% 2
data_even <- temp[,ids==0]
names(data_even)
colnames(temp)[ids==1][-1] <- paste(colnames(data_even), "2", sep="_")
head(temp)


# transpose data so that rows=samples and columns=metabolites
temp2 <- temp %>%
  pivot_longer(!ion,
               names_to="sampleid") %>%
  separate(sampleid,
           into=c("type", "sample_id", "inject_num"),
           sep="_") %>%
  mutate(ion=paste("ion", ion, sep="_"),
         id=paste(type, sample_id, sep="_"),
         inject_num=as.numeric(replace_na(inject_num, "1"))) %>%
  pivot_wider(names_from = ion,
              values_from = value)




####################################################



####################################################
# combine meta with matched IDs
temp <- temp2 %>%
  left_join(pilot_ids,
            by=c("sample_id"="tubeid")) %>%
  left_join(atp_match_ids)



####################################################
pilot_meta <- temp





# pilot_dat <- dat %>%
#   inner_join(pilot_meta)
