####################################################
#  PILOT SAMPLE ANALYSIS  #
####################################################
source("C:/Users/lyhtr/OneDrive - UBC/Thesis/Code/ubc_thesis/00_Setup_both.R")



# load pilot IDs
pilot_ids <- read.csv("C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/ATP/2105_Study2105_02_23_2023_Bhatti_Pilot_Box_Pulls_Plating_2105_Researcher.csv") %>%
  transmute(
    tubeid=as.character(Tube.ID),
    arm=ARM,
    participantkey=participant.Key)

# load metabolomics data
pilot_ions <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Pilot/DATA.xlsx",
                        sheet="ions")
pilot_ion_matrix <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Pilot/DATA.xlsx",
                              sheet="ion_matrix",
                              skip=2)
pilot_annotation <- read_xlsx(path="C:/Users/lyhtr/OneDrive - UBC/Thesis/Data/Pilot/DATA.xlsx",
                        sheet="annotation")




####################################################
## remove blanks and transpose metabolomics data
head(pilot_ions)

# transpose data so each row is one participant
temp <- pilot_ions %>%
  select(ionTopFormula, contains("Sample")) %>%
  pivot_longer(!ionTopFormula,
               names_to="sampleid") %>%
  pivot_wider(names_from = ionTopFormula,
              values_from = value) %>%
  mutate(sampleid=str_replace(sampleid, "Sample_", ""))

pilot_meta <- pilot_ids %>%
  inner_join(temp,
             by=c("tubeid"="sampleid"))



names(pilot_meta)[1:10]
