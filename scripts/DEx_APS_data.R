#### data

# processed data from the core facility
CSF_core_data_processing <- read_excel("data/CSF aPD_individual patient results_after data processing.xlsx", 
                                       sheet = "Quantitative Data")

# raw data
proteinGroups_samples <- read_delim("data/proteinGroups_samples.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

design_APS <- read_delim("data/CSF aPD_individual patient results_after data processing.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

# data with age at lombar puncture
data_age_ctrl <- read_excel("data/updated_miRNA_samples_atypical_cor.xlsx",
                       sheet = "Control")
data_age_PD <- read_excel("data/updated_miRNA_samples_atypical_cor.xlsx", 
                            sheet = "PD")
data_age_aPS <- read_excel("data/updated_miRNA_samples_atypical_cor.xlsx", 
                           sheet = "Atypical")
data_age <- rbind(data_age_ctrl %>% 
                    dplyr::select(Pseudonym,`Alter bei LP / Age at LP`),
                  data_age_PD %>% 
                    dplyr::select(Pseudonym,`Alter bei LP / Age at LP`),
                  data_age_aPS %>% 
                    dplyr::select(Pseudonym,`Alter bei LP / Age at LP`))

# data with NfL info
data_NfL <- read_excel("data/Simoa results all.xlsx")
