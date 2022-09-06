#### check stats of different clinical features in the cohort

# sex as covariate
bio_stats_sex <- table(ExpDesign_APS$Sex,ExpDesign_APS$condition)
fisher.test(bio_stats_sex)

# disease as covariate
bio_stats_disease <- table(ExpDesign_APS$Disease,ExpDesign_APS$condition)
fisher.test(bio_stats_disease, workspace=2e7)

# pathology as covariate
bio_stats_Pathy <- table(ExpDesign_APS$Pathy,ExpDesign_APS$condition)
fisher.test(bio_stats_Pathy)

# Nfl concentration as covariate
bio_stats_Nfl <- table(ExpDesign_APS_Nfl$Nfl_conc,ExpDesign_APS_Nfl$condition)
t.test(bio_stats_Nfl)

# age at lumbar puncture as covariate
bio_stats_age <- table(ExpDesign_APS_age$age,ExpDesign_APS_age$condition)
t.test(bio_stats_age)

## DEP analysis with covariates
# differential analysis with sex as covariate
data_diff_knn_sex <- test_diff(data_proteins_imp_knn,type = "all",
                               design_formula = formula(~0+condition+Sex))
data_diff_MinProb_sex <- test_diff(data_proteins_imp_MinProb,type = "all",
                                   design_formula = formula(~0+condition+Sex))
data_diff_man_sex <- test_diff(data_proteins_imp_man,type = "all",
                               design_formula = formula(~0+condition+Sex))


# differential analysis with Disease as covariate
data_diff_knn_disease <- test_diff(data_proteins_imp_knn,type = "all",
                                   design_formula = formula(~0+condition+Disease))
data_diff_MinProb_disease <- test_diff(data_proteins_imp_MinProb,type = "all",
                                       design_formula = formula(~0+condition+Disease))
data_diff_man_disease <- test_diff(data_proteins_imp_man,type = "all",
                                   design_formula = formula(~0+condition+Disease))

# differential analysis with Nfl as covariate
# data_diff_knn_Nfl<- test_diff_manual(data_proteins_imp_knn_Nfl,type = "manual",
#                               type = "all",
#                               design_formula = formula(~0+condition+Nfl_conc))
# data_diff_MinProb_Nfl <- test_diff(data_proteins_imp_MinProb_Nfl,type = "all",
#                                      design_formula = formula(~0+condition+Nfl_conc))
# data_diff_man_Nfl <- test_diff(data_proteins_imp_man_Nfl,type = "manual",
#                                  design_formula = formula(~0+condition+Nfl_conc))
# 

# differential analysis with age at lumbar puncture as covariate
data_diff_knn_age <- test_diff(data_proteins_imp_knn_age,type = "all",
                               design_formula = formula(~0+condition+age))
data_diff_MinProb_age <- test_diff(data_proteins_imp_MinProb_age,type = "all",
                                   design_formula = formula(~0+condition+age))
data_diff_man_age <- test_diff(data_proteins_imp_man_age,type = "all",
                               design_formula = formula(~0+condition+age))

# add significance 
DEP_knn_sex <- add_rejections(data_diff_knn_sex,alpha = 0.1,lfc = log2(1))
DEP_MinProb_sex <- add_rejections(data_diff_MinProb_sex,alpha = 0.1,lfc = log2(1))
DEP_man_sex <- add_rejections(data_diff_man_sex,alpha = 0.1,lfc = log2(1))

DEP_knn_disease <- add_rejections(data_diff_knn_disease,alpha = 0.1,lfc = log2(1))
DEP_MinProb_disease <- add_rejections(data_diff_MinProb_disease,alpha = 0.1,lfc = log2(1))
DEP_man_disease <- add_rejections(data_diff_man_disease,alpha = 0.1,lfc = log2(1))

DEP_knn_age <- add_rejections(data_diff_knn_age,alpha = 0.1,lfc = log2(1))
DEP_MinProb_age <- add_rejections(data_diff_MinProb_age,alpha = 0.1,lfc = log2(1))
DEP_man_age <- add_rejections(data_diff_man_age,alpha = 0.1,lfc = log2(1))


# PCA with different imputation methods
PCA_DEP_new(DEP_knn_sex, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)
PCA_DEP_new(DEP_MinProb_sex, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)
PCA_DEP_new(DEP_man_sex, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)

PCA_DEP_new(DEP_knn_disease, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)
PCA_DEP_new(DEP_MinProb_disease, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)
PCA_DEP_new(DEP_man_disease, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)

PCA_DEP_new(DEP_knn_age, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)
PCA_DEP_new(DEP_MinProb_age, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)
PCA_DEP_new(DEP_man_age, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition") + geom_text(aes(label = rowname), position = nudge)

# volcano plots
# -> Ctr vs PD
volcano_DEP(DEP_knn_sex, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_sex, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_sex, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_disease, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_disease, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_disease, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_age, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_age, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_age, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Synucleinopathy vs Ctr
volcano_DEP(DEP_knn_sex, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_sex, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_sex, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_disease, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_disease, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_disease, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_age, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_age, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_age, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Synucleinopathy vs PD
volcano_DEP(DEP_knn_sex, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_sex, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_sex, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_disease, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_disease, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_disease, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_age, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_age, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_age, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Synucleinopathy vs Tauopathy
volcano_DEP(DEP_knn_sex, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_sex, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_sex, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_disease, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_disease, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_disease, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_age, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_age, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_age, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Tauopathy vs Ctr
volcano_DEP(DEP_knn_sex, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_sex, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_sex, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_disease, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_disease, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_disease, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_age, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_age, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_age, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Tauopathy vs PD
volcano_DEP(DEP_knn_sex, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_sex, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_sex, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_disease, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_disease, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_disease, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

volcano_DEP(DEP_knn_age, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb_age, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man_age, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)


# get results
data_results_knn_disease <- get_results(DEP_knn_disease)
data_results_knn_disease <- data_results_knn_disease %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_MinProb_disease <- get_results(DEP_MinProb_disease)
data_results_MinProb_disease <- data_results_MinProb_disease %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_man_disease <- get_results(DEP_man_disease)
data_results_man_disease <- data_results_man_disease %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
writexl::write_xlsx(data_results_knn_disease,"data_results_knn_disease.xlsx")
writexl::write_xlsx(data_results_MinProb_disease,"data_results_MinProb_disease.xlsx")
writexl::write_xlsx(data_results_man_disease,"data_results_man_disease.xlsx")
data_results_knn_disease_sign <- data_results_knn_disease %>%
  filter(significant)
data_results_MinProb_disease_sign <- data_results_MinProb_disease %>%
  filter(significant)
data_results_man_disease_sign <- data_results_man_disease %>%
  filter(significant)
writexl::write_xlsx(data_results_knn_disease_sign,"data_results_knn_disease_sign.xlsx")
writexl::write_xlsx(data_results_MinProb_disease_sign,"data_results_MinProb_disease_sign.xlsx")
writexl::write_xlsx(data_results_man_disease_sign,"data_results_man_disease_sign.xlsx")
data_results_knn_sex <- get_results(DEP_knn_sex)
data_results_knn_sex <- data_results_knn_sex %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_MinProb_sex <- get_results(DEP_MinProb_sex)
data_results_MinProb_sex <- data_results_MinProb_sex %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_man_sex <- get_results(DEP_man_sex)
data_results_man_sex <- data_results_man_sex %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
writexl::write_xlsx(data_results_knn_sex,"data_results_knn_sex.xlsx")
writexl::write_xlsx(data_results_MinProb_sex,"data_results_MinProb_sex.xlsx")
writexl::write_xlsx(data_results_man_sex,"data_results_man_sex.xlsx")
data_results_knn_sex_sign <- data_results_knn_sex %>%
  filter(significant)
data_results_MinProb_sex_sign <- data_results_MinProb_sex %>%
  filter(significant)
data_results_man_sex_sign <- data_results_man_sex %>%
  filter(significant)
writexl::write_xlsx(data_results_knn_sex_sign,"data_results_knn_sex_sign.xlsx")
writexl::write_xlsx(data_results_MinProb_sex_sign,"data_results_MinProb_sex_sign.xlsx")
writexl::write_xlsx(data_results_man_sex_sign,"data_results_man_sex_sign.xlsx")
data_results_knn_age <- get_results(DEP_knn_age)
data_results_knn_age <- data_results_knn_age %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_MinProb_age <- get_results(DEP_MinProb_age)
data_results_MinProb_age <- data_results_MinProb_age %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_man_age <- get_results(DEP_man_age)
data_results_man_age <- data_results_man_age %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
writexl::write_xlsx(data_results_knn_age,"data_results_knn_age.xlsx")
writexl::write_xlsx(data_results_MinProb_age,"data_results_MinProb_age.xlsx")
writexl::write_xlsx(data_results_man_age,"data_results_man_age.xlsx")
data_results_knn_age_sign <- data_results_knn_age %>%
  filter(significant)
data_results_MinProb_age_sign <- data_results_MinProb_age %>%
  filter(significant)
data_results_man_age_sign <- data_results_man_age %>%
  filter(significant)
writexl::write_xlsx(data_results_knn_age_sign,"data_results_knn_age_sign.xlsx")
writexl::write_xlsx(data_results_MinProb_age_sign,"data_results_MinProb_age_sign.xlsx")
writexl::write_xlsx(data_results_man_age_sign,"data_results_man_age_sign.xlsx")
