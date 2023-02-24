#### Processing proteomics data with DEP package 

# identify potential contaminant 
proteins_to_remove <- proteinGroups_samples %>% 
  filter(Reverse == "+" | `Potential contaminant` == "+" | `Only identified by site` == "+")

# core facility final list contains potential contaminant and reverse proteins
proteins_to_remove_core <- proteins_to_remove[proteins_to_remove$`Protein IDs` %in% CSF_core_data_processing$ID,]

# remove potential contaminant
data_intensities_proteins <- proteinGroups_samples %>%
  filter(!`Protein IDs` %in% proteins_to_remove$`Protein IDs`)

# Get only the intensities 
data_intensities_proteins <- data_intensities_proteins[,c(1,which(grepl("Intensity",colnames(proteinGroups_samples))))]
names(data_intensities_proteins)[c(2:ncol(data_intensities_proteins))] <- gsub(" ", "_", 
                                                            names(data_intensities_proteins)[c(2:ncol(data_intensities_proteins))])

# check if we have duplicated names of proteins: alles gut
data_intensities_proteins$`Protein IDs` %>% duplicated() %>% any()

# make a summarized object based on the design matrix
ExpDesign_APS <- design_APS %>%
  dplyr::left_join(data_age %>% 
                     dplyr::select(Pseudonym,`Alter bei LP / Age at LP`) %>%
                     dplyr::rename(`Tube ID` = Pseudonym,
                                   age = `Alter bei LP / Age at LP`)) %>%
  dplyr::left_join(data_NfL %>%
                     dplyr::select(`Sample Barcode`,`Replicate Conc.`) %>%
                     dplyr::rename(`Tube ID` = `Sample Barcode`,
                                   Nfl_conc = `Replicate Conc.`)) %>%
  dplyr::select(Name,Condition,Sex,Disease,Pathy,age,Nfl_conc,`Tube ID`) %>%
  dplyr::group_by(Pathy) %>%
  dplyr::mutate(replicate = seq(1:length(Pathy))) %>%
  dplyr::rename(label = Name,
         condition = Pathy) 

ExpDesign_APS[is.na(ExpDesign_APS$age) | is.na(ExpDesign_APS$Nfl_conc),c("age","Nfl_conc"),] <- tibble(age = c(73,NA,75,30,75,47,69,77),
                                                                                                      Nfl_conc = c(9766.004,
                                                                                                                   7033,
                                                                                                                   2358,
                                                                                                                   NA,
                                                                                                                   2019,
                                                                                                                   NA,
                                                                                                                   NA,
                                                                                                                   NA))

# summarised object different for Nfl_conc and age due to missing data in each one
ExpDesign_APS_age <- ExpDesign_APS %>%
  dplyr::filter(!is.na(age)) # 140 IDs

ExpDesign_APS_Nfl <- ExpDesign_APS %>%
  dplyr::filter(!is.na(Nfl_conc)) # 137 IDs

ExpDesign_APS <- ExpDesign_APS %>% 
  dplyr::select(-c("age","Nfl_conc",`Tube ID`))
  
ExpDesign_APS[142,] <- data.table::data.table("Intensity_CSF135", "Ctrl", "f",
                                              "DLB","DLB",13,2,"",45) # 142 IDs

# transform names from data to the same as in the ExpDesign
#colnames(data_intensities_proteins) <- str_replace(colnames(data_intensities_proteins),"LFQ_i","I")
colnames(data_intensities_proteins)[1] <- "Protein.IDs"

# intensity columns
Intensity_columns <- grep("Intensity_",colnames(data_intensities_proteins))

# Summarised object
data_intensities_proteins$ID <- data_intensities_proteins$Protein.IDs
# total of 895 proteins and 142 samples 
data_proteins_SE <- make_se(make_unique(data_intensities_proteins,"Protein.IDs","ID"),
                            Intensity_columns,
                            ExpDesign_APS)

data_proteins_SE_age <-  make_se(make_unique(data_intensities_proteins,"Protein.IDs","ID"),
                                 Intensity_columns,
                                 ExpDesign_APS_age)
data_proteins_SE_Nfl <-  make_se(make_unique(data_intensities_proteins,"Protein.IDs","ID"),
                                 Intensity_columns,
                                 ExpDesign_APS_Nfl)

# check if there is missing values 
plot_frequency(data_proteins_SE)

# looks like the dataset contains missing values which need to be imputed

# lets check the proteins distribution per sample 
df <- assay(data_proteins_SE) %>% data.frame() %>% rownames_to_column() %>% 
  gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin), 
                                                    0, 1))
stat <- df %>% group_by(ID) %>% summarize(n = n(), sum = sum(bin)) %>% 
  left_join(., data.frame(colData(data_proteins_SE)), by = "ID")
plot_numbers(data_proteins_SE)
                                                                                                                                                           miss_val)
# is possible to see one sample with basically no proteins: CTRL 135, only has 3 proteins
# so we decide to remove this sample from the data
data_proteins_filter_new <- data_proteins_SE[,-which(data_proteins_SE$label == "Intensity_CSF135")]
df <- assay(data_proteins_filter_new) %>% data.frame() %>% rownames_to_column() %>% 
  gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin), 
                                                    0, 1))
stat <- df %>% group_by(ID) %>% summarize(n = n(), sum = sum(bin)) %>% 
  left_join(., data.frame(colData(data_proteins_filter_new)), by = "ID")
p <- ggplot(stat, aes(x = ID, y = sum, fill = condition)) + 
  geom_col() + geom_hline(yintercept = unique(stat$n)) + 
  labs(title = "Proteins per sample after Ctrl 135 removal", x = "", y = "Number of proteins") + 
  theme_DEP2()

# check protein coverage after outlier removal
plot_coverage(data_proteins_filter_new)

# first check if we need to filter out proteins with too many missing values
bin_data <- assay(data_proteins_filter_new)
idx <- is.na(assay(data_proteins_filter_new))
bin_data[!idx] <- 1
bin_data[idx] <- 0
keep <- bin_data %>% data.frame() %>% rownames_to_column() %>% 
  gather(ID, value, -rowname) %>% 
  left_join(., data.frame(colData(data_proteins_filter_new)),  by = "ID") %>% 
  group_by(rowname, condition) %>% 
  summarize(miss_val = n() - sum(value),
            total_val = n(),
            thr_val = n()*0.2) %>% 
  filter(miss_val < thr_val) %>% 
  select(rowname,condition,miss_val) %>%
  spread(condition,miss_val) 

# after filtering: 575 proteins and 142 samples
data_proteins_filter <- data_proteins_filter_new[keep$rowname,]
data_proteins_filter_age <- data_proteins_SE_age[keep$rowname,]
data_proteins_filter_Nfl <- data_proteins_SE_Nfl[keep$rowname,]

# normalization 
data_proteins_norm <- normalize_vsn(data_proteins_filter)
data_proteins_norm_age <- normalize_vsn(data_proteins_filter_age)
data_proteins_norm_Nfl <- normalize_vsn(data_proteins_filter_Nfl)

meanSdPlot(data_proteins_norm)
plot_normalization(data_proteins_filter,data_proteins_norm)

## input data for missing values
# first check proteins with missing values
plot_missval(data_proteins_filter)
# -> highly biased by specific samples, and does not look correlated with type of biological condition

# -> check if missing values are biased to lower intense proteins
plot_detect(data_proteins_filter)
# -> conclusion: no, the proteins with missing values have on avg low intensities

# imputation with several methods: MinProb, MAN, KNN
data_proteins_imp_MinProb <- impute(data_proteins_norm, fun = "MinProb", q=0.01)
data_proteins_imp_man <- impute(data_proteins_norm, fun = "man", shift = 1.8, scale = 0.3)
data_proteins_imp_knn <- impute(data_proteins_norm, fun = "knn", rowmax = 0.9)
data_proteins_imp_MinProb_age <- impute(data_proteins_norm_age, fun = "MinProb", q=0.01)
data_proteins_imp_man_age <- impute(data_proteins_norm_age, fun = "man", shift = 1.8, scale = 0.3)
data_proteins_imp_knn_age <- impute(data_proteins_norm_age, fun = "knn", rowmax = 0.9)
data_proteins_imp_MinProb_Nfl <- impute(data_proteins_norm_Nfl, fun = "MinProb", q=0.01)
data_proteins_imp_man_Nfl <- impute(data_proteins_norm_Nfl, fun = "man", shift = 1.8, scale = 0.3)
data_proteins_imp_knn_Nfl <- impute(data_proteins_norm_Nfl, fun = "knn", rowmax = 0.9)

# plot the several imputations
plot_imputation(data_proteins_norm,
                data_proteins_imp_MinProb,
                data_proteins_imp_man,
                data_proteins_imp_knn)

# differential analysis
data_diff_knn <- test_diff(data_proteins_imp_knn,type = "all")
data_diff_MinProb <- test_diff(data_proteins_imp_MinProb,type = "all")
data_diff_man <- test_diff(data_proteins_imp_man,type = "all")

# add significance 
DEP_knn <- add_rejections(data_diff_knn,alpha = 0.1,lfc = log2(1))
DEP_MinProb <- add_rejections(data_diff_MinProb,alpha = 0.1,lfc = log2(1))
DEP_man <- add_rejections(data_diff_man,alpha = 0.1,lfc = log2(1))

# PCA with different imputation methods version 1
pca_knn <- PCA_DEP_new(DEP_knn, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition")
nudge <- position_nudge(y = 1)
pca_knn + geom_text(aes(label = rowname), position = nudge)
pca_MinProb <- PCA_DEP_new(DEP_MinProb, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition")
pca_MinProb + geom_text(aes(label = rowname), position = nudge)
pca_man <- PCA_DEP_new(DEP_man, x = 1, y = 2, n = 500, point_size = 4, indicate = "condition")
pca_man + geom_text(aes(label = rowname), position = nudge)

# Correlation matrix with different imputation methods
plot_cor(DEP_knn, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_cor(DEP_MinProb, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_cor(DEP_man, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# heatmap of top 100 most variant proteins 
heatmap_DEP(DEP_knn,"most variant")
heatmap_DEP(DEP_MinProb,"most variant")
heatmap_DEP(DEP_man,"most variant")

# heatmap of significant proteins 
heatmap_DEP(DEP_knn,"significant proteins")
heatmap_DEP(DEP_MinProb,"significant proteins")
heatmap_DEP(DEP_man,"significant proteins")

# heatmap of top 100 most variant proteins clustered by CTR, PD, Synu
heatmap_DEP(DEP_knn,"most variant",cluster = T)
heatmap_DEP(DEP_MinProb,"most variant", cluster = T)
heatmap_DEP(DEP_man,"most variant", cluster = T)

# heatmap of significant proteins by CTR, PD, Synu
heatmap_DEP(DEP_knn,"significant proteins", cluster = T)
heatmap_DEP(DEP_MinProb,"significant proteins", cluster = T)
heatmap_DEP(DEP_man,"significant proteins", cluster = T)

# -> Ctr vs PD
volcano_DEP(DEP_knn, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man, contrast = "Ctr_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Synucleinopathy vs Ctr
volcano_DEP(DEP_knn, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man, contrast = "Synucleinopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Synucleinopathy vs PD
volcano_DEP(DEP_knn, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man, contrast = "Synucleinopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Synucleinopathy vs Tauopathy
volcano_DEP(DEP_knn, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man, contrast = "Synucleinopathy_vs_Tauopathy", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Tauopathy vs Ctr
volcano_DEP(DEP_knn, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man, contrast = "Tauopathy_vs_Ctr", label_size = 2, add_names = TRUE,adjusted = TRUE)

# -> Tauopathy vs PD
volcano_DEP(DEP_knn, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_MinProb, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)
volcano_DEP(DEP_man, contrast = "Tauopathy_vs_PD", label_size = 2, add_names = TRUE,adjusted = TRUE)

# get results
data_results_knn <- get_results(DEP_knn)
data_results_knn <- data_results_knn %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_MinProb <- get_results(DEP_MinProb)
data_results_MinProb <- data_results_MinProb %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
data_results_man <- get_results(DEP_man)
data_results_man <- data_results_man %>%
  rename(PD_vs_Ctr_ratio = Ctr_vs_PD_ratio,
         PD_vs_Ctr_p.val = Ctr_vs_PD_p.val,
         PD_vs_Ctr_p.adj = Ctr_vs_PD_p.adj,
         PD_vs_Ctr_significant = Ctr_vs_PD_significant) %>%
  mutate(PD_vs_Ctr_ratio = (- PD_vs_Ctr_ratio))
writexl::write_xlsx(data_results_knn,"data_results_knn.xlsx")
writexl::write_xlsx(data_results_MinProb,"data_results_MinProb.xlsx")
writexl::write_xlsx(data_results_man,"data_results_man.xlsx")

# only significant
data_results_knn_sign <- data_results_knn %>%
  filter(significant)
data_results_MinProb_sign <- data_results_MinProb  %>%
  filter(significant)
data_results_man_sign <- data_results_man  %>%
  filter(significant)
writexl::write_xlsx(data_results_knn_sign,"data_results_knn_sign.xlsx")
writexl::write_xlsx(data_results_MinProb_sign,"data_results_MinProb_sign.xlsx")
writexl::write_xlsx(data_results_man_sign,"data_results_man_sign.xlsx")

### auxiliar functions
# PCA with density plots in x- and y-axis 
PCA_DEP_new <- function(dep, x = 1, y = 2, indicate = c("condition", "replicate"), 
                        label = FALSE, 
                        density = TRUE, 
                        n = 500, point_size = 4, label_size = 3, plot = TRUE){
  if (is.integer(x)) 
    x <- as.numeric(x)
  if (is.integer(y)) 
    y <- as.numeric(y)
  if (is.integer(n)) 
    n <- as.numeric(n)
  if (is.integer(point_size)) 
    point_size <- as.numeric(point_size)
  if (is.integer(label_size)) 
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.numeric(x), length(x) == 1, is.numeric(y), length(y) == 
                            1, is.numeric(n), length(n) == 1, is.character(indicate), 
                          is.logical(label), is.numeric(point_size), length(point_size) == 
                            1, is.numeric(label_size), length(label_size) == 
                            1, is.logical(plot), length(plot) == 1)
  if (x > ncol(dep) | y > ncol(dep)) {
    stop(paste0("'x' and/or 'y' arguments are not valid\n", 
                "Run plot_pca() with 'x' and 'y' <= ", ncol(dep), 
                "."), call. = FALSE)
  }
  if (n > nrow(dep)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_pca() with 'n' <= ", 
                nrow(dep), "."), call. = FALSE)
  }
  columns <- colnames(colData(dep))
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_pca() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"), 
                  "' column(s) is/are not present in ", deparse(substitute(dep)), 
                  ".\nValid columns are: '", paste(columns, collapse = "', '"), 
                  "'."), call. = FALSE)
    }
  }
  var <- apply(assay(dep), 1, sd)
  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)], 
  ]
  pca <- prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame() %>% rownames_to_column() %>% 
    left_join(., data.frame(colData(dep)), by = c(rowname = "ID"))
  percent <- round(100 * pca$sdev^2/sum(pca$sdev^2), 1)
  for (feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }
  p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC",  y)))) + 
    labs(title = paste0("PCA plot - top ", n, " variable proteins"), 
         x = paste0("PC", x, ": ", percent[x], "%"), y = paste0("PC", y, ": ", percent[y], "%")) + 
    coord_fixed() + theme_DEP1() +
    scale_color_manual(values = c("#1887ab","#eec76b","#93ae55","#d34467")) +
    #scale_color_brewer(palette="Royal2") +
    theme_minimal()
  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]]), 
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], 
                            shape = pca_df[[indicate[2]]]), size = point_size) + 
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], 
                            shape = pca_df[[indicate[2]]]), size = point_size) + 
      facet_wrap(~pca_df[[indicate[3]]])
    labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (density) {
    p <- ggMarginal(p, groupFill = TRUE)
    
  }
  if (plot) {
    return(p)
  }
  else {
    df <- pca_df %>% select(rowname, paste0("PC", c(x, y)), 
                            match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

# row normalization based on z-score
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

# heatmap of DEP data for top most variant proteins or only significant 
heatmap_DEP <- function(DEP_data,type,cluster = F){
  # default value is to not cluster by CTR, PD and Synu classes 
  if(type == "most variant"){ # top variant proteins
    TOPproteins <- head(order(rowVars(assay(DEP_data)),decreasing = T),100)
    col_title <- "Top 100 most variant proteins"
    row_names = F
  }else{ # only significant 
    TOPproteins_DEP <- which(rowData(DEP_data, use.names = FALSE)$significant,arr.ind = T)
    TOPproteins_Perseus <- which(rowData(DEP_data)$Protein.IDs %in% DE_FDR_0.1_Perseus_sign$`Protein IDs`,arr.ind = T)
    TOPproteins <- unique(c(TOPproteins_DEP,TOPproteins_Perseus))
    col_title <- "Significant proteins"
    row_names = T
  }
  matrix.TOPproteins <- assay(DEP_data)[TOPproteins,]
  matrix.TOPproteins <- matrix.TOPproteins - rowMeans(matrix.TOPproteins)
  matrix.TOPprotein.score <- scale_rows(matrix.TOPproteins)
  annotation.protein <- as.data.frame(colData(DEP_data)[c("condition","Disease","Sex","age","Nfl_conc")])
  protein_names <- rownames(matrix.TOPproteins) 
  rownames(matrix.TOPproteins) <- lapply(strsplit(protein_names,"\\|"), function(x){
    use <- x[3]
    use_aux <- unlist(strsplit(use,"_"))
    use_aux[1]})
  col.protein = list(condition = c("Ctr" = "#1887ab","PD" = "#eec76b", "Synucleinopathy" = "#93ae55", "Tauopathy" = "#d34467"),
                     Disease = c("CBD" = "#B098F2", "Ctr" = "#1887ab", "DLB" = "#ABB8A7", "MSA" = "#EAE29B", "PD" = "#eec76b", "PSP" = "#D69E9A"),
                     Sex = c("f" = "#E0A69C", "m" = "#9CD6E0"),
                     age = circlize::colorRamp2(c(30,90),c("#F7FBFF","#08519C")),
                     Nfl_conc = circlize::colorRamp2(c(714,41120),c("#FCFBFD","#3F007D")))
  if(cluster){
    annotation.protein$condition <- factor(annotation.protein$condition, 
                                           levels = c("Ctr","PD","Synucleinopathy","Tauopathy"))
    ComplexHeatmap::Heatmap(matrix.TOPproteins,name = "z-score",
                            col = circlize::colorRamp2(c(-4,0,4),c("#9C0E47","#F5F5F5","#409E78")),
                            top_annotation = HeatmapAnnotation(df = annotation.protein,
                                                               col = col.protein,
                                                               annotation_name_gp= gpar(fontsize = 17),
                                                               annotation_legend_param = list(labels_gp = gpar(fontsize = 17),
                                                                                              title_gp = gpar(fontsize=18, fontface="bold"))),
                            heatmap_legend_param=list(legend_width=unit(8,"cm"),
                                                      title_gp=gpar(fontsize=18, fontface="bold"),
                                                      labels_gp = gpar(fontsize = 17)),
                            show_column_names = FALSE,show_row_names = row_names,
                            show_column_dend = F, show_row_dend = T,
                            cluster_columns = F,
                            column_names_gp = grid::gpar(fontsize = 17),
                            column_dend_gp = gpar(fontsize = 17),
                            column_title_gp = gpar(fontsize = 20,fontface = "bold"),
                            row_title_gp = gpar(fontsize = 17),
                            row_names_gp = gpar(fontsize = 15),
                            column_order = sort(colnames(matrix.TOPproteins)),
                            column_title = col_title)
  } 
  else{
  ComplexHeatmap::Heatmap(matrix.TOPproteins,name = "z-score",
                          col = circlize::colorRamp2(c(-4,0,4),c("#9C0E47","#F5F5F5","#409E78")),
                          top_annotation = HeatmapAnnotation(df = annotation.protein,
                                                             col = col.protein,
                                                             annotation_name_gp= gpar(fontsize = 17),
                                                             annotation_legend_param = list(labels_gp = gpar(fontsize = 17),
                                                                                            title_gp = gpar(fontsize=18, fontface="bold"))),
                          heatmap_legend_param=list(legend_width=unit(8,"cm"),
                                                    title_gp=gpar(fontsize=18, fontface="bold"),
                                                    labels_gp = gpar(fontsize = 17)),
                          show_column_names = FALSE,show_row_names = row_names,
                          show_column_dend = T, show_row_dend = T,
                          column_names_gp = grid::gpar(fontsize = 17),
                          column_dend_gp = gpar(fontsize = 17),
                          column_title_gp = gpar(fontsize = 20,fontface = "bold"),
                          row_title_gp = gpar(fontsize = 17),
                          row_names_gp = gpar(fontsize = 15),
                          cluster_columns = TRUE,
                          column_title = col_title)
    }
}

# volcano plots (details to save pdf 10x8, hline log10 + 0.15, limits NA or 1.2)
volcano_DEP <- function (dep, contrast, label_size = 3, add_names = TRUE, adjusted = FALSE, 
                         plot = TRUE) 
{
  if (is.integer(label_size)) 
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), 
                          is.character(contrast), length(contrast) == 1, is.numeric(label_size), 
                          length(label_size) == 1, is.logical(add_names), length(add_names) == 
                            1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot), 
                          length(plot) == 1)
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."), 
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), 
         call. = FALSE)
  }
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '", 
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."), 
         call. = FALSE)
  }
  if (length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) == 
      0) {
    valid_cntrsts <- row_data %>% data.frame() %>% select(ends_with("_diff")) %>% 
      colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '", 
                                paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n", 
         valid_cntrsts_msg, call. = FALSE)
  }
  diff <- grep(paste(contrast, "_diff", sep = ""), colnames(row_data))
  if (adjusted) {
    p_values <- grep(paste(contrast, "_p.adj", sep = ""), 
                     colnames(row_data))
  }
  else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""), 
                     colnames(row_data))
  }
  signif <- grep(paste(contrast, "_significant", sep = ""), 
                 colnames(row_data))
  df <- data.frame(x = row_data[, diff], y = -log10(row_data[, p_values]), 
                   significant = row_data[, signif], name = row_data$name) %>% 
    filter(!is.na(significant)) %>% arrange(significant) %>%
    mutate(name = lapply(strsplit(name,"\\|"), function(x){
      use <- x[3]
      use_aux <- unlist(strsplit(use,"_"))
      use_aux[1]
    }))
  if(contrast == "Ctr_vs_PD"){
    contrast = "PD_vs_Ctr"
    df <- df %>%
      mutate(x = (- x))
  }
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  df <- df %>%
    mutate(gene_type = case_when(x >= 0 & significant ~ "up",
                                 x <= (0) & significant ~ "down",
                                 TRUE ~ "ns")) 
  cols <- c("up" = "#d4552b", "down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 0.7, "down" = 0.7, "ns" = 0.5)
  p <- ggplot(data = df, aes(x,y)) + 
    geom_point(aes(colour = gene_type), 
               alpha = 0.5, 
               shape = 16,
               size = 3) + 
    geom_hline(yintercept = log10(10),
               linetype = "dashed") + 
    geom_vline(xintercept = 0,linetype = "dashed") +
    geom_point(data = filter(df, significant),
               aes(colour = gene_type), 
               alpha = 0.5, 
               shape = 16,
               size = 4) + 
    annotate(geom="text", x=-1.9, y=log10(10) + 0.15, label="FDR = 10%",size = 5) +
    geom_label_repel(data = filter(df, significant),
                     aes(label = name),
                     force = 0.3,
                     nudge_x = 0.5,
                     nudge_y = 0.1,
                     size = 5) +
    scale_colour_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(-2.5, 2.5)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-0.2, NA)) +
    labs(title = contrast,
         x = "log2(fold change)",
         y = expression(-log[10] ~ "(adjusted p-value)"),
         colour = "Differential \nExpression") +
    theme_classic() + # Select theme with a white background  
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold",
                                    hjust = 0.5),
          text = element_text(size = 18))
  if (plot) {
    return(p)
  }
  else {
    df <- df %>% select(name, x, y, significant) %>% arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", 
                                  "p_value_-log10")
    if (adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}