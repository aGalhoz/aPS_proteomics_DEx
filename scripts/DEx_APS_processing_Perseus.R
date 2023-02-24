##### Perseus analysis

# data FDR = 10%
DE_FDR0.1_CTR_PD <- read_delim("Perseus analyses/output datasets/DE_FDR0.1_CTR_PD.txt", 
                                delim = "\t", escape_double = FALSE, 
                                locale = locale(decimal_mark = ","), 
                                comment = "#", trim_ws = TRUE)
DE_FDR0.1_CTR_Synu <- read_delim("Perseus analyses/output datasets/DE_FDR0.1_CTR_Synu.txt", 
                               delim = "\t", escape_double = FALSE, 
                               locale = locale(decimal_mark = ","), 
                               comment = "#", trim_ws = TRUE)
DE_FDR0.1_CTR_Tau <- read_delim("Perseus analyses/output datasets/DE_FDR0.1_CTR_Tau.txt", 
                               delim = "\t", escape_double = FALSE, 
                               locale = locale(decimal_mark = ","), 
                               comment = "#", trim_ws = TRUE)
DE_FDR0.1_PD_Synu <- read_delim("Perseus analyses/output datasets/DE_FDR0.1_PD_Synu.txt", 
                                 delim = "\t", escape_double = FALSE, 
                                 locale = locale(decimal_mark = ","), 
                                 comment = "#", trim_ws = TRUE)
DE_FDR0.1_PD_Tau <- read_delim("Perseus analyses/output datasets/DE_FDR0.1_PD_Tau.txt", 
                                delim = "\t", escape_double = FALSE, 
                                locale = locale(decimal_mark = ","), 
                                comment = "#", trim_ws = TRUE)
DE_FDR0.1_Synu_Tau <- read_delim("Perseus analyses/output datasets/DE_FDR0.1_Synu_Tau.txt", 
                               delim = "\t", escape_double = FALSE, 
                               locale = locale(decimal_mark = ","), 
                               comment = "#", trim_ws = TRUE)

# Final lists 
# FDR = 10%
DE_FDR_0.1_Perseus <- DE_FDR0.1_CTR_PD %>%
dplyr::select(-c("Potential contaminant",
                 "Only identified by site","Reverse",
                 "Student's T-test significant")) %>%
  left_join(DE_FDR0.1_CTR_Synu %>%
              dplyr::select(-c("Potential contaminant",
                               "Only identified by site","Reverse",
                               "Student's T-test significant"))) %>%
  left_join(DE_FDR0.1_CTR_Tau %>%
              dplyr::select(-c("Potential contaminant",
                               "Only identified by site","Reverse",
                               "Student's T-test significant"))) %>%
  left_join(DE_FDR0.1_PD_Synu %>%
              dplyr::select(-c("Potential contaminant",
                               "Only identified by site","Reverse",
                               "Student's T-test significant"))) %>%
  left_join(DE_FDR0.1_PD_Tau %>%
              dplyr::select(-c("Potential contaminant",
                               "Only identified by site","Reverse",
                               "Student's T-test significant"))) %>%
  left_join(DE_FDR0.1_Synu_Tau %>%
              dplyr::select(-c("Potential contaminant",
                               "Only identified by site","Reverse",
                               "Student's T-test significant"))) 
DE_FDR_0.1_Perseus <- DE_FDR_0.1_Perseus %>%
  mutate(`Student's T-test Difference Ctr_Synucleinopathy` = (-as.numeric(DE_FDR_0.1_Perseus$`Student's T-test Difference Ctr_Synucleinopathy`)),
         `Student's T-test Difference Ctr_Tauopathy` =  (- as.numeric(DE_FDR_0.1_Perseus$`Student's T-test Difference Ctr_Tauopathy`)),
         `Student's T-test Difference PD_Synucleinopathy` =  (- as.numeric(DE_FDR_0.1_Perseus$`Student's T-test Difference PD_Synucleinopathy`)),
         `Student's T-test Difference PD_Tauopathy` = (- as.numeric(DE_FDR_0.1_Perseus$`Student's T-test Difference PD_Tauopathy`)))

#writexl::write_xlsx(DE_FDR_0.05_Perseus,"DE_FDR_5%_Perseus.xlsx")

DE_FDR_0.1_Perseus_aux <- DE_FDR_0.1_Perseus[,c(1:142,144:147,152:155,157:160,162:165,167:170,172:175)] 
DE_FDR_0.1_Perseus_aux <- apply(DE_FDR_0.1_Perseus_aux,2,as.numeric)
DE_FDR_0.1_Perseus[,c(1:142,144:147,152:155,157:160,162:165,167:170,172:175)] <- DE_FDR_0.1_Perseus_aux
DE_FDR_0.1_Perseus <- DE_FDR_0.1_Perseus %>%
  rename(`Student's T-test Significant PD_Ctr` = `Student's T-test Significant Ctr_PD`,
         `-Log Student's T-test p-value PD_Ctr` = `-Log Student's T-test p-value Ctr_PD`,
         `Student's T-test q-value PD_Ctr` = `Student's T-test q-value Ctr_PD`,
         `Student's T-test Difference PD_Ctr` = `Student's T-test Difference Ctr_PD`,
         `Student's T-test Test statistic PD_Ctr` = `Student's T-test Test statistic Ctr_PD`) %>%
  mutate(`Student's T-test Difference PD_Ctr` = (- `Student's T-test Difference PD_Ctr`))

writexl::write_xlsx(DE_FDR_0.1_Perseus,"DE_FDR_10%_Perseus.xlsx")

# only significant
DE_FDR_0.1_Perseus_temp <- DE_FDR_0.1_Perseus[,143:ncol(DE_FDR_0.1_Perseus)]
DE_FDR_0.1_Perseus_sign <- DE_FDR_0.1_Perseus_temp %>%
  filter(`Student's T-test Significant Tauopathy_PD` == "+" |
           `Student's T-test Significant PD_Ctr` == "+" |
           `Student's T-test Significant Tauopathy_Ctr` == "+" |
           `Student's T-test Significant Synucleinopathy_Tauopathy` == "+" |
           `Student's T-test Significant Synucleinopathy_PD` == "+" |
           `Student's T-test Significant Synucleinopathy_Ctr` == "+")
DE_FDR_0.1_Perseus_PD_Ctr <- DE_FDR_0.1_Perseus_temp %>%
  filter(`Student's T-test Significant PD_Ctr` == "+")
DE_FDR_0.1_Perseus_Synucleinopathy_PD <- DE_FDR_0.1_Perseus_temp %>%
  filter(`Student's T-test Significant Synucleinopathy_PD` == "+")
DE_FDR_0.1_Perseus_Synucleinopathy_Ctr <- DE_FDR_0.1_Perseus_temp %>%
  filter(`Student's T-test Significant Synucleinopathy_Ctr` == "+")

writexl::write_xlsx(DE_FDR_0.1_Perseus_sign,"DE_FDR_10%_Perseus_significant.xlsx")
writexl::write_xlsx(DE_FDR_0.1_Perseus_PD_Ctr,"DE_FDR_10%_Perseus_significant_PD_Ctr.xlsx")
writexl::write_xlsx(DE_FDR_0.1_Perseus_Synucleinopathy_PD,"DE_FDR_10%_Perseus_significant_Synu_PD.xlsx")
writexl::write_xlsx(DE_FDR_0.1_Perseus_Synucleinopathy_Ctr,"DE_FDR_10%_Perseus_significant_Synu_Ctr.xlsx")

# volcano plots 
# -> Ctr vs PD
data_use <- DE_FDR_0.1_Perseus[,c(143:148)]
volcano_plot_perseus(data_use,"PD_Ctr")

# -> Synucleinopathy vs Ctr
data_use <- DE_FDR_0.1_Perseus[,c(148,151:155)]
volcano_plot_perseus(data_use,"Synucleinopathy_Ctr")

# -> Synucleinopathy vs PD
data_use <- DE_FDR_0.1_Perseus[,c(148,161:165)]
volcano_plot_perseus(data_use,"Synucleinopathy_PD")

# -> Synucleinopathy vs Tauopathy
data_use <- DE_FDR_0.1_Perseus[,c(148,171:175)]
volcano_plot_perseus(data_use,"Synucleinopathy_Tauopathy")

# -> Tauopathy vs Ctr
data_use <- DE_FDR_0.1_Perseus[,c(148,156:160)]
volcano_plot_perseus(data_use,"Tauopathy_Ctr")

# -> Tauopathy vs PD
data_use <- DE_FDR_0.1_Perseus[,c(148,166:170)]
volcano_plot_perseus(data_use,"Tauopathy_PD")

## auxiliar functions
volcano_plot_perseus <- function(data,name_title){
  log2fc <- grep("Student's T-test Difference",colnames(data))
  padj <- grep("Student's T-test q-value",colnames(data))
  sign <- grep("Student's T-test Significant",colnames(data))
  print(data[,padj])
  df <- data.frame(x = data[,log2fc], 
                   y = -log10(data[,padj]),
                   significant = data[,sign],
                   name = data$`Protein IDs`) %>%
    mutate(name = lapply(strsplit(name,"\\|"), function(x){
      use <- x[3]
      use_aux <- unlist(strsplit(use,"_"))
      use_aux[1]
    }))
  names(df) <- c("x","y","significant","name")
  df <- df %>%
    mutate(gene_type = case_when(x >= 0 & significant == "+" ~ "up",
                                 x <= (0) & significant == "+" ~ "down",
                                 TRUE ~ "ns")) 
  cols <- c("up" = "#d4552b", "down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 0.7, "down" = 0.7, "ns" = 0.5)
  ggplot(data = df, aes(x,y)) + 
    geom_point(aes(colour = gene_type), 
               alpha = 0.5, 
               shape = 16,
               size = 3) + 
    geom_hline(yintercept = log10(10),
               linetype = "dashed") + 
    geom_vline(xintercept = 0,linetype = "dashed") +
    geom_point(data = filter(df, significant == "+"),
               aes(colour = gene_type), 
               alpha = 0.5, 
               shape = 16,
               size = 4) + 
    annotate(geom="text", x=-1.9, y=log10(10) + 0.15, label="FDR = 10%",size = 5) +
    geom_label_repel(data = filter(df, significant == "+"),
                     aes(label = name),
                     force = 1,
                     nudge_y = 0.2,
                     nudge_x = 0.1,
                     max.overlaps = 30,
                     size = 5) +
    scale_colour_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(-3, 3)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-0.1, NA)) +
    labs(title = name_title,
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
}
