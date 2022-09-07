#### Script to merge and compare DEP and Perseus results

# merge datasets
DE_DEP_Perseus_knn <- comb_DEP_Perseus(data_results_knn)
DE_DEP_Perseus_MinProb <- comb_DEP_Perseus(data_results_MinProb)
DE_DEP_Perseus_man <- comb_DEP_Perseus(data_results_man)

## plot DEP vs Perseus for each comparison
# Ctr vs PD
plot_comparison_methods(DE_DEP_Perseus_knn,"PD_Ctr")
plot_comparison_methods(DE_DEP_Perseus_MinProb,"PD_Ctr")
plot_comparison_methods(DE_DEP_Perseus_man,"PD_Ctr")

# Synucleinopathy vs Ctr
plot_comparison_methods(DE_DEP_Perseus_knn,"Synucleinopathy_Ctr")
plot_comparison_methods(DE_DEP_Perseus_MinProb,"Synucleinopathy_Ctr")
plot_comparison_methods(DE_DEP_Perseus_man,"Synucleinopathy_Ctr")

# Synucleinopathy vs PD
plot_comparison_methods(DE_DEP_Perseus_knn,"Synucleinopathy_PD")
plot_comparison_methods(DE_DEP_Perseus_MinProb,"Synucleinopathy_PD")
plot_comparison_methods(DE_DEP_Perseus_man,"Synucleinopathy_PD")

# Tauopathy vs Ctr
plot_comparison_methods(DE_DEP_Perseus_knn,"Tauopathy_Ctr")
plot_comparison_methods(DE_DEP_Perseus_MinProb,"Tauopathy_Ctr")
plot_comparison_methods(DE_DEP_Perseus_man,"Tauopathy_Ctr")

# Tauopathy vs PD
plot_comparison_methods(DE_DEP_Perseus_knn,"Tauopathy_PD")
plot_comparison_methods(DE_DEP_Perseus_MinProb,"Tauopathy_PD")
plot_comparison_methods(DE_DEP_Perseus_man,"Tauopathy_PD")

# Synucleinopathy vs Tauopathy
plot_comparison_methods(DE_DEP_Perseus_knn,"Synucleinopathy_Tauopathy")
plot_comparison_methods(DE_DEP_Perseus_MinProb,"Synucleinopathy_Tauopathy")
plot_comparison_methods(DE_DEP_Perseus_man,"Synucleinopathy_Tauopathy")

### auxiliar functions
# combination of DEP and Perseus datasets
comb_DEP_Perseus <- function(data_DEP){
  data_DEP_temp <- data_DEP %>%
    dplyr::select(-grep("Intensity",colnames(.))) %>%
    rename("Protein IDs" = name)
  data_DEP_Perseus <- DE_FDR_0.1_Perseus %>%
    full_join(data_DEP_temp) %>%
    dplyr::mutate(PD_Ctr_signlogp_DEP = sign(PD_vs_Ctr_ratio) * (-log10(PD_vs_Ctr_p.adj)),
                  Synucleinopathy_Ctr_signlogp_DEP = sign(Synucleinopathy_vs_Ctr_ratio) * (-log10(Synucleinopathy_vs_Ctr_p.adj)),
                  Synucleinopathy_PD_signlogp_DEP = sign(Synucleinopathy_vs_PD_ratio) * (-log10(Synucleinopathy_vs_PD_p.adj)),
                  Synucleinopathy_Tauopathy_signlogp_DEP = sign(Synucleinopathy_vs_Tauopathy_ratio) * (-log10(Synucleinopathy_vs_Tauopathy_p.adj)),
                  Tauopathy_Ctr_signlogp_DEP = sign(Tauopathy_vs_Ctr_ratio) * (-log10(Tauopathy_vs_Ctr_p.adj)),
                  Tauopathy_PD_signlogp_DEP = sign(Tauopathy_vs_PD_ratio) * (-log10(Tauopathy_vs_PD_p.adj)),
                  PD_Ctr_signlogp_Perseus = sign(`Student's T-test Difference PD_Ctr`) * (-log10(`Student's T-test q-value PD_Ctr`)),
                  Synucleinopathy_Ctr_signlogp_Perseus = sign(`Student's T-test Difference Synucleinopathy_Ctr`) * (-log10(`Student's T-test q-value Synucleinopathy_Ctr`)),
                  Synucleinopathy_PD_signlogp_Perseus = sign(`Student's T-test Difference Synucleinopathy_PD`) * (-log10(`Student's T-test q-value Synucleinopathy_PD`)),
                  Synucleinopathy_Tauopathy_signlogp_Perseus = sign(`Student's T-test Difference Synucleinopathy_Tauopathy`) * (-log10(`Student's T-test q-value Synucleinopathy_Tauopathy`)),
                  Tauopathy_Ctr_signlogp_Perseus = sign(`Student's T-test Difference Tauopathy_Ctr`) * (-log10(`Student's T-test q-value Tauopathy_Ctr`)),
                  Tauopathy_PD_signlogp_Perseus = sign(`Student's T-test Difference Tauopathy_PD`) * (-log10(`Student's T-test q-value Tauopathy_PD`)),
                  PD_Ctr_significant_Perseus = ifelse(`Student's T-test Significant PD_Ctr` == "+",TRUE,FALSE),
                  Synucleinopathy_Ctr_significant_Perseus = ifelse(`Student's T-test Significant Synucleinopathy_Ctr` == "+",TRUE,FALSE),
                  Synucleinopathy_PD_significant_Perseus = ifelse(`Student's T-test Significant Synucleinopathy_PD` == "+",TRUE,FALSE),
                  Synucleinopathy_Tauopathy_significant_Perseus = ifelse(`Student's T-test Significant Synucleinopathy_Tauopathy` == "+",TRUE,FALSE),
                  Tauopathy_Ctr_significant_Perseus = ifelse(`Student's T-test Significant Tauopathy_Ctr` == "+",TRUE,FALSE),
                  Tauopathy_PD_significant_Perseus = ifelse(`Student's T-test Significant Tauopathy_PD` == "+",TRUE,FALSE)) %>%
    dplyr::rename(PD_Ctr_significant_DEP = PD_vs_Ctr_significant,
                  Synucleinopathy_Ctr_significant_DEP = Synucleinopathy_vs_Ctr_significant,
                  Synucleinopathy_PD_significant_DEP = Synucleinopathy_vs_PD_significant,
                  Tauopathy_Ctr_significant_DEP = Tauopathy_vs_Ctr_significant,
                  Tauopathy_PD_significant_DEP = Tauopathy_vs_PD_significant,
                  Synucleinopathy_Tauopathy_significant_DEP = Synucleinopathy_vs_Tauopathy_significant)}

# plot DEP vs Perseus
plot_comparison_methods <- function(data_DEP_Perseus,name_comp){
  data_DEP_Perseus <- as.data.frame(data_DEP_Perseus)
  cols_to_use_DEP <- grep(paste0(name_comp,"_signlogp_DEP"),colnames(data_DEP_Perseus))
  cols_to_use_Perseus <- grep(paste0(name_comp,"_signlogp_Perseus"),colnames(data_DEP_Perseus))
  sign_DEP <- grep(paste0(name_comp,"_significant_DEP"),colnames(data_DEP_Perseus))
  sign_Perseus <- grep(paste0(name_comp,"_significant_Perseus"),colnames(data_DEP_Perseus))
  data_to_plot <- data.frame(x = data_DEP_Perseus[,cols_to_use_DEP],
                             y = data_DEP_Perseus[,cols_to_use_Perseus],
                             sign_DEP = data_DEP_Perseus[,sign_DEP],
                             sign_Perseus = data_DEP_Perseus[,sign_Perseus],
                             name = data_DEP_Perseus$`Protein IDs`) %>%
    mutate(significant = ifelse(sign_DEP == TRUE & is.na(sign_Perseus), "DEP",
                                ifelse((sign_DEP == FALSE | is.na(sign_DEP)) & sign_Perseus == TRUE,"Perseus",
                                       ifelse(sign_DEP == TRUE & sign_Perseus == TRUE,"Both","None"))),
           y = ifelse(!is.na(y),y,0),
           x = ifelse(!is.na(x),x,0))  %>%
    mutate(name = lapply(strsplit(name,"\\|"), function(x){
      use <- x[3]
      use_aux <- unlist(strsplit(use,"_"))
      use_aux[1]
    }))
  cols_significant <- c("DEP" = "seagreen",
                        "Perseus" = "cornflowerblue",
                        "Both" = "tomato2")
  ggplot2::ggplot(data_to_plot,aes(x,y,colour = significant)) +
    geom_point(size = 2.25, color = "grey", alpha = 0.5) +
    geom_point(data = subset(data_to_plot, sign_DEP == TRUE & is.na(sign_Perseus)),size=2.25, alpha=0.6)+
    geom_point(data = subset(data_to_plot,(sign_DEP == FALSE | is.na(sign_DEP)) & sign_Perseus == TRUE), size=2.25,
               alpha=0.6)+
    geom_point(data = subset(data_to_plot, sign_DEP == TRUE & sign_Perseus == TRUE), size=2.25,
               alpha=0.6)+
    theme(panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line()) +
    ggtitle(paste0("Signed p-adjusted value plot comparison between DEP and Perseus for ",name_comp))+
    theme(plot.title = element_text(size = 15, face = "bold"))+
    scale_x_continuous(breaks = seq(-50, 20, by = 5))+
    scale_y_continuous(breaks = seq(-30, 140, by = 10))+
    xlab(paste0("DEP ",name_comp, " (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15, face = "bold"), axis.text=element_text(size=10, face = "bold"))+
    ylab(paste0("Perseus ", name_comp," (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15, face = "bold"), axis.text=element_text(size=10, face = "bold"))+
    geom_vline(xintercept = log10(0.1), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_vline(xintercept = -log10(0.1), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = log10(0.1), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = -log10(0.1), linetype="dashed",color = "grey",alpha = 0.6)+
    #geom_abline(linetype="dashed",color = "grey",alpha = 0.6) +
    annotate(geom = "text", label = "FDR = 10%",x= -log10(0.1)+0.9,y= log10(0.1)-0.15,
             color = "grey") +
    annotate(geom = "text", label = "FDR = 10%",x=log10(0.1)-0.9,y= - log10(0.1)+0.15,
             color = "grey") +
    geom_label_repel(data = subset(data_to_plot, sign_DEP | sign_Perseus), 
                     aes(label=name), size=4,  color="black",
                     nudge_y = 0.3,
                     nudge_x = 0.5, fontface="bold",  
                     min.segment.length = unit(0, "lines"), 
                     max.overlaps = 40) +
    scale_color_manual(values = cols_significant) 
}
