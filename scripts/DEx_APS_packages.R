#### Packages

# installation
install.packages("stringr")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DEP",force = TRUE)
BiocManager::install("gmm",force = TRUE)
install.packages("ggExtra")
install.packages("wesanderson")

# run libraries
library("gmm")
library("DEP")
library(stringr)
library(SummarizedExperiment)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(ggplot2)
library(readxl)
library(readr)
library(DEP)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library("ggExtra")
library(wesanderson)