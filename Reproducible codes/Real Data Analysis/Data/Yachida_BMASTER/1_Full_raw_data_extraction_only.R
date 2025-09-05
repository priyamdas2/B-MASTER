rm(list=ls())
setwd("U:/BMASTER/Data raw/Yachida_BMASTER")
# install.packages("corrplot")
library(corrplot)
library(readxl)
library(BiocManager)
BiocManager::install("microbiome",force = TRUE)
library(devtools)
install_github("Bioconductor/zlibbioc", force = TRUE)
install_github("grimbough/rhdf5filters")

# Read in metadata
meta <- read.table("metadata.tsv", sep = "\t", header = TRUE)

# Read in microbiome data (summarized at genus level)
genera <- read.table("genera.tsv", sep = "\t", header = TRUE)
# 347 samples x 11942 features (first column is "Sample")

# Read in metabolite data
metab <- data.frame(read_excel("mtb.xlsx"))

write.csv(meta,  file = "meta_full.csv",  row.names = FALSE)
write.csv(metab, file = "metab_full.csv", row.names = FALSE)
write.csv(genera, file = "genera_full.csv", row.names = FALSE)