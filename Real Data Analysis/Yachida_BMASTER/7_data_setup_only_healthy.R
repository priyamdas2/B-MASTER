rm(list=ls())

setwd("U:/BMASTER/Data raw/Yachida_BMASTER")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("devtools")

devtools::install_github("microbiome/microbiome")


library(readxl)
library(microbiome)

## ---- Read raw data ----
meta <- read.table("metadata.tsv", sep = "\t", header = TRUE)
genera <- read.table("genera.tsv", sep = "\t", header = TRUE)
metab <- data.frame(read_excel("mtb.xlsx"))

## ---- Read original CRC-analysis feature names ----
X_keep <- as.character(read.table("X_names.csv", header = FALSE, sep = ",")[,1])
Y_keep <- as.character(read.table("Y_names.csv", header = FALSE, sep = ",")[,1])

## ---- Order samples consistently ----
meta <- meta[order(meta$Sample), ]
genera <- genera[order(genera$Sample), ]
metab <- metab[order(metab$Sample), ]

stopifnot(all.equal(meta$Sample, genera$Sample) == TRUE)
stopifnot(all.equal(meta$Sample, metab$Sample) == TRUE)

sample_ids <- meta$Sample

## ---- Drop Sample column ----
genera <- genera[ , names(genera) != "Sample"]
metab <- metab[ , names(metab) != "Sample"]

## ---- Check that original selected features exist ----
missing_X <- setdiff(X_keep, colnames(genera))
missing_Y <- setdiff(Y_keep, colnames(metab))

if (length(missing_X) > 0) {
  stop("These X features are missing from raw genera data: ",
       paste(missing_X, collapse = ", "))
}

if (length(missing_Y) > 0) {
  stop("These Y features are missing from raw metabolite data: ",
       paste(missing_Y, collapse = ", "))
}

## ---- Keep healthy controls only ----
inds_healthy <- which(meta$Study.Group == "Healthy")

meta_healthy <- meta[inds_healthy, ]
genera_healthy <- genera[inds_healthy, X_keep, drop = FALSE]
metab_healthy <- metab[inds_healthy, Y_keep, drop = FALSE]
sample_ids_healthy <- sample_ids[inds_healthy]

## ---- Force numeric matrices ----
genera_healthy <- as.matrix(sapply(genera_healthy, as.numeric))
metab_healthy <- as.matrix(sapply(metab_healthy, as.numeric))

colnames(genera_healthy) <- X_keep
colnames(metab_healthy) <- Y_keep

## ---- CLR transform microbiome ----
genera_healthy <- microbiome::transform(genera_healthy, "clr")

## ---- Log transform metabolites after zero replacement ----
min_value_metab <- min(metab_healthy[metab_healthy > 0], na.rm = TRUE)
metab_healthy[metab_healthy == 0] <- 0.5 * min_value_metab
metab_healthy <- log(metab_healthy)

## ---- Save outputs ----
write.table(genera_healthy, file = "X_healthy.csv",
            row.names = FALSE, col.names = FALSE, sep = ",")

write.table(metab_healthy, file = "Y_healthy.csv",
            row.names = FALSE, col.names = FALSE, sep = ",")

write.csv(meta_healthy, file = "meta_healthy.csv", row.names = FALSE)

write.table(sample_ids_healthy, file = "sample_ids_healthy.csv",
            row.names = FALSE, col.names = FALSE)

cat("Healthy X dimensions:", dim(genera_healthy), "\n")
cat("Healthy Y dimensions:", dim(metab_healthy), "\n")
cat("X columns match original:", identical(colnames(genera_healthy), X_keep), "\n")
cat("Y columns match original:", identical(colnames(metab_healthy), Y_keep), "\n")