## ============================================================
## Create FIS summary table for HEALTHY subjects
## Bioinformatics revision
## ============================================================

rm(list = ls())

library(dplyr)

## ---- Set working directory ----
setwd("U:/BMASTER/Real Data Analysis")

## ---- Load names ----
X.names <- as.matrix(read.table("Data/X_names.csv",
                                header = FALSE, sep = ","))

## ---- Helper to clean genus names ----
extract_name <- function(input_string) {
  
  gen <- sub(".*g__", "", input_string)
  
  if (grepl("^[0-9]+(\\.[0-9]+)?$", gen)) {
    fam <- sub(".*f__", "", input_string)
    fam <- sub("\\.g__.*", "", fam)
    return(paste0(fam, ".g__", gen))
  } else {
    return(gen)
  }
}

xnames_raw <- as.character(X.names[,1])
xnames <- vapply(xnames_raw, extract_name, character(1))

## ---- Helper function to compute FIS ----
compute_FIS <- function(NonZeroLocation) {
  
  h <- colSums(NonZeroLocation)
  w <- ifelse(h > 0, 1 / h, 0)
  
  FIS <- as.numeric(NonZeroLocation %*% w)
  counts <- rowSums(NonZeroLocation)
  
  data.frame(
    FIS = FIS,
    Count = counts
  )
}

## ---- Read 90%, 95%, 99% selected-edge matrices ----
NZ_90 <- as.matrix(
  read.table(
    "Output_HEALTHY_vs_CRC/Output_conf_90_HEALTHY_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
    header = FALSE, sep = ","
  )
)

NZ_95 <- as.matrix(
  read.table(
    "Output_HEALTHY_vs_CRC/Output_conf_95_HEALTHY_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
    header = FALSE, sep = ","
  )
)

NZ_99 <- as.matrix(
  read.table(
    "Output_HEALTHY_vs_CRC/Output_conf_99_HEALTHY_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
    header = FALSE, sep = ","
  )
)

## ---- Read Bayesian p-values ----
PvalMat <- as.matrix(
  read.table(
    "Output_HEALTHY_vs_CRC/Output_HEALTHY_REAL_pValues_NumIter_1000_burnin_100.csv",
    header = FALSE, sep = ","
  )
)

## ---- Compute FIS for each CI threshold ----
FIS90 <- compute_FIS(NZ_90)
FIS95 <- compute_FIS(NZ_95)
FIS99 <- compute_FIS(NZ_99)

## ---- Compute ranks ----
Rank90 <- rank(-FIS90$FIS, ties.method = "min")
Rank95 <- rank(-FIS95$FIS, ties.method = "min")
Rank99 <- rank(-FIS99$FIS, ties.method = "min")

## ---- Median Bayesian p-value among 90%-selected edges ----
P <- nrow(NZ_90)

MedianPval_90_selected <- sapply(seq_len(P), function(p) {
  selected <- NZ_90[p, ] == 1
  if (sum(selected) == 0) return(NA_real_)
  median(PvalMat[p, selected], na.rm = TRUE)
})

## ---- Assemble one summary table ----
FIS_summary_table <- data.frame(
  Rank_90 = Rank90,
  Genus = xnames,
  FIS_90 = round(FIS90$FIS, 4),
  Selected_Metabolites_90 = FIS90$Count,
  Median_Bayesian_pvalue_90_selected_edges = round(MedianPval_90_selected, 5),
  FIS_95 = round(FIS95$FIS, 4),
  Rank_95 = Rank95,
  FIS_99 = round(FIS99$FIS, 4),
  Rank_99 = Rank99,
  stringsAsFactors = FALSE
)

## ---- Order by primary 90% FIS ranking and add cumulative/percentile ranks ----
FIS_summary_table <- FIS_summary_table %>%
  arrange(Rank_90) %>%
  mutate(
    Cumulative_FIS_percent = round(100 * cumsum(FIS_90) / sum(FIS_90), 2),
    Percentile_Rank_90 = round(100 * (1 - (Rank_90 - 1) / (P - 1)), 2)
  ) %>%
  dplyr::select(
    Rank_90,
    Genus,
    FIS_90,
    Cumulative_FIS_percent,
    Selected_Metabolites_90,
    Median_Bayesian_pvalue_90_selected_edges,
    FIS_95,
    Rank_95,
    FIS_99,
    Rank_99,
    Percentile_Rank_90
  )

## ---- Save healthy FIS summary table ----
write.csv(
  FIS_summary_table,
  "Output_HEALTHY_vs_CRC/Output_FIS_HEALTHY_summary_table.csv",
  row.names = FALSE
)



