## ============================================================
## Create one FIS summary table and ranked FIS plot
## for Bioinformatics revision
## ============================================================

rm(list = ls())

library(dplyr)
library(ggplot2)

## ---- Set working directory ----
setwd("U:/BMASTER/Real Data Analysis")

## ---- Load names ----
X.names <- as.matrix(read.table("Data/X_names.csv",
                                header = FALSE, sep = ","))

## ---- Helper to clean genus names ----
## Keeps ordinary genus names unchanged.
## If the genus label is numeric/decimal only, keeps family/genus context.
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
    "Output/Output_conf_90_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
    header = FALSE, sep = ","
  )
)

NZ_95 <- as.matrix(
  read.table(
    "Output/Output_conf_95_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
    header = FALSE, sep = ","
  )
)

NZ_99 <- as.matrix(
  read.table(
    "Output/Output_conf_99_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
    header = FALSE, sep = ","
  )
)

## ---- Read Bayesian p-values ----
PvalMat <- as.matrix(
  read.table(
    "Output/Output_REAL_pValues_NumIter_1000_burnin_100.csv",
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

## ---- Order by primary 90% FIS ranking and add cumulative percentage ----
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

## ---- Save table ----
write.csv(
  FIS_summary_table,
  "Output/Output_FIS_summary_table.csv",
  row.names = FALSE
)

## ============================================================
## Ranked FIS plot across all genera
## ============================================================

p_ranked_fis <- ggplot(FIS_summary_table, aes(x = Rank_90, y = FIS_90)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  geom_vline(xintercept = 50, linetype = "dashed", linewidth = 0.7) +
  annotate(
    "text",
    x = 55,
    y = max(FIS_summary_table$FIS_90, na.rm = TRUE) * 0.95,
    label = "Rank 50",
    hjust = 0,
    size = 4
  ) +
  annotate(
    "text",
    x = 55,
    y = max(FIS_summary_table$FIS_90, na.rm = TRUE) * 0.85,
    label = paste0(
      "Cumulative FIS = ",
      FIS_summary_table$Cumulative_FIS_percent[FIS_summary_table$Rank_90 == 50],
      "%"
    ),
    hjust = 0,
    size = 4
  ) +
  labs(
    x = "Genus rank based on FIS",
    y = "Fractional Influence Score (FIS)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(6, 10, 6, 8)
  )

ggsave(
  "Output/Output_ranked_FIS_all_genera.png",
  plot = p_ranked_fis,
  width = 7.2,
  height = 4.8,
  dpi = 400
)

## ============================================================
## Save polished FIS summary table for Supplement
## ============================================================

FIS_summary_table_polished <- FIS_summary_table %>%
  transmute(
    Genus = Genus,
    `FIS rank (default; 90% CI-based)` = Rank_90,
    `FIS (default; 90% CI-based)` = round(FIS_90, 2),
    `Cumulative FIS (%)` = round(Cumulative_FIS_percent, 2),
    `Num. Selected metabolites` = Selected_Metabolites_90,
    `Median Bayesian p-value (default; 90% CI-based)` =
      round(Median_Bayesian_pvalue_90_selected_edges, 3),
    `FIS (95% CI-based)` = round(FIS_95, 2),
    `Rank (95% CI-based)` = Rank_95,
    `FIS (99% CI-based)` = round(FIS_99, 2),
    `Rank (99% CI-based)` = Rank_99
  )

write.csv(
  FIS_summary_table_polished,
  "Output/Output_FIS_summary_table_polished.csv",
  row.names = FALSE
)