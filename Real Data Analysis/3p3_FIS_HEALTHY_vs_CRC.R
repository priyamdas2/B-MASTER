## ============================================================
## CRC vs HEALTHY FIS comparison: concise table + combined figure
## Bioinformatics revision
## ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

## ---- Input summary tables ----
CRC_FIS <- read.csv(
  "Output/Output_FIS_summary_table.csv",
  stringsAsFactors = FALSE
)

HEALTHY_FIS <- read.csv(
  "Output_HEALTHY_vs_CRC/Output_FIS_HEALTHY_summary_table.csv",
  stringsAsFactors = FALSE
)

## ---- Output folder ----
out_dir <- "Output_HEALTHY_vs_CRC"

## ---- Color palette ----
col_crc     <- "#D55E00"   # red-orange
col_healthy <- "#0072B2"   # blue
col_shared  <- "#009E73"   # green
col_neither <- "grey55"

group_cols <- c(
  "CRC-specific top-50"     = col_crc,
  "Healthy-specific top-50" = col_healthy,
  "Shared top-50"           = col_shared,
  "Neither top-50"          = col_neither
)

## ---- Rename columns for merge ----
CRC_FIS2 <- CRC_FIS %>%
  dplyr::select(
    Genus,
    Rank_CRC = Rank_90,
    FIS_CRC = FIS_90,
    CumFIS_CRC = Cumulative_FIS_percent,
    SelectedMetabolites_CRC = Selected_Metabolites_90,
    MedianPval_CRC = Median_Bayesian_pvalue_90_selected_edges,
    PercentileRank_CRC = Percentile_Rank_90
  )

HEALTHY_FIS2 <- HEALTHY_FIS %>%
  dplyr::select(
    Genus,
    Rank_HEALTHY = Rank_90,
    FIS_HEALTHY = FIS_90,
    CumFIS_HEALTHY = Cumulative_FIS_percent,
    SelectedMetabolites_HEALTHY = Selected_Metabolites_90,
    MedianPval_HEALTHY = Median_Bayesian_pvalue_90_selected_edges,
    PercentileRank_HEALTHY = Percentile_Rank_90
  )

## ---- Full comparison table ----
FIS_compare <- CRC_FIS2 %>%
  full_join(HEALTHY_FIS2, by = "Genus") %>%
  mutate(
    Delta_FIS = FIS_CRC - FIS_HEALTHY,
    Abs_Delta_FIS = abs(Delta_FIS),
    Rank_Shift = Rank_HEALTHY - Rank_CRC,
    Percentile_Shift = PercentileRank_CRC - PercentileRank_HEALTHY,
    Delta_SelectedMetabolites =
      SelectedMetabolites_CRC - SelectedMetabolites_HEALTHY,
    Top50_CRC = Rank_CRC <= 50,
    Top50_HEALTHY = Rank_HEALTHY <= 50,
    Group_Label = case_when(
      Top50_CRC & Top50_HEALTHY ~ "Shared top-50",
      Top50_CRC & !Top50_HEALTHY ~ "CRC-specific top-50",
      !Top50_CRC & Top50_HEALTHY ~ "Healthy-specific top-50",
      TRUE ~ "Neither top-50"
    ),
    Direction = case_when(
      Delta_FIS > 0 ~ "CRC-enriched",
      Delta_FIS < 0 ~ "Healthy-enriched",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(Rank_CRC)

## ============================================================
## One concise manuscript/supplement table
## ============================================================

top_crc_enriched <- FIS_compare %>%
  arrange(desc(Delta_FIS)) %>%
  slice_head(n = 10)

top_healthy_enriched <- FIS_compare %>%
  arrange(Delta_FIS) %>%
  slice_head(n = 10)

top_shared <- FIS_compare %>%
  filter(Group_Label == "Shared top-50") %>%
  arrange(Rank_CRC + Rank_HEALTHY) %>%
  slice_head(n = 10)

Table_diff_master_predictors <- bind_rows(
  top_crc_enriched,
  top_healthy_enriched,
  top_shared
) %>%
  mutate(
    Category = case_when(
      Genus %in% top_crc_enriched$Genus ~ "CRC-enriched",
      Genus %in% top_healthy_enriched$Genus ~ "Healthy-enriched",
      Genus %in% top_shared$Genus ~ "Shared top-50"
    ),
    FIS_CRC = round(FIS_CRC, 2),
    FIS_HEALTHY = round(FIS_HEALTHY, 2),
    Delta_FIS = round(Delta_FIS, 2),
    Percentile_Shift = round(Percentile_Shift, 1)
  ) %>%
  dplyr::select(
    Category,
    Genus,
    Rank_CRC,
    FIS_CRC,
    SelectedMetabolites_CRC,
    Rank_HEALTHY,
    FIS_HEALTHY,
    SelectedMetabolites_HEALTHY,
    Delta_FIS,
    Percentile_Shift
  )

write.csv(
  Table_diff_master_predictors,
  file.path(out_dir, "Table_CRC_HEALTHY_differential_master_predictors.csv"),
  row.names = FALSE
)

## ---- Summary numbers for manuscript text ----
top50_crc <- FIS_compare %>% filter(Top50_CRC) %>% pull(Genus)
top50_healthy <- FIS_compare %>% filter(Top50_HEALTHY) %>% pull(Genus)

overlap_top50 <- intersect(top50_crc, top50_healthy)
crc_only_top50 <- setdiff(top50_crc, top50_healthy)
healthy_only_top50 <- setdiff(top50_healthy, top50_crc)

jaccard_top50 <- length(overlap_top50) / length(union(top50_crc, top50_healthy))

cat("\nTop-50 comparison summary:\n")
cat("CRC top-50:", length(top50_crc), "\n")
cat("Healthy top-50:", length(top50_healthy), "\n")
cat("Shared top-50:", length(overlap_top50), "\n")
cat("CRC-specific top-50:", length(crc_only_top50), "\n")
cat("Healthy-specific top-50:", length(healthy_only_top50), "\n")
cat("Jaccard overlap:", round(jaccard_top50, 4), "\n")
cat("CRC cumulative FIS at top 50:",
    CRC_FIS$Cumulative_FIS_percent[CRC_FIS$Rank_90 == 50], "%\n")
cat("Healthy cumulative FIS at top 50:",
    HEALTHY_FIS$Cumulative_FIS_percent[HEALTHY_FIS$Rank_90 == 50], "%\n")

## ============================================================
## Combined manuscript-ready figure: balanced 3 x 2 panel
## ============================================================

base_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.title.x = element_text(size = 16, margin = margin(t = 4)),
    axis.title.y = element_text(size = 16, margin = margin(r = -6)),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.margin = margin(7, 7, 7, 7)
  )

legend_top_2line <- theme(
  legend.position = "top",
  legend.box = "vertical",
  legend.background = element_rect(color = "grey35", fill = "white", linewidth = 0.35),
  legend.box.background = element_rect(color = "grey35", fill = "white", linewidth = 0.35),
  legend.margin = margin(3, 4, 3, 4),
  legend.key.size = unit(0.42, "cm"),
  legend.text = element_text(size = 12.5)
)

## ---- A. CRC vs Healthy FIS scatterplot ----
p_fis_scatter <- ggplot(
  FIS_compare,
  aes(x = FIS_HEALTHY, y = FIS_CRC)
) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.6, color = "grey35") +
  geom_point(aes(color = Group_Label, shape = Group_Label),
             size = 2.2, alpha = 0.9) +
  scale_color_manual(values = group_cols) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x = "Healthy FIS",
    y = "CRC FIS"
  ) +
  base_theme +
  legend_top_2line

p_fis_scatter <- p_fis_scatter +
  theme(
    axis.title.y = element_text(margin = margin(r = 3)),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---- B. Rank comparison plot ----
p_rank_scatter <- ggplot(
  FIS_compare,
  aes(x = Rank_HEALTHY, y = Rank_CRC)
) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.6, color = "grey35") +
  geom_point(aes(color = Group_Label, shape = Group_Label),
             size = 2.2, alpha = 0.9) +
  scale_color_manual(values = group_cols) +
  scale_x_reverse() +
  scale_y_reverse() +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x = "Healthy FIS rank",
    y = "CRC FIS rank"
  ) +
  base_theme +
  legend_top_2line

p_rank_scatter <- p_rank_scatter +
  theme(
    axis.title.y = element_text(margin = margin(r = 3)),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---- C. Cumulative FIS comparison curve ----
CRC_curve <- CRC_FIS %>%
  arrange(Rank_90) %>%
  mutate(Group = "CRC") %>%
  dplyr::select(Rank_90, Cumulative_FIS_percent, Group)

HEALTHY_curve <- HEALTHY_FIS %>%
  arrange(Rank_90) %>%
  mutate(Group = "Healthy") %>%
  dplyr::select(Rank_90, Cumulative_FIS_percent, Group)

cum_curve <- bind_rows(CRC_curve, HEALTHY_curve)

p_cum_curve <- ggplot(
  cum_curve,
  aes(x = Rank_90, y = Cumulative_FIS_percent, color = Group)
) +
  geom_line(linewidth = 1.0) +
  geom_vline(xintercept = 50,
             linetype = "dashed", linewidth = 0.6, color = "grey35") +
  scale_color_manual(values = c("CRC" = col_crc, "Healthy" = col_healthy)) +
  labs(
    x = "Genus rank based on FIS",
    y = "Cumulative FIS (%)"
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = c(0.20, 0.12),
    legend.background = element_rect(color = "grey40", fill = "white", linewidth = 0.35),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 10)
  )

p_cum_curve <- p_cum_curve +
  theme(
    axis.title.y = element_text(margin = margin(r = 3)),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---- D. Top-50 classification counts ----
class_counts <- FIS_compare %>%
  filter(Top50_CRC | Top50_HEALTHY) %>%
  count(Group_Label) %>%
  mutate(
    Group_Label = factor(
      Group_Label,
      levels = c(
        "Shared top-50",
        "CRC-specific top-50",
        "Healthy-specific top-50"
      )
    ),
    Group_Label_plot = case_when(
      Group_Label == "Shared top-50" ~ "Shared\ntop-50",
      Group_Label == "CRC-specific top-50" ~ "CRC-specific\ntop-50",
      Group_Label == "Healthy-specific top-50" ~ "Healthy-specific\ntop-50"
    )
  )

p_class_counts <- ggplot(
  class_counts,
  aes(x = Group_Label_plot, y = n, fill = Group_Label)
) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = group_cols) +
  labs(
    x = NULL,
    y = "Number of genera"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, lineheight = 0.9)
  )

p_class_counts <- p_class_counts +
  theme(
    axis.title.y = element_text(margin = margin(r = 3)),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---- E. CRC-enriched barplot ----
plot_crc_enriched <- FIS_compare %>%
  arrange(desc(Delta_FIS)) %>%
  slice_head(n = 15) %>%
  mutate(Genus = factor(Genus, levels = rev(Genus)))

p_crc_enriched <- ggplot(
  plot_crc_enriched,
  aes(x = Genus, y = Delta_FIS)
) +
  geom_col(width = 0.75, fill = col_crc) +
  coord_flip() +
  labs(
    x = NULL,
    y = "CRC minus healthy FIS"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10.5)
  )

## ---- F. Healthy-enriched barplot ----
plot_healthy_enriched <- FIS_compare %>%
  arrange(Delta_FIS) %>%
  slice_head(n = 15) %>%
  mutate(
    Delta_FIS_abs = abs(Delta_FIS),
    Genus = factor(Genus, levels = rev(Genus))
  )

p_healthy_enriched <- ggplot(
  plot_healthy_enriched,
  aes(x = Genus, y = Delta_FIS_abs)
) +
  geom_col(width = 0.75, fill = col_healthy) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Healthy minus CRC FIS"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10.5)
  )


## ---- Force equal panel boxes without axis-column alignment ----
combined_fig <-
  (
    (wrap_elements(full = p_fis_scatter) +
       wrap_elements(full = p_rank_scatter)) /
      
      plot_spacer() /
      
      (wrap_elements(full = p_cum_curve) +
         wrap_elements(full = p_class_counts)) /
      
      plot_spacer() /
      
      (wrap_elements(full = p_crc_enriched) +
         wrap_elements(full = p_healthy_enriched)) /
      
      plot_spacer()
  ) +
  plot_layout(
    ncol = 1,
    heights = c(
      1.0,   # row 1
      0.08,  # whitespace
      1.0,   # row 2
      0.08,  # whitespace
      1.0,   # row 3
      0.12   # bottom whitespace
    )
  )

ggsave(
  file.path(out_dir, "Figure_CRC_HEALTHY_master_predictor_comparison.png"),
  combined_fig,
  width = 15.8,
  height = 16.8,
  dpi = 400
)



cat("\nSaved files:\n")
cat("1. Table_CRC_HEALTHY_differential_master_predictors.csv\n")
cat("2. Figure_CRC_HEALTHY_master_predictor_comparison.png\n")