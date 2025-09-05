setwd("U:/BMASTER/Data raw/Yachida_BMASTER")

# ===== Packages =====
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(scales)     # for alpha()
library(ggpattern)

# ===== Load & normalize =====
meta <- readr::read_csv("meta_full.csv", show_col_types = FALSE)

meta <- meta %>%
  mutate(
    .sg = Study.Group %>% str_trim() %>% str_replace_all("\\s+|-", "_") %>% toupper(),
    Study.Group = case_when(
      .sg == "HEALTHY" ~ "Healthy",
      .sg == "MP" ~ "MP",
      .sg %in% c("STAGE_0","STAGE0","S0") ~ "Stage_0",
      .sg %in% c("STAGE_I_II","SI_II","STAGE_I__II") ~ "Stage_I_II",
      .sg %in% c("STAGE_III_IV","SIII_IV","STAGE_III__IV") ~ "Stage_III_IV",
      .sg %in% c("HS","HISTORY_OF_SURGERY","HISTORYOFSURGERY") ~ "HS",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-.sg)

# ===== Count & prepare table (exclude Healthy) =====
groups_keep <- c("MP","Stage_0","Stage_I_II","Stage_III_IV","HS")

tab2 <- meta %>%
  filter(Study.Group %in% groups_keep) %>%
  mutate(Study.Group = factor(Study.Group, levels = groups_keep)) %>%
  count(Study.Group, name = "n") %>%
  right_join(tibble::tibble(Study.Group = factor(groups_keep, levels = groups_keep)),
             by = "Study.Group") %>%
  mutate(n = tidyr::replace_na(n, 0L)) %>%
  mutate(
    group_label = recode(Study.Group,
                         "MP"="MP","Stage_0"="Stage 0","Stage_I_II"="Stage I–II",
                         "Stage_III_IV"="Stage III–IV","HS"="HS"),
    group_label = factor(group_label, levels = c("MP","Stage 0","Stage I–II","Stage III–IV","HS")),
    pct = 100 * n / sum(n),
    lbl = paste0(n, "(", sprintf("%.1f", pct), "%)"),
    is_hs = (Study.Group == "HS")
  )

# ===== Okabe–Ito colors (different pair from genera/metabolites) =====
okabe <- list(green = "#009E73", purple = "#CC79A7")

# ===== Plot =====
p <- ggplot(tab2, aes(y = group_label, x = n)) +
  geom_col(height = 0.70, fill = scales::alpha(okabe$green, 0.35),
           colour = okabe$green, linewidth = 0.6) +
  ggpattern::geom_col_pattern(
    data = dplyr::filter(tab2, is_hs),
    aes(y = group_label, x = n),
    height = 0.70,
    fill = "white",
    colour = okabe$purple,
    pattern = "stripe",
    pattern_fill   = okabe$purple,
    pattern_colour = okabe$purple,
    pattern_angle  = 45,
    pattern_density = 0.35,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.6
  ) +
  geom_text(aes(label = lbl), hjust = -0.1, size = 5, colour = "black") +
  expand_limits(x = max(tab2$n, na.rm = TRUE) * 1.20) +
  labs(x = "Number of subjects", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85"),
    axis.line.x        = element_line(colour = "black"),
    axis.ticks.x       = element_line(colour = "black"),
    axis.title.x       = element_text(size = 15),
    axis.title.y       = element_text(size = 18),
    axis.text          = element_text(size = 15)
  )

ggsave("Fig_StudyGroup_horizontal_color_HS_hatched.png", p,
       width = 6, height = 4.5, units = "in", dpi = 300)


################################################################################

################################################################################
# ---- Packages ----
library(ggplot2)
library(cowplot)   # for side-by-side panel
library(readr)

okabe <- list(blue = "#0072B2")  # single color for both panels

# ---- Read in data ----
X_before <- read_csv("X_before_CLR.csv", col_names = FALSE)
X_after  <- read_csv("X.csv", col_names = FALSE)   # CLR-transformed matrix

# ---- PCA ----
pca_before <- prcomp(as.matrix(X_before), center = TRUE, scale. = FALSE)
pca_after  <- prcomp(as.matrix(X_after),  center = TRUE, scale. = FALSE)

# Variance explained (for axis labels)
ve_before <- round(100 * pca_before$sdev^2 / sum(pca_before$sdev^2), 1)
ve_after  <- round(100 * pca_after$sdev^2  / sum(pca_after$sdev^2), 1)

# Data frames for plotting
df_before <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2])
df_after  <- data.frame(PC1 = pca_after$x[,1],  PC2 = pca_after$x[,2])

# ---- Base theme with panel borders ----
base_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "grey88"),
    panel.grid.minor = element_line(colour = "grey94")
  )

# ---- Scatterplots (same color before & after) ----
g_before <- ggplot(df_before, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.8, alpha = 0.9, colour = okabe$blue) +
  labs(x = paste0("PC1 (", ve_before[1], "%)"),
       y = paste0("PC2 (", ve_before[2], "%)")) +
  base_theme +
  ggtitle("Before CLR")

g_after <- ggplot(df_after, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.8, alpha = 0.9, colour = okabe$blue) +
  labs(x = paste0("PC1 (", ve_after[1], "%)"),
       y = paste0("PC2 (", ve_after[2], "%)")) +
  base_theme +
  ggtitle("After CLR")

# ---- Combine side-by-side with spacer ----
spacer <- ggplot() + theme_void()

p_combined <- cowplot::plot_grid(
  g_before, spacer, g_after,
  nrow = 1,
  rel_widths = c(1, 0.15, 1)   # increase middle value for more space
)

print(p_combined)
ggsave("PCA_before_after_CLR_samecolor.png", p_combined,
       width = 11, height = 4, units = "in", dpi = 300)
################################################################################


################################################################################
# ---- Singular value spectrum ----
library(grid)  # for unit()

# ---- Read data ----
X <- as.matrix(read_csv("X.csv", col_names = FALSE))
Y <- as.matrix(read_csv("Y.csv", col_names = FALSE))
stopifnot(nrow(X) == 220, nrow(Y) == 220)

# ---- Singular values ----
svals_X <- svd(X)$d
svals_Y <- svd(Y)$d

# ---- Data for plotting ----
df_all <- rbind(
  data.frame(Index = seq_along(svals_X), Value = svals_X, Matrix = "Genera (X)"),
  data.frame(Index = seq_along(svals_Y), Value = svals_Y, Matrix = "Metabolites (Y)")
)
df_all$Matrix <- factor(df_all$Matrix, levels = c("Genera (X)", "Metabolites (Y)"))

# ---- Parameters ----
rank_cap <- 220L
tail_start <- 206L
ymin <- 1e-6

# ---- Plot with color lines ----
p_spec <- ggplot(df_all, aes(Index, Value, colour = Matrix)) +
  annotate("rect", xmin = tail_start, xmax = rank_cap,
           ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "grey70") +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("Genera (X)" = okabe$blue,
                                 "Metabolites (Y)" = okabe$vermillion)) +
  scale_y_log10(limits = c(ymin, NA)) +
  labs(x = "Singular value index", y = "Singular value (log scale)", colour = NULL) +
  geom_vline(xintercept = rank_cap, linetype = "dashed", linewidth = 0.6) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position   = "top",
    legend.text       = element_text(size = 14),
    axis.title        = element_text(size = 16),
    axis.text         = element_text(size = 13),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    legend.key.width  = unit(2, "cm"),
    legend.key.height = unit(0.35, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.4)))

ggsave("Design_matrix_spectrum_color.png", p_spec,
       width = 6.0, height = 5.4, units = "in", dpi = 300)
################################################################################
