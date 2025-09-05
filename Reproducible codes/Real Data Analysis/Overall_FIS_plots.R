## ---- Libraries ----
rm(list = ls())
library(ggplot2)
library(dplyr)
library(forcats)

## ---- Load your data (adjust paths as needed) ----
setwd("U:/BMASTER/Real Data Analysis")

X.names <- as.matrix(read.table("Data/X_names.csv", header = FALSE, sep = ",")) 
Y.names <- as.matrix(read.table("Data/Y_names.csv", header = FALSE, sep = ",")) 
NonZeroLocation <- as.matrix(
  read.table("Output_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
             header = FALSE, sep = ",")
)

## ---- Helper to clean genus names ----
extract_name <- function(input_string) {
  sub(".*__", "", input_string)
}

xnames_raw <- as.character(X.names[,1])
xnames <- vapply(xnames_raw, extract_name, character(1))

## ---- Compute FIS ----
P <- nrow(NonZeroLocation)
Q <- ncol(NonZeroLocation)

h <- colSums(NonZeroLocation)              # number of X’s per Y
w <- ifelse(h > 0, 1 / h, 0)               # weights 1/h_q
FIS <- as.numeric(NonZeroLocation %*% w)   # P × 1 vector
counts <- rowSums(NonZeroLocation)         # number of Y influenced

## ---- Assemble data.frame ----
df <- data.frame(
  genus = xnames,
  FIS   = FIS,
  count = counts,
  stringsAsFactors = FALSE
)

## ---- Top K ----
tophowmanynfis <- 50   # choose 20 or 50
df_top <- df %>%
  arrange(desc(FIS), desc(count)) %>%
  slice_head(n = tophowmanynfis) %>%
  mutate(genus = fct_reorder(genus, FIS, .desc = TRUE))

## ---- Label style ----
label_style <- "number"   # options: "number", "mets", "(mets)"
lab <- switch(label_style,
              "number" = sprintf("%d", df_top$count),
              "mets"   = sprintf("%d mets", df_top$count),
              "(mets)" = sprintf("(%d mets)", df_top$count)
)

## ---- Spacing ----
yr <- range(df_top$FIS)
offset <- 0.5 * diff(yr)    # was 0.08
top_extra <- 0.1            # reduce extra headroom

p_fis_x <- ggplot(df_top, aes(x = genus, y = FIS)) +
  geom_segment(aes(xend = genus, y = 0, yend = FIS),
               linewidth = 0.6, color = "gray75") +
  geom_point(size = 3.2, color = "#0072B2") +
  geom_text(aes(y = FIS + offset, label = lab),
            angle = 90, vjust = 0.5, size = 4.2, fontface = "bold",
            color = "#D55E00") +
  scale_y_continuous(expand = expansion(mult = c(0.02, top_extra))) +
  labs(x = NULL, y = "Fractional Influence Score (FIS)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                size = 9, face = "bold"),
    axis.title.y = element_text(size = 13),
    plot.margin  = margin(6, 14, 6, 8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  coord_cartesian(clip = "off")

print(p_fis_x)

ggsave(paste0("FIS_top", tophowmanynfis, "_wide.png"),
       p_fis_x, width = 9, height = 5, dpi = 300)

################################################################################

## ---- Heatmap of contributions, ranked by FIS (Top K = tophowmanynfis) ----

# Needed only if not already loaded in this session
# (signs of coefficients and Y-name cleaner)
if (!exists("Signs")) {
  Signs <- sign(as.matrix(
    read.table("Output_REAL_B_NonSparse_NumIter_1000_burnin_100.csv",
               header = FALSE, sep = ",")
  ))
}
if (!exists("extract_name_y")) {
  extract_name_y <- function(input_string) sub(".*_", "", input_string)
}

# Row indices of the top-K genera by FIS (map names back to row ids)
idx_top_by_FIS <- match(as.character(df_top$genus), xnames)

# Signed selection matrix (−1, 0, +1)
NonZeroLocationWithSigns <- NonZeroLocation * Signs

# Subset rows (top-K by FIS) and keep all metabolites (columns)
Mat.selected <- as.matrix(NonZeroLocationWithSigns[idx_top_by_FIS, ])

# Set row/column names exactly as before
rownames(Mat.selected) <- extract_name(X.names[idx_top_by_FIS])
colnames(Mat.selected) <- extract_name_y(Y.names)

# Save top-K names (optional, same as before)
write.table(rownames(Mat.selected), "TOP_50_master_predictors_FIS.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)

# Long format for ggplot
library(reshape2)
long_data <- melt(Mat.selected)
colnames(long_data) <- c("Row", "Column", "Value")

# Keep the exact same visual order (top FIS at top of heatmap)
long_data$Row <- factor(long_data$Row, levels = rev(rownames(Mat.selected)))

# Color map: −1 (blue), 0 (white), +1 (red)
colors <- c("-1" = "blue", "0" = "white", "1" = "red")

# Heatmap (styling matches your previous plot)
whole_heat <- ggplot(long_data, aes(x = Column, y = Row, fill = as.factor(Value))) +
  geom_tile(color = "gray") +
  scale_fill_manual(
    values = colors,
    labels = c(expression(beta[pq] < 0),
               expression(beta[pq] == 0),
               expression(beta[pq] > 0))
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 3.4, face = "bold", angle = 90, hjust = 1),
    axis.text.y  = element_text(size = 5.4, face = "bold"),
    axis.title.x = element_text(size = 14, vjust = 0),
    axis.title.y = element_text(size = 14),
    panel.grid   = element_blank(),
    plot.title   = element_blank(),
    legend.text  = element_text(size = 13),
    legend.title = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "right",
    legend.box = "horizontal",
    legend.key.width = unit(1.1, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(2, "cm"),
    plot.margin  = margin(5, 5, 5, 2),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  labs(x = "Metabolites", y = "Genera (top 50 master predictors)") +
  coord_cartesian(expand = FALSE)

whole_heat

# Save with same size as before (adjust if needed)
ggsave("Top50genera_heat_FIS.jpeg", plot = whole_heat, width = 10, height = 6, dpi = 400)
