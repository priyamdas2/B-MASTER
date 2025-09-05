setwd("U:/BMASTER/Data raw/Yachida_BMASTER")
# ---- Packages ----
req <- c("ggplot2", "ggrepel")
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
library(ggplot2); library(ggrepel)

# ---- IO helpers ----
read_matrix_no_header <- function(path) {
  df <- read.csv(path, header = FALSE, check.names = FALSE)
  for (j in seq_len(ncol(df))) df[[j]] <- as.numeric(df[[j]])
  for (j in seq_len(ncol(df))) if (anyNA(df[[j]])) df[[j]][is.na(df[[j]])] <- mean(df[[j]], na.rm = TRUE)
  as.data.frame(df)
}
maybe_apply_names <- function(df, names_path) {
  if (!is.null(names_path) && file.exists(names_path)) {
    nm <- read.csv(names_path, header = FALSE, check.names = FALSE)[,1]
    if (length(nm) == ncol(df)) colnames(df) <- nm
  }
  df
}

# ---- Label shorteners ----
shorten_metabolite <- function(v) {
  v2 <- sub("^[^_]*_", "", v)
  v2 <- gsub("_", " ", v2, fixed = TRUE)
  v2 <- gsub("\\.", " ", v2)
  trimws(v2)
}
shorten_microbe <- function(v) {
  parts <- strsplit(v, "__", fixed = TRUE)
  out <- vapply(parts, function(x) x[length(x)], character(1))
  trimws(out)
}

# ---- PCA biplot (with 4-side border + short labels, no title, custom color) ----
pca_biplot <- function(df, k_arrows = 8, file_out = "pca.png",
                       label_shortener = identity,
                       point_color = "#2F76B7",
                       point_size = 2.3) {
  Z <- scale(df, center = TRUE, scale = TRUE)
  pca <- prcomp(Z, center = FALSE, scale. = FALSE)
  scores <- as.data.frame(pca$x[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  
  loads <- as.data.frame(pca$rotation[, 1:2])
  colnames(loads) <- c("L1", "L2")
  mag <- sqrt(loads$L1^2 + loads$L2^2)
  k <- min(k_arrows, nrow(loads))
  idx <- order(mag, decreasing = TRUE)[seq_len(k)]
  loads_top <- loads[idx, , drop = FALSE]
  
  arrow_scale <- 0.9 * mean(apply(scores, 2, sd)) / (mean(apply(loads, 2, sd)) + 1e-8)
  
  arrows_df <- transform(loads_top,
                         x0 = 0, y0 = 0,
                         x1 = L1 * arrow_scale,
                         y1 = L2 * arrow_scale)
  raw_names <- if (!is.null(colnames(df))) rownames(arrows_df) else paste0("V", idx)
  arrows_df$label <- vapply(raw_names, label_shortener, character(1))
  
  pc1_var <- round(100 * summary(pca)$importance[2, 1], 2)
  pc2_var <- round(100 * summary(pca)$importance[2, 2], 2)
  
  p <- ggplot(scores, aes(PC1, PC2)) +
    geom_point(alpha = 0.9, size = point_size, color = point_color) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_segment(data = arrows_df,
                 aes(x = x0, y = y0, xend = x1, yend = y1),
                 arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.4) +
    geom_label_repel(data = arrows_df,
                     aes(x = x1, y = y1, label = label),
                     size = 3, label.size = 0.25,
                     label.padding = unit(0.12, "lines"),
                     fill = "white", min.segment.length = 0) +
    labs(x = paste0("PC1, ", pc1_var, "% variation"),
         y = paste0("PC2, ", pc2_var, "% variation")) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    )
  
  ggsave(file_out, p, width = 11.5, height = 7, dpi = 300)
  invisible(list(file = file_out, labels = arrows_df$label))
}

# ---- Load your data ----
X <- read_matrix_no_header("X.csv")          # microbiome
Y <- read_matrix_no_header("Y.csv")          # metabolites
X <- maybe_apply_names(X, "X_names.csv")
Y <- maybe_apply_names(Y, "Y_names.csv")

if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
if (is.null(colnames(Y))) colnames(Y) <- paste0("M", seq_len(ncol(Y)))

# ---- Make plots ----
res_metab  <- pca_biplot(Y, k_arrows = 8,
                         file_out = "PCA_metabolites.png",
                         label_shortener = shorten_metabolite,
                         point_color = "#D55E00")

res_genera <- pca_biplot(X, k_arrows = 18,
                         file_out = "PCA_genera.png",
                         label_shortener = shorten_microbe,
                         point_color = "#0072B2")

res_metab$labels[1:5]
res_genera$labels[1:5]
