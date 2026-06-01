# ====================== Colors ======================
okabe <- list(
  blue = "#0072B2",
  vermillion = "#D55E00"
)

# ====================== NORMALITY (no ridge) ======================
# "Residuals" = centered Y; standardize per metabolite (divide by column SD)
sdY <- apply(Yc, 2, sd, na.rm = TRUE)
sdY[sdY == 0 | !is.finite(sdY)] <- 1
Z   <- sweep(Yc, 2, sdY, "/")   # standardized residuals per metabolite

# (1) Pooled QQ plot of all standardized residuals (color)
z_vec <- as.numeric(Z); z_vec <- z_vec[is.finite(z_vec)]
nz <- length(z_vec)
df_qq <- tibble(
  theo = sort(qnorm(ppoints(nz))),
  samp = sort(z_vec)
)

p_qq <- ggplot(df_qq, aes(x = theo, y = samp)) +
  geom_point(size = 0.8, alpha = 0.55, colour = okabe$blue) +
  geom_abline(intercept = 0, slope = 1, colour = okabe$vermillion, linewidth = 0.9, linetype = "dashed") +
  labs(x = "Theoretical quantiles (Normal)",
       y = "Standardized residual quantiles") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "grey88"),
    panel.grid.minor = element_line(colour = "grey94")
  )

# (2) Excess kurtosis per metabolite (m4 / var^2 - 3) — colored histogram
kurtosis_excess <- apply(Z, 2, function(v){
  v <- v[is.finite(v)]
  m <- mean(v); s2 <- mean((v - m)^2)
  if (!is.finite(s2) || s2 == 0) return(0)
  m4 <- mean((v - m)^4)
  m4/(s2^2) - 3
})
df_kurt <- tibble(kurtosis = kurtosis_excess)

p_kurt <- ggplot(df_kurt, aes(x = kurtosis)) +
  geom_histogram(bins = 30, colour = "black", fill = okabe$blue, alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = okabe$vermillion, linewidth = 0.9) +
  labs(x = "Excess kurtosis across metabolites", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "grey88"),
    panel.grid.minor = element_line(colour = "grey94")
  )

# ===== Save =====
ggsave("Supp_Normality_QQ_color.png",  p_qq,  width = 6.0, height = 5.0, units = "in", dpi = 300)
ggsave("Supp_Kurtosis_Hist_color.png", p_kurt, width = 6.0, height = 5.0, units = "in", dpi = 300)

# ====================== CORRELATION DIAGNOSTICS ======================

# (A) |corr| heatmap for all metabolites Y
corY  <- cor(Yc, use = "pairwise.complete.obs")
corY[!is.finite(corY)] <- 0
corY_abs <- abs(corY)

# Cluster with distance = 1 - |corr|
if (ncol(Yc) >= 2) {
  ordY <- hclust(as.dist(1 - corY_abs))$order
} else {
  ordY <- 1
}
corY_m <- corY_abs[ordY, ordY]

dfY <- as.data.frame(as.table(corY_m))
names(dfY) <- c("i", "j", "abs_corr")

p_heat_Y <- ggplot(dfY, aes(x = i, y = j, fill = abs_corr)) +
  geom_raster() +
  scale_fill_gradient(name = "|corr|", limits = c(0,1), low = "white", high = "red") +
  labs(x = "Metabolites (clustered)", y = "Metabolites (clustered)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "red", fill = NA, linewidth = 0.6),
    legend.position = "right"
  )

# (B) X–Y block |corr| for top-variance features (coarse view)
top_g <- 50
top_m <- 50
top_g_final <- min(top_g, ncol(Xc))
top_m_final <- min(top_m, ncol(Yc))

# guard if tiny
if (top_g_final < 1 || top_m_final < 1) stop("Not enough features for X–Y block heatmap.")

varX <- apply(Xc, 2, var); varY <- apply(Yc, 2, var)
ix <- order(varX, decreasing = TRUE)[seq_len(top_g_final)]
iy <- order(varY, decreasing = TRUE)[seq_len(top_m_final)]

CX <- scale(Xc[, ix, drop = FALSE], center = TRUE, scale = TRUE)
CY <- scale(Yc[, iy, drop = FALSE], center = TRUE, scale = TRUE)

corr_XY <- abs(cor(CX, CY, use = "pairwise.complete.obs"))  # (g × m)
corr_XY[!is.finite(corr_XY)] <- 0

# cluster rows (X features) & cols (Y features) by within-side |corr|
if (ncol(CX) >= 2) {
  row_ord <- hclust(as.dist(1 - abs(cor(CX, use = "pairwise.complete.obs"))))$order
} else {
  row_ord <- 1
}
if (ncol(CY) >= 2) {
  col_ord <- hclust(as.dist(1 - abs(cor(CY, use = "pairwise.complete.obs"))))$order
} else {
  col_ord <- 1
}

corr_XY_m <- corr_XY[row_ord, col_ord, drop = FALSE]

dfXY <- as.data.frame(as.table(corr_XY_m))
names(dfXY) <- c("Genus", "Metabolite", "abs_corr")

p_heat_XY <- ggplot(dfXY, aes(x = Metabolite, y = Genus, fill = abs_corr)) +
  geom_raster() +
  scale_fill_gradient(name = "|corr|", limits = c(0,1), low = "white", high = "red") +
  labs(x = "Metabolites (top-variance, clustered)",
       y = "Genera (top-variance, clustered)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "red", fill = NA, linewidth = 0.6),
    legend.position = "right"
  )

# ===== Print (or save) =====
print(p_qq)
print(p_kurt)
print(p_heat_Y)
print(p_heat_XY)

# Optional saving:
ggsave("Supp_Heatmap_Y_absCorr.png", p_heat_Y,  width = 6, height = 5, units = "in", dpi = 300)
ggsave("Supp_Heatmap_XY_absCorr_topVar.png", p_heat_XY, width = 6.0, height = 5, units = "in", dpi = 300)
