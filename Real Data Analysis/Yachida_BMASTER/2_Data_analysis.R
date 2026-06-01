rm(list=ls())
library(corrplot)
library(readxl)

library(BiocManager)
BiocManager::install("microbiome",force = TRUE)

setwd("U:/BMASTER/Data raw/Yachida_BMASTER")

meta   <- read.csv("meta_full.csv",  header = TRUE, stringsAsFactors = FALSE)
metab  <- read.csv("metab_full.csv", header = TRUE, stringsAsFactors = FALSE)
genera <- read.csv("genera_full.csv", header = TRUE, stringsAsFactors = FALSE)

# Inspect dimensions
dim(meta)      # 347 x  15
dim(metab)     # 347 x 451
dim(genera)    # 347 x 11943

# Put all three data frames in same order, then drop "Sample" column
meta <- meta[order(meta$Sample), ]



genera <- genera[order(genera$Sample), ]
metab <- metab[order(metab$Sample), ]
all.equal(meta$Sample, genera$Sample)
all.equal(meta$Sample, metab$Sample)
genera <- genera[ , -which(names(genera) == "Sample")]
metab <- metab[ , -which(names(metab) == "Sample")]

# Drop healthy controls (so that data set is just subjects with cancer)
# Another option would be do the opposite (i.e., only analyze healthy controls)
inds_to_drop <- which(meta$Study.Group == "Healthy")
meta <- meta[-inds_to_drop, ]
genera <- genera[-inds_to_drop, ]
metab <- metab[-inds_to_drop, ]

write.csv(meta,  file = "meta_used.csv", row.names = FALSE)

# Dimensions after filtering subjects
dim(genera)  # 220 x 11942
dim(metab)   # 220 x 450

# Need to filter both X and Y to avoid zero-inflated features

# Summarize proportion of non-zero values
genera_prop_nonzero <- colMeans(genera > 0)
metab_prop_nonzero <- colMeans(metab > 0)

# Across features, what is the proportion of values that are non-zero?
# You can see both the microbiome and metabolome data contain some 
# features that are mostly 0s
hist(genera_prop_nonzero)
hist(metab_prop_nonzero)


################################################################################
## === Prevalence cutoff plot (main paper) ===
################################################################################
# Packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
library(grid)  # for unit()

# genera_prop_nonzero, metab_prop_nonzero ∈ [0,1]

# long form
prev_df <- rbind(
  data.frame(prevalence = genera_prop_nonzero, feature = "Genus prevalence"),
  data.frame(prevalence = metab_prop_nonzero,  feature = "Metabolite prevalence")
)

# single cutoff at 0.20
cutoff <- 0.20

# Colorful, publication-friendly palette (Okabe–Ito)
cols <- c(
  "Genus prevalence"      = "#0072B2", # blue
  "Metabolite prevalence" = "#D55E00"  # vermillion
)

p_prev <- ggplot(prev_df, aes(x = prevalence, colour = feature)) +
  geom_density(linewidth = 1.4, adjust = 1.1) +
  geom_vline(xintercept = cutoff, linetype = "dotted", linewidth = 0.9, colour = "grey25") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_colour_manual(values = cols, breaks = names(cols), guide = guide_legend(order = 1)) +
  labs(
    x = "Prevalence (fraction of subjects with nonzero value)",
    y = "Density",
    colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title        = element_blank(),
    legend.position   = "top",
    legend.key.width  = unit(1.6, "cm"),
    legend.key.height = unit(0.35, "cm"),
    panel.grid.major  = element_line(linewidth = 0.4, colour = "grey85"),
    panel.grid.minor  = element_line(linewidth = 0.2, colour = "grey92"),
    axis.line         = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.8, colour = "black")
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.8)))

print(p_prev)

# export
ggsave("prevalence_color_style.png", p_prev, width = 4.5, height = 4.2, units = "in", dpi = 300)

################################################################################

# Microbiome processing: measures are already on relative abundance scale,
# so each row sums to 1
rowSums(genera)

# Filter to features present in at least 20% of samples
genera <- genera[ , -which(genera_prop_nonzero < .2)]
metab <- metab[ , -which(metab_prop_nonzero < .2)]

dim(genera)  # 220 x 3896
dim(metab)   # 220 x 249


# Filter microbiome further to featurs with average abundance at least 0.01%
low_abund_features <- which(colMeans(genera) < 0.0001)
genera <- genera[ , -low_abund_features]


###########
write.table(genera, file= "X_before_CLR.csv", row.names=FALSE, col.names=FALSE, 
            sep = ",")
write.table(metab, file= "Y_before_trans.csv", row.names=FALSE, col.names=FALSE, 
            sep = ",")

###########

# Apply a clr transform to the microbiome data
genera <- microbiome::transform(genera, 'clr')

# For metabolite data, add 1/2 min to zeros and then log transform
min_value_metab <- min(metab[which(metab > 0, arr.ind = TRUE)])
metab[which(metab == 0, arr.ind = TRUE)] <- 0.5 * min_value_metab
metab <- log(metab)

# FINAL DATA
dim(genera) # 220 samples x 287 features as X matrix
dim(metab)  # 220 samples x 249 features as Y matrix

mydata.cor = cor(metab[, 1:20])
corrplot(mydata.cor)

# write.csv(genera, file= "X.csv")
# write.csv(metab, file= "Y.csv")

write.table(genera, file= "X.csv", row.names=FALSE, col.names=FALSE, 
            sep = ",")
write.table(metab, file= "Y.csv", row.names=FALSE, col.names=FALSE, 
            sep = ",")

write.table(colnames(genera), file= "X_names.csv", row.names=FALSE, 
            col.names=FALSE)
write.table(colnames(metab), file= "Y_names.csv", row.names=FALSE,
            col.names=FALSE)

