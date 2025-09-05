rm(list=ls())
setwd("U:/BMASTER/Real Data Analysis")
library(ggplot2)
extract_name <- function(input_string) {
  # Use regular expression to extract everything after the last underscore
  extracted_name <- sub(".*__", "", input_string)
  return(extracted_name)
}

### Number of Master predictors to be plotted ##################################

max.num.X <- 15 

### Names of the variables, X = genera, Y = metabolite #########################

# subset2 = metabolites were identified as differential in cancer vs. control
#           and are hypothesized to be related to intestinal microbiota

X <- read.table("Data/X.csv", header = FALSE, sep = ",")    # p = 287 
Y.set <- read.table("Data/Y.subset2.csv", header = FALSE, sep = ",")    # q = 11
X.names <- read.table("Data/X_names.csv", header = FALSE)    # p = 287 
Y.names.set <- read.table("Data/Y.names.subset2.csv", header = FALSE)   # q = 11 
Y.idx.set <- read.table("Data/Y.indexes.subset2.csv", header = FALSE)

p <- dim(X.names)[1]
q <- dim(Y.names.set)[1]



### Reading analysis outputs ###################################################
posterior.mean.NonSparse.All <- read.table("Output_REAL_B_NonSparse_NumIter_1000_burnin_100.csv",
                                           header = FALSE, sep = ",") 
posterior.mean.All <- read.table("Output_REAL_B_NumIter_1000_burnin_100.csv",
                                 header = FALSE, sep = ",")  # dim = 287 X 10 (p x q) 
non.zero.locations.All <- read.table("Output_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv",
                                     header = FALSE, sep = ",")
p.values.All <- read.table("Output_REAL_pValues_NumIter_1000_burnin_100.csv",
                           header = FALSE, sep = ",")
idx.set <- as.numeric(unlist(Y.idx.set))
posterior.mean.set <- posterior.mean.All[, idx.set]
non.zero.locations.set <- non.zero.locations.All[, idx.set]
p.values.set <- p.values.All[, idx.set]
non.zero.proportion.set <- sum(non.zero.locations.set) / (p * q)  


# subset2 = metabolites were identified as differential in cancer vs. control
#           and are hypothesized to be related to intestinal microbiota

subset2 <- c("X_DCA", "C01921_Glycocholate", "C05122_Taurocholate",
             "C08262_Isovalerate", "C00407_Ile", "C00123_Leu", 
             "C00183_Val", "C00079_Phe", "C00082_Tyr", "C00065_Ser",
             "C00037_Gly")


### Extracting the signs of coefficients of X for important Y subsets ##########

X.coeffs.signs.set <- matrix(0, p, q)

for (i in 1:p) {
  for (j in 1:q) {
    if(posterior.mean.set[i, j] > 0.0001) {
      X.coeffs.signs.set[i, j] <- 1
    }
    if(posterior.mean.set[i, j] < - 0.0001) {
      X.coeffs.signs.set[i, j] <- - 1
    }
  }
}


### Relevance gives number of Y's for which X_i non-zero coeffs ################

X.relevance.subset <- rowSums(abs(X.coeffs.signs.set))

### Calculating mean p-values of the relevant position #########################
num.relevant.pvalues <- array(NA, p)
mean.pvalue.of.nonzero.pos <- array(NA, p)

for (i in 1:p) {
  relevant.pvalues <- p.values.set[i, which(p.values.set[i, ] <= 0.05)]
  num.relevant.pvalues[i] <- length(relevant.pvalues)
  if (length(relevant.pvalues) > 0) {
    mean.pvalue.of.nonzero.pos[i] <- mean(as.numeric(relevant.pvalues))
  }
}

# Selecting top S many X variables based on relevance. If relevances are same, 
# then we select the one with smaller mean of relevance

relevance.minus.pvalue <- num.relevant.pvalues - mean.pvalue.of.nonzero.pos
sorted.relevance.indices <- order(relevance.minus.pvalue, decreasing = TRUE)
selected.X.positions.set <- sorted.relevance.indices[1:max.num.X]

X.selected.names.set <- unlist(X.names)[as.numeric(selected.X.positions.set)]
X.selected.coeffs.set <- posterior.mean.set[selected.X.positions.set, ]
X.selected.coeffs.signs.set <- X.coeffs.signs.set[selected.X.positions.set, ]
X.selected.pvalues.set.temp <- p.values.set[selected.X.positions.set, ]
X.selected.pvalues.set <- X.selected.pvalues.set.temp

for (i in 1:max.num.X) {
  for (j in 1:q) {
    if(as.numeric(X.selected.pvalues.set[i,j]) > 0.05) {
      X.selected.pvalues.set[i,j] <- NA
    }
  }
}


X.selected.names.set
X.selected.names.set.modified <- extract_name(X.selected.names.set[1:max.num.X])

write.table(X.selected.names.set.modified,"TOP_15_Subset_2.csv", sep = ',', 
            row.names = FALSE, col.names = FALSE)

################################################################################
################ Plots #########################################################
################################################################################

subset.names <- c("X_DCA", "Glycocholate", "Taurocholate", "Isovalerate",
                   "L-Isoleucine", "L-Leucine", "L-Valine", "L-Phenylalanine",
                   "L-Tyrosine", "L-Serine", "Glycine")

# To avoid plotting zero coeffs. Zero positions are deleted controlling p-value > 0.05

posterior.mean.NonSparse.set <- posterior.mean.NonSparse.All[, idx.set]

X.selected.coeffs.set 
X.selected.coeffs.NonSparse.set <- posterior.mean.NonSparse.set[selected.X.positions.set, ]



d <- data.frame(
  x = c(rep(subset.names[1], max.num.X), rep(subset.names[2], max.num.X),
        rep(subset.names[3], max.num.X), rep(subset.names[4], max.num.X), 
        rep(subset.names[5], max.num.X), rep(subset.names[6], max.num.X),
        rep(subset.names[7], max.num.X), rep(subset.names[8], max.num.X),
        rep(subset.names[9], max.num.X), rep(subset.names[10], max.num.X),
        rep(subset.names[11], max.num.X)),
  y = rep(X.selected.names.set.modified, q),
  value = as.numeric(unlist(X.selected.pvalues.set)),   # by col
  value2 = factor(unlist(sign(X.selected.coeffs.NonSparse.set)))    # by col
)


library(forcats)
ggp <- ggplot(d, aes(fct_inorder(x), fct_rev(fct_inorder(y)), fill = value2, size = value)) +
  geom_point(shape = 21, stroke = 0) +
  geom_hline(yintercept = seq(.5, max.num.X + 0.5, 1), size = .2) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(10, 2)) +
  scale_fill_manual(breaks = c(1,  -1), values = c("red", "blue"),
                    labels = c("Positive", "Negative")) +
  #  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold",size = 11),
        axis.text.x = element_text(face = "bold", size = 9),  # Make x-axis label bold
        axis.text.y = element_text(face = "bold", size = 10),   # Make y-axis label bold
        axis.title.x = element_text(size = 18), # Make x-axis title bold and larger
        axis.title.y = element_text(size = 16) 
  ) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "top", 
                             order = 1),
         fill = guide_legend(ticks.colour = NA, title.position = "top", 
                             order = 3, override.aes = list(size = 8))) +
  labs(size = "Area = Bayesian p-value", fill = "Sign of coefficient:", 
       x = "Differentially abundant metabolites in cancer vs control",  # Add your X-axis label here
       y = "Genera (Top 15 master predictors)")

ggp

ggsave("Subset_2.png", plot = ggp, width = 12, height = 7, dpi = 400)

write.table(sorted.relevance.indices, "MasterRanks_subset2.csv", sep = ',', 
            row.names = FALSE, col.names = FALSE)


################################################################################
### Saving p-values
################################################################################
## --- Save p-values (rounded to 3 decimals) ---
pval_mat <- as.data.frame(round(X.selected.pvalues.set, 3))
rownames(pval_mat) <- X.selected.names.set.modified
colnames(pval_mat) <- subset.names
write.csv(pval_mat, "Selected_pvalues_subset2.csv", row.names = TRUE)

## --- Save coefficient directions from NON-SPARSE means ---
##     (so colors match the plot). Threshold tiny magnitudes and mask non-significant cells.
thr <- 1e-4
sign_raw <- sign(X.selected.coeffs.NonSparse.set)
sign_raw[abs(X.selected.coeffs.NonSparse.set) <= thr] <- 0

## Mask signs where p-value is NA (i.e., p > 0.05 in your earlier filtering)
sign_raw[is.na(X.selected.pvalues.set)] <- NA

sign_mat <- as.data.frame(sign_raw)
rownames(sign_mat) <- X.selected.names.set.modified
colnames(sign_mat) <- subset.names
write.csv(sign_mat, "Selected_signs_subset2.csv", row.names = TRUE)
