rm(list=ls())
library(ggplot2)
library(GGally)
library(CCA)
library(CCP)

setwd("U:/BMASTER/Real Data Analysis")
extract_name <- function(input_string) {
  # Use regular expression to extract everything after the last underscore
  extracted_name <- sub(".*__", "", input_string)
  return(extracted_name)
}


SubsetNum <- 1   # 1 = Most abundant, 2 = Differential in cancer

X <- as.matrix(read.table("Data/X.csv", header = FALSE, sep = ","))    # p = 287 
Y.set <- as.matrix(read.table(paste0('Data/Y.subset', SubsetNum,'.csv'), header = FALSE, sep = ","))  
MasterRanks.set <- as.matrix(read.table(paste0('MasterRanks_subset', SubsetNum,'.csv'), header = FALSE, sep = ",")) 

P <- dim(X)[2]
N <- dim(X)[1]

Cum.ranked.cca <- rep(NA, P)
Cum.serial.cca <- rep(NA, P)
set.seed(1)
X.shuffled <- X[, sample(1:P)]

for (p in 2:219) {
  print(p)
  x.coors <- MasterRanks.set[1:p]
  cc1 <- cc(Y.set, X[ , x.coors])
  Cum.ranked.cca[p] <- cc1$cor[1]
  
  cc2 <- cc(Y.set, X.shuffled[ , 1:p])
  Cum.serial.cca[p] <- cc2$cor[1]
}
  


data <- data.frame(
  cumulative_top = 2:219,  # X-axis values
  canonical_correlation = Cum.ranked.cca[2:219]  # Y-axis values
)


if(SubsetNum == 1) {
  tt <- "Most abundant metabolites"
}
if(SubsetNum == 2) {
  tt <- "Differential metabolites in cancer"
}



# Create the plot
plot <- ggplot(data, aes(x = cumulative_top, y = canonical_correlation)) +
  geom_line(color = "blue", size = 1) +  # Line plot
  geom_point(color = "red", size = 0.5, alpha = 0.7) +  # Add points for emphasis
  labs(
    title = tt,
    x = "Genera (cumulative ordered master predictors)",
    y = "Canonical correlation"
  ) +
  theme_minimal(base_size = 14) +  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),  # Center title
    axis.title = element_text(size = 16),  # Bold axis titles
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a border
  )

plot
ggsave(paste0('CumCCA_subset',SubsetNum,'.jpeg'), plot = plot, width = 6, height = 4.5, dpi = 400)


