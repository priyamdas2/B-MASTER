# Example dataset: iris (excluding the Species column)
data <- iris[, 1:4]

# Perform PCA
pca_result <- prcomp(data, center = TRUE, scale. = TRUE)

# Print PCA result summary
summary(pca_result)

plot(pca_result, type = "l")

# Scree plot with ggplot2
library(ggplot2)

# Variance explained
var_explained <- (pca_result$sdev)^2 / sum(pca_result$sdev^2)

ggplot(data.frame(PC = 1:length(var_explained), Variance = var_explained), aes(x = PC, y = Variance)) +
  geom_line() +
  geom_point() +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance Explained")



scores <- as.data.frame(pca_result$x)
scores$Species <- iris$Species  # Add species for grouping

ggplot(scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")