rm(list=ls())
library(DescTools)
library(ggplot2)
library(reshape2)
setwd("U:/BMASTER/Simulation study 1/Summary Table and Plots")

set.seed(1)
## Reading data ################################################################
NonZeroLocation <- as.matrix(read.csv('NonZeroLocation_BMASTER_estimated_NumIter_2000_burnin_100.csv',header=FALSE))
NonZeroLocation_est <- as.matrix(read.csv('Output_simReal_NonZeroLocation_BMASTER_rep_1.csv',header=FALSE))


p <- dim(NonZeroLocation)[1]
q <- dim(NonZeroLocation)[2]
# prob <- 0.25
# NonZeroLocation_est <- matrix(rbinom(p * q, size = 1, prob = prob), nrow = p, ncol = q)
A <- t(NonZeroLocation)
B <- t(NonZeroLocation_est)

combined <- matrix(NA, nrow = q, ncol = p)
combined[A == 1 & B == 1] <- "green"   # Both A and B have 1 -> True Positive
combined[A == 0 & B == 0] <- "black"   # Both A and B have 0 -> True Negative
combined[A == 0 & B == 1] <- "red"    # A has 1, B has 0 -> False Negative
combined[A == 1 & B == 0] <- "violet"     # A has 0, B has 1 -> False Positive

TP <- sum(A == 1 & B == 1) 
TN <- sum(A == 0 & B == 0) 
FP <- sum(A == 0 & B == 1) 
FN <- sum(A == 1 & B == 0) 
TPR <- TP/(TP+FN)
FPR <- FP/(FP+TN)
c(TP, TN, FP, FN) / (p * q)

# Convert to a data frame for plotting
plot_data <- melt(combined)
colnames(plot_data) <- c("Row", "Column", "Color")

# Add a factor to map colors to categories
plot_data$Color <- factor(plot_data$Color, levels = c("green", "black","red", "violet"),
                          labels = c("True Positive", "True Negative", "False Positive", "False Negative"))

# Plot the heatmap with x and y labels, and add legend
hmap <- ggplot(plot_data, aes(x = Column, y = Row)) +
  geom_tile(aes(fill = Color), color = "black") +
  scale_fill_manual(values = c("True Positive" = "green", 
                               "True Negative" = "black", 
                               "False Positive" = "red", 
                               "False Negative" = "violet"),
                    labels = c(
                      "True Positive" = expression("True Positive (20.24%) "),
                      "True Negative" = expression("True Negative (75.12%)"),
                      "False Positive" = expression("False Positive (0.68%) "),
                      "False Negative" = expression("False Negative (3.96%)")
                    )) +
  labs(title = "",
       x = "Microbiomes (X)",  # X-axis label
       y = "Metabolites (simulated Y)",  # Y-axis label
       fill = ""
       ) +  # Add title for the legend
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 18), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "top",          
    legend.direction = "horizontal",  
    legend.justification = "center",  
    legend.box = "horizontal"
  ) 

hmap
ggsave("heatmap.jpeg", plot = hmap, width = 12, height = 6, dpi = 400)

################################################################################


MasterRanks <- as.matrix(read.csv('MasterRanks_BMASTER_estimated_NumIter_2000_burnin_100.csv',header=FALSE))
MasterRanks_est <- as.matrix(read.csv('Output_MasterRanks_rep_1.csv',header=FALSE))

array1 <- MasterRanks[,2]  # Replace with your data
array2 <- MasterRanks_est[,2]  # Replace with your data

data <- data.frame(
  Value = c(array1, array2),
  Group = factor(c(rep("True", length(array1)), rep("Estimated", length(array2))),
                 levels = c("True", "Estimated"))
)

MasterFreq <- ggplot(data, aes(x = Value, fill = Group)) +
  geom_histogram(binwidth = 0.8, color = "grey", alpha = 1.7) +
  scale_fill_manual(values = c("True" = "black", "Estimated" = "red")) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  labs(#title = "Frequency of number of Y components influenced",
       x = "Number of Y components influenced",
       y = "Frequency of X components") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),# face = "bold"),
    legend.position = "none",  # Hide legend if not needed
    #legend.position = "top",  # Position the legend at the top
    strip.text = element_text(size = 16, face = "bold"),  # Adjust font size of facet labels (subplot titles)
    strip.background = element_blank(),
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = 1   # Add border around the plot
    )
  ) +
  scale_y_continuous(breaks = seq(0, max(data$Value), by = 2))
MasterFreq
ggsave("MasterFreq.jpeg", plot = MasterFreq, width = 6, height = 5, dpi = 400)