rm(list=ls())
library(ggplot2)
library(GGally)
library(CCA)
library(CCP)
library(reshape2)

setwd("U:/BMASTER/Real Data Analysis")
extract_name <- function(input_string) {
  # Use regular expression to extract everything after the last underscore
  extracted_name <- sub(".*__", "", input_string)
  return(extracted_name)
}

extract_name_y <- function(input_string) {
  # Use regular expression to extract everything after the last underscore
  extracted_name <- sub(".*_", "", input_string)
  return(extracted_name)
}

X.names <- as.matrix(read.table("Data/X_names.csv", header = FALSE, sep = ",")) 
Y.names <- as.matrix(read.table("Data/Y_names.csv", header = FALSE, sep = ",")) 
NonZeroLocation <- as.matrix(read.table("Output_REAL_NonZeroLocation_NumIter_1000_burnin_100.csv", header = FALSE, sep = ",")) 
Signs <- sign(as.matrix(read.table("Output_REAL_B_NonSparse_NumIter_1000_burnin_100.csv", header = FALSE, sep = ",")))
MasterRanks.temp <- read.table("Output_REAL_MasterRanks_NumIter_1000_burnin_100.csv", header = FALSE, sep = ",")
MasterRanks <- MasterRanks.temp[, 1]

TopManyHist <- 50
MasterRanks.heat <- MasterRanks[1:TopManyHist]
NonZeroLocationWithSigns <- NonZeroLocation * Signs
Mat.selected <- as.matrix(NonZeroLocationWithSigns[MasterRanks.heat, ])
rownames(Mat.selected) <- extract_name(X.names[MasterRanks.heat])
colnames(Mat.selected) <- extract_name_y(Y.names)

write.table(rownames(Mat.selected),"TOP_50_master_predictors.csv", sep = ',', 
            row.names = FALSE, col.names = FALSE)

long_data <- melt(Mat.selected)
colnames(long_data) <- c("Row", "Column", "Value")

long_data$Row <- factor(long_data$Row, levels = rev(levels(long_data$Row)))

colors <- c("-1" = "blue", "0" = "white", "1" = "red")

# Plot the heatmap
whole_heat <- ggplot(long_data, aes(x = Column, y = Row, fill = as.factor(Value))) +
  geom_tile(color = "gray") +  # Adds borders between cells
  scale_fill_manual(values = colors, 
                    labels = c(expression(beta[pq] < 0), 
                               c(expression( beta[pq] == 0), expression( beta[pq] > 0)))) +  # Customize legend labels
  theme_minimal() +
  theme(axis.text.x = element_text(size = 2.6, face = "bold", angle = 90, hjust = 1),  # Hide x-axis labels
        axis.text.y = element_text(size = 5.4, face = "bold"),  # Bold row names
        axis.title.x = element_text(size = 14, vjust = 0), #element_text(size = 14, face = "bold"),  # Remove x-axis title
        axis.title.y =  element_text(size = 14),  # Remove y-axis title
        panel.grid = element_blank(),   # Remove gridlines
        plot.title = element_blank(), #element_text(hjust = 0.5, size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_blank(), #element_text(size = 18, face = "bold", color = "white"),
        legend.position = "top",          
        legend.direction = "horizontal",  
        legend.justification = "right",  
        legend.box = "horizontal",
        legend.key.width = unit(1.1, "cm"),  # Increase width of color swatch
        legend.spacing.x = unit(0.5, "cm"),
        legend.spacing.y = unit(2, "cm"),
        plot.margin = margin(5, 5, 5, 2),
        axis.ticks.length = unit(0.2, "cm")) +  # Center title
  labs(title = "Metabolites", fill = "Value", x = "Metabolites",       # X-axis title
       y = "Genera (top 50 master predictors)") +
  #coord_fixed(ratio = 3) +
  coord_cartesian(expand = FALSE)# Keep the aspect ratio square

whole_heat
ggsave("Top50genera_heat.jpeg", plot = whole_heat, width = 9, height = 5, dpi = 400)
####### Histogram plot #########################################################


X.names.ordered <- X.names[MasterRanks]
X.names.ordered.short <- extract_name(X.names.ordered)
Each.X.numYinfluenced <- MasterRanks.temp[, 2]
plot.upto <- 50

data <- data.frame(names = X.names.ordered.short[1:plot.upto], values = Each.X.numYinfluenced[1:plot.upto])
data$names <- factor(data$names, levels = X.names.ordered.short[1:plot.upto])
Top50genera <- ggplot(data, aes(x = names, y = values)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, face = "bold")) +  # Rotate X-axis labels vertically
  labs(x = "Genera (overall top 50 master predictors)", y = "Number of dependent metbolites", title = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title = element_text(size = 14),  # Bold axis titles
        axis.title.x = element_text(vjust = 177),
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add a border)

Top50genera 
ggsave("Top50genera.jpeg", plot = Top50genera, width = 9, height = 5, dpi = 400)

################################################################################

####### Heatmap plot ###########################################################



