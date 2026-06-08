rm(list=ls())
library(DescTools)
library(ggplot2)
setwd("U:/BMASTER/Simulation study 2")

P.array <- c(20, 50, 100, 200, 500, 1000, 2000)
Q.array <- P.array
N.array <- P.array
Num.cases <- length(P.array)

Summary_rhoNonZero <- matrix(NA, 7, 12)
Summary_rhoZero <- matrix(NA, 7, 12)

for (i in 1:Num.cases) {
  P <- P.array[i]
  Summary_rhoNonZero[i, 1:6] <- as.matrix(read.csv(paste0('Output_Scalability_rhoNonZero_summary_P_', P,'_NumIter_100_Nfactor_1.csv'),header=FALSE))
  Summary_rhoNonZero[i, 7] <- as.matrix(read.csv(paste0('Output_Scalability_rhoNonZero_TrueSparsity_P_', P,'_NumIter_100_Nfactor_1.csv'),header=FALSE))
  Summary_rhoNonZero[i, 8] <- P * P
  Summary_rhoNonZero[i, 9] <- (P * P)/(P.array[1] * P.array[1])
  Summary_rhoNonZero[i, 10] <- (1 / 60) * as.matrix(read.csv(paste0('Output_Scalability_rhoNonZero_compTime_P_', P,'_NumIter_100_Nfactor_1.csv'),header=FALSE))
  Summary_rhoNonZero[i, 10] <- round(Summary_rhoNonZero[i, 10], 3)
  Summary_rhoNonZero[i, 11] <- Summary_rhoNonZero[i, 10] / Summary_rhoNonZero[1, 10]
  Summary_rhoNonZero[i, 11] <- round(Summary_rhoNonZero[i, 11], 3)
}
Summary_rhoNonZero[, 12] <- 10 * Summary_rhoNonZero[, 10]

colnames(Summary_rhoNonZero) <- c("TPR", "FPR", "MCC", "AUC", "AUC20",
                                  "Sparsity", "True Sparsity", "PQ",
                                  "PQx", "Comp.time", "Comp.time x",
                                  "Projected time")
write.table(round(Summary_rhoNonZero, 4), "Summary_rhoNonZero.csv", sep = ",", col.names = TRUE, 
            row.names = FALSE)

for (i in 1:Num.cases) {
  P <- P.array[i]
  Summary_rhoZero[i, 1:6] <- as.matrix(read.csv(paste0('Output_Scalability_rhoZero_summary_P_', P,'_NumIter_100_Nfactor_1.csv'),header=FALSE))
  Summary_rhoZero[i, 7] <- as.matrix(read.csv(paste0('Output_Scalability_rhoZero_TrueSparsity_P_', P,'_NumIter_100_Nfactor_1.csv'),header=FALSE))
  Summary_rhoZero[i, 8] <- P * P
  Summary_rhoZero[i, 9] <- (P * P)/(P.array[1] * P.array[1])
  Summary_rhoZero[i, 10] <- (1 / 60) * as.matrix(read.csv(paste0('Output_Scalability_rhoZero_compTime_P_', P,'_NumIter_100_Nfactor_1.csv'),header=FALSE))
  Summary_rhoZero[i, 10] <- round(Summary_rhoZero[i, 10], 3)
  Summary_rhoZero[i, 11] <- Summary_rhoZero[i, 10] / Summary_rhoZero[1, 10]
  Summary_rhoZero[i, 11] <- round(Summary_rhoZero[i, 11], 3)
}
Summary_rhoZero[, 12] <- 10 * Summary_rhoZero[, 10]

colnames(Summary_rhoZero) <- c("TPR", "FPR", "MCC", "AUC", "AUC20",
                                  "Sparsity", "True Sparsity", "PQ",
                                  "PQx", "Comp.time", "Comp.time x",
                                  "Projected time")
write.table(round(Summary_rhoZero, 4), "Summary_rhoZero.csv", sep = ",", col.names = TRUE, 
            row.names = FALSE)


# Example data
data <- data.frame(
  parameters = Summary_rhoNonZero[, 8],
  computation_time = Summary_rhoNonZero[, 10]
)

# Transform data to log10 scale
data$log_parameters <- log10(data$parameters)
data$log_time <- log10(data$computation_time)

slope <- round((data$log_time[7] - data$log_time[1]) / 
                 (data$log_parameters[7] - data$log_parameters[1]), 2)
slope

# Create the plot
gg.rhoNonZero <- ggplot(data, aes(x = log_parameters, y = log_time)) +
  geom_line(color = "blue", size = 1) +  # Line connecting points
  geom_point(color = "red", size = 3) +  # Points on the curve
  labs(
    #title = "B-MASTER",
    x = expression(log[10]("Number of Parameters")),
    y = expression(log[10]("Computation Time"))
  ) +
  annotate(
    "text", 
    x = max(data$log_parameters) * 0.8,  # Position on the x-axis
    y = min(data$log_time) * 0.7,        # Position on the y-axis
    label = "Slope \u2248 1.06",         # Use Unicode for "approximately equals" sign
    size = 6,                            # Text size
    hjust = 0,                           # Horizontal alignment
    vjust = 0.5,                         # Vertical alignment
    color = "black"
  ) + 
  theme_minimal(base_size = 14) +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.title = element_text(face = "bold", size = 16),
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = 1   # Add border around the plot
    )
  )

gg.rhoNonZero
ggsave("CompTime_BMASTER_NonZero.jpeg", plot = gg.rhoNonZero, width = 6, height = 4.5, dpi = 400)



# # Example data
# data <- data.frame(
#   parameters = Summary_rhoZero[, 8],
#   computation_time = Summary_rhoZero[, 10]
# )
# 
# # Transform data to log10 scale
# data$log_parameters <- log10(data$parameters)
# data$log_time <- log10(data$computation_time)
# 
# # Create the plot
# gg.rhoZero <- ggplot(data, aes(x = log_parameters, y = log_time)) +
#   geom_line(color = "blue", size = 1) +  # Line connecting points
#   geom_point(color = "red", size = 3) +  # Points on the curve
#   labs(
#     #title = "B-MASTER",
#     x = expression(log[10]("Number of Parameters")),
#     y = expression(log[10]("Computation Time"))
#   ) +
#   theme_minimal(base_size = 14) +  # Clean theme
#   theme(
#     plot.title = element_text(hjust = 0.5),  # Center the title
#     axis.title = element_text(face = "bold")
#   )


