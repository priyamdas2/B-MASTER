rm(list=ls())
library(DescTools)
setwd("U:/BMASTER/Simulation study 2/remMap scalability")
library(ggplot2)

SampleMultFactorArray <- c(1, 5, 10)
Mat.Remmap_1 <- as.matrix(read.csv(paste0('Time_Remmap_Nfactor_',SampleMultFactorArray[1],'.csv'),header=FALSE))
Mat.Remmap_2 <- as.matrix(read.csv(paste0('Time_Remmap_Nfactor_',SampleMultFactorArray[2],'.csv'),header=FALSE))
Mat.Remmap_3 <- as.matrix(read.csv(paste0('Time_Remmap_Nfactor_',SampleMultFactorArray[3],'.csv'),header=FALSE))

Mat.BMASTER_1 <- as.matrix(read.csv(paste0('Time_BMASTER_Nfactor_',SampleMultFactorArray[1],'.csv'),header=FALSE))
Mat.BMASTER_2 <- as.matrix(read.csv(paste0('Time_BMASTER_Nfactor_',SampleMultFactorArray[2],'.csv'),header=FALSE))
Mat.BMASTER_3 <- as.matrix(read.csv(paste0('Time_BMASTER_Nfactor_',SampleMultFactorArray[3],'.csv'),header=FALSE))

time_points <- as.numeric(Mat.Remmap_1[1, ])

data <- data.frame(
  P = rep(time_points, times = 3),
  Data_Size = rep(c("Small", "Medium", "Large"), each = 6),
  Time = c(Mat.Remmap_1[2, ],   # Small dataset
           Mat.Remmap_2[2, ],   # Medium dataset
           Mat.Remmap_3[2, ])  # Large dataset
)
data$Data_Size <- factor(data$Data_Size, levels = c("Small", "Medium", "Large"))

# Create the plot
plotRemMap <- ggplot(data, aes(x = P, y = Time, color = Data_Size, group = Data_Size)) +
  geom_line(size = 1) +                    
  geom_point(size = 2) +                   
  scale_color_manual(values = c("blue", "green", "red"),
                     labels = c("N = P", "N = 5P", "N = 10P")) + 
  scale_x_continuous(
    breaks = time_points                       # Keep x-axis limits
  ) +
  labs(title = "RemMap",
       x = "P (= Q)",
       y = "Computation Time (seconds)",
       color = "Sample Size") +           # Legend title
  ylim(0, 30) +
  theme_minimal(base_size = 15) +          # Professional theme
  theme(
    plot.title = element_text(
      size = 20, 
      face = "bold", 
      hjust = 0.5  # Center the title
    ),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.position = "top" ,
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = 1   # Add border around the plot
    )
  ) 

plotRemMap

ggsave("CompTime_remMap_sample_size.jpeg", plot = plotRemMap, width = 8, height = 6, dpi = 400)

################################################################################

time_points <- as.numeric(Mat.BMASTER_1[1, ])

data <- data.frame(
  P = rep(time_points, times = 3),
  Data_Size = rep(c("Small", "Medium", "Large"), each = 6),
  Time = c(Mat.BMASTER_1[2, ],   # Small dataset
           Mat.BMASTER_2[2, ],   # Medium dataset
           Mat.BMASTER_3[2, ])  # Large dataset
)
data$Data_Size <- factor(data$Data_Size, levels = c("Small", "Medium", "Large"))

# Create the plot
plotBMASTER <- ggplot(data, aes(x = P, y = Time, color = Data_Size, group = Data_Size)) +
  geom_line(size = 1) +                    
  geom_point(size = 2) +                   
  scale_color_manual(values = c("blue", "green", "red"),
                     labels = c("N = P", "N = 5P", "N = 10P")) + 
  scale_x_continuous(
    breaks = time_points                       # Keep x-axis limits
  ) +
  labs(title = "B-MASTER",
       x = "P (= Q)",
       y = "Computation Time (seconds)",
       color = "Sample Size") +           # Legend title
  ylim(0, 30) +
  theme_minimal(base_size = 15) +          # Professional theme
  theme(
    plot.title = element_text(
      size = 20, 
      face = "bold", 
      hjust = 0.5  # Center the title
    ),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.position = "top" ,
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = 1   # Add border around the plot
    )
  ) 

plotBMASTER

ggsave("CompTime_BMASTER_sample_size.jpeg", plot = plotBMASTER, width = 8, height = 6, dpi = 400)