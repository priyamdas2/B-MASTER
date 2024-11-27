rm(list=ls())
library(DescTools)
setwd("U:/BMASTER/Simulation study 2/remMap scalability")
# install.packages("U:/BMASTER/Simulation artificial/Remmap/remMap_0.2-0.tar.gz", repos = NULL, type = "source")
library(remMap)
library(ggplot2)

set.seed(1)
SampleMultFactor <- 10
P.array <- c(10, 20, 30, 40, 50, 60)
Num_P <- length(P.array)
time.taken <- rep(NA, Num_P)

for (i in 1:Num_P) {
  print(i)
  P <- P.array[i]
  Q <- P
  
  ptm <- proc.time()
  X <- as.matrix(read.csv(paste0('Data_X_P_',P,'_Nfactor_',SampleMultFactor,'.csv'),header=FALSE))
  Y <- as.matrix(read.csv(paste0('Data_Y_P_',P,'_Nfactor_',SampleMultFactor,'.csv'),header=FALSE))
  N <- dim(X)[1]
  ##############################################################################
  
  ####### remMap analysis ######################################################
  
  lamL1.v <- exp(seq(log(1), log(1000), length = 11))
  lamL2.v <- seq(0, 1000, length = 11)
  
  set.seed(1)
  remMap_model <- remMap.CV(X=X, Y=Y,lamL1.v, lamL2.v, C.m = NULL, fold = 5, seed = 1)
  pick <- which.min(as.vector(remMap_model$ols.cv))
  lamL1.pick <- remMap_model$l.index[1,pick]    ##find the optimal (LamL1,LamL2) based on the cv score
  lamL2.pick <- remMap_model$l.index[2,pick]
  
  
  result0 <- remMap(X, Y,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)
  
  time.taken[i] <-  proc.time() - ptm
  print(time.taken[i])
}
time.taken
time.summary <- rbind(P.array, time.taken)
write.table(time.summary, paste0('Time_Remmap_Nfactor_',SampleMultFactor,'.csv'), sep = ",", col.names = FALSE, 
            row.names = FALSE)

# # Example data
# data <- data.frame(
#   parameters = P.array ^ 2,
#   computation_time = time.taken
# )
# 
# # Transform data to log10 scalei = 1
# 
# data$log_parameters <- log10(data$parameters)
# data$log_time <- log10(data$computation_time)
# 
# # Create the plot
# gg.rhoNonZero.remMap <- ggplot(data, aes(x = log_parameters, y = log_time)) +
#   geom_line(color = "blue", size = 1) +  # Line connecting points
#   geom_point(color = "red", size = 3) +  # Points on the curve
#   labs(
#     title = "remMap",
#     x = expression(log[10]("Number of Parameters")),
#     y = expression(log[10]("Computation Time"))
#   ) +
#   theme_minimal(base_size = 14) +  # Clean theme
#   theme(
#     plot.title = element_text(hjust = 0.5),  # Center the title
#     axis.title = element_text()#face = "bold")
#   )
# 
# 
# 
# gg.rhoNonZero.remMap
# ggsave("CompTime_remMap_NonZero.jpeg", plot = gg.rhoNonZero.remMap, width = 6, height = 4.5, dpi = 400)
# 
