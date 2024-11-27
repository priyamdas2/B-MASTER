rm(list=ls())
setwd("U:/BMASTER/Simulation study 1/Summary Table and Plots")
library(ggplot2)
library(splines)
library(stats)
library(fda)
library(zoo)

NumExp <- 10
Summary_BMASTER <- matrix(NA, 10, 6)
Summary_SSLASSO <- matrix(NA, 10, 6)
Summary_mSSL <- matrix(NA, 10, 6)
Summary_remMap <- matrix(NA, 10, 6)



for (i in 1:NumExp) {
  # Construct the file name
  filename <- paste0("Output_simReal_TPR_FPR_MCC_AUC_AUC20_sparsity_BMASTER_rep_", i, ".csv")
  Summary_BMASTER[i, ] <- t(t(read.table(filename, header = FALSE, sep = ",")))
  
  filename <- paste0("Output_simReal_TPR_FPR_MCC_AUC_AUC20_sparsity_SSLasso_rep_", i, ".csv")
  Summary_SSLASSO[i, ] <- t(t(read.table(filename, header = FALSE, sep = ",")))
  
  filename <- paste0("Output_simReal_TPR_FPR_MCC_AUC_AUC20_sparsity_mSSL_dCpe_rep_", i, ".csv")
  Summary_mSSL[i, ] <- t(t(read.table(filename, header = FALSE, sep = ",")))
  
  filename <- paste0("Output_simReal_TPR_FPR_MCC_AUC_AUC20_sparsity_remMapBic_rep_", i, ".csv")
  Summary_remMap[i, ] <- t(t(read.table(filename, header = FALSE, sep = ",")))
}


Table_mean <- matrix(NA, 4, 6)
Table_sd <- matrix(NA, 4, 6)

Table_mean[1, ] <- apply(Summary_BMASTER, 2, mean)
Table_mean[2, ] <- apply(Summary_SSLASSO, 2, mean)
Table_mean[3, ] <- apply(Summary_mSSL, 2, mean)
Table_mean[4, ] <- apply(Summary_remMap, 2, mean)

Table_sd[1, ] <- apply(Summary_BMASTER, 2, sd)
Table_sd[2, ] <- apply(Summary_SSLASSO, 2, sd) 
Table_sd[3, ] <- apply(Summary_mSSL, 2, sd)
Table_sd[4, ] <- apply(Summary_remMap, 2, sd)

colnames(Table_mean) <- c("TPR", "FPR", "MCC", "AUC", "AUC20", "Sparsity")
colnames(Table_sd) <- c("TPR_sd", "FPR_sd", "MCC_sd", "AUC_sd", "AUC20_sd", "Sparsity_sd")

write.table(round(Table_mean, 2), "AA_summary_mean.csv", sep = ",", col.names = TRUE, 
            row.names = FALSE)

write.table(round(Table_sd, 4), "AA_summary_sd.csv", sep = ",", col.names = TRUE, 
            row.names = FALSE)


######### AUC Plot #############################################################


AUC_BMASTER <-  read.table("Output_simReal_AUC_curve_BMASTER_rep_1.csv",
                           header = FALSE, sep = ",") 
AUC_SSLASSO <-  read.table("Output_simReal_AUC_curve_SSLASSO_rep_1.csv",
                           header = FALSE, sep = ",") 
AUC_mSSL <-  read.table("Output_simReal_AUC_curve_mSSL_dCpe_rep_1.csv",
                           header = FALSE, sep = ",") 
AUC_remMap <-  read.table("Output_simReal_AUC_curve_remMapBic_rep_1.csv",
                           header = FALSE, sep = ",") 

# Simulating data: 4 methods, each with FPR and TPR
cutpoint <- 0.2
fpr0 <- AUC_BMASTER[, 1]
tpr0 <- AUC_BMASTER[, 2]
fpr <- fpr0[which(fpr0 <= cutpoint)]
indx <- match(fpr, fpr0)
tpr <- tpr0[indx]
fpr_new <- c(fpr, cutpoint)
tpr_new <- c(tpr, NA)
tpr_new2 <- na.approx(tpr_new, fpr_new, rule = 2) 
method1 <- data.frame(FPR = fpr_new, TPR = tpr_new2)


fpr0 <- AUC_SSLASSO[,1]
tpr0 <- AUC_SSLASSO[,2]
fpr <- fpr0[which(fpr0 <= cutpoint)]
indx <- match(fpr, fpr0)
tpr <- tpr0[indx]
fpr_new <- c(fpr, cutpoint)
tpr_new <- c(tpr, NA)
tpr_new2 <- na.approx(tpr_new, fpr_new, rule = 2) 
method2 <- data.frame(FPR = fpr_new, TPR = tpr_new2)


fpr0 <- AUC_mSSL[,1]
tpr0 <- AUC_mSSL[,2]
fpr <- fpr0[which(fpr0 <= cutpoint)]
indx <- match(fpr, fpr0)
tpr <- tpr0[indx]
fpr_new <- c(fpr, cutpoint)
tpr_new <- c(tpr, NA)
tpr_new2 <- na.approx(tpr_new, fpr_new, rule = 2) 
method3 <- data.frame(FPR = fpr_new, TPR = tpr_new2)

fpr0 <- AUC_remMap[,1]
tpr0 <- AUC_remMap[,2]
fpr <- fpr0[which(fpr0 <= cutpoint)]
indx <- match(fpr, fpr0)
tpr <- tpr0[indx]
fpr_1 <- unique(fpr)
indx2 <- match(fpr_1, fpr)
tpr_1 <- tpr[indx2]
fpr_new <- c(fpr_1, cutpoint)
tpr_new <- c(tpr_1, NA)
tpr_new2 <- na.approx(tpr_new, fpr_new, rule = 2) 
smooth_model <- smooth.spline(fpr_new, tpr_new2)

# Enforce monotonicity
smooth_values <- smooth_model$y
for (i in 2:length(smooth_values)) {
  if (smooth_values[i] < smooth_values[i - 1]) {
    smooth_values[i] <- smooth_values[i - 1]
  }
}
smooth_values[which(smooth_values < 0)] <- 0
tpr_new3 <- smooth_values
method4 <- data.frame(FPR = fpr_new, TPR = tpr_new3)

# fpr <- AUC_remMap[,1]
# tpr <- AUC_remMap[,2]
# method4 <- data.frame(FPR = fpr, TPR = tpr)



# Add method identifiers
method1$Method <- "B-MASTER (0.94)"
method2$Method <- "SSLASSO (0.54)"
method3$Method <- "mSSL (0.44)"
#method4$Method <- "remMap (0.23)"

roc_data <- rbind(method1, method2, method3)
roc_data$Method <- factor(roc_data$Method, levels = c("B-MASTER (0.94)", "SSLASSO (0.54)", "mSSL (0.44)"))


# Plot using ggplot2
AUC_plot <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Method, linetype = Method)) +
  geom_line(size = 1.5, linetype = "solid") +
  labs(title = "",
       x = "False Positive Rate (FPR)",
       y = "True Positive Rate (TPR)",
       color = "Method (AUC20)") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "darkolivegreen4")) +
  theme(legend.text = element_text(size = 9),
        #legend.title = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.direction = "vertical",  
        legend.position = c(.955, 0.044),  # Position the legend inside at bottom-right corner
        legend.justification = c("right", "bottom"),
        #legend.justification = "center",  
        legend.box = "horizontal",
        legend.spacing.x = unit(0.3, "cm"), 
        legend.margin = margin(t = 2.5, r = 2.5, b = 2.5, l = 2.5),
        axis.title.y = element_text(size = 15, margin = margin(l = 5), vjust = 5),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)), 
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        plot.margin = margin(t = 10, r = 7, b = 7, l = 20),
        panel.border = element_rect(
          color = "black", 
          fill = NA, 
          size = 1   # Add border around the plot
        ),
        legend.background = element_rect(fill = "grey", color = "white")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")

AUC_plot
ggsave("AUC_plot.jpeg", plot = AUC_plot, width = 6, height = 5, dpi = 400)