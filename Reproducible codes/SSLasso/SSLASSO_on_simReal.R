rm(list=ls())
library(DescTools)
library(devtools)
setwd("U:/BMASTER/Simulation study 1/SSLasso")
# install.packages("U:/BMASTER/Simulation artificial/SSLasso/SSLASSO_1.1.tar.gz", repos = NULL, type = "source")
# devtools::install_github(repo = "skdeshpande91/multivariate_SSL/mSSL")
library(SSLASSO)
library(mSSL)

set.seed(1)
## Reading data ################################################################
X <- as.matrix(read.csv('X.csv',header=FALSE))
B.True <- read.csv('B_BMASTER_estimated_NumIter_2000_burnin_100.csv', header=FALSE)
nonzero.locations <- read.csv('NonZeroLocation_BMASTER_estimated_NumIter_2000_burnin_100.csv', header=FALSE)
NumExp <- 10

for (dataset_no in 1:NumExp) {
  print(dataset_no)
  ptm <- proc.time()
  Y <- as.matrix(read.csv(paste0('Y_rep_', dataset_no,'.csv'),header=FALSE))
  N <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  ##############################################################################
  
  ####### SSLASSO analysis ######################################################
  
  #fit_sslg <- sep_sslg(X,Y)
  lambda1 <-  1
  nlambda <- 100
  lambda0 <- seq(1, nrow(X), length.out = nlambda)
  L <- length(lambda0)
  B_SSLASSO <- matrix(NA, p, q)
  
  TP_FN_FP_TN <- array(NA, dim = c(nlambda, 4, q))
  
  for (j in 1:q) {
    fit_sslg <- SSLASSO(X, Y[, j], penalty = "adaptive", variance = "fixed", 
                        lambda1, lambda0, nlambda, theta = 0.5, sigma = 1, a = 1, b = ncol(X),  
                        eps = 0.001, max.iter = 500,  counter = 10, warn = FALSE)
    
    NZ_locTrue <- nonzero.locations[, j]
    for (l in 1:L) {
      NZ_locEst <- fit_sslg$select[, l]
      TP_here <- sum(NZ_locEst != 0 & NZ_locTrue == 1) 
      FN_here <- sum(NZ_locEst == 0 & NZ_locTrue == 1)
      FP_here <- sum(NZ_locEst != 0 & NZ_locTrue == 0) 
      TN_here <- sum(NZ_locEst == 0 & NZ_locTrue == 0)
      TP_FN_FP_TN[l, ,j] <- c(TP_here, FN_here, FP_here, TN_here)
    }
    B_SSLASSO[, j] <- fit_sslg$beta[, L]
  }
  time.taken <-  proc.time() - ptm  
  
  TP_FN_FP_TN_sum <- matrix(NA, L, 4)
  FPR.TPR <- matrix(NA, L, 2)
  for (l in 1:L) {
    TP_FN_FP_TN_sum[l, ] <- apply(TP_FN_FP_TN[l, ,1:q], 1, sum)
    FPR.TPR[l, 1] <- TP_FN_FP_TN_sum[l, 3] / (TP_FN_FP_TN_sum[l, 3] + TP_FN_FP_TN_sum[l, 4])
    FPR.TPR[l, 2] <- TP_FN_FP_TN_sum[l, 1] / (TP_FN_FP_TN_sum[l, 1] + TP_FN_FP_TN_sum[l, 2])
  }
  
  

  
  B_est <- B_SSLASSO
  TP <- sum(B_est != 0 & nonzero.locations == 1) 
  FN <- sum(B_est == 0 & nonzero.locations == 1)
  FP <- sum(B_est != 0 & nonzero.locations == 0) 
  TN <- sum(B_est == 0 & nonzero.locations == 0)
  
  
  total.nz <- sum(nonzero.locations)
  total.z <- (p * q) - total.nz
  
  aa <- TP + FP
  bb <- TP + FN
  cc <- TN + FP
  dd <- TN + FN
  
  
  FPR.final <- FP / total.z
  TPR.final <- TP / total.nz
  MCC.final <- (TP * TN - FP * FN) / sqrt(exp(log(aa) + log(bb) + log(cc) + log(dd)))
  Y.pred <- X %*% B_est
  MSE <- sum((Y - Y.pred)^2) / (N*q)
  sparsity.prop <- length(which(B_est == 0)) / (p * q)
  
  ## Finding Master predictors
  
  NumYinfluencedbyEachX <- apply(B_est != 0, 1,sum)
  order.col <- order(NumYinfluencedbyEachX, decreasing = TRUE)
  NumYinfluencedSorted.col <- sort(NumYinfluencedbyEachX, decreasing = TRUE)
  MasterRanks <- cbind(order.col, NumYinfluencedSorted.col)
  
  
  ## AUC calculation ###########################################################
  
  sort.num <- order(FPR.TPR[, 1])
  FPR.TPR.sorted <- rbind(c(0,0), FPR.TPR[sort.num, ])
  unique.FPR.TPR.sorted <- unique(FPR.TPR.sorted)
  unique.FPR.TPR.sorted <- rbind(unique.FPR.TPR.sorted,c(1,1))
  
  
  
  AUC <- AUC(x = unique.FPR.TPR.sorted[, 1], y = unique.FPR.TPR.sorted[, 2], method="trapezoid")
  
  
  #### AUC20 calculation ######################################################
  
  subset.for.AUC20_temp <- unique.FPR.TPR.sorted[which(unique.FPR.TPR.sorted[,1] <= 0.2),]
  last.idx <- dim(subset.for.AUC20_temp)[1]
  AUC20 <- NA
  if (length(subset.for.AUC20_temp) > 2) {
    subset.for.AUC20 <- subset.for.AUC20_temp
    
    if(subset.for.AUC20[last.idx,1] < 0.2) {
      y1 <- unique.FPR.TPR.sorted[last.idx, 2]
      y2 <- unique.FPR.TPR.sorted[last.idx + 1, 2]
      x1 <- unique.FPR.TPR.sorted[last.idx, 1]
      x2 <- unique.FPR.TPR.sorted[last.idx + 1, 1]
      a <- 0.2
      val.at.20 <-  y1 + (a - x1) * (y2 - y1)/(x2-x1)
      subset.for.AUC20 <- rbind(subset.for.AUC20_temp, c(0.2, val.at.20))
    }
    
    AUC20 <- 5 * AUC(x = subset.for.AUC20[, 1], y = subset.for.AUC20[, 2], method="trapezoid")
  }
  
  ###### Final summary #########################################################
  
  output <- c(TPR.final, FPR.final, MCC.final, AUC, AUC20, sparsity.prop)
  
  write.table(unique.FPR.TPR.sorted, paste0('Output_simReal_AUC_curve_SSLasso_rep_',
                                            dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
  
  write.table(t(output), paste0('Output_simReal_TPR_FPR_MCC_AUC_AUC20_sparsity_SSLasso_rep_',
                                dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
  write.table(time.taken[3], paste0('Output_compTime_SSLasso_rep_',
                                    dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
  write.table(MasterRanks, paste0('Output_MasterRanks_SSLasso_rep_',
                                    dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
  
}