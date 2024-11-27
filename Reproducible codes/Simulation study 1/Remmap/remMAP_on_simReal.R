rm(list=ls())
library(DescTools)
setwd("U:/BMASTER/Simulation study 1/Remmap")
# install.packages("U:/BMASTER/Simulation artificial/Remmap/remMap_0.2-0.tar.gz", repos = NULL, type = "source")
library(remMap)

set.seed(1)
## Reading data ################################################################
X <- as.matrix(read.csv('X.csv',header=FALSE))
B.True <- read.csv('B_BMASTER_estimated_NumIter_2000_burnin_100.csv', header=FALSE)
nonzero.locations <- read.csv('NonZeroLocation_BMASTER_estimated_NumIter_2000_burnin_100.csv', header=FALSE)
NumExp <- 10

for (dataset_no in 2:NumExp) {
  ptm <- proc.time()
  Y <- as.matrix(read.csv(paste0('Y_rep_', dataset_no,'.csv'),header=FALSE))
  
  # X <- X[1:50, 1:50]
  # Y <- Y[1:50, 1:50]
  # nonzero.locations <- nonzero.locations[1:50, 1:50]
  
  N <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
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
  
  time.taken <-  proc.time() - ptm
  
  B_est <- result0$phi
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
  
  num.picks <- length(as.vector(remMap_model$ols.cv))
  FPR.TPR <- matrix(NA, num.picks, 2)
  
  
  for (i in 1:num.picks) {
    print(i)
    lamL1.pick <- remMap_model$l.index[1,i]    
    lamL2.pick <- remMap_model$l.index[2,i]
    result <- remMap(X, Y,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)
    TP <- sum(result$phi != 0 & nonzero.locations == 1) 
    FN <- sum(result$phi == 0 & nonzero.locations == 1)
    FP <- sum(result$phi != 0 & nonzero.locations == 0) 
    TN <- sum(result$phi == 0 & nonzero.locations == 0)
    
    FPR <- FP / total.z
    TPR <- TP / total.nz
    FPR.TPR[i, ] <- c(FPR, TPR)
  }
  
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
  
  write.table(unique.FPR.TPR.sorted, paste0('Output_simReal_AUC_curve_remMap_rep_',
                                            dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
  
  write.table(t(output), paste0('Output_simReal_TPR_FPR_MCC_AUC_AUC20_sparsity_remMap_rep_',
                             dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
  write.table(time.taken[3], paste0('Output_compTime_remMap_rep_',
                                 dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
  write.table(MasterRanks, paste0('Output_MasterRanks_remMap_rep_',
                                  dataset_no,'.csv'), sep = ',', 
              row.names = FALSE, col.names = FALSE)
}
