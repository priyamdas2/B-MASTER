rm(list = ls())
setwd("U:/BMASTER/Simulation study 2")
library(DescTools)
#install.packages("U:/BMASTER/Simulation study 2/remMap_0.2-0.tar.gz", repos = NULL, type = "source")
library(remMap)

## ============================================================
## User inputs
## ============================================================

RandSeed <- 1
IsRhoNonZero <- 1        # 0 means rho = 0, 1 means rho = 0.5
SampleMultFactor <- 1

P <- 200
Q <- P
N <- P * SampleMultFactor

## ============================================================
## Folders
## ============================================================

DataFolder <- "Data"
OutputFolder <- "Output other methods"

if (!dir.exists(OutputFolder)) {
  dir.create(OutputFolder)
}

if (IsRhoNonZero == 1) {
  rhoTag <- "rhoNonZero"
  rhoValue <- 0.5
} else {
  rhoTag <- "rhoZero"
  rhoValue <- 0
}

dataPrefix <- file.path(
  DataFolder,
  paste0(
    "SimData_", rhoTag,
    "_P_", P,
    "_Q_", Q,
    "_N_", N,
    "_Nfactor_", SampleMultFactor,
    "_Seed_", RandSeed
  )
)

## ============================================================
## Read simulated data generated from MATLAB
## ============================================================

X.file <- paste0(dataPrefix, "_X.csv")
Y.file <- paste0(dataPrefix, "_Y.csv")
Support.file <- paste0(dataPrefix, "_Support_True.csv")

X <- as.matrix(read.csv(X.file, header = FALSE))
Y <- as.matrix(read.csv(Y.file, header = FALSE))
nonzero.locations <- as.matrix(read.csv(Support.file, header = FALSE))

N <- dim(X)[1]
p <- dim(X)[2]
q <- dim(Y)[2]

cat("Running remMap-BIC for:\n")
cat("rhoTag =", rhoTag, "\n")
cat("P =", p, ", Q =", q, ", N =", N, "\n\n")

## ============================================================
## remMap-BIC analysis
## ============================================================

set.seed(RandSeed)
ptm <- proc.time()

lamL1.v <- exp(seq(log(1), log(1000), length = 11))
lamL2.v <- seq(0, 1000, length = 11)

remMap_model <- remMap.BIC(
  X = X,
  Y = Y,
  lamL1.v = lamL1.v,
  lamL2.v = lamL2.v,
  C.m = NULL
)

pick <- which.min(as.vector(t(remMap_model$BIC)))
result0 <- remMap_model$phi[[pick]]

time.taken <- proc.time() - ptm

B_est <- result0$phi

## ============================================================
## Recovery metrics at BIC-selected model
## ============================================================

TP <- sum(B_est != 0 & nonzero.locations == 1)
FN <- sum(B_est == 0 & nonzero.locations == 1)
FP <- sum(B_est != 0 & nonzero.locations == 0)
TN <- sum(B_est == 0 & nonzero.locations == 0)

total.nz <- sum(nonzero.locations)
total.z <- (p * q) - total.nz

FPR.final <- as.numeric(FP) / as.numeric(total.z)
TPR.final <- as.numeric(TP) / as.numeric(total.nz)

aa <- as.numeric(TP + FP)
bb <- as.numeric(TP + FN)
cc <- as.numeric(TN + FP)
dd <- as.numeric(TN + FN)

# denom <- sqrt(aa * bb * cc * dd)

if (any(c(aa, bb, cc, dd) == 0)) {
  MCC.final <- NA
} else {
  log.denom <- 0.5 * (log(aa) + log(bb) + log(cc) + log(dd))
  MCC.final <- (as.numeric(TP) * as.numeric(TN) - as.numeric(FP) * as.numeric(FN)) / exp(log.denom)
}

Y.pred <- X %*% B_est
MSE <- sum((Y - Y.pred)^2) / (N * q)

sparsity.prop <- sum(B_est == 0) / (p * q)
true.sparsity <- sum(nonzero.locations == 0) / (p * q)

## ============================================================
## Master predictor ranks
## ============================================================

NumYinfluencedbyEachX <- apply(B_est != 0, 1, sum)
order.col <- order(NumYinfluencedbyEachX, decreasing = TRUE)
NumYinfluencedSorted.col <- sort(NumYinfluencedbyEachX, decreasing = TRUE)
MasterRanks <- cbind(order.col, NumYinfluencedSorted.col)

## ============================================================
## AUC and AUC20 over full remMap tuning grid
## ============================================================

num.picks <- length(as.vector(t(remMap_model$BIC)))
FPR.TPR <- matrix(NA, num.picks, 2)

for (i in 1:num.picks) {
  result <- remMap_model$phi[[i]]
  B_tmp <- result$phi
  
  TP <- sum(B_tmp != 0 & nonzero.locations == 1)
  FN <- sum(B_tmp == 0 & nonzero.locations == 1)
  FP <- sum(B_tmp != 0 & nonzero.locations == 0)
  TN <- sum(B_tmp == 0 & nonzero.locations == 0)
  
  FPR <- FP / total.z
  TPR <- TP / total.nz
  
  FPR.TPR[i, ] <- c(FPR, TPR)
}

sort.num <- order(FPR.TPR[, 1], FPR.TPR[, 2])
FPR.TPR.sorted <- rbind(c(0, 0), FPR.TPR[sort.num, ])
unique.FPR.TPR.sorted <- unique(FPR.TPR.sorted)
unique.FPR.TPR.sorted <- rbind(unique.FPR.TPR.sorted, c(1, 1))

AUC <- AUC(
  x = unique.FPR.TPR.sorted[, 1],
  y = unique.FPR.TPR.sorted[, 2],
  method = "trapezoid"
)

## AUC20 calculation
subset.for.AUC20_temp <- unique.FPR.TPR.sorted[
  which(unique.FPR.TPR.sorted[, 1] <= 0.2), 
]

AUC20 <- NA

if (!is.null(dim(subset.for.AUC20_temp)) && nrow(subset.for.AUC20_temp) >= 2) {
  
  last.idx <- nrow(subset.for.AUC20_temp)
  subset.for.AUC20 <- subset.for.AUC20_temp
  
  if (subset.for.AUC20[last.idx, 1] < 0.2 && last.idx < nrow(unique.FPR.TPR.sorted)) {
    
    y1 <- unique.FPR.TPR.sorted[last.idx, 2]
    y2 <- unique.FPR.TPR.sorted[last.idx + 1, 2]
    x1 <- unique.FPR.TPR.sorted[last.idx, 1]
    x2 <- unique.FPR.TPR.sorted[last.idx + 1, 1]
    
    a <- 0.2
    
    if (x2 != x1) {
      val.at.20 <- y1 + (a - x1) * (y2 - y1) / (x2 - x1)
      subset.for.AUC20 <- rbind(subset.for.AUC20_temp, c(0.2, val.at.20))
    }
  }
  
  AUC20 <- 5 * AUC(
    x = subset.for.AUC20[, 1],
    y = subset.for.AUC20[, 2],
    method = "trapezoid"
  )
}

## ============================================================
## Final output
## Same ordering as B-MASTER table:
## TPR, FPR, MCC, AUC, AUC20, estimated sparsity
## ============================================================

output <- c(
  TPR.final,
  FPR.final,
  MCC.final,
  AUC,
  AUC20,
  sparsity.prop
)

output.with.truth <- c(
  TPR.final,
  FPR.final,
  MCC.final,
  AUC,
  AUC20,
  sparsity.prop,
  true.sparsity,
  MSE,
  time.taken[3]
)

outPrefix <- file.path(
  OutputFolder,
  paste0(
    "Output_", rhoTag,
    "_remMapBic",
    "_P_", p,
    "_Q_", q,
    "_N_", N,
    "_Nfactor_", SampleMultFactor,
    "_Seed_", RandSeed
  )
)

write.table(
  t(output),
  paste0(outPrefix, "_TPR_FPR_MCC_AUC_AUC20_sparsity.csv"),
  sep = ",",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  t(output.with.truth),
  paste0(outPrefix, "_summary_with_truth_MSE_time.csv"),
  sep = ",",
  row.names = FALSE,
  col.names = FALSE
)



cat("\nFinished remMap-BIC.\n")
cat("Summary:\n")
cat("TPR      =", TPR.final, "\n")
cat("FPR      =", FPR.final, "\n")
cat("MCC      =", MCC.final, "\n")
cat("AUC      =", AUC, "\n")
cat("AUC20    =", AUC20, "\n")
cat("Sparsity =", sparsity.prop, "\n")
cat("True sparsity =", true.sparsity, "\n")
cat("Time =", time.taken[3], "seconds\n")