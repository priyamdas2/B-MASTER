rm(list = ls())
setwd("U:/BMASTER/Simulation study 2")

library(DescTools)
library(devtools)

## pak::pak(repo = "skdeshpande91/multivariate_SSL/mSSL")

library(mSSL)

## ============================================================
## User inputs
## ============================================================

RandSeed <- 1
IsRhoNonZero <- 0        # 0 means rho = 0, 1 means rho = 0.5
SampleMultFactor <- 1

P <- 20
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

cat("Running mSSL-dCPE for:\n")
cat("rhoTag =", rhoTag, "\n")
cat("P =", p, ", Q =", q, ", N =", N, "\n\n")

## ============================================================
## mSSL-dCPE analysis
## ============================================================

set.seed(RandSeed)
ptm <- proc.time()

fit_mSSL_dCpe <- mSSL_dcpe(X, Y)

B_mSSL_dpe <- fit_mSSL_dCpe$B

time.taken <- proc.time() - ptm

## ============================================================
## Recovery metrics at selected mSSL-dCPE model
## ============================================================

B_est <- B_mSSL_dpe

## Use rounded version for support calculation, consistent with AUC path below
B_est_support <- round(B_est, 5)

TP <- sum(B_est_support != 0 & nonzero.locations == 1)
FN <- sum(B_est_support == 0 & nonzero.locations == 1)
FP <- sum(B_est_support != 0 & nonzero.locations == 0)
TN <- sum(B_est_support == 0 & nonzero.locations == 0)

total.nz <- as.numeric(sum(nonzero.locations))
total.z <- as.numeric((p * q) - total.nz)

FPR.final <- as.numeric(FP) / total.z
TPR.final <- as.numeric(TP) / total.nz

aa <- as.numeric(TP + FP)
bb <- as.numeric(TP + FN)
cc <- as.numeric(TN + FP)
dd <- as.numeric(TN + FN)

if (any(c(aa, bb, cc, dd) == 0)) {
  MCC.final <- NA
} else {
  log.denom <- 0.5 * (log(aa) + log(bb) + log(cc) + log(dd))
  MCC.final <- (
    as.numeric(TP) * as.numeric(TN) -
      as.numeric(FP) * as.numeric(FN)
  ) / exp(log.denom)
}

Y.pred <- X %*% B_est
MSE <- sum((Y - Y.pred)^2) / (N * q)

sparsity.prop <- sum(B_est_support == 0) / (p * q)
true.sparsity <- sum(nonzero.locations == 0) / (p * q)

## ============================================================
## AUC path from mSSL-dCPE B0_path
## ============================================================

B_path <- fit_mSSL_dCpe$B0_path
L <- dim(B_path)[3]

TP_FN_FP_TN <- matrix(NA, L, 4)
FPR.TPR <- matrix(NA, L, 2)

for (l in 1:L) {
  
  B_est_here <- round(B_path[, , l], 5)
  
  TP_here <- sum(B_est_here != 0 & nonzero.locations == 1)
  FN_here <- sum(B_est_here == 0 & nonzero.locations == 1)
  FP_here <- sum(B_est_here != 0 & nonzero.locations == 0)
  TN_here <- sum(B_est_here == 0 & nonzero.locations == 0)
  
  TP_FN_FP_TN[l, ] <- c(TP_here, FN_here, FP_here, TN_here)
  
  TP_l <- as.numeric(TP_here)
  FN_l <- as.numeric(FN_here)
  FP_l <- as.numeric(FP_here)
  TN_l <- as.numeric(TN_here)
  
  FPR.TPR[l, 1] <- FP_l / (FP_l + TN_l)
  FPR.TPR[l, 2] <- TP_l / (TP_l + FN_l)
}

## ============================================================
## AUC calculation
## For duplicate FPR values, keep maximum TPR
## This avoids the "collapsing to unique 'x' values" warning
## ============================================================

sort.num <- order(FPR.TPR[, 1], FPR.TPR[, 2])
FPR.TPR.sorted <- rbind(c(0, 0), FPR.TPR[sort.num, ], c(1, 1))

FPR.TPR.sorted <- FPR.TPR.sorted[
  order(FPR.TPR.sorted[, 1], FPR.TPR.sorted[, 2]),
]

unique.FPR.values <- sort(unique(FPR.TPR.sorted[, 1]))

unique.FPR.TPR.sorted <- cbind(
  unique.FPR.values,
  sapply(unique.FPR.values, function(x) {
    max(FPR.TPR.sorted[FPR.TPR.sorted[, 1] == x, 2], na.rm = TRUE)
  })
)

AUC <- AUC(
  x = unique.FPR.TPR.sorted[, 1],
  y = unique.FPR.TPR.sorted[, 2],
  method = "trapezoid"
)

## ============================================================
## AUC20 calculation
## ============================================================

subset.for.AUC20_temp <- unique.FPR.TPR.sorted[
  which(unique.FPR.TPR.sorted[, 1] <= 0.2),
]

AUC20 <- NA

if (!is.null(dim(subset.for.AUC20_temp)) && nrow(subset.for.AUC20_temp) >= 2) {
  
  last.idx <- nrow(subset.for.AUC20_temp)
  subset.for.AUC20 <- subset.for.AUC20_temp
  
  if (
    subset.for.AUC20[last.idx, 1] < 0.2 &&
    last.idx < nrow(unique.FPR.TPR.sorted)
  ) {
    
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
## Same ordering as B-MASTER:
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
    "_mSSL_dCpe",
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

cat("\nFinished mSSL-dCPE.\n")
cat("Summary:\n")
cat("TPR      =", TPR.final, "\n")
cat("FPR      =", FPR.final, "\n")
cat("MCC      =", MCC.final, "\n")
cat("AUC      =", AUC, "\n")
cat("AUC20    =", AUC20, "\n")
cat("Sparsity =", sparsity.prop, "\n")
cat("True sparsity =", true.sparsity, "\n")
cat("MSE      =", MSE, "\n")
cat("Time     =", time.taken[3], "seconds\n")