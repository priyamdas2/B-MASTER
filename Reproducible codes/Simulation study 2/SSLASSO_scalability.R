rm(list = ls())
setwd("U:/BMASTER/Simulation study 2")

library(DescTools)
library(devtools)

# install.packages("U:/BMASTER/Simulation study 2/SSLASSO_1.1.tar.gz",
#                  repos = NULL, type = "source")
## devtools::install_github(repo = "skdeshpande91/multivariate_SSL/mSSL")

library(SSLASSO)
# library(mSSL)

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

cat("Running SSLASSO for:\n")
cat("rhoTag =", rhoTag, "\n")
cat("P =", p, ", Q =", q, ", N =", N, "\n\n")

## ============================================================
## SSLASSO analysis
## ============================================================

set.seed(RandSeed)
ptm <- proc.time()

lambda1 <- 1
nlambda <- 100
lambda0 <- seq(1, nrow(X), length.out = nlambda)
L <- length(lambda0)

## Store full SSLASSO path across lambda0
B_path_SSLASSO <- array(NA, dim = c(p, q, L))

## Store confusion counts for path-based AUC
TP_FN_FP_TN <- array(NA, dim = c(nlambda, 4, q))

for (j in 1:q) {
  
  cat("Running SSLASSO response j =", j, "of", q, "\n")
  
  fit_sslg <- SSLASSO(
    X,
    Y[, j],
    penalty = "adaptive",
    variance = "fixed",
    lambda1 = lambda1,
    lambda0 = lambda0,
    nlambda = nlambda,
    theta = 0.5,
    sigma = 1,
    a = 1,
    b = ncol(X),
    eps = 0.001,
    max.iter = 500,
    counter = 10,
    warn = FALSE
  )
  
  NZ_locTrue <- nonzero.locations[, j]
  
  for (l in 1:L) {
    
    ## Save coefficient path
    B_path_SSLASSO[, j, l] <- fit_sslg$beta[, l]
    
    ## Save support path for AUC
    NZ_locEst <- fit_sslg$select[, l]
    
    TP_here <- sum(NZ_locEst != 0 & NZ_locTrue == 1)
    FN_here <- sum(NZ_locEst == 0 & NZ_locTrue == 1)
    FP_here <- sum(NZ_locEst != 0 & NZ_locTrue == 0)
    TN_here <- sum(NZ_locEst == 0 & NZ_locTrue == 0)
    
    TP_FN_FP_TN[l, , j] <- c(TP_here, FN_here, FP_here, TN_here)
  }
}

## ============================================================
## Select final SSLASSO model using BIC-type criterion
## This replaces the old all-zero choice: fit_sslg$beta[, L]
## ============================================================

BIC.path <- rep(NA, L)

for (l in 1:L) {
  
  B_l <- B_path_SSLASSO[, , l]
  
  Y.pred.l <- X %*% B_l
  RSS.l <- sum((Y - Y.pred.l)^2)
  df.l <- sum(B_l != 0)
  
  ## Stabilize if RSS is numerically zero
  RSS.l <- max(RSS.l, .Machine$double.eps)
  
  BIC.path[l] <- N * q * log(RSS.l / (N * q)) + log(N * q) * df.l
}

best.l <- which.min(BIC.path)

cat("\nSelected SSLASSO lambda index =", best.l, "out of", L, "\n")
cat("Selected lambda0 =", lambda0[best.l], "\n")
cat("Selected BIC =", BIC.path[best.l], "\n\n")

B_SSLASSO <- B_path_SSLASSO[, , best.l]

time.taken <- proc.time() - ptm

## ============================================================
## AUC path over lambda0 grid
## ============================================================

TP_FN_FP_TN_sum <- matrix(NA, L, 4)
FPR.TPR <- matrix(NA, L, 2)

for (l in 1:L) {
  
  TP_FN_FP_TN_sum[l, ] <- apply(TP_FN_FP_TN[l, , 1:q], 1, sum)
  
  TP_l <- as.numeric(TP_FN_FP_TN_sum[l, 1])
  FN_l <- as.numeric(TP_FN_FP_TN_sum[l, 2])
  FP_l <- as.numeric(TP_FN_FP_TN_sum[l, 3])
  TN_l <- as.numeric(TP_FN_FP_TN_sum[l, 4])
  
  FPR.TPR[l, 1] <- FP_l / (FP_l + TN_l)
  FPR.TPR[l, 2] <- TP_l / (TP_l + FN_l)
}

## ============================================================
## Recovery metrics at selected SSLASSO model
## Current choice: BIC-selected lambda0 solution
## ============================================================

B_est <- B_SSLASSO

TP <- sum(B_est != 0 & nonzero.locations == 1)
FN <- sum(B_est == 0 & nonzero.locations == 1)
FP <- sum(B_est != 0 & nonzero.locations == 0)
TN <- sum(B_est == 0 & nonzero.locations == 0)

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

sparsity.prop <- sum(B_est == 0) / (p * q)
true.sparsity <- sum(nonzero.locations == 0) / (p * q)

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
    "_SSLasso",
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

cat("\nFinished SSLASSO.\n")
cat("Summary:\n")
cat("Selected lambda index =", best.l, "\n")
cat("Selected lambda0      =", lambda0[best.l], "\n")
cat("TPR      =", TPR.final, "\n")
cat("FPR      =", FPR.final, "\n")
cat("MCC      =", MCC.final, "\n")
cat("AUC      =", AUC, "\n")
cat("AUC20    =", AUC20, "\n")
cat("Sparsity =", sparsity.prop, "\n")
cat("True sparsity =", true.sparsity, "\n")
cat("MSE      =", MSE, "\n")
cat("Time     =", time.taken[3], "seconds\n")