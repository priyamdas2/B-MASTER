rm(list = ls())

## ============================================================
## User inputs
## ============================================================

P.values <- c(20, 50, 100, 200)
SampleMultFactor <- 1
RandSeed <- 1

OutputOtherFolder <- "Output other methods"

## ============================================================
## Helper functions
## ============================================================

read_one_row <- function(file) {
  x <- as.numeric(read.csv(file, header = FALSE)[1, ])
  return(x)
}

safe_round3 <- function(x) {
  round(as.numeric(x), 3)
}

find_single_file <- function(pattern, folder = ".") {
  files <- list.files(
    path = folder,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    warning("No file found for pattern: ", pattern)
    return(NA)
  }
  
  if (length(files) > 1) {
    warning("Multiple files found for pattern: ", pattern,
            "\nUsing first one: ", files[1])
  }
  
  return(files[1])
}

## ============================================================
## Function to build one rho-specific table
## ============================================================

build_summary_table <- function(rhoTag) {
  
  all_rows <- data.frame()
  
  for (P in P.values) {
    
    Q <- P
    N <- P * SampleMultFactor
    
    ## --------------------------------------------------------
    ## 1. B-MASTER
    ## Expected file example:
    ## Output_Scalability_rhoZero_summary_P_20_NumIter_100_Nfactor_1.csv
    ## or NumIter_2000 etc.
    ## --------------------------------------------------------
    
    bmaster_pattern <- paste0(
      "^Output_Scalability_", rhoTag,
      "_summary_P_", P,
      "_NumIter_.*_Nfactor_", SampleMultFactor,
      "\\.csv$"
    )
    
    bmaster_file <- find_single_file(bmaster_pattern, folder = ".")
    
    if (!is.na(bmaster_file)) {
      vals <- read_one_row(bmaster_file)
      
      row_here <- data.frame(
        P_Q_N = P,
        Method = "B-MASTER",
        TPR = safe_round3(vals[1]),
        FPR = safe_round3(vals[2]),
        MCC = safe_round3(vals[3]),
        AUC = safe_round3(vals[4]),
        AUC20 = safe_round3(vals[5]),
        Sparsity = safe_round3(vals[6])
      )
      
      all_rows <- rbind(all_rows, row_here)
    }
    
    ## --------------------------------------------------------
    ## 2. SSLASSO
    ## Expected file:
    ## Output other methods/Output_rhoZero_SSLasso_P_20_Q_20_N_20_Nfactor_1_Seed_1_TPR_FPR_MCC_AUC_AUC20_sparsity.csv
    ## --------------------------------------------------------
    
    sslasso_file <- file.path(
      OutputOtherFolder,
      paste0(
        "Output_", rhoTag,
        "_SSLasso",
        "_P_", P,
        "_Q_", Q,
        "_N_", N,
        "_Nfactor_", SampleMultFactor,
        "_Seed_", RandSeed,
        "_TPR_FPR_MCC_AUC_AUC20_sparsity.csv"
      )
    )
    
    if (file.exists(sslasso_file)) {
      vals <- read_one_row(sslasso_file)
      
      row_here <- data.frame(
        P_Q_N = P,
        Method = "SSLASSO",
        TPR = safe_round3(vals[1]),
        FPR = safe_round3(vals[2]),
        MCC = safe_round3(vals[3]),
        AUC = safe_round3(vals[4]),
        AUC20 = safe_round3(vals[5]),
        Sparsity = safe_round3(vals[6])
      )
      
      all_rows <- rbind(all_rows, row_here)
    } else {
      warning("SSLASSO file not found: ", sslasso_file)
    }
    
    ## --------------------------------------------------------
    ## 3. remMap-bic
    ## Expected file:
    ## Output other methods/Output_rhoZero_remMapBic_P_20_Q_20_N_20_Nfactor_1_Seed_1_TPR_FPR_MCC_AUC_AUC20_sparsity.csv
    ## --------------------------------------------------------
    
    remmap_file <- file.path(
      OutputOtherFolder,
      paste0(
        "Output_", rhoTag,
        "_remMapBic",
        "_P_", P,
        "_Q_", Q,
        "_N_", N,
        "_Nfactor_", SampleMultFactor,
        "_Seed_", RandSeed,
        "_TPR_FPR_MCC_AUC_AUC20_sparsity.csv"
      )
    )
    
    if (file.exists(remmap_file)) {
      vals <- read_one_row(remmap_file)
      
      row_here <- data.frame(
        P_Q_N = P,
        Method = "remMap-bic",
        TPR = safe_round3(vals[1]),
        FPR = safe_round3(vals[2]),
        MCC = safe_round3(vals[3]),
        AUC = safe_round3(vals[4]),
        AUC20 = safe_round3(vals[5]),
        Sparsity = safe_round3(vals[6])
      )
      
      all_rows <- rbind(all_rows, row_here)
    } else {
      warning("remMap-bic file not found: ", remmap_file)
    }
  }
  
  return(all_rows)
}

## ============================================================
## Build rho = 0 table
## ============================================================

summary_rhoZero <- build_summary_table("rhoZero")

write.csv(
  summary_rhoZero,
  "Output_Scalability_rhoZero_summary_P_20_50_100_200_ALL_methods.csv",
  row.names = FALSE
)

## ============================================================
## Build rho = 0.5 table
## ============================================================

summary_rhoNonZero <- build_summary_table("rhoNonZero")

write.csv(
  summary_rhoNonZero,
  "Output_Scalability_rhoNonZero_summary_P_20_50_100_200_ALL_methods.csv",
  row.names = FALSE
)

cat("\nFinished creating summary CSV files.\n")
cat("Saved:\n")
cat("Output_Scalability_rhoZero_summary_P_20_50_100_200_ALL_methods.csv\n")
cat("Output_Scalability_rhoNonZero_summary_P_20_50_100_200_ALL_methods.csv\n")