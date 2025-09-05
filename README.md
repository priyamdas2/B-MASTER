# B-MASTER

**Paper Title:**  B-MASTER: Scalable Bayesian Multivariate Regression for Master Predictor Discovery in Colorectal Cancer Microbiome-Metabolite Profiles 

📄 [Read the full paper on arXiv](https://arxiv.org/abs/2412.05998)

This repository provides the necessary code and documentation for reproducing the results presented in the article above. It implements the **B-MASTER** model for Master Predictor Discovery in Colorectal Cancer Microbiome-Metabolite Profiles.

---

## 📄 Paper Summary

The manuscript introduces:

- **B-MASTER** (Bayesian Multivariate regression Analysis for Selecting Targeted Essential Regressors) is a Bayesian multivariate regression framework for identifying *master predictors*.  
- The method is complemented by an efficient Gibbs sampling scheme for posterior inference, with demonstrated scalability in case studies involving up to **4 million parameters**.  
- We establish key theoretical properties:  
  1. Posterior contraction rates under standard design and sparsity conditions  
  2. Robustness under model misspecification  
  3. Row-selection consistency, ensuring reliable recovery of master predictors  
- A comprehensive simulation study highlights the superiority of **B-MASTER** over existing frequentist methods for master predictor selection, as well as over other uni- and multivariate regression techniques.  
- Applied to colorectal cancer data, **B-MASTER** identifies key microbial genera that regulate the overall metabolome profile. In addition, master predictors are investigated for two biologically relevant subsets of metabolites:  
  - The most abundant metabolites  
  - Metabolites differentially abundant in CRC  
---

## 📁 File Overview

| Folder | Description |
|--------|-------------|
| `Reproducible codes/Real Data Analysis/` | Real data analysis scripts (Figures 2, 7-12; Table 2) |
| `Reproducible codes/Simulation Study 1/` | Simulation experiments 1 (Figures 5; Table 1) |
| `Reproducible codes/Simulation Study 2/` | Simulation experiments 2 (Figures 4,6) |
| `Reproducible codes/BMASTER_reproducibility_instructions.pdf` | Full reproduction instructions |
| `B-MASTER Demo` | A demonstration of B-MASTER |
| `images/` | Final versions of key figures used in the manuscript |

---

## 🧮 Figures and Tables in the Paper

To reproduce the tables and figures presented in the paper, please refer to `BMASTER_reproducibility_instructions.pdf`. A brief overview is provided below.


| Output        | Description                                | Script Path                            |
|---------------|--------------------------------------------|----------------------------------------|
| **Figure 1**  | Influence of gut microbiome on metabolomics in CRC    | Concept diagram |
| **Figure 2**  | Data extraction and exploration    | `Reproducible codes/Real Data Analysis/Data/Yachida_BMASTER/`|
| **Figure 3**  | B-MASTER illustration    | Concept diagram |
| **Figure 4, 6**  | B-MASTER scalability  | `Reproducible codes/Simulation Study 2/` |
| **Figure 5**  | Main simulation study  | `Reproducible codes/Simulation Study 1/` |
| **Figure 7**  | Score plots  | `Reproducible codes/Real Data Analysis/Data/Yachida_BMASTER/` |
| **Figure 8-12**  | B-MASTER analysis of microbiome-metabolite interplay  | `Reproducible codes/Real Data Analysis/` |
| **Tables 1** | Simulation results                       | `Reproducible codes/Simulation Study 1/` |
| **Tables 2** | Out-of-sample prediction performance                       | `Reproducible codes/Real Data Analysis/Validation on real data/` |


---

## 📦 Requirements

- **R version**: 4.5.1  
- **MATLAB**: R2021a or later  
- **Required R packages**:
  - `corrplot`, `readxl`, `BiocManager`, `Bioconductor`, `tidyr`, `stringr`, `scales`, `ggpattern`, `cowplot`, `grid`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `patchwork`, `RColorBrewer`, `ggalluvial`, `reshape2`, `purrr`

---
