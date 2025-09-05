# B-MASTER

**Paper Title:**  B-MASTER: Scalable Bayesian Multivariate Regression for Master Predictor Discovery in Colorectal Cancer Microbiome-Metabolite Profiles 

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
