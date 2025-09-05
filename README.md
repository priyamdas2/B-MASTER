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
| **Figure 1**  | Sankey diagrams of treatment sequences     | `Real Data Analysis/Real Data/exploratory_analysis.R` |
| **Figures 5,6a** | Estimated transition probabilities (SMART-MC) | `Real Data Analysis/SMART_MC_Var_effect_plot.R` |
| **Figure 6b** | Odds ratios for across-treatment transitions | `Real Data Analysis/SMART_MC_ODDS_ratio_calculation.m` → `SMART_MC_Odds_ratio_plot.R` |
| **Table 1, S1**, **Figure S1**   | MSCOR benchmark results                    | `MSCOR Benchmark/MSCOR_Benchmark_comparison.m` → `MSCOR_post_evaluation.m` |
| **Tables S2–S5** | Simulation results                       | `Simulation Study/` (see details in `SMART_MC_Reproducibility_and_DEMO_instructions.pdf`)        |
| **Table S6**, **Figure S2** | SMART-MC estimated coefficients; simulated treatment trajectory | `Real Data Analysis/SMART_MC_Real_data.m` |

---
