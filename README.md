# B-MASTER

**Title:**  B-MASTER: Scalable Bayesian Multivariate Regression for Master Predictor Discovery in Colorectal Cancer Microbiome-Metabolite Profiles 

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
### 🕒 Approximate Computation Time

Running B-MASTER for 1,000 iterations (sufficient for posterior sample–based inference) on a dataset with dimensions P = 287, Q = 249, and N = 220 -- matching the case study data -- takes approximately **100–120 minutes** to complete. This benchmark was obtained on a **Windows 10 Enterprise** system with **32 GB RAM** and a **12th Gen Intel(R) Core(TM) i7-12700** processor (12 cores, 20 logical threads, 2.10 GHz base frequency).

---

---

## 📊 Key Visual Summaries

All major figures and visual outputs from the manuscript are available under the `images/` directory. Representative visuals include:

### Microbiome-Metabolite Interplay in Colorectal Cancer
![Microbiome-Metabolite interplay](images/Microbiome_metabolite_diagram.jpg)

This schematic illustrates how gut microbial activity influences metabolite production, which in turn may affect colorectal cancer progression.

### Data Processing and Exploratory Analysis
![Data processing](images/BMASTER_case_study_collage.jpg)

This describes the overview of the colorectal cancer microbiome-metabolome dataset characteristics and processing steps.

### B-MASTER Concept Diagram
![B-MASTER concept](images/BMASTER_concept.jpg)

This diagram visually depicts `master predictors'.

### B-MASTER Scalability
![B-MASTER scalability](images/Comp_times.png)

This diagram highlights the scalability and performance of B-MASTER in higher dimensions compared to remMap.

### Top 50 FIS-based Master Predictors
![Top 50 master predictors](images/FIS_top50_wide.png)

Functional Influence Scores (FIS) for the top 50 genera identified by B-MASTER. Bar height reflects each genus’s FIS; the numbers above the bars indicate, for that genus, the number of metabolites it influences.

### Top 50 FIS-based Master Predictors
![Heatmap](images/Top50genera_heat_FIS.jpeg)

This demonstrates Heatmap illustrating the direction of the impact for each of the identified top 50 master predictors on each metabolite.
---

### Master set of Genera Regulating Most Abundant Metabolites (Subset 1)
![Heatmap](images/Subset_1.png)

The plot demonstrates the direction and statistical significance of the relationships between the most abundant metabolites (Subset 1) on the corresponding top 15 key genera identified in the analysis. Progressively larger circles are associated with smaller Bayesian p-values.
---

### Master set of Genera Regulating Differentially Abundant Metabolites in CRC (Subset 2)
![Heatmap](images/Subset_2.png)

The plot demonstrates the direction and statistical significance of the relationships between the metabolites identified as differentially abundant in the cancer versus control comparison (Subset 2) and the top 15 key genera identified in the analysis.
---

## 🔐 Data Access

Dataset is shared at `Reproducible codes/Real Data Analysis/Data/Yachida_BMASTER/`.

---

## 💬 Contact

For questions, please contact:  
**Priyam Das**  
[dasp4@vcu.edu](mailto:dasp4@vcu.edu)
