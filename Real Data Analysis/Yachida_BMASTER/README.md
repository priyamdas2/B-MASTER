# Raw data extraction and exploratory plots

The analysis can be reproduced by running the following scripts in order:

1. Run `1_Full_raw_data_extraction_only.R` to extract required data of CRC subjects for B-MASTER analysis from raw data.

2. Run `2_Data_analysis.R`, `3_Extra_plots.R`, `4_Extra_supp_plots.R`, `5_patient_summary_plots.R`, `6_Scores_plot.R` in this respective order to perform exploratory analysis and generate exploratory plots.

3. Run `7_data_setup_only_healthy.R` to extract datatset (containing same set of X and Y as considered for the CRC subjects) for healthy subjects.

**P.S.** The raw data file `genera.tsv` and the processed data file `genera_full.csv` are not included in this repository due to their large file sizes.
