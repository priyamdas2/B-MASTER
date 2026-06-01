# Real Data Analysis Workflow

The analysis can be reproduced by running the following scripts in order:

1. Run `0_Extraction_subset_for_analysis.R` to extract data for subset-based analyses. This step can be run before or after the main B-MASTER analyses.
2. Run `BMASTER_Real_Data.m` for the real data analysis of the main cohort, consisting of diseased subjects.
3. Run `BMASTER_Real_Data_HEALTHY.m` for the real data analysis of the healthy cohort.
4. Run `2_FIS_plots.R` to generate FIS score plots.
5. Run `3_FIS_summary_tables.R` to generate FIS summary tables.
6. Run `4_Post_analysis_subset_1_NEW.R` for post-B-MASTER full-fit subset analysis of the first subset, corresponding to the most abundant metabolites.
7. Run `5_Post_analysis_subset_2_NEW.R` for post-B-MASTER full-fit subset analysis of the second subset, corresponding to differentially abundant metabolites in cancer.
8. Run `6_CCA_subset_plots.R` to generate CCA decay plots for both subset 1 and subset 2.

## Codes for raw data extraction and exploratory analyses

See README instructions in `Yachida_BMASTER` folder.
