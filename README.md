# ordinal_regression_STorder
Code and data for 
**"A General Class of Ordinal Regression Models Developed Based on the Notion of Usual Stochastic Order"**

The focus of this repository is the real ordinal data sets analysis in Section 5.

---

## Running order

No strict execution order is required. The three R Markdown files for the empirical analyses correspond to Sections 5.1, 5.2, and 5.3, respectively, and can be run independently from the repository root. These files source `code/common_functions.R`, which contains the main functions for the proposed ordinal regression models, including baseline survival functions, transformation functions, likelihood evaluation, model fitting, and model comparison.

The scripts for Figures 1-3 are also independent and can be run separately from the repository root.

---

## Repository contents

- `data` : data files and descriptions of data sources.

- `R` : R Markdown files and common R functions used to reproduce the empirical analyses.

---

## Reproducing the manuscript results

### Empirical analyses

The empirical analyses in Section 5 can be reproduced by running the R Markdown files in the `code/` folder:

- `5_1_mental_health.Rmd` : Tables 5-7 and Figure 4.
- `5_2_retinopathy.Rmd` : Tables 8-11.
- `5_3_pm25.Rmd` : Tables 12-14.

All files should be run from the root directory of this repository. No separate output directory is required.

### Illustrative figures

- `Figure1_converted_Sf.R` : Figure 1
- `Figure2_baseline_Sfs.R` : Figure 2
- `Figure3_transformation_effects.R` : Figure 3
