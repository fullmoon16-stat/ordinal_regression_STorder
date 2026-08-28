# Code files

This folder contains the R and Python code used to reproduce the figures, simulation studies, and empirical analyses reported in the manuscript and its Supplementary Material.

All files should be run from the root directory of this repository. The required data files are provided in the `data/` folder.

## Common R functions

- `common_functions.R`

  This file contains functions shared by the R analyses, including the baseline survival functions, transformation functions, likelihood evaluation, model fitting, data-driven initial values, prediction, and model-selection utilities for the proposed ordinal regression models.

## Main-manuscript figures

- `Figure1_converted_Sf.R`

  This file reproduces Figure 1, which illustrates the converted survival functions.

- `Figure2_baseline_Sfs.R`

  This file reproduces Figure 2, which compares the baseline survival functions.

- `Figure3_transformation_effects.R`

  This file reproduces Figure 3, which illustrates the effects of the transformation functions.

- `Figure4_simulation_study_GCLM.Rmd`

  This file reproduces Figure 4 from the simulation study for the GCLMs. It also reports the AIC and BIC selection frequencies and the median deviance and MSE for each scenario and competing model.

- `Figure5_different_dispersion.R`

  This file reproduces Figure 5, which illustrates conditional distributions with different dispersion.

- `Figure6_crossing_transformation.R`

  This file reproduces Figure 6, which illustrates crossing conditional distributions under a category-specific transformation.

## Supplementary simulation study

- `FigureS1_simulation_study_GPPOM.Rmd`

  This file reproduces Supplementary Figure S1 from the simulation study for the GPPOMs. It also reports the AIC and BIC selection frequencies and the median deviance and MSE for each simulation scenario and competing model.

## Empirical analyses

- `Table5_APS.Rmd`

  This file reproduces Table 5 from the APS data analysis, including the full-data estimation results and the out-of-sample LogS and RPS values.

- `Table6_retinopathy.Rmd`

  This file reproduces Table 6 from the retinopathy data analysis, including the full-data estimation results and the out-of-sample LogS and RPS values.

## White-wine analyses

- `TableS2_white_wine.py`

  This is the main executable file for Supplementary Table S2, based on the original seven-category white-wine response. It fits the CMNN GCLM, cumulative softmax GCLM, CI model, SI-CS model, and POM under repeated stratified cross-validation. Run it from the repository root using:

  ```bash
  python code/TableS2_white_wine.py
  ```

  The final table is printed directly to the console.

- `TableS3_white_wine.py`

  This is the main executable file for Supplementary Table S3, based on the five-category recoding of the white-wine response. It fits the same five models under repeated stratified cross-validation. Run it from the repository root using:

  ```bash
  python code/TableS3_white_wine.py
  ```

  The final table is printed directly to the console.

### Supporting files for the white-wine analyses

- `white_wine_common.py`

  This file contains the shared Python implementation for Supplementary Tables S2 and S3, including data loading, fold construction, within-fold standardization, neural-network definitions, architecture and baseline selection, early stopping, refitting, prediction, LogS and unnormalized RPS calculation, progress recovery, and final-table construction. It is imported automatically by the two main Python files and should not be run separately.

- `white_wine_POM.R`

  This R helper fits and evaluates the POM using the same cross-validation splits as the neural models. It is called automatically by the two main Python files through `Rscript` and should not normally be run separately.
