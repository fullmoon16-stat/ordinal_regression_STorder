# ordinal_regression_STorder
Code and data for 
**"A General Class of Ordinal Regression Models Developed Based on the Notion of Usual Stochastic Order"**

This repository contains the code and data used to reproduce the empirical analyses in Section 5 and the illustrative figures in the manuscript.


## Repository contents

- `data/`: data files and data documentation. The file `data/data_sources.md` describes the data sources, variable definitions, coding schemes, units, and preprocessing steps.

- `code/`: R Markdown files, R scripts, and common R functions used to reproduce the empirical analyses and illustrative figures. The file `code/code_sources.md` describes the role of each code file.


## Reproducing the manuscript results

All files should be run from the root directory of this repository. The reproduced tables and figures are displayed in the corresponding R Markdown outputs or R plotting windows. No separate output directory is required.

No strict execution order is required. The three R Markdown files for the empirical analyses correspond to Sections 5.1, 5.2, and 5.3, respectively, and can be run independently. These files source `code/common_functions.R`. The scripts for Figures 1--3 are also independent and can be run separately.

### Empirical analyses

The empirical analyses in Section 5 can be reproduced by running the R Markdown files in the `code/` folder:

- `5_1_mental_health.Rmd`: Tables 5--7 and Figure 4.
- `5_2_retinopathy.Rmd`: Tables 8--11.
- `5_3_pm25.Rmd`: Tables 12--14.

### Illustrative figures

The illustrative figures can be reproduced by running the following R scripts in the `code/` folder:

- `Figure1_converted_Sf.R`: Figure 1.
- `Figure2_baseline_Sfs.R`: Figure 2.
- `Figure3_transformation_effects.R`: Figure 3.


