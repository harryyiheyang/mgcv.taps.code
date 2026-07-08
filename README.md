# Companion Analysis Code for the TAPS Manuscript

This repository contains the companion analysis scripts for the manuscript

> Testing Parametric Structure in Genetic Age-Effect Curves: A GAM-Based
> Framework with Application to the UK Biobank

It should be read together with the separate `mgcv.taps` R package repository.
The package implements the TAPS smooth constructors and hypothesis tests. This
repository uses that package to run the simulation studies, UK Biobank analyses,
manuscript figure generation, and the interactive result viewer.

## Repository Roles

| Component | Role | Publicly runnable without restricted data? |
| --- | --- | --- |
| `mgcv.taps` package | Implements TAPS smooths, Wald tests, score tests, package examples, and function documentation | Yes |
| `simulation codes/` | Simulation scripts for type-I error and power experiments | Yes, after installing dependencies |
| `varying coefficient simulation codes/` | Simulation scripts for age-varying coefficient settings | Yes, after installing dependencies |
| `plot codes/` | Downstream scripts that read saved simulation or analysis results and generate manuscript plots | Conditional; requires the corresponding `.rds` or `.csv` outputs |
| `data analysis codes/` | UK Biobank age-varying and retirement-threshold PRS analyses | No; requires restricted UK Biobank individual-level data and locally derived PRS files |
| `online_exhibition/` | Shiny app for browsing precomputed manuscript figures included in this repository | Yes |

## Relationship to the `mgcv.taps` Package

The reproducibility materials are split across two repositories because they
serve different purposes.

The `mgcv.taps` package is the reusable statistical software, with exported
functions, examples, dependency metadata, and roxygen2-generated help files.
This repository is the paper-analysis repository: the scripts here call
`mgcv.taps` to run simulations, UK Biobank analyses, and figure generation.
Accordingly, this README focuses on workflow-level reproducibility and on which
parts can be run without restricted UK Biobank data.

## Installation and Software Dependencies

Install the TAPS package first:

```r
remotes::install_github("harryyiheyang/mgcv.taps")
```

The scripts use several R packages depending on the analysis component:

- Core TAPS and GAM analyses: `mgcv.taps`, `mgcv`, `data.table`, `dplyr`,
  `glue`, `ggplot2`, `scales`, `grid`, `mgcViz`, `gratia`, `glmnet`,
  `rcompanion`.
- Simulation scripts: `mgcv.taps`, `mgcv`, `MASS`, `qgam`, `survival`.
- Plotting scripts: `ggplot2`, `dplyr`, `reshape2`, `RColorBrewer`, `readr`,
  `forcats`, `forestplot`, `FDRestimation`.
- Interactive viewer: `shiny`, `shinydashboard`, `DT`.

The package-level dependencies of `mgcv.taps` are listed in the package
`DESCRIPTION` file.

## What Can Be Run Directly

### Simulation scripts

The simulation scripts do not require UK Biobank data. They generate synthetic
data and then apply TAPS through the `mgcv.taps` package.

Directly runnable simulation scripts include the family/structure scripts under
`simulation codes/` and `varying coefficient simulation codes/`, such as:

- Gaussian, binary, Poisson, ordinal, and survival simulations.
- Linearity, breakpoint, changepoint, and interaction settings.
- Median/quantile GAM variants where `qgam` is used.

Run these scripts from the repository root. The public simulation scripts write
their outputs to repository-local folders created automatically by the scripts:

- `outputs/simulation/` for scripts in `simulation codes/`.
- `outputs/varying_coefficient_simulation/` for scripts in
  `varying coefficient simulation codes/`.

The simulations are stochastic. Users who need bitwise run-to-run
reproducibility should set a seed before running an individual script.

### Interactive result viewer

The Shiny app in `online_exhibition/app.R` can be run directly because it uses
the precomputed PNG figures included under `online_exhibition/Age/` and
`online_exhibition/Retirement/`:

```r
shiny::runApp("online_exhibition")
```

## What Requires Precomputed Outputs

The scripts in `plot codes/` are downstream plotting scripts. They assume that
the corresponding simulation or analysis outputs have already been generated
and saved as `.rds` or `.csv` files. These scripts are useful for documenting
how manuscript figures were assembled, but they are not the first step in the
workflow.

The file `simulation codes/forestplot_plot.R` is also a downstream plotting
script. It depends on precomputed retirement/RDD result files under
`outputs/restricted/RDD/`, which are not included because they are derived from
restricted UK Biobank analyses.

Similarly, `plot codes/manhattenplot.R` depends on precomputed real-data result
files under `outputs/restricted/real_data/`. The trait-category metadata file
used for plotting, `trait_category_mapping.csv`, is included in the repository
because it contains only public trait labels and plotting categories, not
individual-level UK Biobank data. The real-data `.rds` result files are not
generated by the public simulation scripts.

Run order for the public simulation workflow is:

1. Install `mgcv.taps`.
2. Run the relevant scripts in `simulation codes/` or
   `varying coefficient simulation codes/`.
3. Run the corresponding scripts in `plot codes/` to reproduce the figures
   from saved simulation outputs. Public simulation plotting scripts read from
   `outputs/simulation/`.

## What Cannot Be Run Publicly

The UK Biobank real-data scripts cannot be run from this repository alone. They
require restricted individual-level UK Biobank data and locally derived files
that cannot be redistributed under the UK Biobank data-sharing rules.

This applies to:

- `data analysis codes/Age_Interaction_Gaussian.R`
- `data analysis codes/Age_Interaction_Binary.R`
- `data analysis codes/Retirement_Interaction_Gaussian.R`
- `data analysis codes/Retirement_Interaction_Binary.R`

These scripts are provided to document the analysis model specifications,
phenotype/PRS processing steps, covariate adjustment, exclusion criteria, and
output generation. They are intended for authorized users who have access to
the required UK Biobank individual-level data and can prepare the corresponding
local phenotype, disease, PRS, covariate, ancestry, kinship, and withdrawal
files. The local paths in these scripts reflect the authors' analysis
environment and should be remapped by authorized users.

The manuscript reports the analytic sample sizes after exclusions, and the Data
and Code Availability statement describes how qualified researchers may access
UK Biobank individual-level data through the UK Biobank application process.

No UK Biobank individual-level data are included in this repository.

## Notes for Reproducibility Review

The analysis code and the package should be evaluated as complementary but
distinct materials:

1. The `mgcv.taps` package provides the reusable statistical software,
   exported functions, examples, dependency metadata, and roxygen2-generated
   documentation.
2. This repository provides the paper-specific scripts, simulation workflows,
   figure-generation code, and Shiny result viewer.
3. Public users can run the simulation scripts and the Shiny viewer without
   UK Biobank access.
4. The UK Biobank scripts document the full model specifications and
   preprocessing logic, but cannot be made fully self-contained because the
   required individual-level data are access-controlled.
