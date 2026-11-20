# Hyperspectral Cp prediction — reproducibility repository

This repository contains the original fitted model objects and result files supplied with the project, together with portable versions of the six R scripts.

## Important: existing results are preserved

**You do not need to rerun the modelling scripts to use the results already included here.** The existing `.rds` model objects, metadata, cross-validation prediction files, summary tables, and Dataset 26 validation tables have been carried over unchanged from the supplied ZIP archive.

The path-portability update changes file locations in the R scripts only. It does not refit models, regenerate tables, or modify the existing result files.

## Repository contents

- `Code/`: six R scripts, updated to use project-relative paths via `here::here()`.
- `Models/`: existing fitted PLSR, RF, MBL, and SVR model objects and metadata.
- `Results/`: existing nested-CV fold predictions, model comparison tables, dataset-screening table, and Dataset 26 validation metrics.
- `hyperspectral-Cp-prediction.Rproj`: RStudio project file.
- `README_PORTABLE_PATHS.md`: instructions for running the scripts in a portable project layout.

## Open the repository

1. Extract the ZIP archive.
2. Open `hyperspectral-Cp-prediction.Rproj` in RStudio.
3. The existing results can be inspected directly from the `Models/` and `Results/` folders. No model rerun is required.

## If you decide to rerun the workflow in the future

The scripts require the original input file `Spectral_library_06-Oct-2026.csv` in `data/raw/` and the R packages used by the scripts. Install `here` if it is not already available:

```r
install.packages("here")
```

Run the scripts in this order:

1. `Code/01_preprocessing.R`
2. `Code/02_PLSR_nested_CV.R`
3. `Code/03_RF_nested_CV.R`
4. `Code/04_MBL_nested_CV.R`
5. `Code/05_SVR_nested_CV.R`
6. `Code/06_Validation.R`

The raw CSV is not included in this archive. The existing fitted models and results are included, so the above execution sequence is only needed if you intentionally want to reproduce the analysis from raw data.

## Notes on the portable scripts

- Machine-specific `F:/...` and `D:/...` paths were replaced with project-relative paths.
- The modelling procedures, model settings, cross-validation design, preprocessing choices, target definitions, and validation thresholds were not intentionally changed.
- The preprocessing script still changes the working directory to its dated output folder so its relative diagnostic outputs are saved alongside the preprocessing workspace.
- Model scripts save new model objects under `Models/` and new summary outputs under `Results/Model_Summaries/`; the supplied original result tables remain in their original algorithm-specific folders under `Results/`.
- The scripts select the most recently modified preprocessing subfolder. If multiple preprocessing runs exist, verify that the newest folder is the intended one before running model scripts.
