# Portable R scripts: folder setup and execution

These six scripts have had their machine-specific `F:/...` and `D:/...` paths replaced with project-relative paths. The statistical/model-fitting code and parameter settings have not intentionally been changed.

## 1. Open the project

1. Download/extract this repository folder.
2. Open `hyperspectral-Cp-prediction.Rproj` in RStudio.
3. Install the `here` package if needed:

```r
install.packages("here")
```

The scripts use `here::here()` to resolve folders relative to the project root. Run them from the opened RStudio project.

## 2. Required folder layout

```text
hyperspectral-Cp-prediction/
├── hyperspectral-Cp-prediction.Rproj
├── Code/
│   ├── 01_preprocessing.R
│   ├── 02_PLSR_nested_CV.R
│   ├── 03_RF_nested_CV.R
│   ├── 04_MBL_nested_CV.R
│   ├── 05_SVR_nested_CV.R
│   └── 06_Validation.R
├── data/
│   ├── raw/
│   │   └── Spectral_library_06-Oct-2026.csv
│   └── processed/
├── Models/
└── Results/
```

Place the exact input CSV (`Spectral_library_06-Oct-2026.csv`) in `data/raw/`. The preprocessing script creates a dated workspace folder inside `data/processed/`. The model scripts save model objects under `Models/` and summary outputs under `Results/Model_Summaries/`. The validation script writes screening and external-validation outputs under `Results/External_Validation/`.

## 3. Run order

Run the scripts in this order:

1. `Code/01_preprocessing.R`
2. `Code/02_PLSR_nested_CV.R`
3. `Code/03_RF_nested_CV.R`
4. `Code/04_MBL_nested_CV.R`
5. `Code/05_SVR_nested_CV.R`
6. `Code/06_Validation.R` (after the four model scripts have generated their model and metadata files)

Run each script from RStudio using **Source**. Do not run the model scripts before preprocessing has generated a `.RData` workspace.

## 4. Important notes

- The raw spectral CSV is not bundled here; it must be supplied separately and placed in `data/raw/`.
- The preprocessing script still changes the R working directory to its dated output folder so its existing relative diagnostic-output writes are saved alongside the preprocessing workspace. This does not affect the project-relative input/output locations used by the other scripts.
- The scripts continue to select the most recently modified preprocessing subfolder. If you have multiple preprocessing runs, verify that the newest folder is the intended one before running the models.
- Keep the original input CSV filename unchanged, or update `INPUT_CSV`/`Input_CSV` consistently in the scripts.
- No modelling parameters, cross-validation settings, preprocessing choices, target definitions, or validation thresholds were deliberately altered in this path-portability update.
