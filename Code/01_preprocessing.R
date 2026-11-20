#************************************************************************#
# A New Approach to Predict Nutrient Content in Plants                   #
#        Using VNIR-SWIR Spectroscopy and Machine-Learning               #
#                                                                        #
#             PART 1 - PREPROCESSING SPECTRA (v3.5)                      #
#                                                                        #
#              Jhohan Elguera - MATE University - HUN-REN ATK            #
#                 With the Contributions of Dr. Johan PERRET             #
#************************************************************************#

library(hyperSpec)
library(signal)
library(pls)
library(RColorBrewer)
library(rrcov)
library(MASS)

rm(list = ls())

# ============================================================
# 0. PROJECT PATHS AND INPUT LIBRARY
# ============================================================
# Open the accompanying .Rproj file before running this script.
# The 'here' package resolves paths relative to the project root.
PROJECT_ROOT <- here::here()
CSV_DIR      <- here::here("data", "raw")
PREPROC_DIR  <- here::here("data", "processed")

Input_CSV <- "Spectral_library_06-Oct-2026.csv"
INPUT_FILE <- file.path(CSV_DIR, Input_CSV)

if (!file.exists(INPUT_FILE)) {
  stop("Input CSV not found: ", INPUT_FILE,
       "\\nPlace the spectral library in the project's data/raw/ folder.")
}

# ---- Reader for comma-separated / comma-decimal CSVs -------
read_csv_eu <- function(path) {
  df <- read.csv(path, sep = ",", header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE,
                 quote = "\"", comment.char = "")
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      x2  <- gsub(",", ".", x, fixed = TRUE)
      num <- suppressWarnings(as.numeric(x2))
      if (mean(!is.na(num)) > 0.5) return(num)
    }
    x
  })
  df
}

raw_data_library <- read_csv_eu(INPUT_FILE)

# -------- Spectral matrix (columns 11..2111) ----------------
stopifnot(ncol(raw_data_library) >= 2111)
raw_spc <- data.matrix(raw_data_library[, 11:2111])
rownames(raw_spc) <- as.character(raw_data_library[[1]])
colnames(raw_spc) <- gsub("^X", "", colnames(raw_spc))
wl <- suppressWarnings(as.numeric(colnames(raw_spc)))

if (any(is.na(wl))) stop("Column names 11..2111 are not wavelengths.")
cat("Wavelength range :", range(wl), "nm\n")
cat("Wavelength step  :", unique(diff(wl)), "nm\n")
cat("Samples (rows)   :", nrow(raw_spc), "\n")
cat("Bands   (cols)   :", ncol(raw_spc), "\n")

raw_spc <- new("hyperSpec", spc = raw_spc, wavelength = wl)
raw_spc@label$.wavelength <- expression("Wavelength (nm)")
raw_spc@label$spc <- "Reflectance"

# -------- Output directory ----------------------------------
dir_today <- format(Sys.time(), "Preprocessed_spectra_%d_%m_%Y")
dir_today <- file.path(PREPROC_DIR, dir_today)
dir.create(dir_today, showWarnings = FALSE, recursive = TRUE)
setwd(dir_today)

Input_file <- gsub("[.]", "_", Input_CSV)
input_file_name <- paste0("000_Preprocessing_from_", Input_file, ".txt")
cat("Preprocessing source CSV: ", Input_CSV, "\n", file = input_file_name)

# ============================================================
# 0b. NON-FINITE ROW REMOVAL (data integrity)
# ------------------------------------------------------------
# Removes any spectrum containing NA / NaN / Inf. No statistical
# outlier filtering here.
# ============================================================

X_check <- as.matrix(raw_spc$spc)
storage.mode(X_check) <- "double"

finite_row <- apply(X_check, 1, function(x) all(is.finite(x)))
n_drop     <- sum(!finite_row)

cat("\n--- Non-finite value check ---\n")
cat("Total NA   :", sum(is.na(X_check)),  "\n")
cat("Total NaN  :", sum(is.nan(X_check)), "\n")
cat("Total Inf  :", sum(is.infinite(X_check)), "\n")
cat("Rows with any non-finite value:", n_drop, "\n")
cat("Fraction of library removed   :",
    round(100 * n_drop / nrow(X_check), 2), "%\n")

if (n_drop > 0) {
  bad_rows_idx <- which(!finite_row)
  bad_pattern  <- t(sapply(bad_rows_idx, function(i) {
    bc <- which(!is.finite(X_check[i, ]))
    c(n_bad  = length(bc),
      wl_min = if (length(bc)) min(wl[bc]) else NA_real_,
      wl_max = if (length(bc)) max(wl[bc]) else NA_real_)
  }))
  bad_pattern <- data.frame(sample_name = rownames(X_check)[bad_rows_idx],
                            bad_pattern, row.names = NULL)
  write.csv(bad_pattern, "000_dropped_nonfinite_samples.csv", row.names = FALSE)
}

raw_spc <- raw_spc[finite_row, ]
X_check <- NULL
cat("Library after non-finite removal:", nrow(raw_spc), "samples\n\n")

# ============================================================
# 0c. SPECTRAL VALUE SANITY CLEANING (data integrity)
# ------------------------------------------------------------
# Steps:
#   (a) DROP any spectrum that contains a band > 1.5 (physical
#       reflectance ceiling — vegetation/soil reflectance never
#       exceeds ~1.2; detector saturation and CSV corruption
#       produce spikes at 1.5+)
#   (b) null out remaining cells outside [-0.5, 1.5]
#   (c) drop samples where >5% of bands were still affected
#   (d) linearly interpolate the remaining isolated gaps
#   (e) clip remaining negatives / near-zero cells to 1e-4
#       (reflectance is non-negative and log10(1/X) requires
#        X > 0)
# Data-integrity cleaning only — NOT statistical outlier QC.
# ============================================================

# Physical ceiling for reflectance. Vegetation and soil never
# exceed ~1.2; 1.5 is a very permissive upper bound that still
# rejects detector saturation spikes.
ceil_refl <- 1.0

X <- as.matrix(raw_spc$spc)
storage.mode(X) <- "double"

# ---- (a) drop spectra containing any band > ceil_refl ---------
row_max <- apply(X, 1, max, na.rm = TRUE)
bad_row <- row_max > ceil_refl

cat("--- Spectral value sanity check ---\n")
cat(sprintf("Reflectance ceiling used           : %.2f\n", ceil_refl))
cat(sprintf("Samples with any band > %.2f       : %d\n",
            ceil_refl, sum(bad_row)))
cat(sprintf("Fraction of library dropped        : %.2f%%\n",
            100 * sum(bad_row) / nrow(X)))

if (any(bad_row)) {
  drop_report <- data.frame(
    sample_name = rownames(X)[bad_row],
    row_max     = row_max[bad_row]
  )
  write.csv(drop_report, "000_dropped_above_ceiling.csv", row.names = FALSE)
}

X <- X[!bad_row, , drop = FALSE]
cat("Library after ceiling filter:", nrow(X), "samples\n\n")

# ---- (b) cell-level cleaning -------------------------------
lower_bd <- -0.5
upper_bd <-  ceil_refl
bad_cell <- !is.finite(X) | X < lower_bd | X > upper_bd
n_bad_cells <- sum(bad_cell)
cat(sprintf("Cells outside [%.2f, %.2f]: %d  (%.4f%%)\n",
            lower_bd, upper_bd, n_bad_cells,
            100 * n_bad_cells / length(X)))

affected_cells <- which(bad_cell, arr.ind = TRUE)
if (nrow(affected_cells) > 0) {
  cell_report <- data.frame(
    sample  = rownames(X)[affected_cells[, 1]],
    wl_nm   = wl[affected_cells[, 2]],
    value   = X[affected_cells]
  )
  write.csv(cell_report, "000_corrupted_cells.csv", row.names = FALSE)
}
X[bad_cell] <- NA

# ---- (c) row-level filtering (bad-fraction) ----------------
frac_bad <- rowMeans(is.na(X))
keep_row <- frac_bad <= 0.05
cat("Samples with <=5% corrupted bands: kept  :", sum(keep_row),  "\n")
cat("Samples with  >5% corrupted bands: dropped:", sum(!keep_row), "\n")

X <- X[keep_row, , drop = FALSE]
sample_names_kept <- rownames(X)

# ---- (d) interpolate remaining isolated gaps ---------------
n_interp <- 0
for (i in seq_len(nrow(X))) {
  nas <- which(is.na(X[i, ]))
  if (length(nas) == 0) next
  good <- which(!is.na(X[i, ]))
  if (length(good) >= 2) {
    X[i, nas] <- approx(good, X[i, good], xout = nas, rule = 2)$y
    n_interp <- n_interp + length(nas)
  } else {
    X[i, nas] <- median(X[i, ], na.rm = TRUE)
    n_interp <- n_interp + length(nas)
  }
}
cat("Interpolated cells (isolated gaps):", n_interp, "\n")

# ---- (e) clip remaining negatives / near-zero -------------
floor_val <- 1e-4

n_neg <- sum(X < floor_val)
cat(sprintf("Cells below %.4f (negative or near-zero): %d  (%.4f%%)\n",
            floor_val, n_neg, 100 * n_neg / length(X)))

neg_idx <- which(X < floor_val, arr.ind = TRUE)
if (nrow(neg_idx) > 0) {
  wl_of_neg <- wl[neg_idx[, 2]]
  cat("Wavelength range of negative cells:",
      range(wl_of_neg), "nm\n")
  cat("Negatives in the 400-420 nm window :",
      sum(wl_of_neg <= 420), "/", length(wl_of_neg), "\n")
}
X[X < floor_val] <- floor_val

cat("Final spectra value range:",
    round(range(X, na.rm = TRUE), 4), "\n\n")

raw_spc <- new("hyperSpec", spc = X, wavelength = wl)
raw_spc@label$.wavelength <- expression("Wavelength (nm)")
raw_spc@label$spc <- "Reflectance"

# Re-attach the trait table, keeping only the samples that survived
trait_cols <- raw_data_library[, 1:10]
trait_cols <- trait_cols[match(sample_names_kept, as.character(trait_cols[[1]])), ]
rownames(trait_cols) <- NULL

# ============================================================
# 1. OUTLIER QC — REMOVED
# ------------------------------------------------------------
# Statistical outlier detection intentionally disabled so the
# four ML models are compared on the full integrity-cleaned
# library.
# ============================================================

spec_01_refl_raw <- raw_spc
cat("Samples kept for preprocessing (no outlier QC):",
    nrow(spec_01_refl_raw), "\n\n")

set.seed(123)
random_rows <- sample(1:nrow(spec_01_refl_raw),
                      min(500, nrow(spec_01_refl_raw)), replace = FALSE)

title_graph <- c(Input_CSV,
                 paste0(nrow(spec_01_refl_raw),
                        " spectra (integrity-cleaned, no outlier QC)"))

plot(spec_01_refl_raw[random_rows, , ], lines.args = list(pch = 20),
     col = rainbow(300), zeroline = NA,
     plot.args = list(xlim = c(400, 2500)),
     title.args = list(xlab = "Wavelength (nm)"),
     wl.reverse = FALSE, spc.nmax = 500, nxticks = 5)
title(main = title_graph, font.main = 1, cex.main = 0.85)
dev.print(pdf, "spec_01b_refl_raw_clean.pdf")

rm(X, bad_cell, affected_cells, frac_bad, keep_row, cell_report,
   neg_idx, wl_of_neg, bad_row, row_max, drop_report)

# ============================================================
# 2. PREPROCESSING FUNCTIONS (CV-safe, training-agnostic)
# ============================================================

# ---- 2a. Savitzky–Golay ------------------------------------
sg_smooth <- function(X, n = 25, p = 3) {
  Xs <- t(apply(X, 1, sgolayfilt, m = 0, p = p, n = n))
  h  <- (n - 1) / 2
  list(X = Xs[, (h + 1):(ncol(Xs) - h), drop = FALSE], wl_trim = h)
}
sg_deriv1 <- function(X, n = 9, p = 3) {
  Xs <- t(apply(X, 1, sgolayfilt, m = 1, p = p, n = n))
  h  <- (n - 1) / 2
  list(X = Xs[, (h + 1):(ncol(Xs) - h), drop = FALSE], wl_trim = h)
}
sg_deriv2 <- function(X, n = 15, p = 3) {
  Xs <- t(apply(X, 1, sgolayfilt, m = 2, p = p, n = n))
  h  <- (n - 1) / 2
  list(X = Xs[, (h + 1):(ncol(Xs) - h), drop = FALSE], wl_trim = h)
}

# ---- 2b. SNV (per spectrum) --------------------------------
snv <- function(X) {
  t(apply(X, 1, function(v) {
    if (any(!is.finite(v))) return(v)
    s <- sd(v)
    if (!is.finite(s) || s == 0) s <- 1
    (v - mean(v)) / s
  }))
}

# ---- 2c. MSC (reference computed on TRAIN only) ------------
msc_fit   <- function(Xtrain) colMeans(Xtrain, na.rm = TRUE)
msc_apply <- function(X, ref) {
  t(apply(X, 1, function(x) {
    if (any(!is.finite(x)) || any(!is.finite(ref))) return(x)
    fit <- lm(x ~ ref)
    a <- coef(fit)[1]; b <- coef(fit)[2]
    if (!is.finite(a)) a <- 0
    if (!is.finite(b) || b == 0) b <- 1
    (x - a) / b
  }))
}

# ---- 2d. Column centering / z-scoring ----------------------
col_center_fit   <- function(Xtrain) colMeans(Xtrain, na.rm = TRUE)
col_center_apply <- function(X, mu) sweep(X, 2, mu, "-")

col_zscore_fit <- function(Xtrain) {
  list(mu = colMeans(Xtrain, na.rm = TRUE),
       sd = apply(Xtrain, 2, sd, na.rm = TRUE))
}
col_zscore_apply <- function(X, pars) {
  s <- pars$sd; s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(X, 2, pars$mu, "-"), 2, s, "/")
}

# ---- 2e. Absorbance ----------------------------------------
to_absorbance <- function(X, floor = 1e-6, ceil = 1.05) {
  X <- pmin(pmax(X, floor), ceil)
  log10(1 / X)
}

# ============================================================
# 3. APPLY PREPROCESSING
# ============================================================

X_refl  <- as.matrix(spec_01_refl_raw$spc)
storage.mode(X_refl) <- "double"
wl_full <- spec_01_refl_raw@wavelength

stopifnot(all(is.finite(X_refl)))
cat("X_refl finite check: OK\n")

make_hyperspec <- function(X, wl, label) {
  hs <- new("hyperSpec", spc = X, wavelength = wl)
  hs@label$.wavelength <- expression("Wavelength (nm)")
  hs@label$spc <- label
  hs
}

plot_hs <- function(hs, tag, ylim = NULL) {
  if (is.null(ylim)) ylim <- range(hs$spc, na.rm = TRUE)
  plot(hs[random_rows, , ], lines.args = list(pch = 20),
       col = rainbow(300), zeroline = NA, wl.reverse = FALSE,
       spc.nmax = 500, nxticks = 5,
       plot.args = list(xlim = c(400, 2500), ylim = ylim),
       title.args = list(xlab = "Wavelength (nm)", cex.main = 0.85))
  title(main = title_graph, font.main = 1, cex.main = 0.85)
  dev.print(pdf, tag)
}

trim_wl <- function(fit, wl) wl[(fit$wl_trim + 1):(length(wl) - fit$wl_trim)]

# ---- 3a. Reflectance variants ------------------------------
r_sg <- sg_smooth(X_refl, n = 25, p = 3)
r_d1 <- sg_deriv1(X_refl, n = 9,  p = 3)
r_d2 <- sg_deriv2(X_refl, n = 15, p = 3)

spec_02_refl_smooth <- make_hyperspec(r_sg$X, trim_wl(r_sg, wl_full),
                                      "Reflectance (SG smoothed)")
spec_07_refl_deriv1 <- make_hyperspec(r_d1$X, trim_wl(r_d1, wl_full),
                                      "Reflectance (1st derivative)")
spec_08_refl_deriv2 <- make_hyperspec(r_d2$X, trim_wl(r_d2, wl_full),
                                      "Reflectance (2nd derivative)")
spec_05_refl_snv    <- make_hyperspec(snv(X_refl), wl_full, "Reflectance (SNV)")

plot_hs(spec_02_refl_smooth, "spec_02_refl_smooth.pdf")
plot_hs(spec_07_refl_deriv1, "spec_07_refl_deriv1.pdf")
plot_hs(spec_08_refl_deriv2, "spec_08_refl_deriv2.pdf")
plot_hs(spec_05_refl_snv,    "spec_05_refl_snv.pdf")

# ---- 3b. Absorbance variants -------------------------------
X_abs <- to_absorbance(X_refl)

a_sg <- sg_smooth(X_abs, n = 25, p = 3)
a_d1 <- sg_deriv1(X_abs, n = 9,  p = 3)
a_d2 <- sg_deriv2(X_abs, n = 21, p = 3)

spec_09_abs_raw    <- make_hyperspec(X_abs, wl_full, "Absorbance")
spec_10_abs_smooth <- make_hyperspec(a_sg$X, trim_wl(a_sg, wl_full),
                                     "Absorbance (SG smoothed)")
spec_15_abs_deriv1 <- make_hyperspec(a_d1$X, trim_wl(a_d1, wl_full),
                                     "Absorbance (1st derivative)")
spec_16_abs_deriv2 <- make_hyperspec(a_d2$X, trim_wl(a_d2, wl_full),
                                     "Absorbance (2nd derivative)")
spec_13_abs_snv    <- make_hyperspec(snv(X_abs), wl_full, "Absorbance (SNV)")

plot_hs(spec_09_abs_raw,    "spec_09_abs_raw.pdf")
plot_hs(spec_10_abs_smooth, "spec_10_abs_smooth.pdf")
plot_hs(spec_15_abs_deriv1, "spec_15_abs_deriv1.pdf")
plot_hs(spec_16_abs_deriv2, "spec_16_abs_deriv2.pdf")
plot_hs(spec_13_abs_snv,    "spec_13_abs_snv.pdf")

# ---- 3c. Non-finite derivative check -----------------------
for (nm in c("spec_07_refl_deriv1", "spec_08_refl_deriv2",
             "spec_15_abs_deriv1",  "spec_16_abs_deriv2")) {
  hs <- get(nm)
  n_bad <- sum(!is.finite(hs$spc))
  cat(sprintf("%-22s non-finite values before cleanup: %d\n", nm, n_bad))
  if (n_bad > 0) {
    hs$spc[!is.finite(hs$spc)] <- 0
    assign(nm, hs)
  }
}

# ---- 3d. Wavelength range table ----------------------------
wl_table <- data.frame(
  variant = c("Reflectance raw (post-QC)", "Reflectance SG25",
              "Reflectance d1 (SG9)", "Reflectance d2 (SG15)",
              "Reflectance SNV",
              "Absorbance raw", "Absorbance SG25",
              "Absorbance d1 (SG9)", "Absorbance d2 (SG21)",
              "Absorbance SNV"),
  wl_min = c(min(wl_full), min(spec_02_refl_smooth@wavelength),
             min(spec_07_refl_deriv1@wavelength), min(spec_08_refl_deriv2@wavelength),
             min(wl_full),
             min(wl_full), min(spec_10_abs_smooth@wavelength),
             min(spec_15_abs_deriv1@wavelength), min(spec_16_abs_deriv2@wavelength),
             min(wl_full)),
  wl_max = c(max(wl_full), max(spec_02_refl_smooth@wavelength),
             max(spec_07_refl_deriv1@wavelength), max(spec_08_refl_deriv2@wavelength),
             max(wl_full),
             max(wl_full), max(spec_10_abs_smooth@wavelength),
             max(spec_15_abs_deriv1@wavelength), max(spec_16_abs_deriv2@wavelength),
             max(wl_full))
)
write.csv(wl_table, "000_wavelength_ranges.csv", row.names = FALSE)
print(wl_table)

# ============================================================
# 4. SAVE
# ============================================================
save.image(file = file.path(dir_today,
                            format(Sys.time(), "Preprocessed_spectra_%d_%m_%Y.RData")))
cat("\nAll outputs written to:", dir_today, "\n")

