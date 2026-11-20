# ============================================================
# REGENERATE: Dataset screening table for external validation
# Produces Supplementary Table Sx
# ============================================================

suppressPackageStartupMessages(library(hyperSpec))

# --- Portable project paths --------------------------------
# Open the accompanying .Rproj file before running this script.
PROJECT_ROOT <- here::here()
PREPROC_DIR  <- here::here("data", "processed")
CSV_DIR      <- here::here("data", "raw")
INPUT_CSV    <- "Spectral_library_06-Oct-2026.csv"
OUT_DIR      <- here::here("Results", "External_Validation")
TRAIN_GROUP <- "13"

# --- Pre-specified screening thresholds --------------------
MIN_N_SURV     <- 30       # minimum survivors after clean_matrix
CEIL_REFL      <- 1.5      # physical reflectance ceiling
MAX_FRAC_HI    <- 0.05     # max fraction of rows above ceiling
MIN_CP_OVERLAP <- 0.40     # min fraction of Cp inside training range

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load workspace ----------------------------------------
snap <- fileSnapshot(PREPROC_DIR)
latest_dir <- rownames(snap$info[which.max(snap$info$mtime),])
rdata_files <- list.files(file.path(PREPROC_DIR, latest_dir),
                          pattern = "\\.RData$", full.names = TRUE)
if (!length(rdata_files)) stop("No RData in ", latest_dir)
load(rdata_files[which.max(file.mtime(rdata_files))])
cat("Workspace:", latest_dir, "\n")

# --- Load trait table --------------------------------------
read_csv_eu <- function(path) {
  df <- read.csv(path, sep = ",", header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, quote = "\"", comment.char = "")
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      x2 <- gsub(",", ".", x, fixed = TRUE)
      n  <- suppressWarnings(as.numeric(x2))
      if (mean(!is.na(n)) > 0.5) return(n)
    }
    x
  })
  df
}
raw <- read_csv_eu(file.path(CSV_DIR, INPUT_CSV))
hdr <- gsub("^X", "", colnames(raw))
ids_all <- as.character(raw[[1]])
cp_all  <- suppressWarnings(as.numeric(raw[[match("cp", hdr)]]))
grp_all <- as.character(raw[[match("dataset", hdr)]])

# --- Training reference ------------------------------------
tr_rows   <- grp_all == TRAIN_GROUP & is.finite(cp_all)
train_ids <- ids_all[tr_rows]
y_train   <- cp_all[tr_rows]
tr_mean <- mean(y_train); tr_sd <- sd(y_train)
tr_min  <- min(y_train);  tr_max <- max(y_train)
tr_iqr  <- quantile(y_train, c(.25, .75))

cat(sprintf("Training (ds%s): n=%d  mean=%.4e  SD=%.4e  range=[%.4e, %.4e]\n",
            TRAIN_GROUP, length(y_train), tr_mean, tr_sd, tr_min, tr_max))

# --- Helpers -----------------------------------------------
clean_matrix <- function(X) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  X <- X[rowSums(!is.finite(X)) == 0, , drop = FALSE]
  if (nrow(X) > 1) {
    rsd <- apply(X, 1, sd); X <- X[!( !is.finite(rsd) | rsd == 0 ), , drop = FALSE]
  }
  if (ncol(X) > 1) {
    csd <- apply(X, 2, sd); X <- X[, !( !is.finite(csd) | csd == 0 ), drop = FALSE]
  }
  X
}
get_variant <- function(v) {
  hs <- get(v)
  X  <- if (inherits(hs, "hyperSpec")) as.matrix(hs@data$spc) else as.matrix(hs$spc)
  rn <- if (inherits(hs, "hyperSpec")) rownames(hs@data) else rownames(hs)
  rownames(X) <- rn
  X
}

# Post-QC matrix (used to determine sample survival)
X_qc <- clean_matrix(get_variant("spec_07_refl_deriv1"))
clean_ids <- rownames(X_qc)

# Raw matrix for reflectance sanity (pre-QC)
X_raw <- as.matrix(if (inherits(get("raw_spc"), "hyperSpec"))
  get("raw_spc")@data$spc else get("raw_spc")$spc)
rownames(X_raw) <- if (inherits(get("raw_spc"), "hyperSpec"))
  rownames(get("raw_spc")@data) else rownames(get("raw_spc"))
wl_raw <- if (inherits(get("raw_spc"), "hyperSpec"))
  get("raw_spc")@wavelength else get("raw_spc")$wavelength
col500 <- which.min(abs(wl_raw - 500))

tr_500 <- mean(X_raw[train_ids[train_ids %in% rownames(X_raw)], col500], na.rm = TRUE)
cat(sprintf("Training reflectance @ 500 nm: %.4f\n\n", tr_500))

# --- Scan all non-training datasets ------------------------
all_ds <- setdiff(sort(unique(grp_all)), TRAIN_GROUP)
res <- data.frame()

for (ds in all_ds) {
  rows <- which(grp_all == ds & is.finite(cp_all))
  if (!length(rows)) next
  ids_ds <- ids_all[rows]; y_ds <- cp_all[rows]
  
  X_ds <- as.matrix(raw[rows, 11:2111]); storage.mode(X_ds) <- "double"
  
  in_range <- mean(y_ds >= tr_min & y_ds <= tr_max)
  in_iqr   <- mean(y_ds >= tr_iqr[1] & y_ds <= tr_iqr[2])
  mean_500 <- mean(X_ds[, col500], na.rm = TRUE)
  raw_max  <- max(X_ds, na.rm = TRUE)
  raw_min  <- min(X_ds, na.rm = TRUE)
  n_nonfin <- sum(!is.finite(X_ds))
  frac_hi  <- mean(apply(X_ds, 1, function(r) any(r > CEIL_REFL, na.rm = TRUE)))
  n_surv   <- sum(ids_ds %in% clean_ids)
  
  pass_integrity <- (n_surv >= MIN_N_SURV) &
    (raw_max <= CEIL_REFL) &
    (frac_hi <= MAX_FRAC_HI) &
    (n_nonfin < 0.05 * length(X_ds))
  pass_cp <- in_range >= MIN_CP_OVERLAP
  
  reason <- if (pass_integrity && pass_cp) {
    "PASSED - selected for external validation"
  } else if (!pass_integrity && !pass_cp) {
    "FAILED - integrity and Cp overlap"
  } else if (!pass_integrity) {
    if (n_surv < MIN_N_SURV) {
      "FAILED - integrity (<30 samples survive QC)"
    } else if (raw_max > CEIL_REFL) {
      "FAILED - integrity (reflectance > 1.5; DN-scale or corrupted)"
    } else {
      "FAILED - integrity (>5% rows corrupted)"
    }
  } else {
    sprintf("FAILED - Cp overlap (%.0f%% inside range; need >=40%%)",
            100 * in_range)
  }
  
  res <- rbind(res, data.frame(
    dataset              = ds,
    n_raw                = length(ids_ds),
    n_surv               = n_surv,
    cp_mean              = mean(y_ds),
    cp_sd                = sd(y_ds),
    cp_min               = min(y_ds),
    cp_max               = max(y_ds),
    frac_cp_train_range  = round(in_range, 3),
    frac_cp_train_iqr    = round(in_iqr,   3),
    raw_min              = raw_min,
    raw_max              = raw_max,
    n_nonfinite          = n_nonfin,
    frac_rows_above_1.5  = round(frac_hi, 3),
    mean_500             = mean_500,
    brightness_ratio_500 = round(mean_500 / tr_500, 3),
    cp_mean_ratio        = round(mean(y_ds) / tr_mean, 3),
    screen_integrity     = if (pass_integrity) "PASS" else "FAIL",
    screen_cp_overlap    = if (pass_cp)        "PASS" else "FAIL",
    screen_overall       = if (pass_integrity && pass_cp) "PASS" else "FAIL",
    reason               = reason,
    stringsAsFactors     = FALSE
  ))
}

# Sort: passed first, then by Cp overlap descending
res <- res[order(res$screen_overall != "PASS", -res$frac_cp_train_range), ]

# --- Print summary -----------------------------------------
cat("\n===== SCREENING SUMMARY =====\n")
cat("Total non-training datasets :", nrow(res), "\n")
cat("Passed integrity screen     :", sum(res$screen_integrity == "PASS"), "\n")
cat("Passed Cp overlap screen    :", sum(res$screen_cp_overlap == "PASS"), "\n")
cat("Passed both (usable)        :", sum(res$screen_overall == "PASS"), "\n\n")

cat("Datasets passing both criteria:\n")
print(res[res$screen_overall == "PASS",
          c("dataset","n_raw","n_surv","frac_cp_train_range",
            "brightness_ratio_500","cp_mean_ratio")],
      row.names = FALSE)

cat("\nDatasets that FAILED integrity:\n")
print(res[res$screen_integrity == "FAIL",
          c("dataset","n_raw","n_surv","raw_max","n_nonfinite",
            "frac_rows_above_1.5")],
      row.names = FALSE)

cat("\nDatasets that passed integrity but FAILED Cp overlap:\n")
print(res[res$screen_integrity == "PASS" & res$screen_cp_overlap == "FAIL",
          c("dataset","n_surv","frac_cp_train_range","cp_mean_ratio",
            "brightness_ratio_500")],
      row.names = FALSE)

# --- Check the selected campaign ---------------------------
cat("\n===== SELECTED CAMPAIGN (must be Dataset 26) =====\n")
sel <- res[res$screen_overall == "PASS", ]
if (nrow(sel) == 0) {
  cat("No campaign passed both screens.\n")
} else if (nrow(sel) > 1) {
  cat("WARNING: more than one campaign passed both screens.\n")
  print(sel)
} else {
  cat(sprintf("Dataset %s: n_surv = %d, %.1f%% Cp inside range, "
              , sel$dataset, sel$n_surv, 100 * sel$frac_cp_train_range),
      sprintf("brightness ratio %.2f\n", sel$brightness_ratio_500))
  cat(">>> Verify this matches the n (=80) <<<\n")
  if (sel$n_surv != 80) {
    cat(">>> WARNING: n_surv differs from paper's n=80. <<<\n")
    cat(">>> Either the workspace/CSV differs from the run that produced",
        "the paper, or the screening extraction differs. <<<\n")
  }
}

# --- Save --------------------------------------------------
out_csv <- file.path(OUT_DIR, "Supplementary_Table_Sx_dataset_screening.csv")
write.csv(res, out_csv, row.names = FALSE)
cat("\nFull table written to:\n  ", out_csv, "\n", sep = "")


# ============================================================
# EXTERNAL VALIDATION — Dataset 26 as held-out test set
# ============================================================

suppressPackageStartupMessages({
  library(pls); library(ranger); library(e1071); library(resemble)
})

# --- Portable paths ----------------------------------------
PREPROC_DIR     <- here::here("data", "processed")
CSV_DIR         <- here::here("data", "raw")
BEST_MODELS_DIR <- here::here("Models")
INPUT_CSV       <- "Spectral_library_06-Oct-2026.csv"
EXT_GROUP       <- "26"
OUT_DIR         <- file.path(here::here("Results", "External_Validation"),
                             paste0("ds", EXT_GROUP))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TRAIN_GROUP <- "13"
SEED        <- 42
BOOT_B      <- 2000L

EXPECTED_CONFIGS <- list(
  PLSR = list(bin = 50L, variant = "spec_07_refl_deriv1"),
  RF   = list(bin = 50L, variant = "spec_07_refl_deriv1"),
  MBL  = list(bin = 50L, variant = "spec_07_refl_deriv1"),
  SVR  = list(bin =  5L, variant = "spec_02_refl_smooth")
)

# ============================================================
# 1. Workspace, CSV, helpers
# ============================================================
snap <- fileSnapshot(PREPROC_DIR)
latest_dir <- rownames(snap$info[which.max(snap$info$mtime),])
rdata_files <- list.files(file.path(PREPROC_DIR, latest_dir),
                          pattern = "\\.RData$", full.names = TRUE)
if (!length(rdata_files)) stop("No RData in ", latest_dir)
load(rdata_files[which.max(file.mtime(rdata_files))])
cat("Workspace:", latest_dir, "\n")

read_csv_eu <- function(path) {
  df <- read.csv(path, sep = ",", header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, quote = "\"", comment.char = "")
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      x2 <- gsub(",", ".", x, fixed = TRUE)
      n  <- suppressWarnings(as.numeric(x2))
      if (mean(!is.na(n)) > 0.5) return(n)
    }
    x
  })
  df
}
raw <- read_csv_eu(file.path(CSV_DIR, INPUT_CSV))
hdr <- gsub("^X", "", colnames(raw))
ids_all <- as.character(raw[[1]])
cp_all  <- suppressWarnings(as.numeric(raw[[match("cp", hdr)]]))
grp_all <- as.character(raw[[match("dataset", hdr)]])

train_ids <- ids_all[grp_all == TRAIN_GROUP & is.finite(cp_all)]
ext_ids   <- ids_all[grp_all == EXT_GROUP   & is.finite(cp_all)]
if (!length(ext_ids)) stop("No finite-cp observations for dataset ", EXT_GROUP)

cat("Dataset 13 IDs:", length(train_ids),
    "| Dataset", EXT_GROUP, "IDs:", length(ext_ids), "\n")
stopifnot(length(intersect(train_ids, ext_ids)) == 0)

clean_matrix <- function(X) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  X <- X[rowSums(!is.finite(X)) == 0, , drop = FALSE]
  if (nrow(X) > 1) {
    rsd <- apply(X, 1, sd); X <- X[!( !is.finite(rsd) | rsd == 0 ), , drop = FALSE]
  }
  if (ncol(X) > 1) {
    csd <- apply(X, 2, sd); X <- X[, !( !is.finite(csd) | csd == 0 ), drop = FALSE]
  }
  X
}
bin_matrix <- function(X, b) {
  if (b <= 1) return(X)
  p <- ncol(X); grp <- ((seq_len(p) - 1) %/% b) + 1
  M <- matrix(0, p, max(grp)); M[cbind(seq_len(p), grp)] <- 1
  sweep(X %*% M, 2, colSums(M), "/")
}
get_variant <- function(v) {
  hs <- get(v)
  X  <- if (inherits(hs, "hyperSpec")) as.matrix(hs@data$spc) else as.matrix(hs$spc)
  rownames(X) <- if (inherits(hs, "hyperSpec")) rownames(hs@data) else rownames(hs)
  X
}
reproduce_and_apply <- function(variant, bin, train_ids, ext_ids) {
  X_all <- clean_matrix(get_variant(variant))
  tr_ok <- intersect(train_ids, rownames(X_all))
  ex_ok <- intersect(ext_ids,   rownames(X_all))
  if (!length(tr_ok) || !length(ex_ok))
    stop("Empty subset for ", variant, " / bin ", bin)
  Xb_tr <- bin_matrix(X_all[tr_ok, , drop = FALSE], bin)
  Xb_ex <- bin_matrix(X_all[ex_ok, , drop = FALSE], bin)
  if (bin > 1) {
    csd  <- apply(Xb_tr, 2, sd)
    keep <- is.finite(csd) & csd != 0
    Xb_tr <- Xb_tr[, keep, drop = FALSE]
    Xb_ex <- Xb_ex[, keep, drop = FALSE]
  }
  list(train = Xb_tr, ext = Xb_ex,
       y_train = cp_all[match(tr_ok, ids_all)],
       y_ext   = cp_all[match(ex_ok, ids_all)],
       tr_ids  = tr_ok, ex_ids = ex_ok)
}

list_model_files <- function(dir, prefix) {
  f <- list.files(dir, full.names = TRUE)
  f[grepl(paste0("^", prefix, "_cp"), basename(f)) &
      !grepl("_meta",    basename(f)) &
      !grepl("_summary", basename(f))]
}
list_meta_files <- function(dir, prefix) {
  f <- list.files(dir, full.names = TRUE)
  f[grepl(paste0("^", prefix, "_cp.*_meta"), basename(f))]
}
pick_latest <- function(f, label) {
  if (!length(f)) stop("No file for ", label)
  f[which.max(file.mtime(f))]
}

# ============================================================
# 2. Resolve model + meta files, verify configs
# ============================================================
paths <- list()
for (alg in names(EXPECTED_CONFIGS)) {
  d <- file.path(BEST_MODELS_DIR, alg)
  paths[[alg]] <- list(
    model = pick_latest(list_model_files(d, alg), paste(alg, "model")),
    meta  = pick_latest(list_meta_files (d, alg), paste(alg, "meta"))
  )
}
metas <- lapply(paths, function(p) readRDS(p$meta))

cat("\nConfiguration check:\n")
for (alg in names(EXPECTED_CONFIGS)) {
  exp <- EXPECTED_CONFIGS[[alg]]; m <- metas[[alg]]
  ok  <- as.integer(m$bin_size) == exp$bin && m$preprocessing == exp$variant
  cat(sprintf("  %-5s exp bin=%-3d preproc=%-22s | got bin=%-3d preproc=%-22s  %s\n",
              alg, exp$bin, exp$variant, as.integer(m$bin_size),
              m$preprocessing, if (ok) "OK" else "*** MISMATCH ***"))
  if (!ok) stop("Config mismatch for ", alg)
}
configs <- lapply(metas, function(m)
  list(variant = m$preprocessing, bin = as.integer(m$bin_size)))

# ============================================================
# 3. Preprocess + alignment gate
# ============================================================
cache <- list()
for (alg in names(configs)) {
  c_  <- configs[[alg]]
  key <- paste(c_$variant, c_$bin, sep = "|")
  if (is.null(cache[[key]])) {
    cat("\nPreprocessing", c_$variant, "@", c_$bin, "nm\n")
    cache[[key]] <- reproduce_and_apply(c_$variant, c_$bin, train_ids, ext_ids)
    cat(sprintf("  train n=%d  ext n=%d  cols=%d\n",
                nrow(cache[[key]]$train),
                nrow(cache[[key]]$ext),
                ncol(cache[[key]]$ext)))
  }
}

# ---- Alignment gate ----
cat("\nAlignment gate (max |reproduced mean - saved mean|):\n")
fit_plsr <- readRDS(paths$PLSR$model)
fit_svr  <- readRDS(paths$SVR$model)
a_plsr   <- cache[[paste(configs$PLSR$variant, configs$PLSR$bin, sep = "|")]]
a_svr    <- cache[[paste(configs$SVR$variant, configs$SVR$bin, sep = "|")]]
align_plsr <- max(abs(colMeans(a_plsr$train) - as.numeric(fit_plsr$Xmeans)))
align_svr  <- max(abs(colMeans(a_svr$train) -
                        as.numeric(fit_svr$x.scale$`scaled:center`)))
cat(sprintf("  PLSR : %.3e  (tolerance 1e-6)\n", align_plsr))
cat(sprintf("  SVR  : %.3e  (tolerance 1e-6)\n", align_svr))
if (align_plsr > 1e-6 || align_svr > 1e-6)
  warning("Reproduced training matrix does NOT match saved model means. ",
          "Predictions may be misaligned.")

# ============================================================
# 4. Pre-flight: target distribution + brightness check
# ============================================================
y_train <- a_plsr$y_train; y_ext <- a_plsr$y_ext
tr_min <- min(y_train); tr_max <- max(y_train)

cat("\nTarget distribution:\n")
cat(sprintf("  Training (ds13): mean=%.4e  SD=%.4e  range=[%.4e, %.4e]\n",
            mean(y_train), sd(y_train), tr_min, tr_max))
cat(sprintf("  External (ds%s): mean=%.4e  SD=%.4e  range=[%.4e, %.4e]\n",
            EXT_GROUP, mean(y_ext), sd(y_ext), min(y_ext), max(y_ext)))
cat(sprintf("  Cp mean ratio (ext/train): %.2f\n", mean(y_ext) / mean(y_train)))
frac_outside <- mean(y_ext < tr_min | y_ext > tr_max)
cat(sprintf("  Fraction external Cp outside training range: %.1f%%\n",
            100 * frac_outside))

# Brightness check at 500 nm (before any preprocessing)
X_raw_full <- as.matrix(if (inherits(get("raw_spc"), "hyperSpec"))
  get("raw_spc")@data$spc else get("raw_spc")$spc)
rn <- rownames(X_raw_full)
wl_raw <- if (inherits(get("raw_spc"), "hyperSpec"))
  get("raw_spc")@wavelength else get("raw_spc")$wavelength
col500 <- which.min(abs(wl_raw - 500))
br_ratio <- mean(X_raw_full[ext_ids[ext_ids %in% rn], col500], na.rm = TRUE) /
  mean(X_raw_full[train_ids[train_ids %in% rn], col500], na.rm = TRUE)
cat(sprintf("  Brightness ratio (ext/train @ 500 nm): %.2f  (target 0.7-1.4)\n",
            br_ratio))

if (br_ratio < 0.6 || br_ratio > 1.6)
  warning("Brightness ratio is far from 1.0. Expectations: external RPD may be ",
          "< 1. This was the failure mode of Dataset 1 and 27.")
if (mean(y_ext) / mean(y_train) > 3 ||
    frac_outside > 0.5)
  warning("Target distribution is far from training. External predictions may ",
          "regress toward the training mean and inflate RMSEP.")

# ============================================================
# 5. Predict
# ============================================================
preds <- list(); ids_by <- list(); y_by <- list()

a <- cache[[paste(configs$PLSR$variant, configs$PLSR$bin, sep = "|")]]
preds$PLSR <- as.numeric(predict(fit_plsr,
                                 newdata = data.frame(NIR = I(a$ext)),
                                 ncomp   = fit_plsr$ncomp))
ids_by$PLSR <- a$ex_ids; y_by$PLSR <- a$y_ext

a <- cache[[paste(configs$RF$variant, configs$RF$bin, sep = "|")]]
fit_rf <- readRDS(paths$RF$model)
stopifnot(ncol(a$ext) == fit_rf$num.independent.variables)
df_ex <- as.data.frame(a$ext); colnames(df_ex) <- paste0("V", seq_len(ncol(df_ex)))
preds$RF <- predict(fit_rf, data = df_ex)$predictions
ids_by$RF <- a$ex_ids; y_by$RF <- a$y_ext

a <- cache[[paste(configs$MBL$variant, configs$MBL$bin, sep = "|")]]
Xr <- a$train; Xu <- a$ext
colnames(Xr) <- paste0("V", seq_len(ncol(Xr)))
colnames(Xu) <- paste0("V", seq_len(ncol(Xu)))
old <- options(cores = 1L)
fit_mbl_ext <- resemble::mbl(
  Xr = Xr, Yr = a$y_train,
  Xu = Xu, Yu = NULL,
  neighbors   = resemble::neighbors_k(k = metas$MBL$best_k),
  diss_method = resemble::diss_pca(),
  diss_usage  = "weights",
  fit_method  = resemble::fit_wapls(min_ncomp = metas$MBL$ncomp_min,
                                    max_ncomp = metas$MBL$ncomp_max),
  control     = resemble::mbl_control(validation_type = "NNv"),
  seed        = SEED, verbose = FALSE)
options(old)
preds$MBL <- as.numeric(fit_mbl_ext$results[[1]]$pred)
ids_by$MBL <- a$ex_ids; y_by$MBL <- a$y_ext

a <- cache[[paste(configs$SVR$variant, configs$SVR$bin, sep = "|")]]
stopifnot(ncol(a$ext) == length(fit_svr$x.scale$`scaled:center`))
df_ex <- as.data.frame(a$ext); colnames(df_ex) <- paste0("V", seq_len(ncol(df_ex)))
preds$SVR <- as.numeric(predict(fit_svr, newdata = df_ex))
ids_by$SVR <- a$ex_ids; y_by$SVR <- a$y_ext

# ============================================================
# 6. Reconcile, metrics, save
# ============================================================
common_ids <- Reduce(intersect, ids_by)
cat("\nSample reconciliation (common across models):", length(common_ids), "\n")
for (alg in names(preds)) {
  idx <- match(common_ids, ids_by[[alg]])
  preds[[alg]] <- preds[[alg]][idx]
  y_by[[alg]]  <- y_by[[alg]][idx]
}
y_ext <- y_by$PLSR

metrics <- function(y, yhat) {
  e <- yhat - y
  data.frame(n = length(y),
             RMSEP = sqrt(mean(e^2)),
             MAE   = mean(abs(e)),
             bias  = mean(e),
             R2    = cor(y, yhat)^2,
             R2_SS = 1 - sum(e^2) / sum((y - mean(y))^2),
             RPD   = sd(y) / sqrt(mean(e^2)),
             SD_y  = sd(y))
}
boot_rpd <- function(y, yhat, B = BOOT_B) {
  idx <- replicate(B, sample.int(length(y), replace = TRUE))
  v   <- apply(idx, 2, function(i) sd(y[i]) / sqrt(mean((yhat[i] - y[i])^2)))
  quantile(v, c(0.025, 0.975), na.rm = TRUE)
}

ext_metrics <- do.call(rbind, lapply(names(preds), function(m)
  data.frame(Model = m, metrics(y_ext, preds[[m]]))))
rownames(ext_metrics) <- ext_metrics$Model
ci <- t(sapply(preds, function(p) boot_rpd(y_ext, p)))
colnames(ci) <- c("RPD_low", "RPD_high")
ext_metrics  <- cbind(ext_metrics, ci)

internal <- data.frame(
  Model     = names(metas),
  RPD_int   = sapply(metas, function(m) as.numeric(m$cv_summary$RPD)),
  RMSEP_int = sapply(metas, function(m) as.numeric(m$cv_summary$RMSEP)),
  row.names = names(metas)
)
ext_metrics <- merge(ext_metrics, internal, by = "Model")
ext_metrics$RPD_retention   <- ext_metrics$RPD   / ext_metrics$RPD_int
ext_metrics$RMSEP_inflation <- ext_metrics$RMSEP / ext_metrics$RMSEP_int
ext_metrics <- ext_metrics[match(names(preds), ext_metrics$Model), ]

cat("\n===== External validation on Dataset", EXT_GROUP,
    "(n =", length(y_ext), ") =====\n")
print(ext_metrics[, c("Model","n","RPD","RPD_low","RPD_high","RPD_retention",
                      "RMSEP","RMSEP_inflation","R2","R2_SS","bias","SD_y")],
      digits = 4)

write.csv(ext_metrics,
          file.path(OUT_DIR, paste0("external_validation_ds", EXT_GROUP, "_metrics.csv")),
          row.names = FALSE)
saveRDS(list(metrics      = ext_metrics,
             internal     = internal,
             predictions  = preds,
             ids          = common_ids,
             y_ext        = y_ext,
             train_y      = y_train,
             configs      = configs,
             timestamp    = format(Sys.time(), "%Y%m%d_%H%M")),
        file.path(OUT_DIR, paste0("external_validation_ds", EXT_GROUP, "_full.rds")))
cat("\nWritten to:", OUT_DIR, "\n")