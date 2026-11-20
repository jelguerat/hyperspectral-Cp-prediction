#********************************************************************#
# A New Approach to Predict Nutrient Content in Plants               #
#        Using VNIR-SWIR Spectroscopy and Machine-Learning           #
#              Jhohan Elguera - MATE University                      #
#               With contributions of Dr. Johan PERRET               #
#                        PART 4 MEMORY-BASED LEARNING                #
#           FULLY-NESTED CV — ALL CANDIDATES ON EVERY OUTER FOLD     #
#                                    2026                            #
#********************************************************************#
#
#************************************************************************#
# Nested-CV MBL   — Nitrogen prediction (target = cp)                    #
#                                                                        #
# SUBSET VERSION: dataset 13 only                                        #
#   - scans preprocessing variants x bin sizes                           #
#   - saves the winning ranger model to 4-Best_Models/MBL/               #
#                                                                        #
#************************************************************************#

suppressPackageStartupMessages({
  library(resemble)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(tictoc)
})

rm(list = ls())

# ============================================================
# 0. Robust accessor for hyperSpec 0.100+
# ============================================================
get_spc_matrix <- function(obj) {
  m <- tryCatch(obj@data$spc, error = function(e) NULL)
  if (!is.null(m) && is.matrix(m)) return(m)
  m <- tryCatch(slot(obj, "spc"), error = function(e) NULL)
  if (!is.null(m) && is.matrix(m)) return(m)
  m <- tryCatch(obj$spc, error = function(e) NULL)
  if (!is.null(m) && is.matrix(m)) return(m)
  m <- attr(obj, "spc")
  if (!is.null(m) && is.matrix(m)) return(m)
  stop("Could not extract spectral matrix from class ",
       paste(class(obj), collapse = ", "))
}

# ============================================================
# 1. Configuration
# ============================================================
KEEP_DATASETS          <- "13"
GROUP_COL              <- "dataset"
USE_GROUPED_CV         <- FALSE
REMOVE_TEST_OUTLIERS   <- FALSE
USE_SPECTRAL_FILTER    <- FALSE
USE_Y_BOXPLOT          <- FALSE

if (REMOVE_TEST_OUTLIERS || USE_SPECTRAL_FILTER || USE_Y_BOXPLOT) {
  stop("This script requires all outlier-filter options to be FALSE.")
}

TARGET_NAME            <- "cp"
TARGET_COL_FALLBACK    <- 7L

INPUT_CSV <- "Spectral_library_06-Oct-2026.csv"

# Open the accompanying .Rproj file before running this script.
PROJECT_ROOT   <- here::here()
CSV_DIR        <- here::here("data", "raw")
PREPROC_DIR    <- here::here("data", "processed")
SUMMARY_DIR    <- here::here("Results", "Model_Summaries")
BEST_MODELS_DIR <- here::here("Models")
BEST_MODELS_MBL_DIR <- file.path(BEST_MODELS_DIR, "MBL")

BIN_SIZES <- c(5, 25, 50, 75)

N_OUTER <- 5L
N_INNER <- 5L
SEED <- 42L
N_CORES <- max(1L, parallel::detectCores() - 1L)

# MBL hyperparameters
MBL_K_GRID      <- c(10L, 25L, 50L, 100L)
MBL_FIT_METHOD  <- "wapls"
MBL_NCOMP_MIN   <- 1L
MBL_NCOMP_MAX   <- 10L
MBL_DISS_METHOD <- "diss_pca"
MBL_DISS_USAGE  <- "weights"

dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BEST_MODELS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BEST_MODELS_MBL_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Load latest preprocessing workspace
# ============================================================
snap <- fileSnapshot(PREPROC_DIR)
if (is.null(snap$info) || nrow(snap$info) == 0L) {
  stop("No preprocessing directories found in: ", PREPROC_DIR)
}
latest_dir <- rownames(snap$info[which.max(snap$info$mtime), , drop = FALSE])
rdata_files <- list.files(file.path(PREPROC_DIR, latest_dir),
                          pattern = "\\.RData$", full.names = TRUE)
if (length(rdata_files) == 0L) {
  stop("No .RData file found in ", file.path(PREPROC_DIR, latest_dir))
}
load(rdata_files[which.max(file.mtime(rdata_files))])
cat("Loaded preprocessing workspace:", latest_dir, "\n")

# ============================================================
# 3. Identify preprocessing variants
# ============================================================
spec_vars <- sort(ls(pattern = "^spec_\\d{2}_"))
if (length(spec_vars) == 0L) {
  stop("No preprocessing objects matching ^spec_\\d{2}_ were found.")
}
cat("\nAvailable preprocessing variants:\n"); print(spec_vars)

for (v in spec_vars) {
  ok <- tryCatch(is.matrix(get_spc_matrix(get(v))), error = function(e) FALSE)
  cat(sprintf("  %-22s : %s\n", v, if (ok) "OK" else "FAILED"))
  if (!ok) stop("Could not extract spectral matrix from ", v)
}

# ============================================================
# 4. Read raw CSV and target data
# ============================================================
read_csv_eu <- function(path) {
  df <- read.csv(path, sep = ",", header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE,
                 quote = "\"", comment.char = "")
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      x2 <- gsub(",", ".", x, fixed = TRUE)
      num <- suppressWarnings(as.numeric(x2))
      if (mean(!is.na(num)) > 0.5) return(num)
    }
    x
  })
  df
}

raw <- read_csv_eu(file.path(CSV_DIR, INPUT_CSV))
hdr_clean <- gsub("^X", "", colnames(raw))

TARGET_COL <- match(TARGET_NAME, hdr_clean)
if (is.na(TARGET_COL)) {
  cat("Target '", TARGET_NAME, "' not found by name; using column ",
      TARGET_COL_FALLBACK, ".\n", sep = "")
  TARGET_COL <- TARGET_COL_FALLBACK
}
if (TARGET_COL > ncol(raw)) stop("TARGET_COL_FALLBACK exceeds CSV columns.")
cat("Target column:", TARGET_COL, "->", hdr_clean[TARGET_COL], "\n")

if (!GROUP_COL %in% hdr_clean) stop("Grouping column '", GROUP_COL, "' not found.")

ids <- as.character(raw[[1]])
cp_raw <- suppressWarnings(as.numeric(raw[[TARGET_COL]]))
groups <- as.character(raw[[match(GROUP_COL, hdr_clean)]])

ok <- is.finite(cp_raw) & !is.na(groups) & groups != ""
cp_raw <- cp_raw[ok]; groups <- groups[ok]; ids <- ids[ok]
names(cp_raw) <- ids; names(groups) <- ids

if (!is.null(KEEP_DATASETS)) {
  keep <- groups %in% KEEP_DATASETS
  cp_raw <- cp_raw[keep]; groups <- groups[keep]
  cat("\nSubset:", paste(KEEP_DATASETS, collapse = ", "),
      "| n =", length(cp_raw), "\n")
}
if (anyDuplicated(names(cp_raw))) stop("Duplicate sample IDs found in CSV.")
cat("\nTarget summary for analysis subset:\n"); print(summary(cp_raw))
if (length(cp_raw) < 50L) warning("Fewer than 50 target observations available.")

# ============================================================
# 5. Helper functions
# ============================================================
clean_matrix <- function(X, tag = "") {
  X <- as.matrix(X); storage.mode(X) <- "double"
  row_ok <- rowSums(!is.finite(X)) == 0
  if (any(!row_ok)) {
    cat("  [", tag, "] dropping ", sum(!row_ok), " non-finite rows\n", sep = "")
    X <- X[row_ok, , drop = FALSE]
  }
  if (nrow(X) > 1L) {
    rsd <- apply(X, 1, sd); row_bad <- !is.finite(rsd) | rsd == 0
    if (any(row_bad)) {
      cat("  [", tag, "] dropping ", sum(row_bad), " zero-variance rows\n", sep = "")
      X <- X[!row_bad, , drop = FALSE]
    }
  }
  if (ncol(X) > 1L) {
    csd <- apply(X, 2, sd); col_bad <- !is.finite(csd) | csd == 0
    if (any(col_bad)) X <- X[, !col_bad, drop = FALSE]
  }
  X
}

bin_matrix <- function(X, bin_size) {
  if (bin_size <= 1L) return(X)
  p <- ncol(X)
  grp <- ((seq_len(p) - 1L) %/% bin_size) + 1L
  n_bins <- max(grp)
  M <- matrix(0, nrow = p, ncol = n_bins)
  M[cbind(seq_len(p), grp)] <- 1
  counts <- colSums(M)
  Xb <- X %*% M
  sweep(Xb, 2, counts, "/")
}

make_folds <- function(groups, y = NULL, k = 5L, seed = 42L, grouped = TRUE) {
  set.seed(seed)
  n <- length(groups)
  if (!grouped) {
    fold <- integer(n)
    for (g in unique(groups)) {
      idx <- which(groups == g)
      fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
    }
    return(fold)
  }
  grps <- unique(groups); k_use <- min(k, length(grps))
  if (k_use < 2L) stop("Grouped CV requires at least two distinct groups.")
  if (is.null(y)) {
    sizes <- as.numeric(table(groups)[grps])
    ord <- order(sizes, decreasing = TRUE)
    fold_of_group <- integer(length(grps)); names(fold_of_group) <- grps
    fold_load <- integer(k_use)
    for (g in grps[ord]) {
      tgt <- which.min(fold_load)
      fold_of_group[g] <- tgt
      fold_load[tgt] <- fold_load[tgt] + sizes[grps == g]
    }
    return(unname(fold_of_group[groups]))
  }
  g_n <- as.numeric(tapply(y, groups, length)[grps])
  g_sy <- as.numeric(tapply(y, groups, sum)[grps])
  g_n[!is.finite(g_n)] <- 0; g_sy[!is.finite(g_sy)] <- 0
  ord <- order(g_n, decreasing = TRUE)
  grps_ord <- grps[ord]; g_n_ord <- g_n[ord]; g_sy_ord <- g_sy[ord]
  target_n <- sum(g_n_ord) / k_use; target_sy <- sum(g_sy_ord) / k_use
  denom_n <- max(target_n^2, 1e-9); denom_sy <- max(target_sy^2, 1e-9)
  fold_n <- numeric(k_use); fold_sy <- numeric(k_use)
  fold_of_group <- integer(length(grps_ord))
  for (i in seq_along(grps_ord)) {
    new_n <- fold_n + g_n_ord[i]; new_sy <- fold_sy + g_sy_ord[i]
    score <- (new_n - target_n)^2 / denom_n +
      (new_sy - target_sy)^2 / denom_sy
    tgt <- which.min(score)
    fold_of_group[i] <- tgt
    fold_n[tgt] <- new_n[tgt]; fold_sy[tgt] <- new_sy[tgt]
  }
  fold_vec <- integer(n)
  for (i in seq_along(grps_ord))
    fold_vec[groups == grps_ord[i]] <- fold_of_group[i]
  fold_vec
}

make_mbl_fit_method <- function() {
  switch(MBL_FIT_METHOD,
         "pls" = resemble::fit_pls(ncomp = MBL_NCOMP_MAX),
         "wapls" = resemble::fit_wapls(min_ncomp = MBL_NCOMP_MIN,
                                       max_ncomp = MBL_NCOMP_MAX),
         "gpr" = resemble::fit_gpr(),
         stop("Unknown MBL_FIT_METHOD: ", MBL_FIT_METHOD))
}

# ============================================================
# 6. One MBL fit and prediction (unchanged from original)
# ============================================================
fit_mbl_predict <- function(X_train, y_train, X_new, y_new, k, seed) {
  X_train <- as.matrix(X_train); X_new <- as.matrix(X_new)
  colnames(X_train) <- paste0("V", seq_len(ncol(X_train)))
  colnames(X_new) <- paste0("V", seq_len(ncol(X_new)))
  
  if (k >= nrow(X_train)) {
    return(list(pred = rep(NA_real_, nrow(X_new)),
                status = "k_too_large", model = NULL))
  }
  
  fit <- tryCatch(
    resemble::mbl(
      Xr = X_train, Yr = y_train,
      Xu = X_new, Yu = y_new,
      neighbors = resemble::neighbors_k(k = k),
      diss_method = resemble::diss_pca(),
      diss_usage = MBL_DISS_USAGE,
      fit_method = make_mbl_fit_method(),
      control = resemble::mbl_control(validation_type = "NNv"),
      seed = seed, verbose = FALSE),
    error = function(e) NULL)
  
  if (is.null(fit))
    return(list(pred = rep(NA_real_, nrow(X_new)),
                status = "fit_failed", model = NULL))
  
  pred <- tryCatch(fit$results[[1]]$pred, error = function(e) NULL)
  if (is.null(pred) || length(pred) != nrow(X_new))
    return(list(pred = rep(NA_real_, nrow(X_new)),
                status = "prediction_failed", model = NULL))
  
  list(pred = as.numeric(pred), status = "ok", model = fit)
}

# ============================================================
# 7. Select best k for ONE candidate via inner CV
#    Selection: MINIMUM mean fold-level inner RMSEP
# ============================================================
select_k_for_candidate <- function(X, y_train, groups_train,
                                   n_inner, seed, grouped, k_grid) {
  inner_folds <- make_folds(groups_train, y_train, n_inner, seed, grouped)
  fold_ids <- sort(unique(inner_folds))
  
  smallest_train <- min(vapply(fold_ids,
                               function(f) sum(inner_folds != f),
                               integer(1)))
  k_grid_valid <- k_grid[k_grid < smallest_train]
  if (length(k_grid_valid) == 0L) return(NULL)
  
  best_k <- NA_integer_; best_rmsep <- Inf
  
  for (k_value in k_grid_valid) {
    fold_rmsep <- rep(NA_real_, length(fold_ids))
    for (ii in seq_along(fold_ids)) {
      f  <- fold_ids[ii]
      tr <- inner_folds != f
      va <- inner_folds == f
      if (sum(tr) < 20L || sum(va) < 5L) next
      if (k_value >= sum(tr)) next
      
      res <- fit_mbl_predict(
        X_train = X[tr, , drop = FALSE], y_train = y_train[tr],
        X_new   = X[va, , drop = FALSE], y_new   = y_train[va],
        k = k_value, seed = seed + 1000L * ii + k_value)
      
      pred <- res$pred
      if (res$status != "ok" || length(pred) != sum(va) ||
          any(!is.finite(pred))) next
      fold_rmsep[ii] <- sqrt(mean((y_train[va] - pred)^2))
    }
    m <- mean(fold_rmsep, na.rm = TRUE)
    if (is.finite(m) && m < best_rmsep) { best_rmsep <- m; best_k <- k_value }
  }
  
  if (is.na(best_k)) return(NULL)
  list(k = as.integer(best_k), inner_rmsep = best_rmsep)
}

# ============================================================
# 8. Evaluate ALL candidates on ONE outer fold
# ============================================================
evaluate_all_candidates_on_outer <- function(X_train_list, X_test_list,
                                             y_train, y_test,
                                             groups_train, n_inner,
                                             seed, grouped) {
  results <- list()
  
  for (cand in names(X_train_list)) {
    X_tr <- X_train_list[[cand]]
    X_te <- X_test_list[[cand]]
    if (is.null(X_tr) || is.null(X_te)) next
    
    sel <- select_k_for_candidate(X_tr, y_train, groups_train,
                                  n_inner, seed, grouped, MBL_K_GRID)
    if (is.null(sel)) next
    
    # Refit candidate on outer-train with chosen k, predict outer-test
    fit_out <- fit_mbl_predict(
      X_train = X_tr, y_train = y_train,
      X_new   = X_te, y_new   = NULL,
      k = sel$k, seed = seed + 200L)
    
    pred_out <- fit_out$pred
    if (fit_out$status != "ok" || length(pred_out) != nrow(X_te) ||
        any(!is.finite(pred_out))) next
    
    rmsep  <- sqrt(mean((y_test - pred_out)^2))
    mae    <- mean(abs(y_test - pred_out))
    bias   <- mean(pred_out - y_test)
    r2_cor <- if (sd(pred_out) > 0 && sd(y_test) > 0)
      cor(pred_out, y_test)^2 else NA_real_
    denom  <- sum((y_test - mean(y_test))^2)
    r2_ss  <- if (denom > 0) 1 - sum((y_test - pred_out)^2) / denom else NA_real_
    rpd    <- if (is.finite(rmsep) && rmsep > 0) sd(y_test) / rmsep else NA_real_
    
    results[[cand]] <- list(
      metric = data.frame(
        candidate   = cand,
        best_k      = as.integer(sel$k),
        inner_rmsep = sel$inner_rmsep,
        n_train     = nrow(X_tr),
        n_test      = nrow(X_te),
        R2          = r2_cor,
        R2_ss       = r2_ss,
        RMSEP       = rmsep,
        MAE         = mae,
        bias        = bias,
        RPD         = rpd,
        stringsAsFactors = FALSE),
      pred = data.frame(
        sample_id = rownames(X_te),
        y_true    = y_test,
        y_pred    = pred_out,
        best_k    = as.integer(sel$k),
        candidate = cand,
        stringsAsFactors = FALSE))
  }
  
  if (!length(results)) return(NULL)
  list(metrics     = do.call(rbind, lapply(results, `[[`, "metric")),
       predictions = do.call(rbind, lapply(results, `[[`, "pred")))
}

# ============================================================
# 9. Pre-bin spectra and build the common sample set
# ============================================================
cat("\nPre-binning spectra...\n")
binned <- list()

for (v in spec_vars) {
  X <- get_spc_matrix(get(v)); storage.mode(X) <- "double"
  X <- clean_matrix(X, tag = v)
  if (is.null(rownames(X))) { warning("Skipping ", v, ": no sample IDs."); next }
  if (anyDuplicated(rownames(X))) stop("Duplicate spectral sample IDs in ", v)
  
  common <- intersect(rownames(X), names(cp_raw))
  if (length(common) < 30L) {
    cat("  Skipping", v, "- too few overlapping samples\n"); next
  }
  X <- X[common, , drop = FALSE]
  binned[[v]] <- list()
  
  for (bs in BIN_SIZES) {
    Xb <- bin_matrix(X, bs)
    if (ncol(Xb) > 1L) {
      csd <- apply(Xb, 2, sd); col_bad <- !is.finite(csd) | csd == 0
      if (any(col_bad)) Xb <- Xb[, !col_bad, drop = FALSE]
    }
    rownames(Xb) <- common
    binned[[v]][[as.character(bs)]] <- Xb
  }
  n_cols <- vapply(BIN_SIZES,
                   function(b) ncol(binned[[v]][[as.character(b)]]), integer(1))
  cat(sprintf("  %-22s : %d samples x [%s] bands\n",
              v, nrow(X), paste(n_cols, collapse = "/")))
}
if (length(binned) == 0L) stop("No valid preprocessing variants remain.")

# Build candidate keys and matrices
candidate_keys <- expand.grid(
  variant = names(binned), bin_size = BIN_SIZES, stringsAsFactors = FALSE)
candidate_keys$key <- paste(candidate_keys$variant,
                            candidate_keys$bin_size, sep = "__BIN__")

candidate_matrices <- list()
for (i in seq_len(nrow(candidate_keys))) {
  v  <- candidate_keys$variant[i]
  bs <- candidate_keys$bin_size[i]
  if (is.null(binned[[v]][[as.character(bs)]])) next
  candidate_matrices[[candidate_keys$key[i]]] <- binned[[v]][[as.character(bs)]]
}
if (!length(candidate_matrices))
  stop("No candidate preprocessing/binning matrices remain.")

common_ids <- Reduce(intersect, lapply(candidate_matrices, rownames))
common_ids <- common_ids[common_ids %in% names(cp_raw) &
                           common_ids %in% names(groups)]
common_ids <- sort(common_ids)
if (length(common_ids) < 50L)
  stop("Fewer than 50 samples common to all candidates.")

candidate_matrices <- lapply(candidate_matrices,
                             function(X) X[common_ids, , drop = FALSE])
y <- as.numeric(cp_raw[common_ids]); names(y) <- common_ids
g <- as.character(groups[common_ids]); names(g) <- common_ids
if (any(!is.finite(y)))
  stop("Non-finite target values after common-sample matching.")

cat("\nCommon sample set:", length(y), "\n")
cat("Number of candidate configurations:", length(candidate_matrices), "\n")
cat("Candidate preprocessing variants:", length(unique(candidate_keys$variant)), "\n")
cat("Candidate bin sizes:", paste(BIN_SIZES, collapse = ", "), "\n")

# ============================================================
# 10. Fully nested outer CV — all candidates on every fold
# ============================================================
outer_folds <- make_folds(groups = g, y = y, k = N_OUTER, seed = SEED,
                          grouped = USE_GROUPED_CV)
actual_outer_folds <- sort(unique(outer_folds))
if (length(actual_outer_folds) < 2L)
  stop("Outer CV produced fewer than two folds.")

cat("\nStarting fully nested MBL CV (all candidates per fold)...\n")
cat("Outer folds:", length(actual_outer_folds),
    "| Inner folds:", N_INNER, "\n")
cat("Parallel workers:", min(N_CORES, length(actual_outer_folds)), "\n")

tic("Fully nested MBL CV")

n_workers <- min(N_CORES, length(actual_outer_folds))
cl <- parallel::makeCluster(n_workers, type = "PSOCK")
doParallel::registerDoParallel(cl)

clusterExport(cl, varlist = c(
  "candidate_matrices", "candidate_keys", "y", "g", "common_ids",
  "outer_folds", "make_folds", "make_mbl_fit_method",
  "fit_mbl_predict", "select_k_for_candidate",
  "evaluate_all_candidates_on_outer",
  "MBL_K_GRID", "MBL_FIT_METHOD", "MBL_NCOMP_MIN", "MBL_NCOMP_MAX",
  "MBL_DISS_USAGE", "N_INNER", "SEED", "USE_GROUPED_CV"),
  envir = environment())

clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(resemble))
  options(cores = 1L)
  NULL
})

outer_results <- foreach(
  f = actual_outer_folds,
  .packages = "resemble",
  .errorhandling = "pass",
  .export = c(
    "candidate_matrices", "candidate_keys", "y", "g", "common_ids",
    "outer_folds", "make_folds", "make_mbl_fit_method",
    "fit_mbl_predict", "select_k_for_candidate",
    "evaluate_all_candidates_on_outer",
    "MBL_K_GRID", "MBL_FIT_METHOD", "MBL_NCOMP_MIN", "MBL_NCOMP_MAX",
    "MBL_DISS_USAGE", "N_INNER", "SEED", "USE_GROUPED_CV")
) %dopar% {
  tr <- outer_folds != f
  te <- outer_folds == f
  if (sum(tr) < 30L || sum(te) < 5L)
    return(list(error = paste("Outer fold", f, "has too few samples")))
  
  X_tr_list <- lapply(candidate_matrices, function(X) X[tr, , drop = FALSE])
  X_te_list <- lapply(candidate_matrices, function(X) X[te, , drop = FALSE])
  
  res <- tryCatch(
    evaluate_all_candidates_on_outer(
      X_train_list = X_tr_list,
      X_test_list  = X_te_list,
      y_train      = y[tr],
      y_test       = y[te],
      groups_train = g[tr],
      n_inner      = N_INNER,
      seed         = SEED + f,
      grouped      = USE_GROUPED_CV),
    error = function(e) e)
  
  if (inherits(res, "error") || is.null(res))
    return(list(error = paste("Fold", f, ":", conditionMessage(res))))
  
  res$outer_fold <- f
  res
}

toc()
parallel::stopCluster(cl)
foreach::registerDoSEQ()

# ============================================================
# 11. Collect outer-fold results
# ============================================================
cat("\nCollecting outer-fold results...\n")
valid_results <- list(); errors <- character(0)
for (r in outer_results) {
  if (inherits(r, "error") || inherits(r, "condition") ||
      inherits(r, "try-error")) {
    errors <- c(errors, as.character(r)); next
  }
  if (!is.null(r$error)) { errors <- c(errors, r$error); next }
  valid_results[[length(valid_results) + 1L]] <- r
}
if (length(errors)) { cat("\nOuter-fold errors:\n"); print(errors) }
if (!length(valid_results))
  stop("No successful outer-fold evaluations produced.")

all_fold_metrics <- list(); all_fold_preds <- list()

for (r in valid_results) {
  f  <- r$outer_fold
  tr <- which(outer_folds != f); te <- which(outer_folds == f)
  
  # Resolve Preproc / Bin for metrics
  key_row <- match(r$metrics$candidate, candidate_keys$key)
  m_df <- r$metrics
  m_df$Preproc   <- candidate_keys$variant[key_row]
  m_df$Bin       <- as.integer(candidate_keys$bin_size[key_row])
  m_df$OuterFold <- f
  
  # Resolve for predictions
  p_key_row <- match(r$predictions$candidate, candidate_keys$key)
  p_df <- r$predictions
  p_df$Preproc   <- candidate_keys$variant[p_key_row]
  p_df$Bin       <- as.integer(candidate_keys$bin_size[p_key_row])
  p_df$OuterFold <- f
  p_df$group     <- as.character(setNames(g[te], common_ids[te])[p_df$sample_id])
  
  # Order to match attachment
  all_fold_metrics[[length(all_fold_metrics) + 1L]] <- m_df[, c(
    "Bin","Preproc","OuterFold","best_k","inner_rmsep",
    "n_train","n_test","R2","R2_ss","RMSEP","MAE","bias","RPD")]
  all_fold_preds[[length(all_fold_preds) + 1L]] <- p_df[, c(
    "Bin","Preproc","OuterFold","sample_id","group","y_true","y_pred","best_k")]
}

fold_metrics <- do.call(rbind, all_fold_metrics)
fold_preds   <- do.call(rbind, all_fold_preds)
rownames(fold_metrics) <- NULL
rownames(fold_preds)   <- NULL

cat("\nTotal fold-metrics rows   :", nrow(fold_metrics), "\n")
cat("Total fold-prediction rows:", nrow(fold_preds), "\n")

# ============================================================
# 12. Per-fold nested winner (min inner RMSEP)
# ============================================================
selection_history <- do.call(rbind, lapply(
  split(fold_metrics, fold_metrics$OuterFold),
  function(df) df[which.min(df$inner_rmsep), , drop = FALSE]))
rownames(selection_history) <- NULL

cat("\n--- Nested winner per outer fold (min inner RMSEP) ---\n")
print(selection_history[, c(
  "OuterFold","Preproc","Bin","best_k","inner_rmsep","RMSEP","RPD")])

# ============================================================
# 13. Aggregate per candidate (All Models Metrics)
# ============================================================
split_key <- paste(fold_metrics$Preproc, fold_metrics$Bin, sep = "||")

aggregate_per_candidate <- do.call(rbind, lapply(split(fold_metrics, split_key),
                                                 function(df) {
                                                   nf <- length(unique(df$OuterFold))
                                                   data.frame(
                                                     Preproc       = df$Preproc[1],
                                                     Bin_nm        = as.integer(df$Bin[1]),
                                                     n             = length(common_ids),
                                                     n_outer_folds = nf,
                                                     best_k        = mean(df$best_k, na.rm = TRUE),
                                                     R2_mean       = mean(df$R2,    na.rm = TRUE),
                                                     R2_sd         = sd  (df$R2,    na.rm = TRUE),
                                                     R2_SS_mean    = mean(df$R2_ss, na.rm = TRUE),
                                                     R2_SS_sd      = sd  (df$R2_ss, na.rm = TRUE),
                                                     RMSEP_mean    = mean(df$RMSEP, na.rm = TRUE),
                                                     RMSEP_sd      = sd  (df$RMSEP, na.rm = TRUE),
                                                     MAE_mean      = mean(df$MAE,   na.rm = TRUE),
                                                     MAE_sd        = sd  (df$MAE,   na.rm = TRUE),
                                                     bias_mean     = mean(df$bias,  na.rm = TRUE),
                                                     bias_sd       = sd  (df$bias,  na.rm = TRUE),
                                                     RPD_mean      = mean(df$RPD,   na.rm = TRUE),
                                                     RPD_sd        = sd  (df$RPD,   na.rm = TRUE),
                                                     RPD_margin_95 = 1.96 * sd(df$RPD, na.rm = TRUE) / sqrt(nf),
                                                     stringsAsFactors = FALSE)
                                                 }))

aggregate_per_candidate <- aggregate_per_candidate[
  order(-aggregate_per_candidate$R2_mean), ]
aggregate_per_candidate$Rank <- seq_len(nrow(aggregate_per_candidate))
aggregate_per_candidate <- aggregate_per_candidate[, c(
  "Rank","Preproc","Bin_nm","n","n_outer_folds","best_k",
  "R2_mean","R2_sd","R2_SS_mean","R2_SS_sd",
  "RMSEP_mean","RMSEP_sd","MAE_mean","MAE_sd",
  "bias_mean","bias_sd","RPD_mean","RPD_sd","RPD_margin_95")]
rownames(aggregate_per_candidate) <- NULL

cat("\n--- Aggregate outer-CV metrics per candidate (top 10) ---\n")
print(head(aggregate_per_candidate, 10))

# ============================================================
# 14. Winner-config aggregate (Panel D)
# ============================================================
winner_cfg  <- selection_history[1, c("Preproc","Bin")]
winner_rows <- fold_metrics[fold_metrics$Preproc == winner_cfg$Preproc &
                              fold_metrics$Bin   == winner_cfg$Bin, ]
n_folds_win <- nrow(winner_rows)

overall_summary <- data.frame(
  Metric = c("R2","R2_SS","RPD","RPD_margin_95","RMSEP","MAE","bias"),
  Mean_outer_fold = c(mean(winner_rows$R2,    na.rm = TRUE),
                      mean(winner_rows$R2_ss, na.rm = TRUE),
                      mean(winner_rows$RPD,   na.rm = TRUE),
                      NA_real_,
                      mean(winner_rows$RMSEP, na.rm = TRUE),
                      mean(winner_rows$MAE,   na.rm = TRUE),
                      mean(winner_rows$bias,  na.rm = TRUE)),
  SD_outer_fold   = c(sd(winner_rows$R2,    na.rm = TRUE),
                      sd(winner_rows$R2_ss, na.rm = TRUE),
                      sd(winner_rows$RPD,   na.rm = TRUE),
                      NA_real_,
                      sd(winner_rows$RMSEP, na.rm = TRUE),
                      sd(winner_rows$MAE,   na.rm = TRUE),
                      sd(winner_rows$bias,  na.rm = TRUE)),
  # Panel D uses 1.96 (same as All Models for MBL, unlike RF which used t)
  Value = c(NA_real_, NA_real_, NA_real_,
            1.96 * sd(winner_rows$RPD, na.rm = TRUE) / sqrt(n_folds_win),
            NA_real_, NA_real_, NA_real_),
  stringsAsFactors = FALSE)

cat("\n--- Winner-config aggregate (Panel D) ---\n")
print(overall_summary)

# ============================================================
# 15. Final deployment model selection on full data
# ============================================================
cat("\nSelecting deployment configuration on full dataset...\n")

full_inner <- make_folds(g, y, N_INNER, SEED + 100000L, USE_GROUPED_CV)
fold_ids_full <- sort(unique(full_inner))
smallest_full  <- min(vapply(fold_ids_full,
                             function(f) sum(full_inner != f), integer(1)))
k_grid_full <- MBL_K_GRID[MBL_K_GRID < smallest_full]

full_scores <- list(); rc <- 0L
for (cand in names(candidate_matrices)) {
  X <- candidate_matrices[[cand]]
  for (k_value in k_grid_full) {
    fold_rmsep <- rep(NA_real_, length(fold_ids_full))
    for (ii in seq_along(fold_ids_full)) {
      f  <- fold_ids_full[ii]
      tr <- full_inner != f; va <- full_inner == f
      if (sum(tr) < 20L || sum(va) < 5L) next
      res <- fit_mbl_predict(
        X_train = X[tr, , drop = FALSE], y_train = y[tr],
        X_new   = X[va, , drop = FALSE], y_new   = y[va],
        k = k_value, seed = SEED + 500L * ii + k_value)
      pred <- res$pred
      if (res$status != "ok" || length(pred) != sum(va) ||
          any(!is.finite(pred))) next
      fold_rmsep[ii] <- sqrt(mean((y[va] - pred)^2))
    }
    m <- mean(fold_rmsep, na.rm = TRUE)
    if (is.finite(m)) {
      rc <- rc + 1L
      full_scores[[rc]] <- data.frame(
        candidate = cand, k = as.integer(k_value),
        mean_inner_RMSEP = m, stringsAsFactors = FALSE)
    }
  }
}
if (!rc) stop("All MBL candidates failed during full-data inner CV.")

full_scores <- do.call(rbind, full_scores)
full_scores <- full_scores[order(full_scores$mean_inner_RMSEP,
                                 full_scores$k), ]
final_win  <- full_scores[1, , drop = FALSE]
final_key  <- as.character(final_win$candidate)
final_k    <- as.integer(final_win$k)
final_key_row <- match(final_key, candidate_keys$key)
if (is.na(final_key_row)) stop("Final selected candidate not found.")

WINNER_PREPROC <- candidate_keys$variant[final_key_row]
WINNER_BIN     <- as.integer(candidate_keys$bin_size[final_key_row])

cat("\nDeployment selection:\n")
cat("  Preproc:", WINNER_PREPROC, "| Bin:", WINNER_BIN,
    "| k:", final_k, "\n")
cat("  Inner CV RMSEP:", format(final_win$mean_inner_RMSEP, digits = 6), "\n")

# ============================================================
# 16. Fit final deployment model
# ============================================================
cat("\nFitting final MBL model on all common samples...\n")
X_full <- candidate_matrices[[final_key]]
colnames(X_full) <- paste0("V", seq_len(ncol(X_full)))

old_opts <- options(cores = 1L)
fit_deploy <- tryCatch(
  resemble::mbl(
    Xr = X_full, Yr = y,
    Xu = X_full, Yu = y,
    neighbors = resemble::neighbors_k(k = final_k),
    diss_method = resemble::diss_pca(),
    diss_usage = MBL_DISS_USAGE,
    fit_method = make_mbl_fit_method(),
    control = resemble::mbl_control(validation_type = "NNv"),
    seed = SEED + 200000L, verbose = FALSE),
  error = function(e) { warning("Deployment MBL fit failed: ",
                                conditionMessage(e)); NULL })
options(old_opts)

# ============================================================
# 17. Save model, metadata, and outputs
# ============================================================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
tag <- paste(KEEP_DATASETS, collapse = "_")
filter_tag <- paste0("spc", if (USE_SPECTRAL_FILTER) "on" else "off",
                     "_y", if (USE_Y_BOXPLOT) "on" else "off")

model_name <- sprintf("MBL_cp_ds%s_bin%d_%s_%s",
                      tag, WINNER_BIN,
                      gsub("^spec_", "", WINNER_PREPROC), timestamp)

if (!is.null(fit_deploy)) {
  model_file <- file.path(BEST_MODELS_MBL_DIR, paste0(model_name, ".rds"))
  saveRDS(fit_deploy, model_file)
  cat("\nDeployment model saved:", model_file, "\n")
} else {
  model_file <- NA_character_
  warning("Deployment fit failed: no model RDS was saved.")
}

meta <- list(
  algorithm       = "resemble::mbl (Memory-Based Learning)",
  target          = TARGET_NAME, target_col = TARGET_COL,
  subset          = KEEP_DATASETS, n_samples = length(y),
  selection_metric = "mean_inner_RMSEP",
  nested_cv       = TRUE, n_outer = N_OUTER, n_inner = N_INNER,
  seed            = SEED, grouped_cv = USE_GROUPED_CV,
  preprocessing   = WINNER_PREPROC, bin_size = WINNER_BIN,
  best_k          = final_k,
  full_data_inner_mean_RMSEP = final_win$mean_inner_RMSEP,
  fit_method      = MBL_FIT_METHOD,
  ncomp_min       = MBL_NCOMP_MIN, ncomp_max = MBL_NCOMP_MAX,
  diss_method     = MBL_DISS_METHOD, diss_usage = MBL_DISS_USAGE,
  feature_names   = colnames(X_full),
  train_ids       = rownames(X_full),
  y_range         = range(y, na.rm = TRUE),
  nested_cv_aggregate = aggregate_per_candidate,
  winner_cfg_aggregate = overall_summary,
  model_file      = model_file, timestamp = timestamp)
saveRDS(meta, file.path(BEST_MODELS_MBL_DIR, paste0(model_name, "_meta.rds")))

# --- Output prefix (matches attachment naming) ---
out_prefix <- file.path(SUMMARY_DIR,
                        paste0("MBL_cp_ds", tag, "_", filter_tag, "_",
                               if (USE_GROUPED_CV) "Grouped" else "Random",
                               "NestedCV_", timestamp))

# 1) foldmetrics — 40 × 5 = 200 rows
write.csv(fold_metrics, paste0(out_prefix, "_foldmetrics.csv"),
          row.names = FALSE)

# 2) foldpreds — 40 × 5 × ~61 ≈ 12,200 rows
write.csv(fold_preds, paste0(out_prefix, "_foldpreds.csv"),
          row.names = FALSE)

# 3) All Models Metrics — 40 rows
all_models_path <- file.path(SUMMARY_DIR,
                             "Supplementary Material 3 - MBL All Models Metrics.csv")
write.csv(aggregate_per_candidate, all_models_path,
          row.names = FALSE, na = "")

# ============================================================
# Panel A — Search space
# ============================================================
PANEL_COLS <- c("Model","Panel","RowType","OuterFold","Preproc","Bin_nm","best_k",
                "inner_RMSEP","n_train","n_test","R2","R2_SS","RMSEP",
                "MAE","bias","RPD","Metric","Mean","SD","Value","Note")

panel_A <- data.frame(
  Model   = "MBL", Panel = "A", RowType = "SearchSpace",
  Metric  = c("Preprocessing_variants", "Bin_widths_nm", "k_values",
              "Dissimilarity", "Local_model", "WAPLS_components",
              "Outer_CV_folds", "Inner_CV_folds"),
  Mean    = c(length(spec_vars), NA, NA,
              "diss_pca", "WAPLS", NA, N_OUTER, N_INNER),
  Note    = c(paste(length(spec_vars), "spectral preprocessing variants"),
              paste(BIN_SIZES, collapse = ";"),
              paste(MBL_K_GRID, collapse = ";"),
              NA, NA,
              paste0(MBL_NCOMP_MIN, "-", MBL_NCOMP_MAX),
              NA, NA),
  stringsAsFactors = FALSE)

# ============================================================
# Panel B — Nested selected config per outer fold
# ============================================================
panel_B <- data.frame(
  Model       = "MBL", Panel = "B", RowType = "NestedSelected",
  OuterFold   = selection_history$OuterFold,
  Preproc     = selection_history$Preproc,
  Bin_nm      = selection_history$Bin,
  best_k      = selection_history$best_k,
  inner_RMSEP = selection_history$inner_rmsep,
  Note        = paste0("Selected by minimum inner RMSEP in outer fold ",
                       selection_history$OuterFold),
  stringsAsFactors = FALSE)

# ============================================================
# Panel C — Outer metrics of per-fold winner
# ============================================================
panel_C <- data.frame(
  Model     = "MBL", Panel = "C", RowType = "OuterMetrics",
  OuterFold = selection_history$OuterFold,
  Preproc   = selection_history$Preproc,
  Bin_nm    = selection_history$Bin,
  best_k    = selection_history$best_k,
  n_train   = selection_history$n_train,
  n_test    = selection_history$n_test,
  R2        = selection_history$R2,
  R2_SS     = selection_history$R2_ss,
  RMSEP     = selection_history$RMSEP,
  MAE       = selection_history$MAE,
  bias      = selection_history$bias,
  RPD       = selection_history$RPD,
  stringsAsFactors = FALSE)

# ============================================================
# Panel D — Aggregate (winner config)
# ============================================================
panel_D <- data.frame(
  Model   = "MBL", Panel = "D", RowType = "Aggregate",
  Metric  = overall_summary$Metric,
  Mean    = overall_summary$Mean_outer_fold,
  SD      = overall_summary$SD_outer_fold,
  Value   = overall_summary$Value,
  Note    = c("Mean across five nested outer folds",
              NA_character_, NA_character_,
              "Approximate 95% margin across folds",
              NA_character_, NA_character_, NA_character_),
  stringsAsFactors = FALSE)

# ============================================================
# Panel E — Selection frequency (with k)
# ============================================================
selection_history$config_name <- paste0(
  selection_history$Preproc, "_bin", selection_history$Bin,
  "_bestk", selection_history$best_k)

config_folds <- split(selection_history$OuterFold,
                      selection_history$config_name)
config_folds <- config_folds[order(sapply(config_folds, min))]

combine_folds <- function(folds) {
  folds <- sort(unique(folds)); n <- length(folds)
  if (n == 1L) return(as.character(folds))
  if (n == 2L) return(paste(folds, collapse = " and "))
  paste0(paste(folds[-n], collapse = ", "), ", and ", folds[n])
}

panel_E <- data.frame(
  Model   = "MBL", Panel = "E", RowType = "Frequency",
  Metric  = names(config_folds),
  Value   = vapply(config_folds, length, integer(1)),
  Note    = vapply(config_folds, function(x)
    paste0("Selected in fold",
           if (length(x) > 1L) "s" else "", " ",
           combine_folds(x)),
    character(1)),
  stringsAsFactors = FALSE)

# ============================================================
# Combine A-E and write
# ============================================================
fill_cols <- function(df, cols) {
  for (cc in setdiff(cols, names(df))) df[[cc]] <- NA
  df[, cols, drop = FALSE]
}

panels_all <- rbind(
  fill_cols(panel_A, PANEL_COLS),
  fill_cols(panel_B, PANEL_COLS),
  fill_cols(panel_C, PANEL_COLS),
  fill_cols(panel_D, PANEL_COLS),
  fill_cols(panel_E, PANEL_COLS))
rownames(panels_all) <- NULL

panel_path <- file.path(SUMMARY_DIR, "Supplementary Material 3 - MBL.csv")
write.csv(panels_all, panel_path, row.names = FALSE, na = "")

# ============================================================
# Final console summary
# ============================================================
cat("\n=== Saved outputs ===\n")
cat("  ", panel_path,                                    "\n", sep = "")
cat("  ", all_models_path,                                "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldmetrics.csv"),         "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldpreds.csv"),           "\n", sep = "")
if (!is.na(model_file)) cat("  ", model_file,             "\n", sep = "")
cat("  ", file.path(BEST_MODELS_MBL_DIR,
                    paste0(model_name, "_meta.rds")),     "\n", sep = "")
cat("\nDONE.\n")
