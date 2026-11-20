#********************************************************************#
# A New Approach to Predict Nutrient Content in Plants               #
#        Using VNIR-SWIR Spectroscopy and Machine-Learning           #
#              Jhohan Elguera - MATE University                      #
#               With contributions of Dr. Johan PERRET               #
#                        PART 3 RANDOM FOREST MODELING               #
#       FULLY-NESTED CV — ALL CANDIDATES ON EVERY OUTER FOLD         #
#                                    2026                            #
#********************************************************************#

#************************************************************************#
# Nested-CV RF   — Nitrogen prediction (target = cp)                     #
#                                                                        #
# SUBSET VERSION: dataset 13 only                                        #
#   - scans preprocessing variants x bin sizes                           #
#   - saves the winning ranger model to 4-Best_Models/RF/                #
#                                                                        #
#************************************************************************#

suppressPackageStartupMessages({
  library(ranger)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(tictoc)
  library(robustbase)
})

rm(list = ls())

# ============================================================
# CONFIG
# ============================================================
KEEP_DATASETS          <- "13"
GROUP_COL              <- "dataset"
USE_GROUPED_CV         <- FALSE
REMOVE_TEST_OUTLIERS   <- FALSE

USE_SPECTRAL_FILTER    <- FALSE
USE_Y_BOXPLOT          <- FALSE

TARGET_NAME            <- "cp"
TARGET_COL_FALLBACK    <- 7L

INPUT_CSV              <- "Spectral_library_06-Oct-2026.csv"
# Open the accompanying .Rproj file before running this script.
PROJECT_ROOT    <- here::here()
CSV_DIR         <- here::here("data", "raw")
PREPROC_DIR     <- here::here("data", "processed")
SUMMARY_DIR     <- here::here("Results", "Model_Summaries")
BEST_MODELS_DIR <- here::here("Models")
BEST_MODELS_RF_DIR     <- file.path(BEST_MODELS_DIR, "RF")

BIN_SIZES              <- c(5, 25, 50, 75)

N_OUTER                <- 5
N_INNER                <- 5
SEED                   <- 42
N_CORES                <- max(1, detectCores() - 1)

# ---- Random Forest hyperparameters -------------------------
RF_NUM_TREES_INNER     <- 200L
RF_NUM_TREES_FINAL     <- 800L
RF_MIN_NODE_SIZE       <- 5L
RF_SAMPLE_FRACTION     <- 0.632
RF_REPLACE             <- TRUE

SPC_FILTER_N_PC        <- 4
SPC_FILTER_ALPHA       <- 0.975

dir.create(SUMMARY_DIR,        showWarnings = FALSE, recursive = TRUE)
dir.create(BEST_MODELS_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(BEST_MODELS_RF_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load latest preprocessing workspace
# ============================================================
snap <- fileSnapshot(PREPROC_DIR)
latest_dir <- rownames(snap$info[which.max(snap$info$mtime), ])
rdata_files <- list.files(file.path(PREPROC_DIR, latest_dir),
                          pattern = "\\.RData$", full.names = TRUE)
if (length(rdata_files) == 0) stop("No .RData file found in ", latest_dir)
load(rdata_files[which.max(file.mtime(rdata_files))])
cat("Loaded preprocessing workspace:", latest_dir, "\n")

# ============================================================
# 2. Identify available preprocessing variants
# ============================================================
spec_vars <- sort(ls(pattern = "^spec_\\d{2}_"))
cat("Available preprocessing variants:\n"); print(spec_vars)

# ============================================================
# 3. Load raw CSV — cp target + grouping column
# ============================================================
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

raw <- read_csv_eu(file.path(CSV_DIR, INPUT_CSV))
hdr_clean <- gsub("^X", "", colnames(raw))

TARGET_COL <- match(TARGET_NAME, hdr_clean)
if (is.na(TARGET_COL)) {
  cat("'", TARGET_NAME, "' not found by name. Falling back to column ",
      TARGET_COL_FALLBACK, ".\n", sep = "")
  TARGET_COL <- TARGET_COL_FALLBACK
}
cat("Target column :", TARGET_COL, "->", hdr_clean[TARGET_COL], "\n")

if (!GROUP_COL %in% hdr_clean) {
  stop("Grouping column '", GROUP_COL, "' not found.\n",
       "Available columns: ", paste(hdr_clean, collapse = ", "))
}

ids    <- as.character(raw[[1]])
cp_raw <- suppressWarnings(as.numeric(raw[[TARGET_COL]]))
groups <- as.character(raw[[GROUP_COL]])

cat("Full library cp summary:\n"); print(summary(cp_raw))

ok <- is.finite(cp_raw) & !is.na(groups) & groups != ""
cp_raw <- cp_raw[ok]; groups <- groups[ok]; ids <- ids[ok]
names(cp_raw) <- ids; names(groups) <- ids

if (!is.null(KEEP_DATASETS)) {
  keep <- groups %in% KEEP_DATASETS
  cp_raw <- cp_raw[keep]
  groups <- groups[keep]
  cat("\nSubset to dataset(s):", paste(KEEP_DATASETS, collapse = ", "),
      "-> n =", length(cp_raw), "\n")
}
cat("Raw cp summary (subset):\n"); print(summary(cp_raw))

# ============================================================
# 3b. ID overlap diagnostic
# ============================================================
spec_names <- rownames(get(spec_vars[1])$spc)
matched    <- names(cp_raw)[names(cp_raw) %in% spec_names]
unmatched  <- names(cp_raw)[!names(cp_raw) %in% spec_names]
cat("\n--- ID overlap ---\n")
cat("Spectra samples       :", length(spec_names), "\n")
cat("cp samples after QC   :", length(cp_raw), "\n")
cat("cp \u2229 spectra          :", length(matched), "\n")
cat("cp not in spectra     :", length(unmatched), "\n")

# ============================================================
# 3c. Grouping sanity check
# ============================================================
cat("\n--- Grouping structure ---\n")
cat("Number of groups in '", GROUP_COL, "': ",
    length(unique(groups)), "\n", sep = "")
grp_tab <- sort(table(groups), decreasing = TRUE)
cat("Group sizes:\n"); print(grp_tab)

# ============================================================
# 4. Helper functions
# ============================================================
clean_matrix <- function(X, tag = "") {
  X <- as.matrix(X); storage.mode(X) <- "double"
  row_ok <- rowSums(!is.finite(X)) == 0
  if (any(!row_ok)) {
    cat("  [", tag, "] dropping ", sum(!row_ok), " non-finite rows\n", sep = "")
    X <- X[row_ok, , drop = FALSE]
  }
  if (nrow(X) > 1) {
    rsd <- apply(X, 1, sd); row_bad <- !is.finite(rsd) | rsd == 0
    if (any(row_bad)) {
      cat("  [", tag, "] dropping ", sum(row_bad), " zero-variance rows\n", sep = "")
      X <- X[!row_bad, , drop = FALSE]
    }
  }
  if (ncol(X) > 1) {
    csd <- apply(X, 2, sd); col_bad <- !is.finite(csd) | csd == 0
    if (any(col_bad)) {
      cat("  [", tag, "] dropping ", sum(col_bad), " zero-variance columns\n", sep = "")
      X <- X[, !col_bad, drop = FALSE]
    }
  }
  X
}

bin_matrix <- function(X, bin_size) {
  if (bin_size <= 1) return(X)
  p <- ncol(X)
  grp <- ((seq_len(p) - 1) %/% bin_size) + 1
  n_bins <- max(grp)
  M <- matrix(0, nrow = p, ncol = n_bins)
  M[cbind(seq_len(p), grp)] <- 1
  counts <- colSums(M)
  Xb <- X %*% M
  sweep(Xb, 2, counts, "/")
}

make_folds <- function(groups, y = NULL, k = 5, seed = 42, grouped = TRUE) {
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
  if (is.null(y)) {
    sizes <- as.numeric(table(groups)[grps])
    ord   <- order(sizes, decreasing = TRUE)
    fold_of_group <- integer(length(grps)); names(fold_of_group) <- grps
    fold_load <- integer(k_use)
    for (g in grps[ord]) {
      tgt <- which.min(fold_load)
      fold_of_group[g] <- tgt
      fold_load[tgt]   <- fold_load[tgt] + sizes[grps == g]
    }
    return(fold_of_group[groups])
  }
  g_n  <- as.numeric(tapply(y, groups, length)[grps])
  g_sy <- as.numeric(tapply(y, groups, sum)[grps])
  g_n[!is.finite(g_n)]   <- 0; g_sy[!is.finite(g_sy)] <- 0
  ord          <- order(g_n, decreasing = TRUE)
  grps_ord     <- grps[ord]; g_n_ord <- g_n[ord]; g_sy_ord <- g_sy[ord]
  target_n  <- sum(g_n_ord)  / k_use
  target_sy <- sum(g_sy_ord) / k_use
  denom_n   <- max(target_n^2,  1e-9); denom_sy <- max(target_sy^2, 1e-9)
  fold_n   <- numeric(k_use); fold_sy  <- numeric(k_use)
  fold_of_group <- integer(length(grps_ord))
  for (i in seq_along(grps_ord)) {
    new_n  <- fold_n  + g_n_ord[i]; new_sy <- fold_sy + g_sy_ord[i]
    score  <- (new_n  - target_n )^2 / denom_n + (new_sy - target_sy)^2 / denom_sy
    tgt <- which.min(score)
    fold_of_group[i] <- tgt
    fold_n[tgt]  <- new_n[tgt]; fold_sy[tgt] <- new_sy[tgt]
  }
  fold_vec <- integer(n)
  for (i in seq_along(grps_ord))
    fold_vec[groups == grps_ord[i]] <- fold_of_group[i]
  fold_vec
}

# ============================================================
# 5. Pre-bin all spectra matrices
# ============================================================
cat("\nPre-binning spectra...\n")
binned <- list()

for (v in spec_vars) {
  hs <- get(v)
  X  <- as.matrix(hs$spc)
  storage.mode(X) <- "double"
  X <- clean_matrix(X, tag = v)
  if (nrow(X) < 30) { cat("  Skipping", v, "- too few clean rows\n"); next }
  common <- intersect(rownames(X), names(cp_raw))
  if (length(common) < 30) {
    cat("  Skipping", v, "- too few overlapping samples\n"); next
  }
  X <- X[common, , drop = FALSE]
  
  binned[[v]] <- list()
  for (bs in BIN_SIZES) {
    Xb <- bin_matrix(X, bs)
    if (bs > 1) {
      csd <- apply(Xb, 2, sd)
      col_bad <- !is.finite(csd) | csd == 0
      if (any(col_bad)) Xb <- Xb[, !col_bad, drop = FALSE]
    }
    rownames(Xb) <- common
    binned[[v]][[as.character(bs)]] <- Xb
  }
  n_cols <- sapply(BIN_SIZES, function(b) ncol(binned[[v]][[as.character(b)]]))
  cat(sprintf("  %-22s : %d samples x [%s] bands\n",
              v, nrow(X), paste(n_cols, collapse = "/")))
}

# ============================================================
# 6. FULLY NESTED RF CV — evaluate ALL candidates on every outer fold
# ============================================================
cat("\nStarting fully nested RF CV (all candidates per fold)...\n")
cat("Outer folds:", N_OUTER, "| Inner folds:", N_INNER,
    "| Inner trees:", RF_NUM_TREES_INNER,
    "| Final trees:", RF_NUM_TREES_FINAL, "\n")

# ------------------------------------------------------------
# 6.1 Common sample set across candidates
# ------------------------------------------------------------
candidate_keys <- expand.grid(
  variant  = names(binned),
  bin_size = BIN_SIZES,
  stringsAsFactors = FALSE
)
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
  stop("No candidate spectral matrices are available.")

common_ids <- Reduce(intersect, lapply(candidate_matrices, rownames))
common_ids <- common_ids[common_ids %in% names(cp_raw) &
                           common_ids %in% names(groups)]
common_ids <- sort(common_ids)
if (length(common_ids) < 30L)
  stop("Fewer than 30 common observations across candidates.")
candidate_matrices <- lapply(candidate_matrices,
                             function(X) X[common_ids, , drop = FALSE])
y_all <- unname(cp_raw[common_ids]); names(y_all) <- common_ids
g_all <- unname(groups[common_ids]); names(g_all) <- common_ids

cat("Common observations:", length(common_ids),
    "| Candidate configurations:", length(candidate_matrices), "\n")

# ------------------------------------------------------------
# 6.2 Metric helper
# ------------------------------------------------------------
calc_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  if (length(obs) < 2L)
    return(c(R2 = NA_real_, R2_ss = NA_real_, RMSEP = NA_real_,
             MAE = NA_real_, bias = NA_real_, RPD = NA_real_))
  rmsep  <- sqrt(mean((obs - pred)^2)); sd_obs <- sd(obs)
  c(R2    = if (sd(obs) > 0 && sd(pred) > 0) cor(obs, pred)^2 else NA_real_,
    R2_ss = if (sum((obs - mean(obs))^2) > 0)
      1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2) else NA_real_,
    RMSEP = rmsep,
    MAE   = mean(abs(obs - pred)),
    bias  = mean(pred - obs),
    RPD   = if (is.finite(rmsep) && rmsep > 0) sd_obs / rmsep else NA_real_)
}

# ------------------------------------------------------------
# 6.3 RF fitting helpers
# ------------------------------------------------------------
fit_rf_model <- function(X, y, mtry, num_trees, seed,
                         importance = "none", num_threads = 1L) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  df <- data.frame(y = as.numeric(y), X, check.names = FALSE)
  ranger::ranger(
    y ~ ., data = df,
    num.trees       = as.integer(num_trees),
    mtry            = as.integer(mtry),
    min.node.size   = RF_MIN_NODE_SIZE,
    sample.fraction = RF_SAMPLE_FRACTION,
    replace         = RF_REPLACE,
    respect.unordered.factors = "ignore",
    importance      = importance,
    seed            = as.integer(seed),
    num.threads     = as.integer(num_threads),
    verbose         = FALSE)
}

predict_rf_model <- function(model, X) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  df <- data.frame(X, check.names = FALSE)
  as.numeric(predict(model, data = df)$predictions)
}

# ------------------------------------------------------------
# 6.4 Evaluate ALL candidates on ONE outer fold
#     For each candidate:
#       - Inner CV over mtry grid -> pick mtry with MIN mean inner RMSEP
#       - Refit candidate on outer-train with best mtry
#       - Predict outer-test, compute outer metrics
# ------------------------------------------------------------
evaluate_all_candidates_on_outer <- function(
    X_train_list, X_test_list, y_train, y_test,
    groups_train, n_inner, seed, grouped) {
  
  inner_fold <- make_folds(groups_train, y_train, n_inner, seed, grouped)
  k_inner    <- length(unique(inner_fold))
  res_list   <- list()
  
  for (cand in names(X_train_list)) {
    X_tr <- X_train_list[[cand]]
    X_te <- X_test_list[[cand]]
    if (is.null(X_tr) || is.null(X_te)) next
    p <- ncol(X_tr)
    if (p < 2L) next
    
    mtry_grid <- unique(as.integer(round(c(sqrt(p), p / 3, p / 6))))
    mtry_grid <- sort(unique(mtry_grid[
      is.finite(mtry_grid) & mtry_grid >= 1L & mtry_grid <= p]))
    if (!length(mtry_grid)) next
    
    # --- inner CV per mtry ---
    rmsep_by_mtry <- rep(NA_real_, length(mtry_grid))
    for (g in seq_along(mtry_grid)) {
      fold_rmsep <- rep(NA_real_, k_inner)
      for (j in seq_len(k_inner)) {
        idx_tr <- which(inner_fold != j)
        idx_va <- which(inner_fold == j)
        if (length(idx_tr) < 10L || length(idx_va) < 2L) next
        
        fit <- tryCatch(
          fit_rf_model(X_tr[idx_tr, , drop = FALSE], y_train[idx_tr],
                       mtry_grid[g], RF_NUM_TREES_INNER,
                       seed + 100L * j, num_threads = 1L),
          error = function(e) NULL)
        if (is.null(fit)) next
        
        pred <- tryCatch(
          predict_rf_model(fit, X_tr[idx_va, , drop = FALSE]),
          error = function(e) NULL)
        if (is.null(pred) || length(pred) != length(idx_va)) next
        
        fold_rmsep[j] <- sqrt(mean((y_train[idx_va] - pred)^2))
      }
      rmsep_by_mtry[g] <- mean(fold_rmsep, na.rm = TRUE)
    }
    if (all(!is.finite(rmsep_by_mtry))) next
    
    best_g      <- which.min(rmsep_by_mtry)
    best_mtry   <- mtry_grid[best_g]
    inner_rmsep <- rmsep_by_mtry[best_g]
    
    # --- refit candidate on outer-train, predict outer-test ---
    fit_out <- tryCatch(
      fit_rf_model(X_tr, y_train, best_mtry, RF_NUM_TREES_FINAL,
                   seed + 200L, num_threads = 1L),
      error = function(e) NULL)
    if (is.null(fit_out)) next
    
    pred_out <- tryCatch(
      predict_rf_model(fit_out, X_te),
      error = function(e) NULL)
    if (is.null(pred_out) || length(pred_out) != nrow(X_te) ||
        any(!is.finite(pred_out))) next
    
    m <- calc_metrics(y_test, pred_out)
    
    res_list[[cand]] <- list(
      metric = data.frame(
        candidate   = cand,
        mtry        = as.integer(best_mtry),
        inner_rmsep = inner_rmsep,
        n_train     = nrow(X_tr),
        n_test      = nrow(X_te),
        R2          = unname(m["R2"]),
        R2_ss       = unname(m["R2_ss"]),
        RMSEP       = unname(m["RMSEP"]),
        MAE         = unname(m["MAE"]),
        bias        = unname(m["bias"]),
        RPD         = unname(m["RPD"]),
        stringsAsFactors = FALSE),
      pred = data.frame(
        sample_id = rownames(X_te),
        y_true    = y_test,
        y_pred    = pred_out,
        candidate = cand,
        stringsAsFactors = FALSE))
  }
  
  if (!length(res_list)) return(NULL)
  list(metrics     = do.call(rbind, lapply(res_list, `[[`, "metric")),
       predictions = do.call(rbind, lapply(res_list, `[[`, "pred")))
}

# ------------------------------------------------------------
# 6.5 Parallel outer CV loop
# ------------------------------------------------------------
outer_fold <- make_folds(groups = g_all, y = y_all, k = N_OUTER,
                         seed = SEED, grouped = USE_GROUPED_CV)
outer_levels <- sort(unique(outer_fold))
k_outer <- length(outer_levels)

cl <- parallel::makeCluster(min(N_CORES, k_outer), type = "PSOCK")
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(ranger))
})
parallel::clusterExport(
  cl,
  varlist = c(
    "candidate_matrices", "candidate_keys",
    "y_all", "g_all", "common_ids",
    "outer_fold", "outer_levels",
    "N_INNER", "SEED", "USE_GROUPED_CV",
    "RF_NUM_TREES_INNER", "RF_NUM_TREES_FINAL",
    "RF_MIN_NODE_SIZE", "RF_SAMPLE_FRACTION", "RF_REPLACE",
    "make_folds", "calc_metrics",
    "fit_rf_model", "predict_rf_model",
    "evaluate_all_candidates_on_outer"
  ),
  envir = environment()
)

tic("Fully nested RF CV (all candidates per outer fold)")

outer_results <- foreach::foreach(
  f = outer_levels,
  .packages = "ranger",
  .errorhandling = "pass"
) %dopar% {
  idx_tr <- which(outer_fold != f)
  idx_te <- which(outer_fold == f)
  if (length(idx_tr) < 30L || length(idx_te) < 5L) {
    return(list(error = paste("Fold", f, "has insufficient samples")))
  }
  
  X_tr_list <- lapply(candidate_matrices, function(X) X[idx_tr, , drop = FALSE])
  X_te_list <- lapply(candidate_matrices, function(X) X[idx_te, , drop = FALSE])
  
  res <- tryCatch(
    evaluate_all_candidates_on_outer(
      X_train_list = X_tr_list,
      X_test_list  = X_te_list,
      y_train      = y_all[idx_tr],
      y_test       = y_all[idx_te],
      groups_train = g_all[idx_tr],
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

# ------------------------------------------------------------
# 6.6 Assemble outer-fold outputs
# ------------------------------------------------------------
good <- list(); errors <- character(0)
for (r in outer_results) {
  if (inherits(r, "error") || inherits(r, "condition")) {
    errors <- c(errors, conditionMessage(r)); next
  }
  if (!is.null(r$error)) { errors <- c(errors, r$error); next }
  good[[length(good) + 1L]] <- r
}
if (length(errors)) { cat("\nOuter-fold errors:\n"); print(errors) }
if (!length(good)) stop("No successful outer-fold evaluations produced.")

all_fold_metrics <- list(); all_fold_preds <- list()

for (r in good) {
  f <- r$outer_fold
  idx_tr <- which(outer_fold != f); idx_te <- which(outer_fold == f)
  
  # Resolve Preproc / Bin for metrics
  key_row <- match(r$metrics$candidate, candidate_keys$key)
  m_df <- r$metrics
  m_df$Preproc       <- candidate_keys$variant[key_row]
  m_df$Bin           <- as.integer(candidate_keys$bin_size[key_row])
  m_df$OuterFold     <- f
  m_df$min.node.size <- RF_MIN_NODE_SIZE
  
  # Resolve for predictions
  pred_key_row <- match(r$predictions$candidate, candidate_keys$key)
  p_df <- r$predictions
  p_df$Preproc   <- candidate_keys$variant[pred_key_row]
  p_df$Bin       <- as.integer(candidate_keys$bin_size[pred_key_row])
  p_df$OuterFold <- f
  p_df$group     <- as.character(setNames(g_all[idx_te], common_ids[idx_te])[p_df$sample_id])
  
  # Column order matching the attached files
  all_fold_metrics[[length(all_fold_metrics) + 1L]] <- m_df[, c(
    "Bin","Preproc","OuterFold","mtry","min.node.size","inner_rmsep",
    "n_train","n_test","R2","R2_ss","RMSEP","MAE","bias","RPD")]
  all_fold_preds[[length(all_fold_preds) + 1L]] <- p_df[, c(
    "Bin","Preproc","OuterFold","sample_id","group","y_true","y_pred")]
}

fold_metrics_all <- do.call(rbind, all_fold_metrics)
fold_preds_all   <- do.call(rbind, all_fold_preds)
rownames(fold_metrics_all) <- NULL
rownames(fold_preds_all)   <- NULL

cat("\nTotal fold-metrics rows   :", nrow(fold_metrics_all), "\n")
cat("Total fold-prediction rows:", nrow(fold_preds_all), "\n")

# ------------------------------------------------------------
# 6.7 Per-fold nested winner (min inner RMSEP within fold)
# ------------------------------------------------------------
selection_history <- do.call(rbind, lapply(
  split(fold_metrics_all, fold_metrics_all$OuterFold),
  function(df) df[which.min(df$inner_rmsep), , drop = FALSE]))
rownames(selection_history) <- NULL

cat("\n--- Nested winner per outer fold (min inner RMSEP) ---\n")
print(selection_history)

# ------------------------------------------------------------
# 6.8 Aggregate per candidate (All Models Metrics)
# ------------------------------------------------------------
split_key <- paste(fold_metrics_all$Preproc, fold_metrics_all$Bin, sep = "||")

aggregate_per_candidate <- do.call(rbind, lapply(split(fold_metrics_all, split_key), function(df) {
  nf <- length(unique(df$OuterFold))
  data.frame(
    Preproc       = df$Preproc[1],
    Bin_nm        = as.integer(df$Bin[1]),
    n             = length(common_ids),
    n_outer_folds = nf,
    mtry          = mean(df$mtry, na.rm = TRUE),   # mean across folds (numeric)
    R2_mean       = mean(df$R2,     na.rm = TRUE),
    R2_sd         = sd  (df$R2,     na.rm = TRUE),
    R2_SS_mean    = mean(df$R2_ss,  na.rm = TRUE),
    R2_SS_sd      = sd  (df$R2_ss,  na.rm = TRUE),
    RMSEP_mean    = mean(df$RMSEP,  na.rm = TRUE),
    RMSEP_sd      = sd  (df$RMSEP,  na.rm = TRUE),
    MAE_mean      = mean(df$MAE,    na.rm = TRUE),
    MAE_sd        = sd  (df$MAE,    na.rm = TRUE),
    bias_mean     = mean(df$bias,   na.rm = TRUE),
    bias_sd       = sd  (df$bias,   na.rm = TRUE),
    RPD_mean      = mean(df$RPD,    na.rm = TRUE),
    RPD_sd        = sd  (df$RPD,    na.rm = TRUE),
    RPD_margin_95 = 1.96 * sd(df$RPD, na.rm = TRUE) / sqrt(nf),
    stringsAsFactors = FALSE)
}))

aggregate_per_candidate <- aggregate_per_candidate[
  order(-aggregate_per_candidate$R2_mean), ]
aggregate_per_candidate$Rank <- seq_len(nrow(aggregate_per_candidate))
aggregate_per_candidate <- aggregate_per_candidate[, c(
  "Rank","Preproc","Bin_nm","n","n_outer_folds","mtry",
  "R2_mean","R2_sd","R2_SS_mean","R2_SS_sd",
  "RMSEP_mean","RMSEP_sd","MAE_mean","MAE_sd",
  "bias_mean","bias_sd","RPD_mean","RPD_sd","RPD_margin_95")]
rownames(aggregate_per_candidate) <- NULL

cat("\n--- Aggregate outer-CV metrics per candidate (top 10) ---\n")
print(head(aggregate_per_candidate, 10))

# ------------------------------------------------------------
# 6.9 Winner-config aggregate (Panel D)
# ------------------------------------------------------------
winner_cfg  <- selection_history[1, c("Preproc","Bin")]
winner_rows <- fold_metrics_all[
  fold_metrics_all$Preproc == winner_cfg$Preproc &
    fold_metrics_all$Bin   == winner_cfg$Bin, ]
n_folds_win <- nrow(winner_rows)
t_val <- qt(0.975, df = n_folds_win - 1)

overall_summary <- data.frame(
  Metric = c("R2","R2_SS","RPD","RPD_margin_95","RMSEP","MAE","bias"),
  Mean_outer_fold = c(mean(winner_rows$R2,     na.rm = TRUE),
                      mean(winner_rows$R2_ss,  na.rm = TRUE),
                      mean(winner_rows$RPD,    na.rm = TRUE),
                      NA_real_,
                      mean(winner_rows$RMSEP,  na.rm = TRUE),
                      mean(winner_rows$MAE,    na.rm = TRUE),
                      mean(winner_rows$bias,   na.rm = TRUE)),
  SD_outer_fold   = c(sd(winner_rows$R2,     na.rm = TRUE),
                      sd(winner_rows$R2_ss,  na.rm = TRUE),
                      sd(winner_rows$RPD,    na.rm = TRUE),
                      NA_real_,
                      sd(winner_rows$RMSEP,  na.rm = TRUE),
                      sd(winner_rows$MAE,    na.rm = TRUE),
                      sd(winner_rows$bias,   na.rm = TRUE)),
  Value = c(NA_real_, NA_real_, NA_real_,
            t_val * sd(winner_rows$RPD, na.rm = TRUE) / sqrt(n_folds_win),
            NA_real_, NA_real_, NA_real_),
  stringsAsFactors = FALSE
)
cat("\n--- Winner-config aggregate (Panel D) ---\n")
print(overall_summary)

# ============================================================
# 7. Final deployment model — full-data inner CV
# ============================================================
cat("\nSelecting final deployment configuration on full data...\n")

full_inner <- make_folds(g_all, y_all, N_INNER, SEED + 100000L, USE_GROUPED_CV)
k_full <- length(unique(full_inner))

final_scores <- list(); rc <- 0L

for (cand in names(candidate_matrices)) {
  X  <- candidate_matrices[[cand]]
  p  <- ncol(X)
  if (p < 2L) next
  
  mtry_grid <- unique(as.integer(round(c(sqrt(p), p / 3, p / 6))))
  mtry_grid <- sort(unique(mtry_grid[
    is.finite(mtry_grid) & mtry_grid >= 1L & mtry_grid <= p]))
  if (!length(mtry_grid)) next
  
  for (mt in mtry_grid) {
    fold_rmsep <- rep(NA_real_, k_full)
    for (j in seq_len(k_full)) {
      idx_tr <- which(full_inner != j)
      idx_va <- which(full_inner == j)
      if (length(idx_tr) < 10L || length(idx_va) < 2L) next
      
      fit <- tryCatch(
        fit_rf_model(X[idx_tr, , drop = FALSE], y_all[idx_tr],
                     mt, RF_NUM_TREES_INNER,
                     SEED + 500L * j, num_threads = 1L),
        error = function(e) NULL)
      if (is.null(fit)) next
      
      pred <- tryCatch(
        predict_rf_model(fit, X[idx_va, , drop = FALSE]),
        error = function(e) NULL)
      if (is.null(pred) || length(pred) != length(idx_va)) next
      
      fold_rmsep[j] <- sqrt(mean((y_all[idx_va] - pred)^2))
    }
    mr <- mean(fold_rmsep, na.rm = TRUE)
    if (is.finite(mr)) {
      rc <- rc + 1L
      final_scores[[rc]] <- data.frame(
        candidate        = cand,
        mtry             = as.integer(mt),
        inner_mean_RMSEP = mr,
        stringsAsFactors = FALSE)
    }
  }
}
if (!rc) stop("All RF candidates failed during full-data inner CV.")
final_scores <- do.call(rbind, final_scores)
final_scores <- final_scores[order(final_scores$inner_mean_RMSEP,
                                   final_scores$mtry), ]
final_win     <- final_scores[1, , drop = FALSE]
final_key     <- as.character(final_win$candidate)
final_mtry    <- as.integer(final_win$mtry)
final_key_row <- match(final_key, candidate_keys$key)
if (is.na(final_key_row)) stop("Final selected candidate not found.")

WINNER_PREPROC <- candidate_keys$variant[final_key_row]
WINNER_BIN     <- as.integer(candidate_keys$bin_size[final_key_row])

X_full <- candidate_matrices[[final_key]]
y_full <- y_all

fit_deploy <- fit_rf_model(
  X = X_full, y = y_full,
  mtry        = final_mtry,
  num_trees   = RF_NUM_TREES_FINAL,
  seed        = SEED,
  importance  = "permutation",
  num_threads = N_CORES)

cat("\nFinal deployment RF model:\n")
cat("  Preprocessing :", WINNER_PREPROC, "\n")
cat("  Bin width     :", WINNER_BIN, "\n")
cat("  mtry          :", final_mtry, "\n")
cat("  Trees         :", RF_NUM_TREES_FINAL, "\n")
cat("  Inner-CV RMSEP:", format(final_win$inner_mean_RMSEP, digits = 6), "\n")
cat("  Full n        :", length(y_full), "\n")

# ============================================================
# 8. Save model, metadata, and the four attached CSV outputs
# ============================================================
timestamp  <- format(Sys.time(), "%Y%m%d_%H%M")
model_name <- sprintf("RF_cp_ds%s_bin%d_%s_%s",
                      paste(KEEP_DATASETS, collapse = "_"),
                      WINNER_BIN,
                      gsub("^spec_", "", WINNER_PREPROC),
                      timestamp)
model_file <- file.path(BEST_MODELS_RF_DIR, paste0(model_name, ".rds"))
saveRDS(fit_deploy, model_file)

meta <- list(
  algorithm       = "ranger::ranger (Random Forest)",
  target          = TARGET_NAME,
  target_col      = TARGET_COL,
  subset          = KEEP_DATASETS,
  bin_size        = WINNER_BIN,
  preprocessing   = WINNER_PREPROC,
  n_samples       = length(y_full),
  mtry            = final_mtry,
  min.node.size   = RF_MIN_NODE_SIZE,
  num.trees       = RF_NUM_TREES_FINAL,
  sample.fraction = RF_SAMPLE_FRACTION,
  replace         = RF_REPLACE,
  seed            = SEED,
  feature_names   = colnames(X_full),
  train_ids       = rownames(X_full),
  y_range         = range(y_full, na.rm = TRUE),
  final_inner_cv        = as.list(final_win),
  final_inner_cv_scores = final_scores,
  outer_cv_metrics      = fold_metrics_all,
  outer_cv_selections   = selection_history,
  aggregate_per_candidate = aggregate_per_candidate,
  timestamp       = timestamp)
saveRDS(meta, file.path(BEST_MODELS_RF_DIR, paste0(model_name, "_meta.rds")))

# ------------------------------------------------------------
# Output prefix (matches attachment naming convention)
# ------------------------------------------------------------
tag <- paste(KEEP_DATASETS, collapse = "_")
filter_tag <- paste0("spc", if (USE_SPECTRAL_FILTER) "on" else "off",
                     "_y", if (USE_Y_BOXPLOT) "on" else "off")
out_prefix <- file.path(
  SUMMARY_DIR,
  paste0("RF_cp_ds", tag, "_", filter_tag, "_",
         if (USE_GROUPED_CV) "Grouped" else "Random",
         "NestedCV_", timestamp))

# 1) foldmetrics — 40 × 5 = 200 rows
write.csv(fold_metrics_all,
          paste0(out_prefix, "_foldmetrics.csv"),
          row.names = FALSE)

# 2) foldpreds — 40 × 5 × n_test ≈ 12,200 rows
write.csv(fold_preds_all,
          paste0(out_prefix, "_foldpreds.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# Panel A — Search space
# ------------------------------------------------------------
PANEL_COLS <- c("Model","Panel","RowType","OuterFold","Preproc","Bin_nm","mtry",
                "inner_RMSEP","n_train","n_test","R2","R2_SS","RMSEP",
                "MAE","bias","RPD","Metric","Mean","SD","Value","Note")

panel_A <- data.frame(
  Model   = "RF",
  Panel   = "A",
  RowType = "SearchSpace",
  Metric  = c("Preprocessing_variants", "Bin_widths_nm", "mtry_candidates",
              "min_node_size", "Inner_trees", "Final_trees",
              "sample_fraction", "Outer_CV_folds", "Inner_CV_folds"),
  Mean    = c(length(spec_vars), NA, NA,
              RF_MIN_NODE_SIZE, RF_NUM_TREES_INNER, RF_NUM_TREES_FINAL,
              RF_SAMPLE_FRACTION, N_OUTER, N_INNER),
  Note    = c(paste(length(spec_vars), "spectral preprocessing variants"),
              paste(BIN_SIZES, collapse = ";"),
              "sqrt(p);p/3;p/6",
              NA, NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Panel B — Nested selected config per outer fold
# ------------------------------------------------------------
panel_B <- data.frame(
  Model       = "RF",
  Panel       = "B",
  RowType     = "NestedSelected",
  OuterFold   = selection_history$OuterFold,
  Preproc     = selection_history$Preproc,
  Bin_nm      = selection_history$Bin,
  mtry        = selection_history$mtry,
  inner_RMSEP = selection_history$inner_rmsep,
  Note        = paste0("Selected by minimum inner RMSEP in outer fold ",
                       selection_history$OuterFold),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Panel C — Outer metrics of per-fold winner
# ------------------------------------------------------------
panel_C <- data.frame(
  Model     = "RF",
  Panel     = "C",
  RowType   = "OuterMetrics",
  OuterFold = selection_history$OuterFold,
  Preproc   = selection_history$Preproc,
  Bin_nm    = selection_history$Bin,
  mtry      = selection_history$mtry,
  n_train   = selection_history$n_train,
  n_test    = selection_history$n_test,
  R2        = selection_history$R2,
  R2_SS     = selection_history$R2_ss,
  RMSEP     = selection_history$RMSEP,
  MAE       = selection_history$MAE,
  bias      = selection_history$bias,
  RPD       = selection_history$RPD,
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Panel D — Aggregate (winner-config nested-CV summary)
# RPD_margin_95 uses t(0.975, df = n_folds - 1)
# (All Models Metrics uses 1.96 instead — see section 6.8)
# ------------------------------------------------------------
panel_D <- data.frame(
  Model   = "RF",
  Panel   = "D",
  RowType = "Aggregate",
  Metric  = overall_summary$Metric,
  Mean    = overall_summary$Mean_outer_fold,
  SD      = overall_summary$SD_outer_fold,
  Value   = overall_summary$Value,
  Note    = c("Mean across five nested outer folds",
              NA_character_, NA_character_,
              "Approximate 95% margin across folds",
              NA_character_, NA_character_, NA_character_),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Panel E — Frequency of winning configurations (with mtry)
# Sorted by first appearance (i.e. min fold where selected)
# ------------------------------------------------------------
selection_history$config_name <- paste0(
  selection_history$Preproc,
  "_bin",  selection_history$Bin,
  "_mtry", selection_history$mtry)

config_folds <- split(selection_history$OuterFold,
                      selection_history$config_name)
config_folds <- config_folds[order(sapply(config_folds, min))]

combine_folds <- function(folds) {
  folds <- sort(unique(folds))
  n <- length(folds)
  if (n == 1L) return(as.character(folds))
  if (n == 2L) return(paste(folds, collapse = " and "))
  paste0(paste(folds[-n], collapse = ", "), ", and ", folds[n])
}

panel_E <- data.frame(
  Model   = "RF",
  Panel   = "E",
  RowType = "Frequency",
  Metric  = names(config_folds),
  Value   = vapply(config_folds, length, integer(1)),
  Note    = vapply(config_folds, function(x)
    paste0("Selected in fold",
           if (length(x) > 1L) "s" else "", " ",
           combine_folds(x)),
    character(1)),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Combine Panels A–E and write to SUMMARY_DIR
# ------------------------------------------------------------
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

panel_rf_path <- file.path(SUMMARY_DIR,
                           "Supplementary Material 2 - RF.csv")
write.csv(panels_all, panel_rf_path,
          row.names = FALSE, na = "")

# ------------------------------------------------------------
# All Models Metrics CSV
# ------------------------------------------------------------
all_models_path <- file.path(SUMMARY_DIR,
                             "Supplementary Material 2 - RF All Models Metrics.csv")
write.csv(aggregate_per_candidate, all_models_path,
          row.names = FALSE, na = "")

# ------------------------------------------------------------
# Final console summary
# ------------------------------------------------------------
cat("\n=== Saved outputs ===\n")
cat("  ", panel_rf_path,                           "\n", sep = "")
cat("  ", all_models_path,                         "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldmetrics.csv"),  "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldpreds.csv"),    "\n", sep = "")
cat("  ", model_file,                              "\n", sep = "")
cat("  ", file.path(BEST_MODELS_RF_DIR,
                    paste0(model_name, "_meta.rds")), "\n", sep = "")
cat("\nDone.\n")
