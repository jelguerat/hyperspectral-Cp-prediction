#********************************************************************#
# A New Approach to Predict Nutrient Content in Plants               #
#        Using VNIR-SWIR Spectroscopy and Machine-Learning           #
#              Jhohan Elguera - MATE University                      #
#               With contributions of Dr. Johan PERRET               #
#                        PART 5  SVR MODELING                        #
#         FULLY-NESTED CV — ALL CANDIDATES ON EVERY OUTER FOLD       #
#                                    2026                            #
#********************************************************************#

#************************************************************************#
# Nested-CV SVR   — Nitrogen prediction (target = cp)                    #
#                                                                        #
# SUBSET VERSION: dataset 13 only                                        #
#   - scans preprocessing variants x bin sizes                           #
#   - saves the winning ranger model to 4-Best_Models/SVR/               #
#                                                                        #
#************************************************************************#

suppressPackageStartupMessages({
  library(e1071)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(tictoc)
})

rm(list = ls())

# ============================================================
# 0. CONFIGURATION
# ============================================================
KEEP_DATASETS <- "13"
GROUP_COL     <- "dataset"
USE_GROUPED_CV <- FALSE

REMOVE_TEST_OUTLIERS <- FALSE
USE_SPECTRAL_FILTER  <- FALSE
USE_Y_BOXPLOT        <- FALSE

TARGET_NAME          <- "cp"
TARGET_COL_FALLBACK  <- 7L

INPUT_CSV   <- "Spectral_library_06-Oct-2026.csv"
# Open the accompanying .Rproj file before running this script.
PROJECT_ROOT    <- here::here()
CSV_DIR         <- here::here("data", "raw")
PREPROC_DIR     <- here::here("data", "processed")
SUMMARY_DIR     <- here::here("Results", "Model_Summaries")
BEST_MODELS_DIR <- here::here("Models")
BEST_MODELS_SVR_DIR <- file.path(BEST_MODELS_DIR, "SVR")

BIN_SIZES <- c(5, 25, 50, 75)

N_OUTER <- 5L
N_INNER <- 5L
SEED    <- 42L
N_CORES <- max(1L, parallel::detectCores() - 1L)

SVR_KERNEL      <- "radial"
SVR_EPSILON     <- 0.1
SVR_COST_GRID   <- c(50, 100, 200, 500, 1000, 2000)
SVR_GAMMA_GRID  <- c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
SVR_COST_GRID_STR  <- "50;100;200;500;1000;2000"
SVR_GAMMA_GRID_STR <- "1e-4;3e-4;1e-3;3e-3;1e-2"
SVR_MAX_TRAIN   <- 4000L

dir.create(SUMMARY_DIR,         recursive = TRUE, showWarnings = FALSE)
dir.create(BEST_MODELS_DIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(BEST_MODELS_SVR_DIR, recursive = TRUE, showWarnings = FALSE)

if (USE_SPECTRAL_FILTER || USE_Y_BOXPLOT || REMOVE_TEST_OUTLIERS) {
  stop("This nested-selection script assumes all outlier-filter toggles are FALSE.")
}

# ============================================================
# 1. LOAD LATEST PREPROCESSING WORKSPACE
# ============================================================
snap <- fileSnapshot(PREPROC_DIR)
latest_dir <- rownames(snap$info[which.max(snap$info$mtime), ])
if (length(latest_dir) == 0L || is.na(latest_dir)) stop("No preprocessing directory.")

latest_path <- file.path(PREPROC_DIR, latest_dir)
rdata_files <- list.files(latest_path, pattern = "\\.RData$", full.names = TRUE)
if (length(rdata_files) == 0L) stop("No .RData file found in ", latest_path)
load(rdata_files[which.max(file.mtime(rdata_files))])
cat("Loaded preprocessing workspace:", latest_path, "\n")

# ============================================================
# 2. IDENTIFY PREPROCESSING VARIANTS
# ============================================================
spec_vars <- sort(ls(pattern = "^spec_\\d{2}_"))
if (length(spec_vars) == 0L) stop("No preprocessing objects found.")
cat("\nAvailable preprocessing variants:\n"); print(spec_vars)

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================
get_spc <- function(obj, tag = "") {
  if (!isS4(obj) && (is.list(obj) || is.environment(obj)) && !is.null(obj$spc))
    return(as.matrix(obj$spc))
  if (isS4(obj)) {
    sn <- slotNames(obj)
    if (inherits(obj, "hyperSpec")) return(as.matrix(obj@data$spc))
    if ("spc" %in% sn) return(as.matrix(obj@spc))
    if (inherits(obj, "Spectra") || "peaksData" %in% sn) {
      m <- tryCatch(peaksData(obj)[[1]], error = function(e) NULL)
      if (!is.null(m)) return(as.matrix(m))
    }
    if ("data" %in% sn && is.data.frame(obj@data) && "spc" %in% names(obj@data))
      return(as.matrix(obj@data$spc))
    stop("Cannot identify spectral data in S4 object: ", tag)
  }
  if (is.matrix(obj) || is.data.frame(obj)) return(as.matrix(obj))
  stop("Cannot extract spectra from class ", paste(class(obj), collapse = "/"), ". Tag: ", tag)
}

read_csv_eu <- function(path) {
  df <- read.csv(path, sep = ",", header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, quote = "\"", comment.char = "")
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

clean_matrix <- function(X, tag = "") {
  X <- as.matrix(X); storage.mode(X) <- "double"
  if (is.null(rownames(X))) stop("Spectral matrix has no sample IDs: ", tag)
  row_ok <- rowSums(!is.finite(X)) == 0
  if (any(!row_ok)) { cat("  [", tag, "] dropping ", sum(!row_ok), " non-finite rows\n", sep = ""); X <- X[row_ok, , drop = FALSE] }
  if (nrow(X) > 1L) {
    rsd <- apply(X, 1, sd); row_bad <- !is.finite(rsd) | rsd == 0
    if (any(row_bad)) { cat("  [", tag, "] dropping ", sum(row_bad), " zero-variance rows\n", sep = ""); X <- X[!row_bad, , drop = FALSE] }
  }
  if (ncol(X) > 1L) {
    csd <- apply(X, 2, sd); col_bad <- !is.finite(csd) | csd == 0
    if (any(col_bad)) X <- X[, !col_bad, drop = FALSE]
  }
  X
}

bin_matrix <- function(X, bin_size) {
  if (bin_size <= 1L) return(X)
  p <- ncol(X); grp <- ((seq_len(p) - 1L) %/% bin_size) + 1L
  M <- matrix(0, nrow = p, ncol = max(grp))
  M[cbind(seq_len(p), grp)] <- 1
  sweep(X %*% M, 2, colSums(M), "/")
}

make_folds <- function(groups, y = NULL, k = 5L, seed = 42L, grouped = TRUE) {
  set.seed(seed); n <- length(groups)
  if (!grouped) {
    fold <- integer(n)
    for (g in unique(groups)) {
      idx <- which(groups == g)
      fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
    }
    return(fold)
  }
  grps <- unique(groups); k_use <- min(k, length(grps))
  if (k_use < 2L) stop("Grouped CV requires at least two groups.")
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
    return(as.integer(fold_of_group[groups]))
  }
  g_n  <- as.numeric(tapply(y, groups, length)[grps])
  g_sy <- as.numeric(tapply(y, groups, sum)[grps])
  g_n[!is.finite(g_n)] <- 0; g_sy[!is.finite(g_sy)] <- 0
  ord <- order(g_n, decreasing = TRUE)
  grps_ord <- grps[ord]; g_n_ord <- g_n[ord]; g_sy_ord <- g_sy[ord]
  target_n <- sum(g_n_ord)/k_use; target_sy <- sum(g_sy_ord)/k_use
  denom_n <- max(target_n^2, 1e-9); denom_sy <- max(target_sy^2, 1e-9)
  fold_n <- numeric(k_use); fold_sy <- numeric(k_use)
  fold_of_group <- integer(length(grps_ord))
  for (i in seq_along(grps_ord)) {
    new_n <- fold_n + g_n_ord[i]; new_sy <- fold_sy + g_sy_ord[i]
    score <- (new_n - target_n)^2/denom_n + (new_sy - target_sy)^2/denom_sy
    tgt <- which.min(score)
    fold_of_group[i] <- tgt
    fold_n[tgt] <- new_n[tgt]; fold_sy[tgt] <- new_sy[tgt]
  }
  fold_vec <- integer(n)
  for (i in seq_along(grps_ord)) fold_vec[groups == grps_ord[i]] <- fold_of_group[i]
  fold_vec
}

# ============================================================
# 4. SVR FITTING
# ============================================================
fit_svr <- function(X_train, y_train, X_test,
                    cost, gamma, epsilon = 0.1) {
  fit <- tryCatch(
    e1071::svm(x = as.matrix(X_train), y = as.numeric(y_train),
               type = "eps-regression", kernel = "radial",
               cost = cost, gamma = gamma, epsilon = epsilon, scale = TRUE),
    error = function(e) NULL)
  if (is.null(fit))
    return(list(fit = NULL, pred = rep(NA_real_, nrow(X_test))))
  pred <- tryCatch(
    as.numeric(predict(fit, newdata = as.matrix(X_test))),
    error = function(e) rep(NA_real_, nrow(X_test)))
  list(fit = fit, pred = pred)
}

# ============================================================
# 5. INNER CV — per candidate: search (cost, gamma)
#    Selection: MIN mean fold-level inner RMSEP
# ============================================================
select_svr_config <- function(X_list, y_train, groups_train,
                              n_inner, seed, grouped,
                              cost_grid, gamma_grid, epsilon,
                              max_train = Inf) {
  if (length(X_list) == 0L) stop("No candidates supplied.")
  if (length(y_train) != length(groups_train)) stop("y / groups length mismatch.")
  
  if (is.finite(max_train) && length(y_train) > max_train) {
    set.seed(seed)
    take <- sort(sample.int(length(y_train), max_train))
    y_train <- y_train[take]; groups_train <- groups_train[take]
    X_list <- lapply(X_list, function(X) X[take, , drop = FALSE])
  }
  
  if (length(y_train) < n_inner * 5L)
    stop("Too few training observations for inner CV.")
  
  inner_folds <- make_folds(groups_train, y_train, n_inner, seed, grouped)
  fold_values <- sort(unique(inner_folds))
  if (length(fold_values) < 2L) stop("Inner CV < 2 folds.")
  
  audit_rows <- list(); audit_i <- 0L
  
  for (key in names(X_list)) {
    X <- X_list[[key]]
    if (nrow(X) != length(y_train)) stop("Candidate ", key, " mismatch.")
    
    grid <- expand.grid(cost = cost_grid, gamma = gamma_grid,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    
    for (gi in seq_len(nrow(grid))) {
      cost_i <- grid$cost[gi]; gamma_i <- grid$gamma[gi]
      fold_rpd <- rep(NA_real_, length(fold_values))
      fold_rmsep <- rep(NA_real_, length(fold_values))
      
      for (fi in seq_along(fold_values)) {
        f <- fold_values[fi]; tr <- inner_folds != f; va <- inner_folds == f
        if (sum(tr) < 10L || sum(va) < 5L) next
        
        res <- fit_svr(X[tr, , drop = FALSE], y_train[tr],
                       X[va, , drop = FALSE],
                       cost = cost_i, gamma = gamma_i, epsilon = epsilon)
        pred <- res$pred
        if (length(pred) != sum(va) || any(!is.finite(pred))) next
        rmsep_i <- sqrt(mean((y_train[va] - pred)^2))
        if (!is.finite(rmsep_i) || rmsep_i <= 0) next
        fold_rmsep[fi] <- rmsep_i
        fold_rpd[fi]   <- sd(y_train[va]) / rmsep_i
      }
      
      mean_rpd   <- if (all(!is.finite(fold_rpd)))   NA_real_ else mean(fold_rpd,   na.rm = TRUE)
      mean_rmsep <- if (all(!is.finite(fold_rmsep))) NA_real_ else mean(fold_rmsep, na.rm = TRUE)
      
      audit_i <- audit_i + 1L
      audit_rows[[audit_i]] <- data.frame(
        candidate = key, cost = cost_i, gamma = gamma_i,
        mean_inner_RMSEP = mean_rmsep,
        mean_inner_RPD   = mean_rpd,
        n_inner_folds_valid = sum(is.finite(fold_rmsep)),
        stringsAsFactors = FALSE)
    }
  }
  
  audit <- do.call(rbind, audit_rows)
  valid <- which(is.finite(audit$mean_inner_RMSEP))
  if (length(valid) == 0L) stop("All inner-CV SVR configs failed.")
  
  best <- audit[valid[which.min(audit$mean_inner_RMSEP[valid])], , drop = FALSE]
  
  list(best_key      = as.character(best$candidate),
       best_cost     = as.numeric(best$cost),
       best_gamma    = as.numeric(best$gamma),
       best_mean_RPD = as.numeric(best$mean_inner_RPD),
       best_mean_RMSEP = as.numeric(best$mean_inner_RMSEP),
       audit = audit)
}

# ============================================================
# 6. Per-candidate best (cost, gamma)
#    Selection: MIN mean inner RMSEP for each candidate
# ============================================================
pick_best_per_candidate <- function(audit, use_rmsep = TRUE) {
  candidates <- unique(audit$candidate)
  res <- lapply(candidates, function(cand) {
    sub <- audit[audit$candidate == cand, , drop = FALSE]
    if (use_rmsep) {
      valid <- which(is.finite(sub$mean_inner_RMSEP))
      if (!length(valid)) return(NULL)
      idx <- valid[which.min(sub$mean_inner_RMSEP[valid])]
    } else {
      valid <- which(is.finite(sub$mean_inner_RPD))
      if (!length(valid)) return(NULL)
      idx <- valid[which.max(sub$mean_inner_RPD[valid])]
    }
    b <- sub[idx, , drop = FALSE]
    data.frame(candidate  = as.character(b$candidate),
               cost       = as.numeric(b$cost),
               gamma      = as.numeric(b$gamma),
               inner_rmse = as.numeric(b$mean_inner_RMSEP),
               inner_RPD  = as.numeric(b$mean_inner_RPD),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

# ============================================================
# 7. LOAD RAW CSV, TARGET
# ============================================================
raw <- read_csv_eu(file.path(CSV_DIR, INPUT_CSV))
hdr_clean <- gsub("^X", "", colnames(raw))

TARGET_COL <- match(TARGET_NAME, hdr_clean)
if (is.na(TARGET_COL)) {
  cat("Target not found, using column", TARGET_COL_FALLBACK, "\n")
  TARGET_COL <- TARGET_COL_FALLBACK
}
cat("Target column:", TARGET_COL, "->", hdr_clean[TARGET_COL], "\n")
if (!GROUP_COL %in% hdr_clean) stop("Group column not found.")

ids     <- as.character(raw[[1]])
cp_raw  <- suppressWarnings(as.numeric(raw[[TARGET_COL]]))
groups  <- as.character(raw[[GROUP_COL]])
ok      <- is.finite(cp_raw) & !is.na(groups) & groups != ""
ids <- ids[ok]; cp_raw <- cp_raw[ok]; groups <- groups[ok]
names(cp_raw) <- ids; names(groups) <- ids

if (!is.null(KEEP_DATASETS)) {
  keep <- groups %in% KEEP_DATASETS
  cp_raw <- cp_raw[keep]; groups <- groups[keep]
  cat("\nSubset:", paste(KEEP_DATASETS, collapse = ", "), "-> n =", length(cp_raw), "\n")
}
cat("\nTarget summary (g/cm^2):\n"); print(summary(cp_raw))

# ============================================================
# 8. PRE-BIN ALL SPECTRA
# ============================================================
cat("\nPre-binning spectra...\n")
binned <- list()

for (v in spec_vars) {
  X <- get_spc(get(v), tag = v); storage.mode(X) <- "double"
  X <- clean_matrix(X, tag = v)
  if (nrow(X) < 30L) { cat("  Skipping", v, "- too few rows\n"); next }
  common <- intersect(rownames(X), names(cp_raw))
  if (length(common) < 30L) { cat("  Skipping", v, "- too few overlaps\n"); next }
  common <- names(cp_raw)[names(cp_raw) %in% common]
  X <- X[common, , drop = FALSE]
  binned[[v]] <- list()
  
  for (bs in BIN_SIZES) {
    Xb <- bin_matrix(X, bs)
    if (ncol(Xb) > 1L) {
      csd <- apply(Xb, 2, sd); col_bad <- !is.finite(csd) | csd == 0
      if (any(col_bad)) Xb <- Xb[, !col_bad, drop = FALSE]
    }
    rownames(Xb) <- common
    colnames(Xb) <- paste0("V", seq_len(ncol(Xb)))
    if (ncol(Xb) < 1L) next
    binned[[v]][[as.character(bs)]] <- Xb
  }
  n_cols <- vapply(BIN_SIZES, function(bs) {
    Xb <- binned[[v]][[as.character(bs)]]
    if (is.null(Xb)) return(NA_integer_); ncol(Xb)
  }, integer(1))
  cat(sprintf("  %-24s : %d samples; bands [%s]\n", v, nrow(X), paste(n_cols, collapse = "/")))
}
if (length(binned) == 0L) stop("No usable variants.")

# ============================================================
# 9. COMMON SAMPLE SET
# ============================================================
valid_variants <- names(binned)[vapply(binned,
                                       function(x) all(as.character(BIN_SIZES) %in% names(x)), logical(1))]
if (!length(valid_variants)) stop("No variant has all bins.")
binned <- binned[valid_variants]

candidate_keys <- unlist(lapply(valid_variants, function(v)
  paste0(v, "__bin", BIN_SIZES)), use.names = FALSE)

X_candidates <- list()
for (v in valid_variants) for (bs in BIN_SIZES)
  X_candidates[[paste0(v, "__bin", bs)]] <- binned[[v]][[as.character(bs)]]

common_ids <- Reduce(intersect, lapply(X_candidates, rownames))
common_ids <- names(cp_raw)[names(cp_raw) %in% common_ids]
if (length(common_ids) < max(30L, N_OUTER * 5L))
  stop("Too few common observations: n = ", length(common_ids))

X_candidates <- lapply(X_candidates, function(X) X[common_ids, , drop = FALSE])
y <- as.numeric(cp_raw[common_ids]); names(y) <- common_ids
g <- as.character(groups[common_ids]); names(g) <- common_ids

candidate_info <- do.call(rbind, lapply(names(X_candidates), function(key) {
  bm <- regmatches(key, regexpr("__bin[0-9]+$", key))
  data.frame(candidate = key,
             Preproc   = sub("__bin[0-9]+$", "", key),
             Bin       = as.integer(sub("__bin", "", bm)),
             stringsAsFactors = FALSE)
}))
rownames(candidate_info) <- NULL

cat("\n--- Common analysis dataset ---\n")
cat("Samples:", length(y), "\n")
cat("Candidates:", length(X_candidates), "\n")
cat("Target summary:\n"); print(summary(y))

# ============================================================
# 10. OUTER NESTED CV — evaluate ALL candidates on every fold
# ============================================================
cat("\nStarting fully nested SVR CV (all candidates per fold)...\n")
cat("Outer folds:", N_OUTER, "| Inner folds:", N_INNER,
    "| Candidates:", length(X_candidates), "\n")

outer_folds <- make_folds(g, y, N_OUTER, SEED, USE_GROUPED_CV)
k_outer_actual <- length(unique(outer_folds))
if (k_outer_actual < 2L) stop("Outer CV produced < 2 folds.")

N_WORKERS <- min(N_CORES, k_outer_actual)
cl <- makeCluster(N_WORKERS, type = "PSOCK")
registerDoParallel(cl)
clusterEvalQ(cl, { suppressPackageStartupMessages(library(e1071)); NULL })

tic("Fully nested SVR CV")

outer_results <- foreach(
  f = sort(unique(outer_folds)),
  .packages = "e1071",
  .export = c("X_candidates", "candidate_info", "y", "g", "outer_folds",
              "make_folds", "fit_svr", "select_svr_config",
              "pick_best_per_candidate",
              "SVR_COST_GRID", "SVR_GAMMA_GRID", "SVR_EPSILON",
              "SVR_MAX_TRAIN", "N_INNER", "SEED", "USE_GROUPED_CV"),
  .errorhandling = "stop"
) %dopar% {
  
  te <- which(outer_folds == f); tr <- which(outer_folds != f)
  if (length(tr) < 30L || length(te) < 5L)
    stop("Insufficient samples in outer fold ", f)
  
  X_train_list <- lapply(X_candidates, function(X) X[tr, , drop = FALSE])
  y_train <- y[tr]; g_train <- g[tr]; y_te <- y[te]
  
  # --- Step 1: inner CV over all candidates × (cost, gamma) ---
  sel <- select_svr_config(
    X_list = X_train_list, y_train = y_train, groups_train = g_train,
    n_inner = N_INNER, seed = SEED + 100L * f,
    grouped = USE_GROUPED_CV,
    cost_grid = SVR_COST_GRID, gamma_grid = SVR_GAMMA_GRID,
    epsilon = SVR_EPSILON, max_train = SVR_MAX_TRAIN)
  
  # --- Step 2: per-candidate best (cost, gamma) by min inner RMSEP ---
  per_cand <- pick_best_per_candidate(sel$audit, use_rmsep = TRUE)
  
  # --- Step 3: refit each candidate, predict outer-test ---
  fold_metrics_list <- list(); fold_preds_list <- list()
  
  for (i in seq_len(nrow(per_cand))) {
    cand <- per_cand$candidate[i]
    cost_i <- per_cand$cost[i]
    gamma_i <- per_cand$gamma[i]
    inner_rmse_i <- per_cand$inner_rmse[i]
    
    info <- candidate_info[match(cand, candidate_info$candidate), , drop = FALSE]
    
    X_sel <- X_candidates[[cand]]
    X_tr <- X_sel[tr, , drop = FALSE]
    X_te <- X_sel[te, , drop = FALSE]
    y_tr <- y[tr]
    
    # Apply training cap consistently (uses same seed as original winner refit)
    if (is.finite(SVR_MAX_TRAIN) && length(y_tr) > SVR_MAX_TRAIN) {
      set.seed(SEED + 1000L + f)
      take <- sort(sample.int(length(y_tr), SVR_MAX_TRAIN))
      X_tr <- X_tr[take, , drop = FALSE]
      y_tr <- y_tr[take]
    }
    
    final_fit <- fit_svr(X_tr, y_tr, X_te,
                         cost = cost_i, gamma = gamma_i, epsilon = SVR_EPSILON)
    pred <- final_fit$pred
    if (is.null(final_fit$fit) || length(pred) != length(y_te) ||
        any(!is.finite(pred))) next
    
    rmsep  <- sqrt(mean((y_te - pred)^2))
    mae    <- mean(abs(y_te - pred))
    bias   <- mean(pred - y_te)
    r2_cor <- if (sd(y_te) > 0 && sd(pred) > 0) cor(pred, y_te)^2 else NA_real_
    denom  <- sum((y_te - mean(y_te))^2)
    r2_ss  <- if (denom > 0) 1 - sum((y_te - pred)^2)/denom else NA_real_
    rpd    <- if (is.finite(rmsep) && rmsep > 0) sd(y_te)/rmsep else NA_real_
    
    fold_metrics_list[[length(fold_metrics_list) + 1L]] <- data.frame(
      Bin = info$Bin, Preproc = info$Preproc, OuterFold = f,
      kernel = "radial", cost = cost_i, gamma = gamma_i,
      epsilon = SVR_EPSILON, inner_rmse = inner_rmse_i,
      n_train = nrow(X_tr), n_test = length(y_te),
      R2 = r2_cor, R2_ss = r2_ss, RMSEP = rmsep,
      MAE = mae, bias = bias, RPD = rpd,
      stringsAsFactors = FALSE)
    
    fold_preds_list[[length(fold_preds_list) + 1L]] <- data.frame(
      Bin = info$Bin, Preproc = info$Preproc, OuterFold = f,
      sample_id = common_ids[te], group = g[te],
      y_true = y_te, y_pred = pred,
      kernel = "radial", cost = cost_i, gamma = gamma_i,
      stringsAsFactors = FALSE)
  }
  
  list(metrics = do.call(rbind, fold_metrics_list),
       predictions = do.call(rbind, fold_preds_list))
}

toc()
stopCluster(cl); registerDoSEQ()

# ============================================================
# 11. AGGREGATE OUTER-CV RESULTS
# ============================================================
fold_metrics <- do.call(rbind, lapply(outer_results, function(x) x$metrics))
fold_preds   <- do.call(rbind, lapply(outer_results, function(x) x$predictions))

# Sort: Bin, Preproc, OuterFold (matches attachment)
fold_metrics <- fold_metrics[order(fold_metrics$Bin, fold_metrics$Preproc,
                                   fold_metrics$OuterFold), ]
fold_preds   <- fold_preds[order(fold_preds$Bin, fold_preds$Preproc,
                                 fold_preds$OuterFold), ]
rownames(fold_metrics) <- NULL; rownames(fold_preds) <- NULL

cat("\nTotal fold-metrics rows   :", nrow(fold_metrics), "(expected 200)\n")
cat("Total fold-prediction rows:", nrow(fold_preds),   "(expected ~12200)\n")

# ============================================================
# 12. Per-fold nested winner (MIN inner_rmse across candidates)
# ============================================================
selection_history <- do.call(rbind, lapply(
  split(fold_metrics, fold_metrics$OuterFold),
  function(df) df[which.min(df$inner_rmse), , drop = FALSE]))
rownames(selection_history) <- NULL

cat("\n--- Nested winner per outer fold (min inner_rmse) ---\n")
print(selection_history[, c("OuterFold", "Preproc", "Bin", "cost", "gamma",
                            "inner_rmse", "RMSEP", "RPD")])

# ============================================================
# 13. Aggregate per candidate (All Models Metrics)
# ============================================================
split_key <- paste(fold_metrics$Preproc, fold_metrics$Bin, sep = "||")

aggregate_per_candidate <- do.call(rbind, lapply(split(fold_metrics, split_key),
                                                 function(df) {
                                                   nf <- length(unique(df$OuterFold))
                                                   data.frame(
                                                     Preproc    = df$Preproc[1],
                                                     Bin_nm     = as.integer(df$Bin[1]),
                                                     n          = length(common_ids),
                                                     n_outer_folds = nf,
                                                     cost       = median(df$cost,  na.rm = TRUE),
                                                     gamma      = median(df$gamma, na.rm = TRUE),
                                                     R2_mean    = mean(df$R2,      na.rm = TRUE),
                                                     R2_sd      = sd  (df$R2,      na.rm = TRUE),
                                                     R2_SS_mean = mean(df$R2_ss,   na.rm = TRUE),
                                                     R2_SS_sd   = sd  (df$R2_ss,   na.rm = TRUE),
                                                     RMSEP_mean = mean(df$RMSEP,   na.rm = TRUE),
                                                     RMSEP_sd   = sd  (df$RMSEP,   na.rm = TRUE),
                                                     MAE_mean   = mean(df$MAE,     na.rm = TRUE),
                                                     MAE_sd     = sd  (df$MAE,     na.rm = TRUE),
                                                     bias_mean  = mean(df$bias,    na.rm = TRUE),
                                                     bias_sd    = sd  (df$bias,    na.rm = TRUE),
                                                     RPD_mean   = mean(df$RPD,     na.rm = TRUE),
                                                     RPD_sd     = sd  (df$RPD,     na.rm = TRUE),
                                                     RPD_margin_95 = 1.96 * sd(df$RPD, na.rm = TRUE) / sqrt(nf),
                                                     stringsAsFactors = FALSE)
                                                 }))

aggregate_per_candidate <- aggregate_per_candidate[
  order(-aggregate_per_candidate$R2_mean), ]
aggregate_per_candidate$Rank <- seq_len(nrow(aggregate_per_candidate))
aggregate_per_candidate <- aggregate_per_candidate[, c(
  "Rank", "Preproc", "Bin_nm", "n", "n_outer_folds", "cost", "gamma",
  "R2_mean", "R2_sd", "R2_SS_mean", "R2_SS_sd",
  "RMSEP_mean", "RMSEP_sd", "MAE_mean", "MAE_sd",
  "bias_mean", "bias_sd", "RPD_mean", "RPD_sd", "RPD_margin_95")]
rownames(aggregate_per_candidate) <- NULL

cat("\n--- Aggregate outer-CV metrics per candidate (top 10) ---\n")
print(head(aggregate_per_candidate, 10))

# ============================================================
# 14. Winner-config aggregate (Panel D) — uses t-margin, not 1.96
# ============================================================
win_rows <- fold_metrics[match(
  paste(selection_history$Preproc, selection_history$Bin, selection_history$OuterFold),
  paste(fold_metrics$Preproc, fold_metrics$Bin, fold_metrics$OuterFold)), ]
n_folds_win <- nrow(win_rows)

t_crit <- qt(0.975, df = n_folds_win - 1)

overall_summary <- data.frame(
  Metric = c("R2", "R2_SS", "RPD", "RPD_margin_95", "RMSEP", "MAE", "bias"),
  Mean_outer_fold = c(mean(win_rows$R2,    na.rm = TRUE),
                      mean(win_rows$R2_ss, na.rm = TRUE),
                      mean(win_rows$RPD,   na.rm = TRUE),
                      NA_real_,
                      mean(win_rows$RMSEP, na.rm = TRUE),
                      mean(win_rows$MAE,   na.rm = TRUE),
                      mean(win_rows$bias,  na.rm = TRUE)),
  SD_outer_fold = c(sd(win_rows$R2,    na.rm = TRUE),
                    sd(win_rows$R2_ss, na.rm = TRUE),
                    sd(win_rows$RPD,   na.rm = TRUE),
                    NA_real_,
                    sd(win_rows$RMSEP, na.rm = TRUE),
                    sd(win_rows$MAE,   na.rm = TRUE),
                    sd(win_rows$bias,  na.rm = TRUE)),
  Value = c(NA_real_, NA_real_, NA_real_,
            t_crit * sd(win_rows$RPD, na.rm = TRUE) / sqrt(n_folds_win),
            NA_real_, NA_real_, NA_real_),
  stringsAsFactors = FALSE)

cat("\n--- Winner-config aggregate (Panel D, t-margin) ---\n")
print(overall_summary)

# ============================================================
# 15. FINAL DEPLOYMENT MODEL
# ============================================================
cat("\nSelecting deployment config on full data...\n")
full_selection <- select_svr_config(
  X_list = X_candidates, y_train = y, groups_train = g,
  n_inner = N_INNER, seed = SEED + 10000L, grouped = USE_GROUPED_CV,
  cost_grid = SVR_COST_GRID, gamma_grid = SVR_GAMMA_GRID,
  epsilon = SVR_EPSILON, max_train = SVR_MAX_TRAIN)

WINNER_KEY <- full_selection$best_key
winner_info <- candidate_info[match(WINNER_KEY, candidate_info$candidate), , drop = FALSE]
WINNER_PREPROC <- as.character(winner_info$Preproc)
WINNER_BIN     <- as.integer(winner_info$Bin)
WINNER_COST    <- full_selection$best_cost
WINNER_GAMMA   <- full_selection$best_gamma

cat("\nDeployment selection:\n")
cat("  Preproc:", WINNER_PREPROC, "| Bin:", WINNER_BIN,
    "| cost:", WINNER_COST, "| gamma:", WINNER_GAMMA, "\n")

X_full <- X_candidates[[WINNER_KEY]]; y_full <- y
if (is.finite(SVR_MAX_TRAIN) && length(y_full) > SVR_MAX_TRAIN) {
  set.seed(SEED + 20000L)
  take <- sort(sample.int(length(y_full), SVR_MAX_TRAIN))
  X_deploy <- X_full[take, , drop = FALSE]; y_deploy <- y_full[take]
} else { X_deploy <- X_full; y_deploy <- y_full }

fit_deploy <- e1071::svm(x = as.matrix(X_deploy), y = as.numeric(y_deploy),
                         type = "eps-regression", kernel = "radial",
                         cost = WINNER_COST, gamma = WINNER_GAMMA,
                         epsilon = SVR_EPSILON, scale = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
tag <- paste(KEEP_DATASETS, collapse = "_")
filter_tag <- paste0("spc", if (USE_SPECTRAL_FILTER) "on" else "off",
                     "_y", if (USE_Y_BOXPLOT) "on" else "off")

model_name <- sprintf("SVR_cp_ds%s_bin%d_%s_%s",
                      tag, WINNER_BIN, gsub("^spec_", "", WINNER_PREPROC), timestamp)
model_file <- file.path(BEST_MODELS_SVR_DIR, paste0(model_name, ".rds"))
saveRDS(fit_deploy, model_file)

meta <- list(
  algorithm = "e1071::svm, epsilon regression, radial kernel",
  target = TARGET_NAME, target_units = "g/cm^2", target_column = TARGET_COL,
  subset = KEEP_DATASETS,
  selection_criterion = "Minimum mean fold-level inner-CV RMSEP",
  cv_outer_folds = N_OUTER, cv_inner_folds = N_INNER,
  seed = SEED, grouped_cv = USE_GROUPED_CV,
  preprocessing = WINNER_PREPROC, bin_size = WINNER_BIN,
  kernel = SVR_KERNEL, cost = WINNER_COST, gamma = WINNER_GAMMA,
  epsilon = SVR_EPSILON,
  inner_mean_RMSEP = full_selection$best_mean_RMSEP,
  inner_mean_RPD   = full_selection$best_mean_RPD,
  n_samples_available = length(y_full),
  n_samples_model_fit = length(y_deploy),
  n_features = ncol(X_deploy),
  feature_names = colnames(X_deploy),
  train_ids = rownames(X_deploy),
  y_range = range(y_deploy, na.rm = TRUE),
  n_support_vectors = length(fit_deploy$index),
  timestamp = timestamp)
saveRDS(meta, file.path(BEST_MODELS_SVR_DIR, paste0(model_name, "_meta.rds")))

# ============================================================
# 16. WRITE OUTPUT FILES
# ============================================================
out_prefix <- file.path(SUMMARY_DIR,
                        paste0("SVR_cp_ds", tag, "_", filter_tag, "_",
                               if (USE_GROUPED_CV) "Grouped" else "Random",
                               "NestedCV_", timestamp))

write.csv(fold_metrics, paste0(out_prefix, "_foldmetrics.csv"), row.names = FALSE)
write.csv(fold_preds,   paste0(out_prefix, "_foldpreds.csv"),   row.names = FALSE)

all_models_path <- file.path(SUMMARY_DIR,
                             "Supplementary Material 4 - SVR All Models Metrics.csv")
write.csv(aggregate_per_candidate, all_models_path, row.names = FALSE, na = "")

# ============================================================
# 17. PANELS A-E
# ============================================================
PANEL_COLS <- c("Model", "Panel", "RowType", "OuterFold", "Preproc", "Bin_nm",
                "cost", "gamma", "inner_RMSEP", "n_train", "n_test",
                "R2", "R2_SS", "RMSEP", "MAE", "bias", "RPD",
                "Metric", "Mean", "SD", "Value", "Note")

# --- Panel A: Search space (9 rows) ---
panel_A <- data.frame(
  Model   = "SVR", Panel = "A", RowType = "SearchSpace",
  Metric  = c("Preprocessing_variants", "Bin_widths_nm",
              "cost_values", "gamma_values",
              "epsilon", "kernel", "scale",
              "Outer_CV_folds", "Inner_CV_folds"),
  Mean    = c(length(spec_vars), NA, NA, NA,
              SVR_EPSILON, SVR_KERNEL, "TRUE",
              N_OUTER, N_INNER),
  Note    = c(paste(length(spec_vars), "spectral preprocessing variants"),
              paste(BIN_SIZES,       collapse = ";"),
              SVR_COST_GRID_STR,
              SVR_GAMMA_GRID_STR,
              NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE)

# --- Panel B: Nested selected config per outer fold ---
panel_B <- data.frame(
  Model       = "SVR", Panel = "B", RowType = "NestedSelected",
  OuterFold   = selection_history$OuterFold,
  Preproc     = selection_history$Preproc,
  Bin_nm      = selection_history$Bin,
  cost        = selection_history$cost,
  gamma       = selection_history$gamma,
  inner_RMSEP = selection_history$inner_rmse,
  Note        = paste0("Selected by minimum inner RMSEP in outer fold ",
                       selection_history$OuterFold),
  stringsAsFactors = FALSE)

# --- Panel C: Outer metrics of per-fold winner ---
panel_C <- data.frame(
  Model     = "SVR", Panel = "C", RowType = "OuterMetrics",
  OuterFold = selection_history$OuterFold,
  Preproc   = selection_history$Preproc,
  Bin_nm    = selection_history$Bin,
  cost      = selection_history$cost,
  gamma     = selection_history$gamma,
  n_train   = selection_history$n_train,
  n_test    = selection_history$n_test,
  R2        = selection_history$R2,
  R2_SS     = selection_history$R2_ss,
  RMSEP     = selection_history$RMSEP,
  MAE       = selection_history$MAE,
  bias      = selection_history$bias,
  RPD       = selection_history$RPD,
  stringsAsFactors = FALSE)

# --- Panel D: Aggregate (winner config, t-margin) ---
panel_D <- data.frame(
  Model   = "SVR", Panel = "D", RowType = "Aggregate",
  Metric  = overall_summary$Metric,
  Mean    = overall_summary$Mean_outer_fold,
  SD      = overall_summary$SD_outer_fold,
  Value   = overall_summary$Value,
  Note    = c("Mean across five nested outer folds",
              "Mean across five nested outer folds",
              "Mean and SD across five nested outer folds",
              "Approximate 95% margin across folds",
              "Mean and SD across five nested outer folds",
              "Mean and SD across five nested outer folds",
              "Mean and SD across five nested outer folds"),
  stringsAsFactors = FALSE)

# --- Panel E: Selection frequency ---
selection_history$config_name <- paste0(
  selection_history$Preproc, "_bin", selection_history$Bin,
  "_cost", selection_history$cost, "_gamma", selection_history$gamma)

config_folds <- split(selection_history$OuterFold, selection_history$config_name)
config_folds <- config_folds[order(sapply(config_folds, min))]

combine_folds <- function(folds) {
  folds <- sort(unique(folds)); n <- length(folds)
  if (n == 1L) return(as.character(folds))
  if (n == 2L) return(paste(folds, collapse = " and "))
  paste0(paste(folds[-n], collapse = ", "), ", and ", folds[n])
}

panel_E <- data.frame(
  Model   = "SVR", Panel = "E", RowType = "Frequency",
  Metric  = names(config_folds),
  Value   = vapply(config_folds, length, integer(1)),
  Note    = vapply(config_folds, function(x)
    paste0("Selected in fold",
           if (length(x) > 1L) "s" else "", " ",
           combine_folds(x)),
    character(1)),
  stringsAsFactors = FALSE)

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

panel_path <- file.path(SUMMARY_DIR, "Supplementary Material 4 - SVR.csv")
write.csv(panels_all, panel_path, row.names = FALSE, na = "")

# ============================================================
# 18. CONSOLE SUMMARY
# ============================================================
cat("\n=== Saved outputs ===\n")
cat("  ", panel_path,                            "\n", sep = "")
cat("  ", all_models_path,                        "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldmetrics.csv"), "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldpreds.csv"),   "\n", sep = "")
cat("  ", model_file,                             "\n", sep = "")
cat("  ", file.path(BEST_MODELS_SVR_DIR,
                    paste0(model_name, "_meta.rds")), "\n", sep = "")

cat("\nDONE.\n")