#********************************************************************#
# A New Approach to Predict Nutrient Content in Plants               #
#        Using VNIR-SWIR Spectroscopy and Machine-Learning           #
#              Jhohan Elguera - MATE University                      #
#               With contributions of Dr. Johan PERRET               #
#                        PART 2 PLSR MODELING                        #
#         FULLY-NESTED CV — ALL CANDIDATES ON EVERY OUTER FOLD       #
#                                    2026                            #
#********************************************************************#

#************************************************************************#
# Nested-CV PLSR — Nitrogen prediction (target = cp)                     #
#                                                                        #
# SUBSET VERSION: dataset 13 only                                        #
#   - scans preprocessing variants x bin sizes                           #
#   - saves the winning ranger model to 4-Best_Models/PLSR/              #
#                                                                        #
#************************************************************************#

suppressPackageStartupMessages({
  library(pls)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(tictoc)
  library(robustbase)
  library(dplyr)
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
BEST_MODELS_PLSR_DIR   <- file.path(BEST_MODELS_DIR, "PLSR")

BIN_SIZES              <- c(5, 25, 50, 75)

N_OUTER                <- 5
N_INNER                <- 5
MAX_NCOMP              <- 25
SEED                   <- 42
N_CORES                <- max(1, detectCores() - 1)

SPC_FILTER_N_PC        <- 4
SPC_FILTER_ALPHA       <- 0.975

dir.create(SUMMARY_DIR,           showWarnings = FALSE, recursive = TRUE)
dir.create(BEST_MODELS_DIR,       showWarnings = FALSE, recursive = TRUE)
dir.create(BEST_MODELS_PLSR_DIR,  showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load latest preprocessing workspace
# ============================================================
snap <- fileSnapshot(PREPROC_DIR)
latest_dir <- rownames(snap$info[which.max(snap$info$mtime),])
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

extract_spectra <- function(obj, tag, target_ids) {
  data_part <- NULL
  if (isS4(obj)) {
    slots <- methods::slotNames(obj)
    if ("data" %in% slots) data_part <- methods::slot(obj, "data")
  } else if (is.list(obj) || is.data.frame(obj)) {
    data_part <- obj
  }
  if (is.null(data_part)) stop("Cannot extract spectral data from ", tag)
  if (is.data.frame(data_part) && "spc" %in% names(data_part)) {
    X <- as.matrix(data_part[["spc"]])
  } else if (is.list(data_part) && !is.null(data_part[["spc"]])) {
    X <- as.matrix(data_part[["spc"]])
  } else if (is.matrix(obj) || is.data.frame(obj)) {
    X <- as.matrix(obj)
  } else stop("No 'spc' spectral field found in ", tag)
  storage.mode(X) <- "double"
  if (!nrow(X) || !ncol(X)) stop("Empty spectral matrix in ", tag)
  
  id_candidates <- list()
  if (is.data.frame(data_part) || is.list(data_part)) {
    metadata_names <- setdiff(names(data_part), "spc")
    id_like <- grep("(^id$|sample|specimen|accession|record)", metadata_names,
                    ignore.case = TRUE, value = TRUE)
    metadata_names <- c(id_like, setdiff(metadata_names, id_like))
    for (nm in metadata_names) {
      value <- data_part[[nm]]
      if (is.atomic(value) && length(value) == nrow(X))
        id_candidates[[paste0("metadata column '", nm, "'")]] <- as.character(value)
    }
  }
  if (!is.null(rownames(X))) id_candidates[["spectral row names"]] <- rownames(X)
  if (!is.null(rownames(data_part)) && length(rownames(data_part)) == nrow(X))
    id_candidates[["data row names"]] <- rownames(data_part)
  if (!length(id_candidates)) stop("No sample-ID candidates found for ", tag)
  
  scores <- vapply(id_candidates, function(z) {
    z <- trimws(as.character(z))
    sum(!is.na(z) & z != "" & z %in% target_ids)
  }, numeric(1))
  best_name <- names(which.max(scores))[1]
  best_ids  <- trimws(as.character(id_candidates[[best_name]]))
  if (max(scores) < 1) stop("Could not match sample IDs for ", tag)
  rownames(X) <- best_ids
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
  denom_n   <- max(target_n^2,  1e-9); denom_sy  <- max(target_sy^2, 1e-9)
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

compute_outlier_filter <- function(X_train,
                                   n_pc = SPC_FILTER_N_PC,
                                   alpha = SPC_FILTER_ALPHA) {
  if (nrow(X_train) < 20 || ncol(X_train) < 2) return(NULL)
  pca <- tryCatch(prcomp(X_train, center = TRUE, scale = FALSE),
                  error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  n_pc <- min(n_pc, ncol(pca$x)); if (n_pc < 2) return(NULL)
  pairs <- combn(n_pc, 2)
  filters <- vector("list", ncol(pairs))
  for (j in seq_len(ncol(pairs))) {
    p1 <- pairs[1, j]; p2 <- pairs[2, j]
    sub <- pca$x[, c(p1, p2), drop = FALSE]
    mcd <- tryCatch(covMcd(sub), error = function(e) NULL)
    if (is.null(mcd)) { filters[[j]] <- NULL; next }
    filters[[j]] <- list(pair = c(p1, p2), center = mcd$center, cov = mcd$cov)
  }
  filters <- Filter(Negate(is.null), filters)
  if (length(filters) == 0) return(NULL)
  list(pca = pca, n_pc = n_pc, filters = filters, cutoff = qchisq(alpha, df = 2))
}

apply_outlier_filter <- function(filt, X) {
  if (is.null(filt)) return(rep(TRUE, nrow(X)))
  sc <- predict(filt$pca, newdata = X)[, 1:filt$n_pc, drop = FALSE]
  out <- logical(nrow(sc))
  for (f in filt$filters) {
    d <- mahalanobis(sc[, f$pair, drop = FALSE], f$center, f$cov)
    d[!is.finite(d)] <- Inf
    out <- out | (d > filt$cutoff)
  }
  !out
}

# ============================================================
# 3b. ID overlap diagnostic
# ============================================================
spec_names <- rownames(extract_spectra(get(spec_vars[1]), spec_vars[1], names(cp_raw)))
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
# 5. Pre-bin all spectra matrices
# ============================================================
cat("\nPre-binning spectra...\n")
binned <- list()

for (v in spec_vars) {
  hs <- get(v)
  X <- tryCatch(
    extract_spectra(hs, tag = v, target_ids = names(cp_raw)),
    error = function(e) { message("Skipping ", v, ": ", conditionMessage(e)); NULL })
  if (is.null(X)) next
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
# 6. FULLY NESTED CV — evaluate ALL candidates on EVERY outer fold
# ============================================================
cat("\nStarting fully nested PLSR CV (all candidates per fold)...\n")
cat("Outer folds:", N_OUTER, "| Inner folds:", N_INNER,
    "| Maximum components:", MAX_NCOMP, "\n")

# ------------------------------------------------------------
# 6.1 Build the common sample set across every candidate
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
  candidate_matrices[[candidate_keys$key[i]]] <-
    binned[[v]][[as.character(bs)]]
}
if (length(candidate_matrices) == 0L)
  stop("No candidate spectral matrices are available.")

common_ids <- Reduce(intersect, lapply(candidate_matrices, rownames))
common_ids <- common_ids[common_ids %in% names(cp_raw) &
                           common_ids %in% names(groups)]
common_ids <- sort(common_ids)
if (length(common_ids) < 30L)
  stop("Fewer than 30 common observations across all candidates.")

candidate_matrices <- lapply(candidate_matrices,
                             function(X) X[common_ids, , drop = FALSE])
y_all <- unname(cp_raw[common_ids]); names(y_all) <- common_ids
g_all <- unname(groups[common_ids]); names(g_all) <- common_ids

cat("Common observations:", length(common_ids),
    "| Candidate configurations:", length(candidate_matrices), "\n")

# ------------------------------------------------------------
# 6.2 Metrics helper
# ------------------------------------------------------------
calc_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  if (length(obs) < 2L)
    return(c(R2=NA_real_, R2_ss=NA_real_, RMSEP=NA_real_,
             MAE=NA_real_, bias=NA_real_, RPD=NA_real_))
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
# 6.3 Evaluate ALL candidates on ONE outer fold
#     • inner CV per candidate → per-candidate best ncomp
#       (selection = MIN mean inner RMSEP)
#     • refit each candidate on outer-train, score on outer-test
# ------------------------------------------------------------
evaluate_all_candidates_on_outer <- function(
    X_train_list, X_test_list, y_train, y_test,
    groups_train, n_inner, max_ncomp, seed, grouped) {
  
  inner_fold <- make_folds(groups_train, y_train, n_inner, seed, grouped)
  k_inner    <- length(unique(inner_fold))
  res_list   <- list()
  
  for (cand in names(X_train_list)) {
    X_tr <- X_train_list[[cand]]
    X_te <- X_test_list[[cand]]
    if (is.null(X_tr) || is.null(X_te)) next
    
    max_comp <- min(max_ncomp, ncol(X_tr), nrow(X_tr) - 1L)
    if (max_comp < 1L) next
    
    # --- inner CV: observation × component prediction matrix ---
    pred_by_comp <- matrix(NA_real_, nrow = nrow(X_tr), ncol = max_comp)
    for (f in seq_len(k_inner)) {
      idx_tr <- which(inner_fold != f)
      idx_va <- which(inner_fold == f)
      if (length(idx_tr) < 10L || length(idx_va) < 2L) next
      
      X_tr_f <- X_tr[idx_tr, , drop = FALSE]
      X_va_f <- X_tr[idx_va, , drop = FALSE]
      mc <- min(max_comp, ncol(X_tr_f), nrow(X_tr_f) - 1L)
      if (mc < 1L) next
      
      cal_df <- data.frame(y = y_train[idx_tr], NIR = I(X_tr_f))
      fit <- tryCatch(
        pls::plsr(y ~ NIR, ncomp = mc, data = cal_df,
                  scale = FALSE, method = "oscorespls"),
        error = function(e) NULL)
      if (is.null(fit)) next
      
      pr <- tryCatch(
        predict(fit,
                newdata = data.frame(NIR = I(X_va_f)),
                ncomp = seq_len(mc)),
        error = function(e) NULL)
      if (is.null(pr)) next
      pr <- as.array(pr)
      
      for (nc in seq_len(mc)) {
        p <- tryCatch(as.numeric(pr[, 1L, nc]),
                      error = function(e) rep(NA_real_, length(idx_va)))
        if (length(p) == length(idx_va)) pred_by_comp[idx_va, nc] <- p
      }
    }
    
    # --- per-comp mean inner RMSEP across inner folds ---
    rmsep_by_comp <- rep(NA_real_, ncol(pred_by_comp))
    for (nc in seq_len(ncol(pred_by_comp))) {
      fold_rmsep <- rep(NA_real_, k_inner)
      for (f in seq_len(k_inner)) {
        idx_va <- which(inner_fold == f)
        p   <- pred_by_comp[idx_va, nc]
        obs <- y_train[idx_va]
        ok  <- is.finite(p) & is.finite(obs)
        if (sum(ok) < 2L) next
        fold_rmsep[f] <- sqrt(mean((obs[ok] - p[ok])^2))
      }
      rmsep_by_comp[nc] <- mean(fold_rmsep, na.rm = TRUE)
    }
    if (all(!is.finite(rmsep_by_comp))) next
    
    best_ncomp  <- which.min(rmsep_by_comp)
    inner_rmsep <- rmsep_by_comp[best_ncomp]
    
    # --- refit on outer-train, predict outer-test ---
    cal_df_fin <- data.frame(y = y_train, NIR = I(X_tr))
    fit_outer <- tryCatch(
      pls::plsr(y ~ NIR, ncomp = best_ncomp, data = cal_df_fin,
                scale = FALSE, method = "oscorespls"),
      error = function(e) NULL)
    if (is.null(fit_outer)) next
    
    pred_outer <- tryCatch(
      as.numeric(predict(fit_outer,
                         newdata = data.frame(NIR = I(X_te)),
                         ncomp = best_ncomp)[, 1, 1]),
      error = function(e) rep(NA_real_, nrow(X_te)))
    if (any(!is.finite(pred_outer))) next
    
    m <- calc_metrics(y_test, pred_outer)
    
    res_list[[cand]] <- list(
      metric = data.frame(
        candidate   = cand,
        ncomp       = best_ncomp,
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
        y_pred    = pred_outer,
        ncomp     = best_ncomp,
        candidate = cand,
        stringsAsFactors = FALSE))
  }
  
  if (!length(res_list)) return(NULL)
  list(metrics     = do.call(rbind, lapply(res_list, `[[`, "metric")),
       predictions = do.call(rbind, lapply(res_list, `[[`, "pred")))
}

# ------------------------------------------------------------
# 6.4 Outer CV loop — evaluate ALL candidates in every fold
# ------------------------------------------------------------
outer_fold <- make_folds(groups = g_all, y = y_all, k = N_OUTER,
                         seed = SEED, grouped = USE_GROUPED_CV)
k_outer <- length(unique(outer_fold))

all_fold_metrics <- list()
all_fold_preds   <- list()

for (f in seq_len(k_outer)) {
  idx_tr <- which(outer_fold != f)
  idx_te <- which(outer_fold == f)
  if (length(idx_tr) < 30L || length(idx_te) < 5L) {
    warning("Skipping outer fold ", f, " - insufficient samples.")
    next
  }
  
  cat("\nOuter fold", f, "of", k_outer,
      " (n_tr =", length(idx_tr), "| n_te =", length(idx_te), ")\n")
  
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
      max_ncomp    = MAX_NCOMP,
      seed         = SEED + f,
      grouped      = USE_GROUPED_CV),
    error = function(e) {
      warning("Fold ", f, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) next
  
  # -- resolve Preproc / Bin for the metrics table --
  key_row <- match(res$metrics$candidate, candidate_keys$key)
  res$metrics$Preproc   <- candidate_keys$variant[key_row]
  res$metrics$Bin       <- as.integer(candidate_keys$bin_size[key_row])
  res$metrics$OuterFold <- f
  
  # -- resolve Preproc / Bin / group for the predictions table --
  pred_key_row <- match(res$predictions$candidate, candidate_keys$key)
  res$predictions$Preproc   <- candidate_keys$variant[pred_key_row]
  res$predictions$Bin       <- as.integer(candidate_keys$bin_size[pred_key_row])
  res$predictions$OuterFold <- f
  res$predictions$group     <- as.character(
    setNames(g_all[idx_te], common_ids[idx_te])[res$predictions$sample_id])
  
  # Column order matching the attached files
  all_fold_metrics[[f]] <- res$metrics[, c(
    "Bin","Preproc","OuterFold","ncomp","inner_rmsep",
    "n_train","n_test","R2","R2_ss","RMSEP","MAE","bias","RPD")]
  all_fold_preds[[f]] <- res$predictions[, c(
    "Bin","Preproc","OuterFold","sample_id","group",
    "y_true","y_pred","ncomp")]
  
  cat("  Candidates evaluated: ", nrow(res$metrics), "\n", sep = "")
}

if (!length(all_fold_metrics))
  stop("No outer folds produced any results.")

fold_metrics_all <- do.call(rbind, all_fold_metrics)
fold_preds_all   <- do.call(rbind, all_fold_preds)
rownames(fold_metrics_all) <- NULL
rownames(fold_preds_all)   <- NULL

cat("\nTotal fold-metrics rows   :", nrow(fold_metrics_all), "\n")
cat("Total fold-predictions rows:", nrow(fold_preds_all), "\n")

# ------------------------------------------------------------
# 6.5 Per-fold nested winner (Panel B / C / E of PLSR.csv)
#     Selection rule = MINIMUM inner_mean_RMSEP within each fold
# ------------------------------------------------------------
selection_history <- do.call(rbind, lapply(
  split(fold_metrics_all, fold_metrics_all$OuterFold),
  function(df) df[which.min(df$inner_rmsep), , drop = FALSE]))
rownames(selection_history) <- NULL

# Per-fold winner's predictions (one row per held-out sample)
sel_keys <- paste(selection_history$OuterFold,
                  selection_history$Preproc,
                  selection_history$Bin, sep = "||")
pred_keys <- paste(fold_preds_all$OuterFold,
                   fold_preds_all$Preproc,
                   fold_preds_all$Bin, sep = "||")
fold_preds <- fold_preds_all[pred_keys %in% sel_keys, ]
rownames(fold_preds) <- NULL

cat("\n--- Nested winner per outer fold (min inner RMSEP) ---\n")
print(selection_history)

cat("\n--- Selection frequency ---\n")
freq_df <- as.data.frame(table(paste(selection_history$Preproc,
                                     selection_history$Bin, sep = "_bin")))
colnames(freq_df) <- c("Config", "Count")
print(freq_df)

# ------------------------------------------------------------
# 6.6 Aggregate metrics per candidate across outer folds
#     (All Models Metrics.csv, one row per candidate)
# ------------------------------------------------------------
split_key <- paste(fold_metrics_all$Preproc, fold_metrics_all$Bin, sep = "||")

aggregate_per_candidate <- do.call(rbind, lapply(split(fold_metrics_all, split_key), function(df) {
  nf <- length(unique(df$OuterFold))
  data.frame(
    Preproc       = df$Preproc[1],
    Bin_nm        = as.integer(df$Bin[1]),
    n             = length(common_ids),
    n_outer_folds = nf,
    ncomp_median  = median(df$ncomp, na.rm = TRUE),
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
  "Rank","Preproc","Bin_nm","n","n_outer_folds","ncomp_median",
  "R2_mean","R2_sd","R2_SS_mean","R2_SS_sd",
  "RMSEP_mean","RMSEP_sd","MAE_mean","MAE_sd",
  "bias_mean","bias_sd","RPD_mean","RPD_sd","RPD_margin_95")]
rownames(aggregate_per_candidate) <- NULL

cat("\n--- Aggregate outer-CV metrics per candidate (top 10) ---\n")
print(head(aggregate_per_candidate, 10))

# ------------------------------------------------------------
# 6.7 Winner-config aggregate (Panel-D equivalent)
# ------------------------------------------------------------
winner_cfg  <- selection_history[1, c("Preproc","Bin")]
winner_rows <- fold_metrics_all[
  fold_metrics_all$Preproc == winner_cfg$Preproc &
    fold_metrics_all$Bin   == winner_cfg$Bin, ]

overall_summary <- data.frame(
  Metric          = c("R2","R2_ss","RMSEP","MAE","bias","RPD"),
  Mean_outer_fold = c(mean(winner_rows$R2,    na.rm = TRUE),
                      mean(winner_rows$R2_ss, na.rm = TRUE),
                      mean(winner_rows$RMSEP, na.rm = TRUE),
                      mean(winner_rows$MAE,   na.rm = TRUE),
                      mean(winner_rows$bias,  na.rm = TRUE),
                      mean(winner_rows$RPD,   na.rm = TRUE)),
  SD_outer_fold   = c(sd(winner_rows$R2,    na.rm = TRUE),
                      sd(winner_rows$R2_ss, na.rm = TRUE),
                      sd(winner_rows$RMSEP, na.rm = TRUE),
                      sd(winner_rows$MAE,   na.rm = TRUE),
                      sd(winner_rows$bias,  na.rm = TRUE),
                      sd(winner_rows$RPD,   na.rm = TRUE)))
cat("\n--- Winner-config aggregate (Panel D) ---\n"); print(overall_summary)

# ============================================================
# 7. Final deployment model — full-data inner CV
# ============================================================
cat("\nSelecting final deployment configuration on full data...\n")

select_full <- function(X_list, y, groups, n_inner, max_ncomp, seed, grouped) {
  inner_fold <- make_folds(groups, y, n_inner, seed, grouped)
  k_inner    <- length(unique(inner_fold))
  rows <- list(); rc <- 0L
  
  for (cand in names(X_list)) {
    X  <- X_list[[cand]]
    mc <- min(max_ncomp, ncol(X), nrow(X) - 1L)
    if (mc < 1L) next
    
    pred_by_comp <- matrix(NA_real_, nrow(X), mc)
    for (f in seq_len(k_inner)) {
      idx_tr <- which(inner_fold != f)
      idx_va <- which(inner_fold == f)
      if (length(idx_tr) < 10L || length(idx_va) < 2L) next
      
      X_tr <- X[idx_tr, , drop = FALSE]
      X_va <- X[idx_va, , drop = FALSE]
      mcf  <- min(mc, ncol(X_tr), nrow(X_tr) - 1L)
      if (mcf < 1L) next
      
      cal_df <- data.frame(y = y[idx_tr], NIR = I(X_tr))
      fit <- tryCatch(
        pls::plsr(y ~ NIR, ncomp = mcf, data = cal_df,
                  scale = FALSE, method = "oscorespls"),
        error = function(e) NULL)
      if (is.null(fit)) next
      
      pr <- tryCatch(
        predict(fit, newdata = data.frame(NIR = I(X_va)),
                ncomp = seq_len(mcf)),
        error = function(e) NULL)
      if (is.null(pr)) next
      pr <- as.array(pr)
      
      for (nc in seq_len(mcf)) {
        p <- tryCatch(as.numeric(pr[, 1, nc]),
                      error = function(e) rep(NA_real_, length(idx_va)))
        if (length(p) == length(idx_va)) pred_by_comp[idx_va, nc] <- p
      }
    }
    
    for (nc in seq_len(mc)) {
      fold_rmsep <- rep(NA_real_, k_inner)
      for (f in seq_len(k_inner)) {
        idx_va <- which(inner_fold == f)
        p   <- pred_by_comp[idx_va, nc]
        obs <- y[idx_va]
        ok  <- is.finite(p) & is.finite(obs)
        if (sum(ok) < 2L) next
        fold_rmsep[f] <- sqrt(mean((obs[ok] - p[ok])^2))
      }
      mr <- mean(fold_rmsep, na.rm = TRUE)
      if (is.finite(mr)) {
        rc <- rc + 1L
        rows[[rc]] <- data.frame(candidate        = cand,
                                 ncomp            = nc,
                                 inner_mean_RMSEP = mr,
                                 stringsAsFactors = FALSE)
      }
    }
  }
  if (!rc) stop("All PLSR candidates failed during full-data CV.")
  sc <- do.call(rbind, rows)
  sc[which.min(sc$inner_mean_RMSEP), , drop = FALSE]
}

final_win <- select_full(candidate_matrices, y_all, g_all,
                         N_INNER, MAX_NCOMP, SEED + 10000L, USE_GROUPED_CV)
final_key     <- as.character(final_win$candidate)
final_ncomp   <- as.integer(final_win$ncomp)
final_key_row <- match(final_key, candidate_keys$key)
if (is.na(final_key_row)) stop("Final selected candidate not found.")

WINNER_PREPROC <- candidate_keys$variant[final_key_row]
WINNER_BIN     <- as.integer(candidate_keys$bin_size[final_key_row])

X_full <- candidate_matrices[[final_key]]
y_full <- y_all

fit_deploy <- pls::plsr(y_full ~ X_full, ncomp = final_ncomp,
                        scale = FALSE, method = "oscorespls")

cat("\nFinal deployment model:\n")
cat("  Preprocessing:", WINNER_PREPROC, "\n")
cat("  Bin width    :", WINNER_BIN, "\n")
cat("  Components   :", final_ncomp, "\n")
cat("  Full n       :", length(y_full), "\n")

# ============================================================
# 8. Save model, metadata, and the four attached CSV outputs
# ============================================================
timestamp  <- format(Sys.time(), "%Y%m%d_%H%M")
model_name <- sprintf("PLSR_cp_ds%s_bin%d_%s_%s",
                      paste(KEEP_DATASETS, collapse = "_"),
                      WINNER_BIN,
                      gsub("^spec_", "", WINNER_PREPROC),
                      timestamp)
model_file <- file.path(BEST_MODELS_PLSR_DIR, paste0(model_name, ".rds"))
saveRDS(fit_deploy, model_file)

meta <- list(
  algorithm           = "pls::plsr (oscorespls, scale = FALSE)",
  target              = TARGET_NAME,
  target_col          = TARGET_COL,
  subset              = KEEP_DATASETS,
  bin_size            = WINNER_BIN,
  preprocessing       = WINNER_PREPROC,
  n_samples           = length(y_full),
  ncomp               = final_ncomp,
  feature_names       = colnames(X_full),
  train_ids           = rownames(X_full),
  y_range             = range(y_full, na.rm = TRUE),
  selection_history   = selection_history,
  overall_summary     = overall_summary,
  all_candidate_aggregate = aggregate_per_candidate,
  timestamp           = timestamp)
saveRDS(meta, file.path(BEST_MODELS_PLSR_DIR, paste0(model_name, "_meta.rds")))

# ------------------------------------------------------------------
# Output prefix
# ------------------------------------------------------------------
tag <- paste(KEEP_DATASETS, collapse = "_")
filter_tag <- paste0("spc", if (USE_SPECTRAL_FILTER) "on" else "off",
                     "_y", if (USE_Y_BOXPLOT) "on" else "off")
out_prefix <- file.path(
  SUMMARY_DIR,
  paste0("PLSR_cp_ds", tag, "_", filter_tag, "_",
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

# 3) All Models Metrics — 40 rows, ranked by R2_mean
write.csv(aggregate_per_candidate,
          file.path(SUMMARY_DIR,
                    paste0("PLSR_All_Models_Metrics_ds", tag, "_", timestamp, ".csv")),
          row.names = FALSE)

# 4) PLSR Panels A–E — one long-format CSV
PANEL_COLS <- c("Panel","RowType","OuterFold","Preproc","Bin_nm","ncomp",
                "inner_RMSEP","n_train","n_test","R2","R2_SS","RMSEP",
                "MAE","bias","RPD","Metric","Mean","SD","Value","Note")

panel_A <- data.frame(
  Panel      = "A",
  RowType    = "SearchSpace",
  Metric     = c("Preprocessing_variants", "Bin_widths_nm",
                 "PLSR_components", "Outer_CV_folds", "Inner_CV_folds"),
  Value      = c(length(spec_vars),
                 paste(BIN_SIZES, collapse = ";"),
                 paste0("1-", MAX_NCOMP),
                 N_OUTER, N_INNER),
  Note       = c(paste(length(spec_vars), "spectral preprocessing variants"),
                 "", "", "", ""),
  stringsAsFactors = FALSE)

panel_B <- data.frame(
  Panel      = "B",
  RowType    = "NestedSelected",
  OuterFold  = selection_history$OuterFold,
  Preproc    = selection_history$Preproc,
  Bin_nm     = selection_history$Bin,
  ncomp      = selection_history$ncomp,
  inner_RMSEP = selection_history$inner_rmsep,
  Note       = paste0("Selected by minimum inner RMSEP within outer fold ",
                      selection_history$OuterFold),
  stringsAsFactors = FALSE)

panel_C <- data.frame(
  Panel      = "C",
  RowType    = "OuterMetrics",
  OuterFold  = selection_history$OuterFold,
  Preproc    = selection_history$Preproc,
  Bin_nm     = selection_history$Bin,
  ncomp      = selection_history$ncomp,
  n_train    = selection_history$n_train,
  n_test     = selection_history$n_test,
  R2         = selection_history$R2,
  R2_SS      = selection_history$R2_ss,
  RMSEP      = selection_history$RMSEP,
  MAE        = selection_history$MAE,
  bias       = selection_history$bias,
  RPD        = selection_history$RPD,
  stringsAsFactors = FALSE)

panel_D <- data.frame(
  Panel      = "D",
  RowType    = "Aggregate",
  Metric     = overall_summary$Metric,
  Mean       = overall_summary$Mean_outer_fold,
  SD         = overall_summary$SD_outer_fold,
  Note       = "Mean across five outer folds for per-fold nested selections",
  stringsAsFactors = FALSE)

panel_E <- data.frame(
  Panel      = "E",
  RowType    = "Frequency",
  Metric     = freq_df$Config,
  Value      = freq_df$Count,
  Note       = paste0("Selected in ", freq_df$Count,
                      "/", length(unique(selection_history$OuterFold)),
                      " outer folds by inner CV"),
  stringsAsFactors = FALSE)

fill_cols <- function(df) {
  for (cc in setdiff(PANEL_COLS, colnames(df))) df[[cc]] <- NA
  df[, PANEL_COLS, drop = FALSE]
}

plsr_panels <- dplyr::bind_rows(
  fill_cols(panel_A),
  fill_cols(panel_B),
  fill_cols(panel_C),
  fill_cols(panel_D),
  fill_cols(panel_E))

write.csv(plsr_panels,
          file.path(SUMMARY_DIR,
                    paste0("PLSR_Panels_ds", tag, "_", timestamp, ".csv")),
          row.names = FALSE, na = "")

# Target descriptive statistics
target_stats <- data.frame(
  Statistic = c("N","Mean","SD","Min","Q1 (25%)","Median",
                "Q3 (75%)","Max","Range"),
  Value = c(length(y_full), mean(y_full), sd(y_full), min(y_full),
            as.numeric(quantile(y_full, 0.25)), median(y_full),
            as.numeric(quantile(y_full, 0.75)), max(y_full),
            max(y_full) - min(y_full)))
write.csv(target_stats,
          file.path(SUMMARY_DIR,
                    sprintf("Target_Statistics_ds%s_n%d_%s.csv",
                            tag, length(y_full), timestamp)),
          row.names = FALSE)

cat("\nFiles written:\n")
cat("  ", model_file, "\n", sep = "")
cat("  ", file.path(BEST_MODELS_PLSR_DIR, paste0(model_name, "_meta.rds")), "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldmetrics.csv"), "\n", sep = "")
cat("  ", paste0(out_prefix, "_foldpreds.csv"), "\n", sep = "")
cat("  ", file.path(SUMMARY_DIR,
                    paste0("PLSR_All_Models_Metrics_ds", tag, "_", timestamp, ".csv")), "\n", sep = "")
cat("  ", file.path(SUMMARY_DIR,
                    paste0("PLSR_Panels_ds", tag, "_", timestamp, ".csv")), "\n", sep = "")
cat("  ", file.path(SUMMARY_DIR,
                    sprintf("Target_Statistics_ds%s_n%d_%s.csv",
                            tag, length(y_full), timestamp)), "\n", sep = "")
cat("DONE.\n")