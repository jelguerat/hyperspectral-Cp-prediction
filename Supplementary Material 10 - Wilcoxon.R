# ======================================================================
# Paired statistical comparison — NESTED SELECTED models only
#   Paired Wilcoxon signed-rank tests on absolute prediction errors
#   for the six pairwise model comparisons.
#   + matched-pairs rank-biserial correlation (effect size)
#   + four-panel boxplot with free y-axis per model
#
#   Jhohan Elguera - MATE University - 2026
# ======================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(utils)
})

# ----------------------------------------------------------------------
# 0. Paths
# ----------------------------------------------------------------------
SUMMARY_DIR <- "F:/Nitrogen/Machine_Learning_Rstudio_2026/3-Machine_Learning_All_Summary"
SUPP_DIR    <- "F:/Nitrogen/Machine_Learning_Rstudio_2026/7-Supplementary"
OUT_DIR     <- "F:/Nitrogen/Machine_Learning_Rstudio_2026/6-Statistics"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

models <- c("PLSR", "RF", "MBL", "SVR")

# ----------------------------------------------------------------------
# 1. Load helper
# ----------------------------------------------------------------------
latest_file <- function(dir, prefix, suffix) {
  files <- list.files(dir,
                      pattern = paste0("^", prefix, ".*", suffix, "$"),
                      full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files[which.max(file.mtime(files))]
}

# ----------------------------------------------------------------------
# 2. Read REAL nested selections from S1-S4
# ----------------------------------------------------------------------
SELECTION_FILES <- c(
  PLSR = file.path(SUPP_DIR, "PLSR.csv"),
  RF   = file.path(SUPP_DIR, "RF.csv"),
  MBL  = file.path(SUPP_DIR, "MLB.csv"),
  SVR  = file.path(SUPP_DIR, "SVR.csv")
)

stopifnot(all(file.exists(SELECTION_FILES)))

supp_all <- dplyr::bind_rows(
  lapply(names(SELECTION_FILES), function(m) {
    df <- read.csv(SELECTION_FILES[[m]], stringsAsFactors = FALSE, 
                   check.names = FALSE, colClasses = "character")
    df$Model <- m
    df
  })
)

supp_all <- supp_all %>%
  filter(RowType == "NestedSelected") %>%
  mutate(
    OuterFold = as.integer(OuterFold),
    Bin       = as.numeric(Bin_nm),
    ncomp     = as.numeric(ncomp),
    mtry      = as.numeric(mtry),
    best_k    = as.numeric(best_k),
    cost      = as.numeric(cost),
    gamma     = as.numeric(gamma)
  ) %>%
  select(Model, OuterFold, Preproc, Bin, ncomp, mtry, best_k, cost, gamma)

# ----------------------------------------------------------------------
# 3. Per-model loader for NESTED selected predictions
# ----------------------------------------------------------------------
read_preds_nested <- function(m) {
  fp <- latest_file(SUMMARY_DIR, paste0(m, "_cp_ds13"), "_foldpreds.csv")
  if (is.na(fp)) stop("No foldpreds file for ", m)
  
  p <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  s <- supp_all %>% filter(Model == m)
  
  # Dynamically join on hyperparameters to get the exact selected predictions
  join_cols <- c("OuterFold", "Preproc", "Bin")
  for (col in c("ncomp", "mtry", "best_k", "cost", "gamma")) {
    if (col %in% names(p) && col %in% names(s)) {
      if (any(!is.na(p[[col]])) && any(!is.na(s[[col]]))) {
        join_cols <- c(join_cols, col)
      }
    }
  }
  
  p_sel <- p %>% inner_join(s, by = join_cols)
  
  if (nrow(p_sel) == 0) {
    stop("Filtered to zero rows for ", m, 
         " — check that foldpreds contains predictions for the nested selected configs.")
  }
  
  p_sel$Model <- m
  cat(sprintf("  %-5s selected predictions: %d rows\n", m, nrow(p_sel)))
  
  p_sel[, c("sample_id", "OuterFold", "y_true", "y_pred", "Model")]
}

# ----------------------------------------------------------------------
# 4. Load data
# ----------------------------------------------------------------------
cat("Loading nested selected predictions per model...\n")
pred_list <- lapply(models, read_preds_nested)
names(pred_list) <- models

sizes <- sapply(pred_list, nrow)
cat("Rows per model (should be ~305 = n samples):\n")
print(sizes)
if (any(sizes > 1000))
  warning("Some models still have >1000 rows after filtering. ",
          "Check the Preproc/Bin/hyperparameter strings for whitespace or case differences.")

# ----------------------------------------------------------------------
# 5. Build the paired wide table
# ----------------------------------------------------------------------
wide <- pred_list[[1]][, c("sample_id", "OuterFold", "y_true")]
for (m in models) {
  pm <- pred_list[[m]][, c("sample_id", "OuterFold", "y_pred")]
  names(pm)[names(pm) == "y_pred"] <- paste0("y_pred.", m)
  wide <- merge(wide, pm, by = c("sample_id", "OuterFold"), all = FALSE)
  gc(verbose = FALSE)
}

pred_cols <- paste0("y_pred.", models)
wide <- wide[complete.cases(wide[, c("y_true", pred_cols)]), , drop = FALSE]

n_pairs <- nrow(wide)
cat(sprintf("\nPaired observations: %d\n", n_pairs))
stopifnot(n_pairs >= 5)

# ----------------------------------------------------------------------
# 6. Absolute errors per model
# ----------------------------------------------------------------------
for (m in models)
  wide[[paste0("abs_err.", m)]] <- abs(wide$y_true - wide[[paste0("y_pred.", m)]])

# ----------------------------------------------------------------------
# 7. Paired Wilcoxon tests — six pairs
# ----------------------------------------------------------------------
pairs <- t(combn(models, 2))
colnames(pairs) <- c("Model_A", "Model_B")

res <- data.frame(
  Model_A      = character(nrow(pairs)),
  Model_B      = character(nrow(pairs)),
  median_A     = numeric(nrow(pairs)),
  median_B     = numeric(nrow(pairs)),
  median_diff  = numeric(nrow(pairs)),
  mean_diff    = numeric(nrow(pairs)),
  n_pairs      = integer(nrow(pairs)),
  W_statistic  = numeric(nrow(pairs)),
  p_value      = numeric(nrow(pairs)),
  p_adj_BH     = numeric(nrow(pairs)),
  direction    = character(nrow(pairs)),
  stringsAsFactors = FALSE
)

better <- function(a, b, ma, mb) {
  if (ma < mb) a else if (mb < ma) b else "tie"
}

for (i in seq_len(nrow(pairs))) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  ea <- wide[[paste0("abs_err.", a)]]
  eb <- wide[[paste0("abs_err.", b)]]
  wt <- wilcox.test(ea, eb, paired = TRUE,
                    alternative = "two.sided",
                    exact = FALSE, correct = TRUE)
  d       <- ea - eb
  ma      <- median(ea); mb <- median(eb)
  res$Model_A[i]     <- a
  res$Model_B[i]     <- b
  res$median_A[i]    <- ma
  res$median_B[i]    <- mb
  res$median_diff[i] <- median(d)
  res$mean_diff[i]   <- mean(d)
  res$n_pairs[i]     <- length(d)
  res$W_statistic[i] <- unname(wt$statistic)
  res$p_value[i]     <- wt$p.value
  res$direction[i]   <- better(a, b, ma, mb)
}
res$p_adj_BH <- p.adjust(res$p_value, method = "BH")

cat("\n--- Paired Wilcoxon on |error| (nested selected models) ---\n")
print(res[, c("Model_A","Model_B","median_A","median_B",
              "median_diff","p_value","p_adj_BH","direction")],
      row.names = FALSE, digits = 4)

# ----------------------------------------------------------------------
# 7b. Matched-pairs rank-biserial correlation (effect size)
# ----------------------------------------------------------------------
rank_biserial_paired <- function(ea, eb) {
  d <- ea - eb
  d <- d[!is.na(d) & d != 0]        # Wilcoxon drops zero differences
  if (length(d) == 0) return(NA_real_)
  r     <- rank(abs(d))
  r_pos <- sum(r[d > 0])            # rank mass where A has larger error
  r_neg <- sum(r[d < 0])            # rank mass where B has larger error
  (r_pos - r_neg) / (r_pos + r_neg)
}

res$r_rb <- mapply(
  function(a, b) rank_biserial_paired(wide[[paste0("abs_err.", a)]],
                                      wide[[paste0("abs_err.", b)]]),
  res$Model_A, res$Model_B
)

cat("\n--- Rank-biserial effect sizes (|error|) ---\n")
print(res[, c("Model_A","Model_B","median_diff","p_value","p_adj_BH","r_rb")],
      row.names = FALSE, digits = 3)

write.csv(res, file.path(OUT_DIR, "paired_wilcoxon_abs_error.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------
# 8. Same tests on squared errors
# ----------------------------------------------------------------------
res_sq <- res
res_sq[, c("median_A","median_B","median_diff","mean_diff",
           "p_value","p_adj_BH","W_statistic","r_rb")] <- NA_real_
res_sq$direction <- NA_character_

for (i in seq_len(nrow(pairs))) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  ea <- (wide$y_true - wide[[paste0("y_pred.", a)]])^2
  eb <- (wide$y_true - wide[[paste0("y_pred.", b)]])^2
  wt <- wilcox.test(ea, eb, paired = TRUE,
                    alternative = "two.sided",
                    exact = FALSE, correct = TRUE)
  d <- ea - eb
  res_sq$median_A[i]    <- median(ea)
  res_sq$median_B[i]    <- median(eb)
  res_sq$median_diff[i] <- median(d)
  res_sq$mean_diff[i]   <- mean(d)
  res_sq$W_statistic[i] <- unname(wt$statistic)
  res_sq$p_value[i]     <- wt$p.value
  res_sq$direction[i]   <- better(a, b, median(ea), median(eb))
}
res_sq$p_adj_BH <- p.adjust(res_sq$p_value, method = "BH")

res_sq$r_rb <- mapply(
  function(a, b) rank_biserial_paired(
    (wide$y_true - wide[[paste0("y_pred.", a)]])^2,
    (wide$y_true - wide[[paste0("y_pred.", b)]])^2),
  res_sq$Model_A, res_sq$Model_B
)

cat("\n--- Paired Wilcoxon on squared error (nested selected models) ---\n")
print(res_sq[, c("Model_A","Model_B","median_diff",
                 "p_value","p_adj_BH","r_rb","direction")],
      row.names = FALSE, digits = 4)

write.csv(res_sq, file.path(OUT_DIR, "paired_wilcoxon_sq_error.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------
# 9. Paper-ready table (with effect size)
# ----------------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 1e-4, sprintf("%.2e", p), sprintf("%.3f", p))

interp_rb <- function(r) {
  a <- abs(r)
  ifelse(is.na(r), "NA",
         ifelse(a < 0.10, "negligible",
                ifelse(a < 0.30, "small",
                       ifelse(a < 0.50, "moderate", "large"))))
}

paper <- data.frame(
  Comparison        = paste0(res$Model_A, " vs ", res$Model_B),
  Median_err_A      = sprintf("%.2e", res$median_A),
  Median_err_B      = sprintf("%.2e", res$median_B),
  Delta_median_pct  = sprintf("%+.1f%%",
                              100 * res$median_diff /
                                ((res$median_A + res$median_B) / 2)),
  W                 = formatC(res$W_statistic, format = "d", big.mark = ","),
  p                 = fmt_p(res$p_value),
  p_adj_BH          = fmt_p(res$p_adj_BH),
  r_rb              = sprintf("%+.2f", res$r_rb),
  Effect            = interp_rb(res$r_rb),
  Direction         = res$direction,
  stringsAsFactors  = FALSE
)

cat("\n--- Table for the manuscript (with effect size) ---\n")
print(paper, row.names = FALSE)

write.csv(paper, file.path(OUT_DIR, "paired_wilcoxon_table_paper.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------
# 10. Figure — Option B, four-panel boxplot, free y-axis per model
# ----------------------------------------------------------------------
cols <- c("#8C8C8C", "#E69F00", "#009E73", "#D55E00")

draw_panel <- function(i, y, col_i, show_ylabel) {
  y <- y[is.finite(y) & y > 0]
  if (length(y) == 0) return(invisible(NULL))
  
  rng <- range(y)
  pad <- diff(log10(rng)) * 0.05
  ylim_i <- 10^c(log10(rng[1]) - pad, log10(rng[2]) + pad)
  
  boxplot(y,
          col       = col_i,
          border    = "grey20",
          boxwex    = 0.55,
          staplewex = 0.4,
          whisklty  = 2,
          log       = "y",
          ylim      = ylim_i,
          outline   = FALSE,
          xaxt      = "n",
          ylab      = if (show_ylabel)
            expression("|observed - predicted| (g cm"^-2*")")
          else "",
          xlab      = "")
  
  set.seed(42)
  yj <- y
  if (length(yj) > 200) yj <- sample(yj, 200)
  points(jitter(rep(1, length(yj)), amount = 0.08),
         yj, pch = 21, cex = 0.5,
         col = "grey40", bg = "white", lwd = 0.3)
  
  med <- median(y)
  text(1, ylim_i[2], sprintf("med = %.2e", med),
       adj = c(0.5, 1.2), cex = 0.8, col = "grey20")
  
  mtext(models[i], side = 1, line = 1.4, font = 2, cex = 1.1)
  box()
}

pdf(file.path(OUT_DIR, "paired_wilcoxon_abs_error_boxplot.pdf"),
    width = 8, height = 3.2)
par(mfrow = c(1, 4),
    mar   = c(3.5, 3.8, 2.5, 0.5),
    las   = 1,
    xaxs  = "i", yaxs = "i")
for (i in seq_along(models)) {
  draw_panel(i,
             wide[[paste0("abs_err.", models[i])]],
             cols[i],
             show_ylabel = (i == 1))
}
dev.off()

png(file.path(OUT_DIR, "paired_wilcoxon_abs_error_boxplot.png"),
    width = 2400, height = 960, res = 300)
par(mfrow = c(1, 4),
    mar   = c(3.5, 3.8, 2.5, 0.5),
    las   = 1,
    xaxs  = "i", yaxs = "i")
for (i in seq_along(models)) {
  draw_panel(i,
             wide[[paste0("abs_err.", models[i])]],
             cols[i],
             show_ylabel = (i == 1))
}
dev.off()

cat("\nFiles written to:", OUT_DIR, "\n")
cat("  paired_wilcoxon_abs_error.csv\n")
cat("  paired_wilcoxon_sq_error.csv\n")
cat("  paired_wilcoxon_table_paper.csv\n")
cat("  paired_wilcoxon_abs_error_boxplot.pdf\n")
cat("  paired_wilcoxon_abs_error_boxplot.png\n")
