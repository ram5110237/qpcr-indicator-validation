#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Synthetic example data for every module, with known truth:
#   reference_ct.tsv      8 candidate reference genes + 1 target, 3 conditions x 8 samples
#   standard_curve.tsv    10-fold dilution series, 3 runs x 6 replicates, Poisson sampling at low copy
#   marker_panel.tsv      120 isolates: sequencing genotype + three candidate marker scores
#   dose_response.tsv     8 isolates x 9 doses x 3 replicates (mycelial growth, % of control)
#   multi_env.tsv         30 entries x 4 environments x 3 replicates of a quantitative indicator
# All values are artificial.
# ---------------------------------------------------------------------------
set.seed(20261016)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "data/example"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(df, f) utils::write.table(df, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# 1. Reference genes ----------------------------------------------------------
cond <- rep(c("control", "fungicide", "heat"), each = 8)
n <- length(cond)
input <- stats::rnorm(n, 0, 0.8)                         # shared template-input effect
genes <- data.frame(
  gene  = c("EF1a", "UBQ", "H3", "18S", "TUB", "ACT", "CYP", "GAPDH"),
  base  = c(20, 22, 23, 9, 24, 21, 26, 19),
  noise = c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.55),
  fung  = c(0, 0, 0, 0.8, 0, 0.3, 0, 1.5),              # condition-dependent shifts
  heat  = c(0, 0, 0.1, 0.2, 0.5, 0, 0.4, -0.8)
)
ct <- sapply(seq_len(nrow(genes)), function(i) {
  g <- genes[i, ]
  g$base + input + stats::rnorm(n, 0, g$noise) + g$fung * (cond == "fungicide") + g$heat * (cond == "heat")
})
colnames(ct) <- genes$gene
target <- 27 + input + stats::rnorm(n, 0, 0.3) - log2(3) * (cond == "fungicide")   # 3-fold induction
w(data.frame(sample = sprintf("S%02d", 1:n), condition = cond, round(ct, 3), CYP51 = round(target, 3), check.names = FALSE),
  "reference_ct.tsv")

# 2. Standard curve -------------------------------------------------------------
levels_q <- c(1e6, 1e5, 1e4, 1e3, 1e2, 10, 5, 2, 1)
sc <- expand.grid(rep = 1:6, quantity = levels_q, run = paste0("run", 1:3))
run_eff <- c(run1 = 0, run2 = 0.25, run3 = -0.15)
copies <- stats::rpois(nrow(sc), sc$quantity * 0.95)        # sampling of template molecules
sc$cq <- ifelse(copies == 0, NA,
                38.2 - 3.40 * log10(pmax(copies, 1)) + run_eff[sc$run] + stats::rnorm(nrow(sc), 0, 0.12))
sc$detected <- !is.na(sc$cq)
w(sc[, c("run", "quantity", "rep", "cq", "detected")], "standard_curve.tsv")

# 3. Marker panel -----------------------------------------------------------------
m <- 120
resistant <- stats::rbinom(m, 1, 0.4) == 1
w(data.frame(
  isolate = sprintf("I%03d", 1:m),
  sequencing_resistant = resistant,
  asqpcr_dcq = round(ifelse(resistant, stats::rnorm(m, -2, 1.3), stats::rnorm(m, 2, 1.3)), 3),
  hrm_score  = round(ifelse(resistant, stats::rnorm(m, 0.65, 0.2), stats::rnorm(m, 0.4, 0.2)), 3),
  expr_ratio = round(ifelse(resistant, stats::rnorm(m, 1.3, 0.8), stats::rnorm(m, 1.0, 0.8)), 3)
), "marker_panel.tsv")

# 4. Dose response ----------------------------------------------------------------
ec50 <- c(REF = 0.05, A = 0.04, B = 0.08, C = 0.30, D = 0.60, E = 2.5, F = 5, G = 12)
doses <- c(0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)
dr <- expand.grid(rep = 1:3, dose = doses, isolate = names(ec50), stringsAsFactors = FALSE)
true_resp <- 100 / (1 + exp(1.2 * (log(pmax(dr$dose, 1e-12)) - log(ec50[dr$isolate]))))
dr$response <- round(pmax(0, true_resp + stats::rnorm(nrow(dr), 0, 4)), 2)
w(dr[, c("isolate", "dose", "rep", "response")], "dose_response.tsv")
w(data.frame(isolate = names(ec50), true_ec50 = ec50), "dose_response_truth.tsv")

# 5. Multi-environment indicator -----------------------------------------------------
entries <- sprintf("G%02d", 1:30); envs <- paste0("E", 1:4)
me <- expand.grid(rep = 1:3, env = envs, entry = entries, stringsAsFactors = FALSE)
g_eff <- stats::setNames(stats::rnorm(30, 0, 1), entries)
e_eff <- stats::setNames(c(-1, 0, 0.5, 1.5), envs)
ge <- stats::rnorm(30 * 4, 0, sqrt(0.3)); names(ge) <- paste(rep(entries, each = 4), envs)
me$value <- round(10 + g_eff[me$entry] + e_eff[me$env] + ge[paste(me$entry, me$env)] +
                    stats::rnorm(nrow(me), 0, sqrt(0.5)), 3)
w(me[, c("entry", "env", "rep", "value")], "multi_env.tsv")
message("Example data written to ", normalizePath(out))
