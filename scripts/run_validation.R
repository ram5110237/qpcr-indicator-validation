#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
#   Rscript scripts/run_validation.R config/config.yaml
# Runs each module whose input is given in the config; missing inputs are skipped.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
get_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a))) else getwd() }
repo <- normalizePath(file.path(get_dir(), ".."))
for (f in list.files(file.path(repo, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args)) args[1] else file.path(repo, "config", "config.yaml")
cfg <- yaml::read_yaml(cfg_path)
base <- dirname(normalizePath(cfg_path))
p <- function(x) if (is.null(x)) NULL else normalizePath(file.path(base, x), mustWork = FALSE)
out <- p(cfg$outdir); dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
rd <- function(f) utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
wt <- function(df, f) utils::write.table(df, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)
sv <- function(pl, f, w = 6, h = 4.5) ggsave(file.path(out, "figures", f), pl, width = w, height = h, dpi = 300)
th <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
fmt <- function(x, d = 3) formatC(x, digits = d, format = "fg")
report <- c("# Indicator validation report", "", paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
md_table <- function(df) {
  df[] <- lapply(df, function(x) if (is.numeric(x)) fmt(x) else as.character(x))
  c(paste("|", paste(names(df), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

# 1. Reference genes ------------------------------------------------------------------
if (!is.null(cfg$reference_genes)) {
  log_msg("Reference-gene stability")
  rc <- cfg$reference_genes
  d <- rd(p(rc$file))
  ct <- as.matrix(d[, unlist(rc$candidates)])
  res <- consensus_ranking(ct, groups = d[[rc$group_column]])
  wt(res$consensus, "refgenes_consensus_ranking.tsv")
  wt(res$genorm$M_steps, "refgenes_genorm_M_steps.tsv"); wt(res$genorm$pairwise_V, "refgenes_genorm_V.tsv")
  wt(res$normfinder, "refgenes_normfinder.tsv"); wt(res$bestkeeper, "refgenes_bestkeeper.tsv")
  wt(res$delta_ct, "refgenes_delta_ct.tsv")
  gm <- res$genorm$ranking
  sv(ggplot(gm, aes(stats::reorder(gene, -M_final), M_final)) + geom_point(size = 3) +
       geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
       labs(x = NULL, y = "geNorm M (at exclusion step)", title = "geNorm stability (lower = more stable)") + th,
     "genorm_M.png")
  sv(ggplot(res$genorm$pairwise_V, aes(comparison, V)) + geom_col(fill = "#4393c3") +
       geom_hline(yintercept = 0.15, linetype = 2) + labs(x = NULL, y = "Pairwise variation V",
       title = paste("Optimal number of reference genes:", res$genorm$optimal_n)) + th, "genorm_V.png")
  best <- res$consensus$gene[seq_len(res$genorm$optimal_n)]
  report <- c(report, "## Reference-gene stability", "", md_table(res$consensus), "",
              paste0("geNorm suggests ", res$genorm$optimal_n, " reference genes (V < 0.15). ",
                     "Recommended set: ", paste(best, collapse = ", "), "."), "")

  if (!is.null(rc$target)) {
    re <- relative_expression(d[[rc$target]], ct[, best, drop = FALSE], d[[rc$group_column]], rc$control)
    worst <- utils::tail(res$consensus$gene, 1)
    re_bad <- relative_expression(d[[rc$target]], ct[, worst, drop = FALSE], d[[rc$group_column]], rc$control)
    cmp <- rbind(cbind(normaliser = paste(best, collapse = "+"), re$contrasts),
                 cbind(normaliser = paste(worst, "(least stable)"), re_bad$contrasts))
    wt(cmp, "relative_expression_target.tsv")
    report <- c(report, paste0("### Relative expression of ", rc$target, ": effect of the normaliser"), "",
                md_table(cmp[, c("normaliser", "comparison", "fold_change", "ci_low", "ci_high", "pvalue")]), "")
  }
}

# 2. Assay validation ----------------------------------------------------------------------
if (!is.null(cfg$standard_curve)) {
  log_msg("Assay validation")
  sc <- rd(p(cfg$standard_curve$file))
  quant_levels <- sc$quantity >= cfg$standard_curve$min_quantity_for_curve
  curve <- standard_curve(sc$quantity[quant_levels], sc$cq[quant_levels])
  lod <- lod_probit(sc$quantity, sc$detected)
  loq <- loq_cv(sc$quantity, sc$cq, curve, cfg$standard_curve$max_cv_loq)
  prec <- precision_components(sc)
  wt(curve, "assay_standard_curve.tsv"); wt(lod, "assay_lod.tsv"); wt(loq$table, "assay_loq_table.tsv")
  wt(detection_table(sc$quantity, sc$detected), "assay_detection_rates.tsv"); wt(prec, "assay_precision.tsv")
  sv(ggplot(sc[!is.na(sc$cq), ], aes(quantity, cq, colour = run)) + geom_point(alpha = 0.7) +
       scale_x_log10() + geom_abline(intercept = curve$intercept, slope = curve$slope, colour = "black") +
       labs(x = "Template copies (log scale)", y = "Cq",
            title = sprintf("Standard curve: E = %.1f%%, R² = %.3f", curve$efficiency_pct, curve$r_squared)) + th,
     "standard_curve.png")
  report <- c(report, "## Assay performance", "",
              sprintf("Efficiency %.1f%% (95%% CI %.1f-%.1f), slope %.3f, R² %.4f, MIQE criteria met: %s.",
                      curve$efficiency_pct, curve$efficiency_ci_low, curve$efficiency_ci_high,
                      curve$slope, curve$r_squared, curve$passes_miqe),
              sprintf("LOD95 (probit): %s copies (95%% CI %s-%s). LOQ (CV <= %s%%): %s copies.",
                      fmt(lod$lod), fmt(lod$lod_ci_low), fmt(lod$lod_ci_high),
                      cfg$standard_curve$max_cv_loq, fmt(loq$loq)), "",
              md_table(prec[order(-prec$quantity), ]), "")
}

# 3. Marker validation ----------------------------------------------------------------------
if (!is.null(cfg$markers)) {
  log_msg("Marker validation")
  mc <- cfg$markers
  mk <- rd(p(mc$file))
  ev <- marker_panel_evaluation(mk, unlist(mc$scores), mc$reference, k = mc$cv_folds, repeats = mc$cv_repeats)
  wt(ev$single, "markers_roc_single.tsv"); wt(ev$combined, "markers_cv_combined.tsv")
  call <- if (mc$call_direction == "below") mk[[mc$call_marker]] < mc$call_threshold else mk[[mc$call_marker]] > mc$call_threshold
  dm <- diagnostic_metrics(call, mk[[mc$reference]])
  wt(dm, "markers_diagnostic_metrics.tsv")
  rocs <- do.call(rbind, lapply(unlist(mc$scores), function(s) {
    r <- roc_marker(mk[[s]], mk[[mc$reference]])$roc
    data.frame(marker = s, fpr = 1 - r$specificities, tpr = r$sensitivities)
  }))
  sv(ggplot(rocs[order(rocs$marker, rocs$fpr, rocs$tpr), ], aes(fpr, tpr, colour = marker)) + geom_step() +
       geom_abline(linetype = 2, colour = "grey60") + coord_equal() +
       labs(x = "1 - specificity", y = "Sensitivity", title = "Marker ROC curves") + th, "marker_roc.png", 6, 5)
  report <- c(report, "## Marker validation", "",
              md_table(ev$single[, c("marker", "auc", "auc_ci_low", "auc_ci_high", "youden_threshold", "cv_auc_mean")]), "",
              sprintf("Combined model (%s): cross-validated AUC %.3f ± %.3f.", ev$combined$model,
                      ev$combined$cv_auc_mean, ev$combined$cv_auc_sd), "",
              sprintf("Diagnostic accuracy of %s %s %s against %s:", mc$call_marker,
                      ifelse(mc$call_direction == "below", "<", ">"), mc$call_threshold, mc$reference), "",
              md_table(dm), "")
}

# 4. Dose response -------------------------------------------------------------------------
if (!is.null(cfg$dose_response)) {
  log_msg("Dose-response (EC50)")
  dc <- cfg$dose_response
  dr <- rd(p(dc$file))
  fits <- fit_isolates(dr, dc$reference_isolate, fix_lower = dc$fix_lower)
  wt(fits, "ec50_resistance_factors.tsv")
  curves <- predict_curves(dr, fits)
  pts <- dr[dr$dose > 0, ]
  sv(ggplot(pts, aes(dose, response, colour = isolate)) + geom_point(alpha = 0.6) +
       geom_line(data = curves) + scale_x_log10() +
       labs(x = "Fungicide concentration (mg/L, log scale)", y = "Growth (% of control)",
            title = "Four-parameter log-logistic fits") + th, "dose_response.png", 7, 4.5)
  report <- c(report, "## Fungicide sensitivity", "",
              md_table(fits[, c("isolate", "ec50", "ec50_ci_low", "ec50_ci_high", "rf", "sensitivity_class")]), "")
}

# 5. Multi-environment reproducibility -------------------------------------------------------
if (!is.null(cfg$multi_env)) {
  log_msg("Multi-environment reproducibility")
  ec <- cfg$multi_env
  me <- rd(p(ec$file))
  mod <- fit_indicator_model(me, ec$response, ec$entry, ec$env)
  boot <- bootstrap_indicator(mod, nsim = ec$bootstrap_n)
  bl <- entry_blups(mod); fw <- finlay_wilkinson(me, ec$response, ec$entry, ec$env)
  wt(mod$summary, "multienv_variance_components.tsv"); wt(boot, "multienv_h2_icc_bootstrap.tsv")
  wt(bl, "multienv_entry_blups.tsv"); wt(fw, "multienv_finlay_wilkinson.tsv")
  report <- c(report, "## Reproducibility across environments", "", md_table(mod$summary), "", md_table(boot), "")
}

report <- c(report, "## Session", "", "```", utils::capture.output(utils::sessionInfo())[1:3], "```")
writeLines(report, file.path(out, "REPORT.md"))
file.copy(cfg_path, file.path(out, "config_used.yaml"), overwrite = TRUE)
log_msg("Done. Results in ", out)
