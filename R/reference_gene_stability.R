# ---------------------------------------------------------------------------
# reference_gene_stability.R
# Expression-stability ranking of candidate qPCR reference genes with four
# established algorithms plus a consensus ranking.
#
# Input everywhere: `ct`, a samples x genes matrix of Cq values (no missing).
# Relative quantities assume the per-gene amplification efficiency `eff`
# (fold per cycle, 2 = 100%), which defaults to 2 for every gene.
# ---------------------------------------------------------------------------

#' log2 relative quantities from Cq, efficiency-aware (higher = more template)
log2_rq <- function(ct, eff = NULL) {
  if (is.null(eff)) eff <- rep(2, ncol(ct))
  stopifnot(length(eff) == ncol(ct))
  sweep(-ct, 2, log2(eff), "*")
}

check_ct <- function(ct) {
  if (!is.matrix(ct)) ct <- as.matrix(ct)
  if (anyNA(ct)) stop("Cq matrix contains missing values; impute or drop before ranking")
  if (ncol(ct) < 3) stop("Need at least 3 candidate genes")
  if (is.null(colnames(ct))) stop("Cq matrix needs gene names as column names")
  ct
}

#' geNorm (Vandesompele et al. 2002)
#'
#' M = mean pairwise variation of a gene with all other remaining genes; the
#' least stable gene is excluded stepwise. Pairwise variation V(n/n+1) tells how
#' many genes are needed (V < 0.15 is the usual cut-off).
genorm <- function(ct, eff = NULL) {
  ct <- check_ct(ct)
  lq <- log2_rq(ct, eff)
  m_values <- function(genes) {
    vapply(genes, function(j) mean(vapply(setdiff(genes, j),
                                          function(k) stats::sd(lq[, j] - lq[, k]), 0)), 0)
  }
  remaining <- colnames(ct); steps <- list(); order_out <- character()
  while (length(remaining) > 2) {
    m <- m_values(remaining)
    steps[[length(steps) + 1]] <- data.frame(n_genes = length(remaining), gene = names(m), M = m)
    worst <- names(which.max(m))
    order_out <- c(worst, order_out)
    remaining <- setdiff(remaining, worst)
  }
  m <- m_values(remaining)
  steps[[length(steps) + 1]] <- data.frame(n_genes = 2, gene = names(m), M = m)
  ranking <- c(remaining, order_out)    # most stable first (last two tie)

  # pairwise variation V(n/n+1) using normalisation factors from the ranked genes
  nf <- function(n) rowMeans(lq[, ranking[seq_len(n)], drop = FALSE])   # log2 geometric mean
  v <- vapply(2:(length(ranking) - 1), function(n) stats::sd(nf(n) - nf(n + 1)), 0)
  list(
    ranking = data.frame(rank = seq_along(ranking), gene = ranking,
                         M_final = vapply(ranking, function(g) {
                           s <- do.call(rbind, steps); utils::tail(s$M[s$gene == g], 1)
                         }, 0), row.names = NULL),
    M_steps = do.call(rbind, steps),
    pairwise_V = data.frame(comparison = paste0("V", 2:(length(ranking) - 1), "/", 3:length(ranking)),
                            V = v),
    optimal_n = { k <- which(v < 0.15); if (length(k)) k[1] + 1 else length(ranking) }
  )
}

#' NormFinder (Andersen, Jensen & Orntoft 2004), model-based
#'
#' Estimates intra-group variance per gene corrected for the sample-mean bias,
#' and, when groups are supplied, the inter-group variation shrunk towards zero.
#' Lower stability value = more stable.
normfinder <- function(ct, groups = NULL, eff = NULL) {
  ct <- check_ct(ct)
  y <- log2_rq(ct, eff)                         # samples x genes
  k <- ncol(y)
  if (is.null(groups)) groups <- rep("all", nrow(y))
  groups <- factor(groups)
  G <- nlevels(groups)

  sig2 <- matrix(NA, k, G, dimnames = list(colnames(y), levels(groups)))
  ng <- table(groups)
  for (g in levels(groups)) {
    yg <- y[groups == g, , drop = FALSE]
    if (nrow(yg) < 2) stop("Each group needs at least 2 samples")
    r <- sweep(sweep(yg, 2, colMeans(yg)), 1, rowMeans(yg)) + mean(yg)
    s2 <- colSums(r^2) / ((nrow(yg) - 1) * (1 - 2 / k))
    sig2[, g] <- s2 - sum(s2) / (k * (k - 1))
  }
  sig2 <- pmax(sig2, 0)

  if (G == 1) {
    stab <- sqrt(sig2[, 1])
    return(data.frame(gene = names(stab), stability = stab, intragroup_var = sig2[, 1],
                      row.names = NULL)[order(stab), ])
  }
  grp_mean <- sapply(levels(groups), function(g) colMeans(y[groups == g, , drop = FALSE]))  # k x G
  d <- grp_mean - rowMeans(grp_mean) - matrix(colMeans(grp_mean), k, G, byrow = TRUE) + mean(grp_mean)
  var_mean <- sweep(sig2, 2, as.numeric(ng), "/")           # sigma^2 / n_g
  gamma2 <- max(0, sum(d^2) / ((k - 1) * (G - 1)) - sum(var_mean) / (k * G))
  shrink <- if (gamma2 > 0) gamma2 / (gamma2 + var_mean) else matrix(0, k, G)
  d_tilde <- shrink * d
  # posterior variance of the shrunk group difference (falls back to sigma^2/n_g
  # when no inter-group variation is detected, so genes remain distinguishable)
  var_tilde <- if (gamma2 > 0) gamma2 * var_mean / (gamma2 + var_mean) else var_mean
  rho <- rowMeans(abs(d_tilde) + sqrt(var_tilde))
  out <- data.frame(gene = rownames(sig2), stability = rho,
                    intragroup_var = rowMeans(sig2), intergroup_diff = rowMeans(abs(d)),
                    row.names = NULL)
  out[order(out$stability), ]
}

#' BestKeeper (Pfaffl et al. 2004): descriptive Cq variability and correlation
#' with the BestKeeper index (geometric mean Cq of all candidates).
bestkeeper <- function(ct) {
  ct <- check_ct(ct)
  idx <- exp(rowMeans(log(ct)))
  out <- data.frame(
    gene = colnames(ct),
    mean_ct = colMeans(ct),
    sd_ct = apply(ct, 2, function(x) mean(abs(x - mean(x)))),     # BestKeeper uses mean absolute deviation
    cv_pct = 100 * apply(ct, 2, stats::sd) / colMeans(ct),
    r_index = apply(ct, 2, function(x) stats::cor(x, idx)),
    p_index = apply(ct, 2, function(x) stats::cor.test(x, idx)$p.value),
    row.names = NULL
  )
  out$acceptable_sd <- out$sd_ct <= 1
  out[order(out$sd_ct), ]
}

#' Comparative delta-Ct method (Silver et al. 2006)
delta_ct <- function(ct) {
  ct <- check_ct(ct)
  g <- colnames(ct)
  sds <- vapply(g, function(j) mean(vapply(setdiff(g, j), function(k) stats::sd(ct[, j] - ct[, k]), 0)), 0)
  out <- data.frame(gene = g, mean_sd_delta_ct = sds, row.names = NULL)
  out[order(out$mean_sd_delta_ct), ]
}

#' Consensus ranking: geometric mean of the ranks from each method
#' (the aggregation used by RefFinder).
consensus_ranking <- function(ct, groups = NULL, eff = NULL) {
  gn <- genorm(ct, eff); nf <- normfinder(ct, groups, eff); bk <- bestkeeper(ct); dc <- delta_ct(ct)
  genes <- colnames(ct)
  rk <- data.frame(
    gene = genes,
    geNorm = match(genes, gn$ranking$gene),
    NormFinder = match(genes, nf$gene),
    BestKeeper = match(genes, bk$gene),
    deltaCt = match(genes, dc$gene)
  )
  rk$geomean_rank <- exp(rowMeans(log(rk[, -1])))
  rk <- rk[order(rk$geomean_rank), ]
  rk$consensus_rank <- seq_len(nrow(rk))
  list(consensus = rk, genorm = gn, normfinder = nf, bestkeeper = bk, delta_ct = dc)
}
