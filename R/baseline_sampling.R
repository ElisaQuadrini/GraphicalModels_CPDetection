# =============================================================================
# BASELINE SAMPLING
# Baseline graph measure and predictive probability for a new cluster.
# Used by: sim2_changepoint.R, sim3_clustering.R
# Requires: model_core.R (init_stats, add_point_stats, log_marginal_G)
# =============================================================================

# sample_baseline_graphs: sample S decomposable graphs from the baseline
# measure p(G), approximated by taking k_mix random steps with
# move_decomposable starting from the empty graph. Under pi_edge = 0.5 this
# targets the uniform distribution over decomposable graphs.
# k_mix defaults to 3*p*(p-1)/2 (proportional to the number of possible
# edges), which is a reasonable default to decorrelate from the empty graph.
sample_baseline_graphs <- function(S, p, k_mix = NULL) {

  graphs_list <- vector("list", S)

  A <- matrix(0, p, p)

  for (s in 1:S) {
    for (k in 1:(ifelse(is.null(k_mix), 3 * p * (p - 1) / 2, k_mix))) {
      move_res <- move_decomposable(A)
      A <- move_res$A_new
    }
    graphs_list[[s]] <- A
  }

  return(graphs_list)
}

# normalize_weights: normalize a vector of log-weights (with -Inf allowed)
# into probabilities summing to 1, using a numerically stable log-sum-exp
normalize_weights <- function(prob) {

  n      <- length(prob)
  weight <- numeric(n)
  finite <- prob != -Inf

  if (any(finite)) {
    const          <- mean(prob[finite])
    exp_vals       <- exp(prob[finite] - const)
    weight[finite] <- exp_vals / sum(exp_vals)
  }

  return(weight)
}

# log_pred_new: log predictive probability for a single new observation x
# under an empty cluster (n = 0). G is a baseline graph sampled once by the
# caller and reused for the cluster-assignment step. Since log_marginal_G
# on empty statistics is 0, only the numerator (stats including x) is needed.
log_pred_new <- function(x, G, mu0, n0, delta0, D0) {

  p <- length(x)

  stats_empty <- init_stats(p)
  stats_x     <- add_point_stats(stats_empty, x, mu0, n0)

  cs <- mpd(G)
  cl <- lapply(cs$cliques, as.integer)
  se <- lapply(cs$separators, as.integer)

  log_marginal_G(cl, se, stats_x, delta0, D0)
}
