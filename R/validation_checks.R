# =============================================================================
# VALIDATION CHECKS
# Shared diagnostic and validation functions for the thesis Illustrations
# chapters. Used by: sim1_base_model.R, sim2_changepoint.R, sim3_clustering.R
# Requires: mclust (adjustedRandIndex), coda (for ESS/Geweke, used directly
#           in the simulation scripts)
# =============================================================================

# graph_recovery_metrics: standard binary classification metrics for edge
# recovery, comparing an estimated adjacency/PIP matrix to the true one.
# If G_hat contains values other than 0/1, it is thresholded first.
# Only the upper triangle is used (undirected graph, no self-loops).
graph_recovery_metrics <- function(G_hat, G_true, threshold = 0.5) {

  if (any(G_hat != 0 & G_hat != 1)) {
    G_hat <- (G_hat > threshold) * 1
  }

  ut <- upper.tri(G_true)

  est  <- G_hat[ut]
  true <- G_true[ut]

  TP <- sum(est == 1 & true == 1)
  FP <- sum(est == 1 & true == 0)
  FN <- sum(est == 0 & true == 1)
  TN <- sum(est == 0 & true == 0)

  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall      <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1          <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else NA
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  hamming     <- FP + FN

  data.frame(
    TP = TP, FP = FP, FN = FN, TN = TN,
    precision = precision, recall = recall,
    specificity = specificity, F1 = f1,
    hamming_distance = hamming
  )
}

# bayes_p_value: posterior predictive Bayesian p-value, defined as the
# proportion of replicated statistics t_rep at least as extreme as the
# observed statistic t_obs. Values close to 0 or 1 signal model
# misspecification; values close to 0.5 indicate good fit.
bayes_p_value <- function(t_obs, t_rep) {
  mean(t_rep >= t_obs, na.rm = TRUE)
}

# compute_block_pip: posterior inclusion probability matrix for a
# contiguous block of observations (change-point model), computed as the
# average, over post-burnin iterations, of the graph assigned to the
# cluster label that is modal within the block at that iteration.
compute_block_pip <- function(start_obs, end_obs, xi_chain, G_chain, post_idx, p) {

  pip_matrix <- matrix(0, p, p)
  counter    <- 0

  for (t in post_idx) {

    labels_in_block <- xi_chain[t, start_obs:end_obs]
    modal_label      <- as.character(names(which.max(table(labels_in_block))))

    G_t <- G_chain[[t]][[modal_label]]
    if (!is.null(G_t)) {
      pip_matrix <- pip_matrix + G_t
      counter    <- counter + 1
    }
  }

  pip_matrix / counter
}

# plot_matrix: grayscale heatmap of an adjacency or PIP matrix, with axes
# labeled by variable index
plot_matrix <- function(mat, title) {
  p <- nrow(mat)
  image(1:p, 1:p, mat[, p:1], col = gray((32:0) / 32), main = title,
        xaxt = "n", yaxt = "n", xlab = "Variables", ylab = "Variables")
  axis(1, at = 1:p)
  axis(2, at = 1:p, labels = p:1)
  box()
}

# compute_ari: Adjusted Rand Index between an estimated and a true
# partition vector. Thin wrapper around mclust::adjustedRandIndex, kept
# here so that all three scripts call the same entry point.
compute_ari <- function(xi_estimated, xi_true) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for compute_ari().")
  }
  mclust::adjustedRandIndex(xi_estimated, xi_true)
}

# compute_psm: posterior similarity (co-clustering) matrix, psm[i, j] is
# the posterior probability that observations i and j are assigned to the
# same cluster, estimated from the post-burnin partition chain.
compute_psm <- function(xi_chain, post_idx) {

  n   <- ncol(xi_chain)
  psm <- matrix(0, n, n)

  for (t in post_idx) {
    psm <- psm + (outer(xi_chain[t, ], xi_chain[t, ], "==") * 1)
  }

  psm / length(post_idx)
}
