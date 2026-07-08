# =============================================================================
# COLLAPSED SAMPLER FOR MIXTURES OF DECOMPOSABLE GAUSSIAN GRAPHICAL MODELS
# The collapsed sampler marginalizes out {mu_l*} and {K_l*}, operating only
# on the space of cluster assignments {xi_j} and graphs {G_l*}.

source("R/graph_moves.R")
source("R/model_core.R")
source("R/baseline_sampling.R")
source("R/mcmc_updates.R")
source("R/validation_checks.R")

library(gRbase)
library(BDgraph)
library(mvtnorm)
library(igraph)
library(coda)
library(clue)
library(mclust)

set.seed(42)

# =============================================================================
# SECTION 1: DATA SIMULATION
# We simulate data from L_true = 3 clusters, each with a different sparse
# precision matrix defined on a decomposable graph.

p      <- 5       # number of variables
n      <- 150     # total number of observations
n1     <- 50      # obs per cluster (balanced for simulation)
L_true <- 3        # true number of clusters

# cluster 1 graph: two cliques {1,2,3} and {3,4,5}
Adj1 <- matrix(0, p, p)
Adj1[1:3, 1:3] <- 1
Adj1[3:5, 3:5] <- 1
diag(Adj1) <- 0

# cluster 2 graph: cliques {1,2} and {2,3,4,5}
Adj2 <- matrix(0, p, p)
Adj2[1:2, 1:2] <- 1
Adj2[2:5, 2:5] <- 1
diag(Adj2) <- 0

# cluster 3 graph: complete graph (all edges)
Adj3 <- matrix(1, p, p)
diag(Adj3) <- 0

b_gen <- 3
D_gen <- diag(p)

Omega1_true <- rgwish(n = 1, adj = Adj1, b = b_gen, D = D_gen)
Omega2_true <- rgwish(n = 1, adj = Adj2, b = b_gen, D = D_gen)
Omega3_true <- rgwish(n = 1, adj = Adj3, b = b_gen, D = D_gen)

mu1_true <- c( 2,  2,  0, -1, -1)
mu2_true <- c(-2,  0,  1,  1,  0)
mu3_true <- c( 0, -2,  2,  0,  1)

X1 <- rmvnorm(n1, mean = mu1_true, sigma = solve(Omega1_true))
X2 <- rmvnorm(n1, mean = mu2_true, sigma = solve(Omega2_true))
X3 <- rmvnorm(n1, mean = mu3_true, sigma = solve(Omega3_true))

X       <- rbind(X1, X2, X3)        # n x p
xi_true <- rep(1:L_true, each = n1) # true assignment vector

# Shuffle observations (remove block structure)
perm      <- sample(n)
X         <- X[perm, ]
xi_true   <- xi_true[perm]

pairs(X)


# =============================================================================
# SECTION 2: PRIOR HYPERPARAMETERS

pi_edge <- 0.5   # non-informative graph prior

mu0    <- rep(0, p)
n0     <- 1
delta0 <- p + 2         # G-Wishart degrees of freedom (must be > p-1)
D0     <- diag(p)

c_alpha <- 1
d_alpha <- 1

# =============================================================================
# SECTION 3: MODEL-SPECIFIC HELPER FUNCTIONS
# All shared primitives (log_prior_graph, init_stats, add/remove_point_stats,
# marg_S, log_marginal_G, sample_baseline_graphs, normalize_weights,
# log_pred_new, update_graph_MH, update_alpha0_gibbs) are sourced from R/.
# This model has no additional helper functions beyond the shared ones.


# =============================================================================
# SECTION 4: MCMC INITIALIZATION

n_iter  <- 3000
n_graph <- 7

S_baseline      <- 200
baseline_graphs <- sample_baseline_graphs(S_baseline, p, k_mix = NULL)

xi_curr <- rep(1, n)
L_curr  <- 1

G_list <- list(matrix(0, p, p))

alpha0 <- 1

stats_list      <- vector("list", 1)
cliques_list    <- vector("list", 1)
separators_list <- vector("list", 1)

stats_list[[1]] <- init_stats(p)

for (i in 1:n) {
  stats_list[[1]] <- add_point_stats(
    stats_list[[1]], X[i, ], mu0, n0
  )
}

cs_init <- mpd(G_list[[1]])

cliques_list[[1]]    <- lapply(cs_init$cliques, as.integer)
separators_list[[1]] <- lapply(cs_init$separators, as.integer)


xi_chain     <- matrix(NA, n_iter, n)
L_chain      <- numeric(n_iter)
alpha0_chain <- numeric(n_iter)
G_chain      <- vector("list", n_iter)  # stores active cluster graphs at each iteration

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)


# =============================================================================
# MAIN MCMC LOOP

for (t in 1:n_iter) {

  # ====================================================================
  # STEP 1: UPDATE CLUSTER ASSIGNMENTS {xi_j}

  for (j in 1:n) {

    x_j   <- X[j, ]
    old_l <- xi_curr[j]

    # ------------------------------------------------------------------
    # REMOVE observation j

    stats_list[[old_l]] <- remove_point_stats(stats_list[[old_l]], x_j, mu0, n0)

    xi_curr[j] <- NA

    active_clusters <- which(sapply(stats_list, function(s) s$n_l > 0))
    L_minus_j       <- length(active_clusters)

    log_weights <- numeric(L_minus_j + 1)


    # ------------------------------------------------------------------
    # EXISTING CLUSTERS

    for (idx in seq_along(active_clusters)) {

      l <- active_clusters[idx]

      stats_without <- stats_list[[l]]
      stats_with    <- add_point_stats(stats_without, x_j, mu0, n0)

      cl <- cliques_list[[l]]
      se <- separators_list[[l]]

      log_marg_with    <- log_marginal_G(cl, se, stats_with, delta0, D0)
      log_marg_without <- log_marginal_G(cl, se, stats_without, delta0, D0)

      # posterior predictive for x_j given cluster l (ratio of marginals)
      log_weights[idx] <- log(stats_without$n_l) + (log_marg_with - log_marg_without)
    }


    # ------------------------------------------------------------------
    # NEW CLUSTER

    G_new <- baseline_graphs[[sample(length(baseline_graphs), 1)]]

    log_weights[L_minus_j + 1] <- log(alpha0) + log_pred_new(x_j, G_new, mu0, n0, delta0, D0)


    # ------------------------------------------------------------------
    # SAMPLE ASSIGNMENT

    weights     <- normalize_weights(log_weights)
    sampled_idx <- sample(length(weights), 1, prob = weights)

    if (sampled_idx <= L_minus_j) {

      l_new <- active_clusters[sampled_idx]

    } else {

      l_new <- length(G_list) + 1

      G_list[[l_new]] <- G_new

      cs_new <- mpd(G_new)

      cliques_list[[l_new]]    <- lapply(cs_new$cliques, as.integer)
      separators_list[[l_new]] <- lapply(cs_new$separators, as.integer)

      stats_list[[l_new]] <- init_stats(p)
    }

    # ------------------------------------------------------------------
    # ADD observation j to new cluster

    xi_curr[j] <- l_new

    stats_list[[l_new]] <- add_point_stats(
      stats_list[[l_new]], x_j, mu0, n0
    )
  }


  # ====================================================================
  # STEP 2: UPDATE GRAPHS {G_l} VIA MH

  active_clusters <- unique(xi_curr)

  for (l in active_clusters) {

    updated <- update_graph_MH(
      G_curr  = G_list[[l]],
      stats_l = stats_list[[l]],
      n_graph = n_graph,
      pi_edge = pi_edge,
      p       = p,
      delta0  = delta0,
      D0      = D0
    )

    G_list[[l]]          <- updated$G
    cliques_list[[l]]    <- updated$cliques
    separators_list[[l]] <- updated$separators
  }

  # ====================================================================
  # STEP 3: UPDATE alpha0 VIA AUXILIARY VARIABLE GIBBS STEP

  K_curr <- length(unique(xi_curr))
  alpha0 <- update_alpha0_gibbs(alpha0, K_curr, n, c_alpha, d_alpha)


  # ====================================================================
  # STORE RESULTS

  xi_chain[t, ] <- xi_curr
  L_chain[t]    <- length(unique(xi_curr))
  alpha0_chain[t] <- alpha0

  active_t     <- unique(xi_curr)
  G_chain[[t]] <- setNames(
    lapply(active_t, function(l) G_list[[l]]),
    as.character(active_t)
  )

  setTxtProgressBar(pb, t)
}

close(pb)

# =============================================================================
# SECTION 5: POSTERIOR INFERENCE

burnin <- 500

L_post <- L_chain[(burnin + 1):n_iter]

cat("\nPosterior summary of number of clusters L:\n")
print(table(L_post))
cat("Posterior mode of L:", as.numeric(names(which.max(table(L_post)))), "\n")
cat("True L:", L_true, "\n")

alpha0_post <- alpha0_chain[(burnin + 1):n_iter]

cat("\nPosterior summary of alpha0:\n")
print(summary(alpha0_post))

# Quantitative convergence diagnostics (consistent with sim2_changepoint.R)
ess_L     <- effectiveSize(as.mcmc(L_post))
ess_alpha <- effectiveSize(as.mcmc(alpha0_post))
geweke_L  <- geweke.diag(as.mcmc(L_post))$z

cat("\n==================================================\n")
cat("            MCMC PERFORMANCE METRICS              \n")
cat("==================================================\n")
cat(sprintf("Effective Sample Size for L:      %.2f\n", ess_L))
cat(sprintf("Effective Sample Size for alpha0: %.2f\n", ess_alpha))
cat(sprintf("Geweke Z-score for L:             %.2f\n", geweke_L))
cat("==================================================\n")


# 5.1: Co-clustering (posterior similarity) matrix
# compute_psm is sourced from R/validation_checks.R
xi_post      <- xi_chain[(burnin + 1):n_iter, ]
post_idx     <- (burnin + 1):n_iter
coclustering <- compute_psm(xi_chain, post_idx)


# 5.2: Point estimate of cluster assignment (max co-clustering)
dist_coclustering <- as.dist(1 - coclustering)

hc <- hclust(dist_coclustering, method = "average")
plot(hc)
L_hat <- as.numeric(names(which.max(table(L_post))))
xi_estimated <- cutree(hc, k = L_hat)

cat("\nEstimated cluster sizes:\n")
print(table(xi_estimated))
cat("True cluster sizes:\n")
print(table(xi_true[order(perm)]))



# =============================================================================
# SECTION 6: DIAGNOSTIC PLOTS (saved to PDF)
# =============================================================================

# 6.1: Trace plots of L and alpha0 (saved to PDF)
tryCatch({
  pdf("figures/trace_L_alpha0_clustering.pdf", width = 10, height = 8)
  
  par(mfrow = c(2, 1))
  
  # Trace plot for the number of active clusters
  plot(L_chain, type = "l", col = "steelblue",
       main = "Trace Plot: Number of Active Clusters L",
       xlab = "Iteration", ylab = "L")
  abline(h = L_true, col = "red", lty = 2, lwd = 2)
  abline(v = burnin, col = "gray50", lty = 3)
  legend("topright",
         legend = c("Sampled L", "True L", "Burn-in"),
         col    = c("steelblue", "red", "gray50"),
         lty    = c(1, 2, 3))
  
  # Trace plot for the concentration parameter
  plot(alpha0_chain, type = "l", col = "darkorange",
       main = "Trace Plot: Concentration Parameter alpha0",
       xlab = "Iteration", ylab = "alpha0")
  abline(v = burnin, col = "gray50", lty = 3)
  
}, finally = dev.off())

# 6.2: Co-clustering matrix heatmap (saved to PDF, square aspect ratio)
tryCatch({
  pdf("figures/coclustering_matrix_clustering.pdf", width = 7, height = 7)
  
  layout(1)
  ord <- order(xi_true)
  image(1:n, 1:n,
        coclustering[ord, ord],
        col  = heat.colors(100, rev = TRUE),
        main = "Posterior Co-Clustering Matrix (Ordered by True Labels)",
        xlab = "Observation Index",
        ylab = "Observation Index")
  
}, finally = dev.off())


# =============================================================================
# SECTION 7: POSTERIOR PREDICTIVE CHECKS (PPC)
# =============================================================================

S_ppc <- 200
keep  <- sample((burnin + 1):n_iter, S_ppc)

X_rep <- array(NA, c(S_ppc, n, p))

for (s in 1:S_ppc) {
  t_s  <- keep[s]
  xi_s <- xi_chain[t_s, ]
  L_s  <- max(xi_s)
  
  mu_list_s <- vector("list", L_s)
  K_list_s  <- vector("list", L_s)
  
  for (l in 1:L_s) {
    in_l <- which(xi_s == l)
    X_l  <- X[in_l, , drop = FALSE]
    n_l  <- length(in_l)
    
    if (n_l > 0) {
      Xbar_l <- colMeans(X_l)
      S_l    <- matrix(0, p, p)
      for (i in 1:n_l) {
        d   <- matrix(X_l[i,] - Xbar_l, ncol = 1)
        S_l <- S_l + d %*% t(d)
      }
      dbar <- matrix(Xbar_l - mu0, ncol = 1)
      A_l  <- (n_l * n0) / (n_l + n0) * (dbar %*% t(dbar))
      
      delta_post_l <- delta0 + n_l
      D_post_l     <- D0 + S_l + A_l
      mu_post_l    <- (n_l * Xbar_l + n0 * mu0) / (n_l + n0)
      
      G_l_s <- G_list[[min(l, length(G_list))]]
      
      K_s <- tryCatch(
        rgwish(n = 1, adj = G_l_s, b = delta_post_l, D = D_post_l),
        error = function(e) diag(p)
      )
      K_list_s[[l]] <- K_s
      
      mu_cov_l       <- solve((n_l + n0) * K_s)
      mu_list_s[[l]] <- as.numeric(rmvnorm(1, mu_post_l, mu_cov_l))
      
    } else {
      K_list_s[[l]]  <- diag(p)
      mu_list_s[[l]] <- mu0
    }
  }
  
  for (j in 1:n) {
    l_j           <- xi_s[j]
    mu_j          <- mu_list_s[[l_j]]
    Sigma_j       <- solve(K_list_s[[l_j]])
    X_rep[s, j, ] <- as.numeric(rmvnorm(1, mu_j, Sigma_j))
  }
}

# 7.1: Compute Observed vs Replicated Cluster-Specific Statistics
t_obs_mean_clust <- matrix(NA, L_true, p)
t_obs_var_clust  <- matrix(NA, L_true, p)

for(l in 1:L_true) {
  X_cluster_real <- X[xi_estimated == l, , drop = FALSE]
  t_obs_mean_clust[l, ] <- colMeans(X_cluster_real)
  t_obs_var_clust[l, ]  <- apply(X_cluster_real, 2, var)
}

t_rep_mean_clust <- array(NA, c(S_ppc, L_true, p))
t_rep_var_clust  <- array(NA, c(S_ppc, L_true, p))

for (s in 1:S_ppc) {
  for (l in 1:L_true) {
    in_l_rep <- which(xi_estimated == l)
    
    if (length(in_l_rep) > 1) {
      t_rep_mean_clust[s, l, ] <- colMeans(X_rep[s, in_l_rep, ])
      t_rep_var_clust[s, l, ]  <- apply(X_rep[s, in_l_rep, ], 2, var)
    } else {
      t_rep_mean_clust[s, l, ] <- rep(NA, p)
      t_rep_var_clust[s, l, ]  <- rep(NA, p)
    }
  }
}

# 7.2: Save Cluster Means Histograms to PDF
tryCatch({
  pdf("figures/ppc_means_per_cluster.pdf", width = 3 * L_true, height = 2.5 * p)
  
  par(mfrow = c(p, L_true), mar = c(3, 3, 2.5, 1), oma = c(1, 1, 4, 0))
  for (j in 1:p) {
    for (l in 1:L_true) {
      rep_vals <- t_rep_mean_clust[, l, j]
      rep_vals <- rep_vals[!is.na(rep_vals)]
      
      hist(rep_vals,
           breaks = 20, col = "lightblue", border = "white",
           main = paste("Cluster", l, "- Var", j),
           xlab = "", ylab = "", cex.main = 0.95)
      
      abline(v = t_obs_mean_clust[l, j], col = "red", lwd = 2)
    }
  }
  title("Posterior Predictive Checks: MEANS per Cluster and Variable", outer = TRUE, cex.main = 1.2)
  
}, finally = dev.off())

# 7.3: Save Cluster Variances Histograms to PDF
tryCatch({
  pdf("figures/ppc_variances_per_cluster.pdf", width = 3 * L_true, height = 2.5 * p)
  
  par(mfrow = c(p, L_true), mar = c(3, 3, 2.5, 1), oma = c(1, 1, 4, 0))
  for (j in 1:p) {
    for (l in 1:L_true) {
      rep_vals <- t_rep_var_clust[, l, j]
      rep_vals <- rep_vals[!is.na(rep_vals)]
      
      hist(rep_vals,
           breaks = 20, col = "lightgreen", border = "white",
           main = paste("Cluster", l, "- Var", j),
           xlab = "", ylab = "", cex.main = 0.95)
      
      abline(v = t_obs_var_clust[l, j], col = "red", lwd = 2)
    }
  }
  title("Posterior Predictive Checks: VARIANCES per Cluster and Variable", outer = TRUE, cex.main = 1.2)
  
}, finally = dev.off())

# 7.4: Print Bayesian p-values summary table to console
p_mean_tbl <- sapply(1:L_true, function(l) sapply(1:p, function(j) {
  rep_vals <- t_rep_mean_clust[, l, j]
  bayes_p_value(t_obs_mean_clust[l, j], rep_vals[!is.na(rep_vals)])
}))
p_var_tbl <- sapply(1:L_true, function(l) sapply(1:p, function(j) {
  rep_vals <- t_rep_var_clust[, l, j]
  bayes_p_value(t_obs_var_clust[l, j], rep_vals[!is.na(rep_vals)])
}))

cat("\nBayesian p-values (rows = variables, columns = clusters)\n")
cat("Means:\n"); print(round(p_mean_tbl, 3))
cat("Variances:\n"); print(round(p_var_tbl, 3))


# =============================================================================
# SECTION 8: GRAPH RECOVERY EVALUATION
# =============================================================================

conf_mat <- table(Estimated = xi_estimated, True = xi_true)
cost_mat <- max(conf_mat) - conf_mat
perm_lsap <- solve_LSAP(t(cost_mat))

cat("\nOptimal cluster alignment (true cluster k -> estimated label perm[k]):\n")
for (k in seq_len(ncol(conf_mat))) {
  cat(sprintf("  True cluster %d  ->  Estimated label %d\n", k, perm_lsap[k]))
}

n_true      <- ncol(conf_mat)
est_to_true <- integer(nrow(conf_mat))
for (k in seq_len(n_true)) {
  est_to_true[perm_lsap[k]] <- k
}

# 8.1: Compute Per-Cluster Posterior Inclusion Probabilities (PIP)
L_hat_val   <- as.numeric(names(which.max(table(L_post))))
pip_cluster <- array(0, c(p, p, L_hat_val))
pip_count   <- numeric(L_hat_val)

for (t_idx in (burnin + 1):n_iter) {
  xi_s <- xi_chain[t_idx, ]
  G_t  <- G_chain[[t_idx]]
  
  for (k in seq_len(L_hat_val)) {
    obs_in_k <- which(xi_estimated == k)
    if (length(obs_in_k) == 0) next
    
    mcmc_labs_k <- xi_s[obs_in_k]
    l_match     <- as.integer(names(which.max(table(mcmc_labs_k))))
    l_char      <- as.character(l_match)
    
    if (!is.na(l_match) && l_char %in% names(G_t)) {
      pip_cluster[, , k] <- pip_cluster[, , k] + G_t[[l_char]]
      pip_count[k]       <- pip_count[k] + 1
    }
  }
}

for (k in seq_len(L_hat_val)) {
  if (pip_count[k] > 0)
    pip_cluster[, , k] <- pip_cluster[, , k] / pip_count[k]
}

Adj_list <- list(Adj1, Adj2, Adj3)

# 8.2: Save PIP vs True Adjacency Heatmaps to PDF
tryCatch({
  pdf("figures/graph_recovery_pip_heatmaps_clustering.pdf", width = 8, height = 3.5 * L_hat_val)
  
  par(mfrow = c(L_hat_val, 2), mar = c(3, 3, 3, 1))
  
  for (k in seq_len(L_hat_val)) {
    # Estimated Cluster PIP Heatmap
    image(1:p, 1:p, t(pip_cluster[, , k]),
          zlim = c(0, 1),
          col  = heat.colors(100, rev = TRUE),
          main = paste("Estimated Cluster", k, "- PIP Matrix"),
          xlab = "Node Index", ylab = "Node Index")
    
    true_k <- est_to_true[k]
    
    # Matched True Adjacency Heatmap
    if (true_k >= 1 && true_k <= length(Adj_list)) {
      image(1:p, 1:p, t(Adj_list[[true_k]]),
            zlim = c(0, 1),
            col  = c("white", "steelblue"),
            main = paste("True Cluster", true_k, "- True Adjacency"),
            xlab = "Node Index", ylab = "Node Index")
    } else {
      plot.new()
      title("No True Cluster Matched")
    }
  }
  
}, finally = dev.off())

# 8.3: Compute and Display Recovery Metrics
cat("\nGraph recovery metrics per cluster:\n")
for (k in seq_len(L_hat_val)) {
  true_k <- est_to_true[k]
  if (true_k >= 1 && true_k <= length(Adj_list)) {
    metrics_k <- graph_recovery_metrics(pip_cluster[, , k], Adj_list[[true_k]])
    cat(sprintf("Estimated Cluster %d <-> True Cluster %d:\n", k, true_k))
    print(metrics_k)
  }
}
