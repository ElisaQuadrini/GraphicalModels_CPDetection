# ================================================================================================
# SPLIT AND MERGE COLLAPSED SAMPLER FOR DECOMPOSABLE GRAPH WITH MULTIVARIATE GAUSSIAN ORDERED DATA

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

set.seed(42)


# =============================================================================
# SECTION 1: DATA SIMULATION (FOR ORDERED CHANGE-POINT PARTITIONS)
# We simulate data from K_true = 3 ordered blocks. Observations within each block
# follow a distinct sparse precision matrix defined on a decomposable graph.
# NO SHUFFLING is performed to preserve the temporal/sequential structure.

p      <- 5        # number of variables (dimensions)
nk     <- 50       # number of observations per ordered block
K_true <- 3        # true number of blocks (change-points)
n      <- nk * K_true # total number of observations

# --- Define Decomposable Graphs for each block ---

# Block 1 graph: two cliques {1,2,3} and {3,4,5}
Adj1 <- matrix(0, p, p)
Adj1[1:3, 1:3] <- 1
Adj1[3:5, 3:5] <- 1
diag(Adj1) <- 0

# Block 2 graph: cliques {1,2} and {2,3,4,5}
Adj2 <- matrix(0, p, p)
Adj2[1:2, 1:2] <- 1
Adj2[2:5, 2:5] <- 1
diag(Adj2) <- 0

# Block 3 graph: complete graph (all edges)
Adj3 <- matrix(1, p, p)
diag(Adj3) <- 0

# --- G-Wishart hyperparameters for data generation ---
b_gen <- 3
D_gen <- diag(p)

# Sample true precision matrices (Omega) from G-Wishart given the graphs
Omega1_true <- rgwish(n = 1, adj = Adj1, b = b_gen, D = D_gen)
Omega2_true <- rgwish(n = 1, adj = Adj2, b = b_gen, D = D_gen)
Omega3_true <- rgwish(n = 1, adj = Adj3, b = b_gen, D = D_gen)

# True means (distinct per ordered block)
mu1_true <- c( 2,  2,  0, -1, -1)
mu2_true <- c(-2,  0,  1,  1,  0)
mu3_true <- c( 0, -2,  2,  0,  1)

# --- Sample ordered observations from each block ---
X1 <- rmvnorm(nk, mean = mu1_true, sigma = solve(Omega1_true))
X2 <- rmvnorm(nk, mean = mu2_true, sigma = solve(Omega2_true))
X3 <- rmvnorm(nk, mean = mu3_true, sigma = solve(Omega3_true))

# Combine data sequentially (preserving time-series block structure)
X       <- rbind(X1, X2, X3)        # n x p matrix
xi_true <- rep(1:K_true, each = nk) # true sequential block labels (1,1..., 2,2..., 3,3...)

# --- Diagnostic plot for sequential structural changes ---
plot(rowMeans(X), type = "b", col = xi_true, pch = 20,
     main = "Sequential Data Matrix (Colors represent true structural blocks)",
     xlab = "Observation Index (Time/Sequence)", ylab = "Row Means")
abline(v = c(nk, 2*nk), col = "darkgray", lty = 2, lwd = 2)

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
# SECTION 3: HELPER FUNCTIONS SPECIFIC TO THIS MODEL
# All shared primitives (log_prior_graph, init_stats, add/remove_point_stats,
# marg_S, log_marginal_G, sample_baseline_graphs, normalize_weights,
# log_pred_new, update_graph_MH, update_alpha0_gibbs) are sourced from R/.

# log_eppf_dp_restricted: log EPPF of a restricted Dirichlet Process for
# ordered partitions (change-point models). This is specific to the
# change-point model and has no counterpart in the exchangeable mixture.
log_eppf_dp_restricted <- function(alpha, n_all) {
  k <- length(n_all)
  n <- sum(n_all)

  return(lgamma(n + 1) - lgamma(k + 1) +
           (k) * log(alpha) +
           lgamma(alpha) - lgamma(alpha + n) -
           sum(log(n_all)))
}


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
  # STEP 1: UPDATE PARTITION {xi_curr} VIA SPLIT & MERGE MH MOVES

  n_all_s <- as.vector(table(xi_curr))
  K_curr  <- length(n_all_s)

  p_split <- 0.5
  if (K_curr == 1) {
    move <- "split"
  } else if (K_curr == n) {
    move <- "merge"
  } else {
    move <- sample(c("split", "merge"), 1, prob = c(p_split, 1 - p_split))
  }

  xi_star         <- xi_curr
  stats_list_star <- stats_list
  G_list_star     <- G_list
  cliques_star    <- cliques_list
  separators_star <- separators_list

  if (move == "split") {

    splittable <- which(n_all_s > 1)
    n_gt1      <- length(splittable)
    if (n_gt1 == 0) next

    j   <- sample(splittable, 1)
    n_j <- n_all_s[j]
    cut <- sample(1:(n_j - 1), 1)

    start_j   <- if (j == 1) 1 else sum(n_all_s[1:(j - 1)]) + 1
    split_pos <- start_j + cut - 1

    xi_star[(split_pos + 1):n] <- xi_star[(split_pos + 1):n] + 1
    K_star <- K_curr + 1

    stats_j_v1 <- init_stats(p)
    for (idx in start_j:split_pos) {
      stats_j_v1 <- add_point_stats(stats_j_v1, X[idx, ], mu0, n0)
    }
    stats_j_v2 <- init_stats(p)
    for (idx in (split_pos + 1):(start_j + n_j - 1)) {
      stats_j_v2 <- add_point_stats(stats_j_v2, X[idx, ], mu0, n0)
    }

    stats_list_star      <- append(stats_list_star, list(stats_j_v2), after = j)
    stats_list_star[[j]] <- stats_j_v1

    # Left block inherits the current graph G_j; right block receives a
    # baseline sample
    G_new  <- baseline_graphs[[sample(length(baseline_graphs), 1)]]
    cs_new <- mpd(G_new)

    G_list_star     <- append(G_list_star, list(G_new), after = j)
    cliques_star    <- append(cliques_star, list(lapply(cs_new$cliques, as.integer)), after = j)
    separators_star <- append(separators_star, list(lapply(cs_new$separators, as.integer)), after = j)

    log_marg_curr <- log_marginal_G(cliques_list[[j]], separators_list[[j]], stats_list[[j]], delta0, D0)

    log_marg_prop <- log_marginal_G(cliques_star[[j]], separators_star[[j]], stats_list_star[[j]], delta0, D0) +
      log_marginal_G(cliques_star[[j+1]], separators_star[[j+1]], stats_list_star[[j+1]], delta0, D0)

    if (K_curr == 1) {
      log_q_ratio <- log(1 - p_split) + log(n - 1)
    } else {
      log_q_ratio <- log(1 - p_split) - log(p_split) + log(n_gt1 * (n_j - 1)) - log(K_curr)
    }

  } else { # move == "merge"

    j <- sample(1:(K_curr - 1), 1)

    boundary <- sum(n_all_s[1:j])
    xi_star[(boundary + 1):n] <- xi_star[(boundary + 1):n] - 1
    K_star <- K_curr - 1

    n_all_star_tmp <- as.vector(table(xi_star))
    n_gt1_star     <- sum(n_all_star_tmp > 1)
    merged_size    <- n_all_s[j] + n_all_s[j + 1]

    stats_merged <- stats_list[[j]]
    start_j1     <- boundary + 1
    end_j1       <- boundary + n_all_s[j + 1]
    for (idx in start_j1:end_j1) {
      stats_merged <- add_point_stats(stats_merged, X[idx, ], mu0, n0)
    }

    stats_list_star[[j]] <- stats_merged
    stats_list_star      <- stats_list_star[-(j + 1)]

    G_list_star     <- G_list_star[-(j + 1)]
    cliques_star    <- cliques_star[-(j + 1)]
    separators_star <- separators_star[-(j + 1)]

    log_marg_curr <- log_marginal_G(cliques_list[[j]], separators_list[[j]], stats_list[[j]], delta0, D0) +
      log_marginal_G(cliques_list[[j+1]], separators_list[[j+1]], stats_list[[j+1]], delta0, D0)

    log_marg_prop <- log_marginal_G(cliques_star[[j]], separators_star[[j]], stats_list_star[[j]], delta0, D0)

    if (K_curr == n) {
      log_q_ratio <- if (K_star == 1) log(n - 1) else log(p_split) + log(n - 1)
    } else {
      log_q_ratio <- log(p_split) - log(1 - p_split) + log(K_curr - 1) - log(n_gt1_star * (merged_size - 1))
    }
  }

  # ====================================================================
  # TARGET RATIO CALCULATION & METROPOLIS-HASTINGS ACCEPTANCE

  n_all_star    <- as.vector(table(xi_star))
  log_eppf_curr <- log_eppf_dp_restricted(alpha0, n_all_s)
  log_eppf_prop <- log_eppf_dp_restricted(alpha0, n_all_star)

  log_ratio <- (log_eppf_prop - log_eppf_curr) +
    (log_marg_prop - log_marg_curr) +
    log_q_ratio

  if (log(runif(1)) < log_ratio) {
    xi_curr         <- xi_star
    stats_list      <- stats_list_star
    G_list          <- G_list_star
    cliques_list    <- cliques_star
    separators_list <- separators_star
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
# SECTION 5: PERFORMANCE CHECK AND DIAGNOSTICS (SAVED TO PDF)
# =============================================================================

burn_in <- 500
post_idx <- (burn_in + 1):n_iter

# 5.1: Save MCMC chain trace plots and convergence graphics to PDF
tryCatch({
  pdf("figures/changepoint_mcmc_diagnostics.pdf", width = 10, height = 8)
  
  par(mfrow = c(2, 2))
  
  # Trace plot for the number of active clusters (K)
  plot(L_chain, type = "l", col = "darkblue", main = "Trace Plot of K",
       xlab = "Iteration", ylab = "Number of Active Clusters (K)")
  abline(h = K_true, col = "red", lty = 2, lwd = 2)
  
  # Posterior barplot of K
  barplot(table(L_chain[post_idx]) / length(post_idx),
          main = "Posterior Distribution of K",
          col = "lightblue", xlab = "K", ylab = "Posterior Probability")
  
  # Trace plot for the concentration parameter alpha0
  plot(alpha0_chain, type = "l", col = "darkgreen", main = "Trace Plot of alpha0",
       xlab = "Iteration", ylab = "alpha0")
  
  # Autocorrelation function plot for K
  acf(L_chain[post_idx], main = "ACF of K (Post Burn-in)", lag.max = 50)
  
}, finally = dev.off())

# Compute quantitative convergence diagnostics for K and alpha0
ess_K     <- effectiveSize(as.mcmc(L_chain[post_idx]))
ess_alpha <- effectiveSize(as.mcmc(alpha0_chain[post_idx]))
geweke_K  <- geweke.diag(as.mcmc(L_chain[post_idx]))$z

cat("\n==================================================\n")
cat("            MCMC PERFORMANCE METRICS              \n")
cat("==================================================\n")
cat(sprintf("Effective Sample Size for K:      %.2f\n", ess_K))
cat(sprintf("Effective Sample Size for alpha0: %.2f\n", ess_alpha))
cat(sprintf("Geweke Z-score for K:             %.2f\n", geweke_K))
cat("==================================================\n")


# 5.2: Save Posterior Similarity Matrix (PSM) Heatmap to PDF
tryCatch({
  pdf("figures/changepoint_psm_matrix.pdf", width = 7, height = 7)
  
  par(mfrow = c(1, 1))
  psm <- compute_psm(xi_chain, post_idx)
  
  image(1:n, 1:n, psm, col = heat.colors(32, rev = TRUE),
        main = "Posterior Similarity Matrix (PSM) with True Change-Points",
        xlab = "Observation Index (Time)", ylab = "Observation Index (Time)")
  abline(v = c(nk, 2*nk), col = "blue", lty = 2, lwd = 2)
  abline(h = c(nk, 2*nk), col = "blue", lty = 2, lwd = 2)
  
}, finally = dev.off())


# 5.3: Save Graph Posterior Inclusion Probabilities (PIP) Heatmaps to PDF
pip_block1 <- compute_block_pip(1, 50, xi_chain, G_chain, post_idx, p)
pip_block2 <- compute_block_pip(51, 100, xi_chain, G_chain, post_idx, p)
pip_block3 <- compute_block_pip(101, 150, xi_chain, G_chain, post_idx, p)

tryCatch({
  pdf("figures/changepoint_graph_recovery_pip.pdf", width = 10, height = 7)
  
  par(mfrow = c(2, 3))
  plot_matrix(Adj1, "True Graph: Block 1")
  plot_matrix(Adj2, "True Graph: Block 2")
  plot_matrix(Adj3, "True Graph: Block 3")
  plot_matrix(pip_block1, "Estimated PIP: Block 1")
  plot_matrix(pip_block2, "Estimated PIP: Block 2")
  plot_matrix(pip_block3, "Estimated PIP: Block 3")
  
}, finally = dev.off())

# Display graph recovery metrics per block on the console
cat("\nGraph recovery metrics per block:\n")
cat("Block 1:\n"); print(graph_recovery_metrics(pip_block1, Adj1))
cat("Block 2:\n"); print(graph_recovery_metrics(pip_block2, Adj2))
cat("Block 3:\n"); print(graph_recovery_metrics(pip_block3, Adj3))


# =============================================================================
# DETECTING AND ESTIMATING CHANGE-POINTS (SAVED TO PDF)
# =============================================================================

# Compute marginal change-point probabilities
cp_probabilities <- numeric(n - 1)
for (t in post_idx) {
  xi_t <- xi_chain[t, ]
  cps_at_iteration <- (xi_t[-1] != xi_t[-n]) * 1
  cp_probabilities <- cp_probabilities + cps_at_iteration
}
cp_probabilities <- cp_probabilities / length(post_idx)

# Save Change-Point Probabilities Plot to PDF
tryCatch({
  pdf("figures/changepoint_probabilities_plot.pdf", width = 9, height = 5)
  
  par(mfrow = c(1, 1))
  plot(1:(n-1), cp_probabilities, type = "h", lwd = 2, col = "darkred",
       main = "Posterior Marginal Change-Point Probabilities",
       xlab = "Observation Index (Time/Sequence)", ylab = "Probability of a Change-Point")
  abline(v = c(nk, 2*nk), col = "blue", lty = 2, lwd = 1.5)
  
}, finally = dev.off())

# Extract estimated change-point locations via threshold
threshold <- 0.5
estimated_cps <- which(cp_probabilities > threshold)

cat("\n==================================================\n")
cat("            CHANGE-POINT ESTIMATION REPORT          \n")
cat("==================================================\n")
cat("True Change-Points simulated at observations:", paste(c(nk, 2*nk), collapse = ", "), "\n")
if(length(estimated_cps) > 0) {
  cat("Estimated Change-Points detected at observations:", paste(estimated_cps, collapse = ", "), "\n")
} else {
  cat("No clear change-points detected with threshold =", threshold, "\n")
}

# Compute the MAP (Maximum A Posteriori) partition estimator
partition_strings <- apply(xi_chain[post_idx, ], 1, paste, collapse = "-")
map_partition_string <- names(which.max(table(partition_strings)))
map_partition <- as.numeric(strsplit(map_partition_string, "-")[[1]])
map_cps <- which(map_partition[-1] != map_partition[-n])

cat("Estimated Change-Points from MAP partition:      ", paste(map_cps, collapse = ", "), "\n")
cat("==================================================\n")



# =============================================================================
# SECTION 6: POSTERIOR PREDICTIVE CHECKS (PPC) (SAVED TO PDF)
# =============================================================================

S_ppc <- 200
keep  <- sample(post_idx, S_ppc)
K_map <- length(unique(map_partition))
X_rep <- array(NA, c(S_ppc, n, p))

# 1. Generate replicated datasets from the posterior predictive distribution
for (s in 1:S_ppc) {
  t_s  <- keep[s]
  xi_s <- xi_chain[t_s, ]
  L_s  <- max(xi_s)
  
  for (l in 1:L_s) {
    in_l <- which(xi_s == l)
    if (length(in_l) == 0) next
    
    X_l <- X[in_l, , drop = FALSE]
    n_l <- length(in_l)
    
    # Compute cluster-specific empirical means and sum of squares
    Xbar_l <- colMeans(X_l)
    S_l    <- matrix(0, p, p)
    for (i in 1:n_l) {
      d   <- matrix(X_l[i,] - Xbar_l, ncol = 1)
      S_l <- S_l + d %*% t(d)
    }
    dbar <- matrix(Xbar_l - mu0, ncol = 1)
    A_l  <- (n_l * n0) / (n_l + n0) * (dbar %*% t(dbar))
    
    # Update hyperparameter values for the full conditionals
    delta_post_l <- delta0 + n_l
    D_post_l     <- D0 + S_l + A_l
    mu_post_l    <- (n_l * Xbar_l + n0 * mu0) / (n_l + n0)
    
    # Extract graph structure for current iteration and component
    G_l_s <- G_chain[[t_s]][[as.character(l)]]
    if (is.null(G_l_s)) G_l_s <- diag(p) * 0
    
    # Draw precision matrix from the G-Wishart full conditional distribution
    K_s <- tryCatch(
      rgwish(n = 1, adj = G_l_s, b = delta_post_l, D = D_post_l),
      error = function(e) diag(p)
    )
    
    # Draw mean vector conditional on the sampled precision matrix
    mu_cov_l  <- solve((n_l + n0) * K_s)
    mu_s_l    <- as.numeric(rmvnorm(1, mu_post_l, mu_cov_l))
    Sigma_s_l <- solve(K_s)
    
    # Simulate replicated data points for the active observations
    X_rep[s, in_l, ] <- rmvnorm(n_l, mu_s_l, Sigma_s_l)
  }
}

# 2. Compute observed test statistics per true block/regime
t_obs_mean <- matrix(NA, K_true, p)
t_obs_var  <- matrix(NA, K_true, p)

for (k in 1:K_true) {
  X_k <- X[xi_true == k, , drop = FALSE]
  t_obs_mean[k, ] <- colMeans(X_k)
  t_obs_var[k, ]  <- apply(X_k, 2, var)
}

# 3. Compute replicated test statistics per true block (Corrected 2D Mapping)
t_rep_mean <- array(NA, c(S_ppc, K_true, p))
t_rep_var  <- array(NA, c(S_ppc, K_true, p))

for (s in 1:S_ppc) {
  for (k in 1:K_true) {
    # Extract the 3D slice corresponding to sample 's' and true block 'k'
    X_rep_k <- X_rep[s, xi_true == k, , drop = FALSE]
    
    # Drop the redundant first dimension (sample index 's') to obtain a 2D matrix
    X_rep_k <- array(X_rep_k, dim = dim(X_rep_k)[-1])
    
    # Validate block cardinality to prevent undefined sample variances
    if (length(which(xi_true == k)) > 1) {
      # Compute column-wise means and variances across the variables (Margin = 2)
      t_rep_mean[s, k, ] <- colMeans(X_rep_k)
      t_rep_var[s, k, ]  <- apply(X_rep_k, 2, var)
    } else {
      # Fallback for single observation blocks where variance is undefined
      t_rep_mean[s, k, ] <- X_rep_k
      t_rep_var[s, k, ]  <- rep(NA, p)
    }
  }
}

# Save PPC Means Histograms to PDF
tryCatch({
  pdf("figures/changepoint_ppc_means.pdf", width = 3 * K_true, height = 2.5 * p)
  
  par(mfrow = c(p, K_true), mar = c(3, 3, 2.5, 1), oma = c(1, 1, 4, 0))
  for (j in 1:p) {
    for (k in 1:K_true) {
      hist(t_rep_mean[, k, j], breaks = 20, col = "lightblue", border = "white",
           main = paste("Block", k, "- Var", j), xlab = "", ylab = "", cex.main = 0.95)
      abline(v = t_obs_mean[k, j], col = "red", lwd = 2)
    }
  }
  title("Posterior Predictive Checks: MEANS per Block and Variable", outer = TRUE, cex.main = 1.2)
  
}, finally = dev.off())

# Save PPC Variances Histograms to PDF
tryCatch({
  pdf("figures/changepoint_ppc_variances.pdf", width = 3 * K_true, height = 2.5 * p)
  
  par(mfrow = c(p, K_true), mar = c(3, 3, 2.5, 1), oma = c(1, 1, 4, 0))
  for (j in 1:p) {
    for (k in 1:K_true) {
      hist(t_rep_var[, k, j], breaks = 20, col = "lightgreen", border = "white",
           main = paste("Block", k, "- Var", j), xlab = "", ylab = "", cex.main = 0.95)
      abline(v = t_obs_var[k, j], col = "red", lwd = 2)
    }
  }
  title("Posterior Predictive Checks: VARIANCES per Block and Variable", outer = TRUE, cex.main = 1.2)
  
}, finally = dev.off())

# Print Bayesian p-values summary table
p_mean <- sapply(1:K_true, function(k) sapply(1:p, function(j) bayes_p_value(t_obs_mean[k, j], t_rep_mean[, k, j])))
p_var  <- sapply(1:K_true, function(k) sapply(1:p, function(j) bayes_p_value(t_obs_var[k, j],  t_rep_var[, k, j])))

cat("\nBayesian p-values (rows = variables, columns = blocks)\n")
cat("Means:\n"); print(round(p_mean, 3))
cat("Variances:\n"); print(round(p_var, 3))
