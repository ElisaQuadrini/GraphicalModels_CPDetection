# =============================================================================
# SPARSE GRAPH SIMULATION FOR BASE MODEL ILLUSTRATION
# We simulate data from a SPARSE decomposable graph: a chain structure
# 1 - 2 - 3 - 4 - 5, i.e. cliques {1,2}, {2,3}, {3,4}, {4,5} with singleton
# separators {2}, {3}, {4}. Only 4 edges out of 10 possible pairs are present,
# making this considerably sparser than the base model's default graph.

source("R/graph_moves.R")
source("R/model_core.R")
source("R/validation_checks.R")

library(gRbase)
library(BDgraph)
library(mvtnorm)
library(igraph)
library(ggraph)
library(coda)

set.seed(123)
q <- 5
n <- 200

# True sparse decomposable graph: chain 1-2-3-4-5
Adj_true <- matrix(0, q, q)
Adj_true[1, 2] <- Adj_true[2, 1] <- 1
Adj_true[2, 3] <- Adj_true[3, 2] <- 1
Adj_true[3, 4] <- Adj_true[4, 3] <- 1
Adj_true[4, 5] <- Adj_true[5, 4] <- 1

# G-Wishart hyperparameters (true)
b_true <- 3
D_true <- diag(q)
Omega_true <- rgwish(n = 1, adj = Adj_true, b = b_true, D = D_true)

Sigma_true <- solve(Omega_true)
mu_true <- rnorm(q)
X <- rmvnorm(n, mean = mu_true, sigma = Sigma_true)

# =============================================================================
# FIGURE 1: TRUE GRAPH PLOT (SAVED TO PDF)
# =============================================================================
Adj_Omega <- (abs(Omega_true) > 1e-6) * 1
diag(Adj_Omega) <- 0

G_omega <- graph_from_adjacency_matrix(
  Adj_Omega,
  mode = "undirected",
  diag = FALSE
)

tryCatch({
  pdf("figures/graph_true_basemodel.pdf", width = 6, height = 6)

  layout(1)
  plot(
    G_omega,
    layout = layout_in_circle(G_omega),
    vertex.color = "lightblue",
    vertex.size = 25,
    vertex.label.color = "black",
    main = "True Sparse Decomposable Graph (Chain Structure)"
  )

}, finally = dev.off())

# =============================================================================
# FIGURE 2: PAIRS PLOT OF SIMULATED DATA (SAVED TO PDF)
# =============================================================================
tryCatch({
  pdf("figures/pairs_data_basemodel.pdf", width = 8, height = 8)

  layout(1)
  pairs(X, main = "Simulated Data from Sparse Decomposable Graph")

}, finally = dev.off())

###########################
## Prior Hyperparameters ##

# on G
pi_edge <- 0.5 # non informative
M <- q*(q-1)/2

# on mu and Omega
mu1 <- rep(0, q)
n1 <- 1
delta_1 <- q + 2
D1 <- diag(q)

# sufficient statistics needed for the posterior parameters update
Xbar <- colMeans(X)

S_data <- matrix(0, q, q)
for (i in 1:n) {
  d <- matrix(X[i, ] - Xbar, ncol = 1)
  S_data <- S_data + d %*% t(d)
}

dbar <- matrix(Xbar - mu1, ncol = 1)
A_data <- (n * n1) / (n + n1) * (dbar %*% t(dbar))

# log_prior_graph is sourced from R/model_core.R

# log_marginal_full_graph: log marginal likelihood for the whole graph,
# computed directly via the G-Wishart normalizing constant (gnorm). This is
# the single-cluster, non-factorized counterpart of log_marginal_G used in
# the clustering/change-point models, kept local since it is specific to
# this base model.
log_marginal_full_graph <- function(A, delta, D, n, S){

  delta_post <- delta + n
  D_post     <- D + S

  log_norm_prior <- gnorm(adj=A, b=delta, D=D)
  log_norm_post  <- gnorm(adj=A, b=delta_post, D=D_post)

  return(log_norm_post - log_norm_prior)
}



###########################
## Gibbs-MetropolisHast MCMC ##

# initialization
n_iter <- 20000
burnin   <- 1000

Adj_curr <- matrix(0, q, q)   # start from empty graph
Omega_curr <- diag(q)
mu_curr <- Xbar

Adj_chain <- array(NA, c(q, q, n_iter))
Omega_chain <- array(NA, c(q, q, n_iter))
mu_chain    <- matrix(NA, n_iter, q)

accept <- rep(0, n_iter)


for(t in 1:n_iter){

  # =========================
  # 1) PROPOSE NEW GRAPH

  move_curr <- move_decomposable(Adj_curr)
  Adj_prop  <- move_curr$A_new
  n_moves_curr <- move_curr$n_moves

  move_prop_tmp <- move_decomposable(Adj_prop)
  n_moves_prop  <- move_prop_tmp$n_moves

  # -------------------------
  # PRIOR

  log_prior_curr <- log_prior_graph(Adj_curr, pi_edge, q)
  log_prior_prop <- log_prior_graph(Adj_prop, pi_edge, q)

  # -------------------------
  # MARGINAL LIKELIHOOD

  log_like_curr <- log_marginal_full_graph(
    Adj_curr,
    delta_1,
    D1,
    n,
    S_data)

  log_like_prop <- log_marginal_full_graph(
    Adj_prop,
    delta_1,
    D1,
    n,
    S_data)

  # -------------------------
  # PROPOSAL CORRECTION

  log_q_corr <- log(n_moves_curr) - log(n_moves_prop)

  # -------------------------
  # MH RATIO

  log_alpha <- (log_prior_prop + log_like_prop) -
    (log_prior_curr + log_like_curr) +
    log_q_corr

  if(log(runif(1)) < log_alpha){
    Adj_curr <- Adj_prop
    accept[t] <- 1
  }

  # =========================
  # 2) SAMPLE Omega | G, X

  delta_post <- delta_1 + n
  D_post     <- D1 + S_data

  Omega_curr <- rgwish(
    n = 1,
    adj = Adj_curr,
    b = delta_post,
    D = D_post)

  # =========================
  # 3) SAMPLE mu | Omega, X

  mu_mean <- (n * Xbar + n1 * mu1) / (n + n1)
  mu_cov  <- solve((n + n1) * Omega_curr)

  mu_curr <- mvtnorm::rmvnorm(1, mu_mean, mu_cov)

  # =========================
  # SAVE

  Adj_chain[,,t]   <- Adj_curr
  Omega_chain[,,t] <- Omega_curr
  mu_chain[t,]     <- mu_curr
}


#################
## Diagnostic ##

# =========================
# graph structure

n_edges_chain <- apply(Adj_chain, 3, function(A) sum(A[lower.tri(A)]))

layout(1)
plot(n_edges_chain, type="l",
     main="Trace plot: number of edges",
     ylab="Edges",
     xlab="Iteration")

# Extract post-burn-in adjacency matrices using the MCMC burnin parameter 
Adj_post <- Adj_chain[, , (burnin + 1):n_iter]

# Compute Posterior Inclusion Probabilities (PIP)
pip <- apply(Adj_post, c(1, 2), mean)
round(pip,2)

# maximum a posteriori
Adj_map <- (pip > 0.5) * 1

# graph recovery metrics (shared across the three Illustrations chapters)
recovery_metrics <- graph_recovery_metrics(pip, Adj_true)
print(recovery_metrics)
# only the edge (3,5) is not included -> we have 2 false negatives because
# the graph is undirected (we consider both edges (3,5), (5,3))
# this can be explained by the correlation between 3 and 5, which is
# weak (see cor(X)) -> the data does not strongly support that edge
cor(X)

# acceptance rate
mean(accept)
acf(n_edges_chain) # reasonable, since we move in the neighborhood of the previous decomposable graph


# ==============================================================================
# CONFIGURATION AND FIGURES DIRECTORY SETUP
# ==============================================================================
if (!dir.exists("figures")) {
  dir.create("figures")
}

# ------------------------------------------------------------------------------
# 1. TRACEPLOTS & ACF FOR MU (BASE MODEL)
# Traceplots are saved to PDF, ACF plots are shown on screen only.
# ------------------------------------------------------------------------------

# --- Save traceplots only ---
pdf("figures/mu_plots_basemodel.pdf", width = 5, height = 2.2 * q)
par(mfrow = c(q, 1), mar = c(2, 4, 4, 2))

for (j in 1:q) {
  plot(
    mu_chain[, j],
    type = "l",
    col = "black",
    lwd = 0.5,
    main = paste("Trace plot: mu[", j, "]", sep = ""),
    xlab = "Iteration",
    ylab = paste("mu[", j, "]", sep = "")
  )
}
dev.off()

# --- Show ACF plots on screen (not saved) ---
par(mfrow = c(q, 1), mar = c(2, 4, 4, 2))
for (j in 1:q) {
  acf(mu_chain[, j],
      main = paste("ACF plot: mu[", j, "]", sep = ""),
      xlab = "Lag",
      ylab = paste("mu[", j, "]", sep = "")
  )
}


# ------------------------------------------------------------------------------
# 2. TRACEPLOTS & ACF FOR OMEGA (BASE MODEL)
# Traceplots are saved to PDF, ACF plots are shown on screen only.
# ------------------------------------------------------------------------------
idx <- which(upper.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
n_elements <- nrow(idx)

# --- Save traceplots only ---
pdf("figures/omega_diagnostics_basemodel.pdf", width = 10, height = 12)
par(mfrow = c(5, 3), mar = c(4, 4, 3, 1))
for (k in 1:n_elements) {
  i <- idx[k, 1]
  j <- idx[k, 2]
  
  plot(Omega_chain[i, j, ], type = "l", lwd = 0.5,
       main = paste("Trace: Omega[", i, ",", j, "]", sep = ""),
       xlab = "Iteration", ylab = "Value")
}
dev.off()

# --- Show ACF plots on screen (not saved) ---
par(mfrow = c(5, 3), mar = c(4, 4, 3, 1))
for (k in 1:n_elements) {
  i <- idx[k, 1]
  j <- idx[k, 2]
  
  x <- Omega_chain[i, j, ]
  if (sd(x) > 0) {
    acf(x,
        main = paste("ACF: Omega[", i, ",", j, "]", sep = ""),
        xlab = "Lag")
  } else {
    plot.new()
    title(main = paste("ACF: Omega[", i, ",", j, "] (constant)"), cex.main = 0.8)
  }
}

# ==============================================================================
# EFFECTIVE SAMPLE SIZE (ESS) COMPUTATION
# ==============================================================================
Omega_post <- Omega_chain[, , (burnin + 1):n_iter]
n_post <- dim(Omega_post)[3]

ESS_Omega <- matrix(NA, q, q)
for (i in 1:q) {
  for (j in 1:q) {
    chain_ij <- mcmc(Omega_post[i, j, ])
    ESS_Omega[i, j] <- effectiveSize(chain_ij)
  }
}
cat("\nEffective Sample Size for Omega (Rounded):\n")
print(round(ESS_Omega, 1))


# ==============================================================================
# POSTERIOR PREDICTIVE CHECKS (PPC)
# ==============================================================================
# Sample MCMC iterations strictly from the post-burn-in period (1,001 to 20,000)
S_ppc <- 1000                      
keep <- sample((burnin + 1):n_iter, S_ppc)
X_rep <- array(NA, c(S_ppc, n, q))

# Generate posterior predictive replicated data
for (s in 1:S_ppc) {
  mu_s    <- mu_chain[keep[s], ]
  Omega_s <- Omega_chain[, , keep[s]]
  Sigma_s <- solve(Omega_s)
  
  X_rep[s, , ] <- mvtnorm::rmvnorm(
    n = n,
    mean = mu_s,
    sigma = Sigma_s)
}

# Compute test statistics
t_rep_mean <- apply(X_rep, c(1, 3), mean)
t_rep_var  <- apply(X_rep, c(1, 3), var)
t_obs_mean <- colMeans(X)
t_obs_var  <- apply(X, 2, var)

# Compute marginal Bayesian p-values
p_mean <- sapply(1:q, function(j) bayes_p_value(t_obs_mean[j], t_rep_mean[, j]))
p_var  <- sapply(1:q, function(j) bayes_p_value(t_obs_var[j],  t_rep_var[, j]))

cat("\nBayesian p-values (mean):", round(p_mean, 3), "\n")
cat("Bayesian p-values (var):",  round(p_var, 3),  "\n")


# ------------------------------------------------------------------------------
# 3. PDF FOR POSTERIOR PREDICTIVE CHECKS ON MARGINAL STATS (BASE MODEL)
# ------------------------------------------------------------------------------
pdf("figures/ppc_marginal_stats_basemodel.pdf", width = 10, height = 7)

# Histograms for Means (Flexible grid based on dimension q)
n_cols_grid <- min(3, q)
n_rows_grid <- ceiling(q / n_cols_grid)
par(mfrow = c(n_rows_grid, n_cols_grid), mar = c(4, 4, 3, 1))

for (j in 1:q) {
  hist(t_rep_mean[, j],
       breaks = 40,
       col = "lightgrey",
       main = paste("PPC Mean: X", j, sep = ""),
       xlab = "Value")
  abline(v = t_obs_mean[j], col = "red", lwd = 2)
}

# Histograms for Variances (Outputs to a new page within the same PDF)
par(mfrow = c(n_rows_grid, n_cols_grid), mar = c(4, 4, 3, 1))
for (j in 1:q) {
  hist(t_rep_var[, j],
       breaks = 40,
       col = "lightgrey",
       main = paste("PPC Var: X", j, sep = ""),
       xlab = "Value")
  abline(v = t_obs_var[j], col = "red", lwd = 2)
}
dev.off()


# ==============================================================================
# NETWORK CORRELATIONS & STRUCTURAL PPC (EDGES VS NON-EDGES)
# ==============================================================================
pairs_all <- which(upper.tri(matrix(0, q, q)), arr.ind = TRUE)
K <- nrow(pairs_all)

t_corr_obs <- numeric(K)
t_corr_rep <- matrix(NA, S_ppc, K)

for (k in 1:K) {
  i <- pairs_all[k, 1]
  j <- pairs_all[k, 2]
  
  t_corr_obs[k] <- cor(X[, i], X[, j])
  for (s in 1:S_ppc) {
    t_corr_rep[s, k] <- cor(X_rep[s, , i], X_rep[s, , j])
  }
}

# Identify graph structure using the MAP graph (PIP > 0.5 threshold)
pip <- apply(Adj_chain[, , (burnin + 1):n_iter], c(1, 2), mean)
Adj_map <- (pip > 0.5) * 1

is_edge <- logical(K)
for (k in 1:K) {
  i <- pairs_all[k, 1]
  j <- pairs_all[k, 2]
  is_edge[k] <- Adj_map[i, j] == 1
}

idx_edge   <- which(is_edge)
idx_noedge <- which(!is_edge)


# ------------------------------------------------------------------------------
# 4. PDF FOR PPC NETWORK CORRELATIONS AND GLOBAL DISCREPANZA (BASE MODEL)
# ------------------------------------------------------------------------------
pdf("figures/ppc_network_correlations_basemodel.pdf", width = 8, height = 6)

# Predictive correlations for structural non-edges
if(length(idx_noedge) > 0) {
  layout(1)
  hist(t_corr_rep[, idx_noedge],
       breaks = 40,
       col = "lightblue",
       main = "Predictive correlations (non-edges)",
       xlab = "Correlation")
  abline(v = t_corr_obs[idx_noedge], col = rgb(1, 0, 0, 0.4), lwd = 2)
}

# Predictive correlations for structural edges
if(length(idx_edge) > 0) {
  layout(1)
  hist(t_corr_rep[, idx_edge],
       breaks = 40,
       col = "lightgreen",
       main = "Predictive correlations (edges)",
       xlab = "Correlation")
  abline(v = t_corr_obs[idx_edge], col = rgb(1, 0, 0, 0.4), lwd = 2)
}

# Global structural discrepancy metric over non-edges
T_obs <- sum(t_corr_obs[idx_noedge]^2)
T_rep <- numeric(S_ppc)

for (s in 1:S_ppc) {
  T_rep[s] <- sum(t_corr_rep[s, idx_noedge]^2)
}

layout(1)
hist(T_rep,
     breaks = 40,
     col = "lightblue",
     main = "Global non-edge correlation discrepancy",
     xlab = "Discrepancy Value")
abline(v = T_obs, col = "red", lwd = 2)

dev.off()

cat("Bayesian p-value (global non-edge correlation discrepancy):",
    round(bayes_p_value(T_obs, T_rep), 3), "\n")
