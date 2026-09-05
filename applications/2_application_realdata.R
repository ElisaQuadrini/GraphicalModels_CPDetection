# ================================================================================================
# BAYESIAN CHANGE-POINT DETECTION VIA COLLAPSED GIBBS SAMPLER WITH DECOMPOSABLE GRAPHS
# Empirical application: 9 US Industry Portfolios (Weekly Log-Returns: 2005 - 2025)
# Run from the repository root: source("application/2_application_realdata.R")
# ================================================================================================

# Shared model functions (see R/ for definitions, reused across chapters)
source("R/graph_moves.R")     # move_decomposable
source("R/model_core.R")      # init_stats, add/remove_point_stats, marg_S,
                               # log_marginal_G, lgamma_mv, log_prior_graph
source("R/baseline_sampling.R") # sample_baseline_graphs

library(gRbase)
library(BDgraph)
library(mvtnorm)
library(igraph)
library(coda)

set.seed(42)

# Ensure output directory exists (figures are not tracked as empty folders by git)
dir.create("figures/application", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SECTION 1: EMPIRICAL DATA LOADING
# =============================================================================

# Load the financial multivariate time-series data matrix produced by
# application/1_eda_realdata.R
Y <- readRDS("data/industry_portfolios_2005_2025.rds")

# Reassign Y to X to maintain naming consistency with the original sampler's core logic
X <- Y           
p <- ncol(X)     # p = 9 covariate sectors (NoDur, Durbl, Manuf, Enrgy, HiTec, Telcm, Shops, Hlth, Utils)
n <- nrow(X)     # n ~ 1096 weekly observations spanning 21 years

cat("--- REAL DATA SHAPE INITIALIZED ---\n")
cat("Dimensions (n weeks x p sectors):", n, "x", p, "\n\n")

date_index <- as.Date(rownames(X))
cat("Date range recovered from X rownames:", format(range(date_index)), "\n")

# =============================================================================
# SECTION 2: CALIBRATED PRIOR HYPERPARAMETERS
# All parameters are adjusted to match empirical returns and allow multiple regimes

# 1. Graph Prior: Uniform non-informative Bernoulli prior over edge inclusion (50% chance per link)
pi_edge <- 0.5   

# 2. Normal-G-Wishart Location/Scale Priors
mu0    <- rep(0, p)     # Weekly log-returns are naturally centered around 0 in the long run
n0     <- 1           # Weak shrinkage (low prior sample size) to let data override the prior mean quickly
delta0 <- p + 2         # Degrees of freedom (9 + 2 = 11), minimal value ensuring finite expected covariance matrix

# CRITICAL FIX: Extremely tight prior scale to force sensitivity to volatility shocks
D0     <- diag(0.01, p) 

# 3. Dirichlet Process Concentration Parameter Prior (Gamma distribution)
c_alpha <- 1            # Reset concentration to let the data rule
d_alpha <- 1

# =============================================================================
# SECTION 3: APPLICATION-SPECIFIC MODEL FUNCTIONS
#
# log_prior_graph, init_stats, add_point_stats, remove_point_stats, lgamma_mv,
# marg_S, log_marginal_G (from R/model_core.R), move_decomposable (from
# R/graph_moves.R), and sample_baseline_graphs (from R/baseline_sampling.R)
# are sourced above and intentionally not redefined here.
#
# NOTE: log_eppf_dp_restricted below is currently duplicated from
# simulations/sim3_changepoint.R. If it is not already centralized in R/,
# it is a good candidate for promotion to a shared file (e.g.
# R/model_core.R), since both the Chapter 3 illustration and this
# application script rely on the same restricted EPPF.
# =============================================================================

# Restricted EPPF for ordered partitions (forces sequential clusters instead of free mixtures)
log_eppf_dp_restricted <- function(alpha, n_all) {
  k <- length(n_all)
  n <- sum(n_all)
  return(lgamma(n + 1) - lgamma(k + 1) + (k) * log(alpha) + lgamma(alpha) - lgamma(alpha + n) - sum(log(n_all)))
}

# =============================================================================
# SECTION 4: MCMC INITIALIZATION (UNBIASED SEQUENTIAL INITIALIZATION)
# =============================================================================

n_iter  <- 5000         # 5000 iterations are sufficient now that D0 is calibrated
n_graph <- 3            # reduced in order to speed up the code execution

# Speed optimization for baseline graph dictionaries
S_baseline      <- 50
baseline_graphs <- sample_baseline_graphs(S_baseline, p, k_mix = 25)

# UNBIASED INITIALIZATION: Split the dataset exactly in half (K = 2)
# This avoids injecting prior knowledge about 2008 or 2020, testing the true sampler power
xi_curr <- c(rep(1, floor(n/2)), rep(2, n - floor(n/2)))
K_start <- length(table(xi_curr))

# Initialize empty graphs for the starting clusters
G_list  <- list(matrix(0, p, p), matrix(0, p, p))
alpha0  <- 1

# Dynamic stats and clique cache list creation based on initial unbiased K
stats_list      <- vector("list", K_start)
cliques_list    <- vector("list", K_start)
separators_list <- vector("list", K_start) 

for(l in 1:K_start) {
  stats_list[[l]] <- init_stats(p)
  obs_in_l <- which(xi_curr == l)
  for (idx in obs_in_l) {
    stats_list[[l]] <- add_point_stats(stats_list[[l]], X[idx, ], mu0, n0)
  }
  cs_init <- mpd(G_list[[l]])
  cliques_list[[l]]    <- lapply(cs_init$cliques, as.integer)
  separators_list[[l]] <- lapply(cs_init$separators, as.integer)
}

# Storage allocation for chains
xi_chain     <- matrix(NA, n_iter, n)
L_chain      <- numeric(n_iter)
alpha0_chain <- numeric(n_iter)
G_chain      <- vector("list", n_iter)

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

# =============================================================================
# MAIN MCMC INTERACTION LOOP
# =============================================================================

for (t in 1:n_iter) {
  
  # --- STEP 1: METROPOLIS-HASTINGS PARTITION UPDATE (SPLIT & MERGE) ---
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
  
  n_all_star    <- as.vector(table(xi_star))
  log_eppf_curr <- log_eppf_dp_restricted(alpha0, n_all_s)
  log_eppf_prop <- log_eppf_dp_restricted(alpha0, n_all_star)
  
  log_ratio <- (log_eppf_prop - log_eppf_curr) + (log_marg_prop - log_marg_curr) + log_q_ratio
  
  if (log(runif(1)) < log_ratio) {
    xi_curr         <- xi_star
    stats_list      <- stats_list_star
    G_list          <- G_list_star
    cliques_list    <- cliques_star
    separators_list <- separators_star
  }
  
  # --- STEP 2: METROPOLIS-HASTINGS CONDITIONAL GRAPH STRUCTURE UPDATES ---
  # NOTE: this block re-implements a per-cluster graph MH step inline. If it
  # matches update_graph_MH() in R/mcmc_updates.R, replace this loop with a
  # call to that function to keep the sampler logic in a single place.
  active_clusters <- unique(xi_curr)
  
  for (l in active_clusters) {
    stats_l <- stats_list[[l]] 
    G_curr  <- G_list[[l]]
    
    for (mh_step in 1:n_graph) {
      move_curr    <- move_decomposable(G_curr)
      G_prop       <- move_curr$A_new
      n_moves_curr <- move_curr$n_moves
      nodes_star   <- move_curr$nodes
      type_move    <- move_curr$type
      
      n_moves_prop <- move_decomposable(G_prop)$n_moves
      
      if (type_move == "add") {
        cliques_prop <- mpd(G_prop)$cliques
        C.star <- as.integer(unlist(cliques_prop[sapply(lapply(cliques_prop, intersect, nodes_star), length) == 2]))
      } else {
        cliques_curr <- mpd(G_curr)$cliques
        C.star <- as.integer(unlist(cliques_curr[sapply(lapply(cliques_curr, intersect, nodes_star), length) == 2]))
      }
      
      C.u <- setdiff(C.star, nodes_star[1])
      C.v <- setdiff(C.star, nodes_star[2])
      C.0 <- setdiff(C.star, nodes_star)
      
      if (type_move == "add") {
        log_marg_ratio <- marg_S(C.star, stats_l, delta0, D0) + marg_S(C.0, stats_l, delta0, D0) -
          marg_S(C.u, stats_l, delta0, D0) - marg_S(C.v, stats_l, delta0, D0)
      } else {
        log_marg_ratio <- marg_S(C.u, stats_l, delta0, D0) + marg_S(C.v, stats_l, delta0, D0) -
          marg_S(C.star, stats_l, delta0, D0) - marg_S(C.0, stats_l, delta0, D0)
      }
      
      log_prior_ratio <- log_prior_graph(G_prop, pi_edge, p) - log_prior_graph(G_curr, pi_edge, p)
      log_alpha       <- log_marg_ratio + log_prior_ratio + log(n_moves_curr) - log(n_moves_prop)
      
      if (!is.finite(log_alpha)) {
        log_alpha <- -Inf 
      }
      
      if (log(runif(1)) < log_alpha) {
        G_curr <- G_prop
        cs_new               <- mpd(G_curr)
        cliques_list[[l]]    <- lapply(cs_new$cliques,    as.integer)
        separators_list[[l]] <- lapply(cs_new$separators, as.integer)
      }
    }
    G_list[[l]] <- G_curr
  }
  
  # --- STEP 3: CONCENTRATION PARAMETER UPDATE (AUXILIARY GIBBS MOVE) ---
  # NOTE: if this matches update_alpha0_gibbs() in R/mcmc_updates.R, replace
  # this block with a call to that function instead of the inline version.
  K_curr <- length(unique(xi_curr))
  eta    <- rbeta(1, alpha0 + 1, n)
  ratio  <- (c_alpha + K_curr - 1) / (n * (-log(eta)))
  g      <- ratio / (1 + ratio)
  
  if (runif(1) < g) {
    alpha0 <- rgamma(1, shape = c_alpha + K_curr,     rate = d_alpha - log(eta))
  } else {
    alpha0 <- rgamma(1, shape = c_alpha + K_curr - 1, rate = d_alpha - log(eta))
  }
  
  # --- LOGGING MCMC CHAINS ---
  xi_chain[t, ]   <- xi_curr
  L_chain[t]      <- length(unique(xi_curr))
  alpha0_chain[t] <- alpha0
  
  active_t     <- unique(xi_curr)
  G_chain[[t]] <- setNames(lapply(active_t, function(l) G_list[[l]]), as.character(active_t))
  
  setTxtProgressBar(pb, t)
}
close(pb)

# =============================================================================
# SECTION 5: OUTPUT DIAGNOSTICS & POSTERIOR PLOTS (EXTENDED BURN-IN)
# =============================================================================

burn_in  <- 2000
post_idx <- (burn_in + 1):n_iter

# --- 1. MCMC Diagnostics Plot ---
pdf("figures/application/realdata_mcmc_diagnostics.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
plot(L_chain, type = "l", col = "darkblue", main = "Trace Plot of K", xlab = "Iteration", ylab = "Active Clusters (K)")
barplot(table(L_chain[post_idx]) / length(post_idx), main = "Posterior of K", col = "lightblue", xlab = "K")
plot(alpha0_chain, type = "l", col = "darkgreen", main = "Trace Plot of alpha0", xlab = "Iteration", ylab = "alpha0")

if (sd(L_chain[post_idx]) > 0) {
  acf(L_chain[post_idx], main = "ACF of K", lag.max = 50)
} else {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "ACF of K (Constant)")
  text(1, 1, "K is static post-burn-in", cex = 1)
}
par(mfrow = c(1, 1))
dev.off()

# --- Calculation of Posterior Similarity Matrix & Change-Point Probabilities ---
psm <- matrix(0, nrow = n, ncol = n)
for (t in post_idx) {
  psm <- psm + (outer(xi_chain[t, ], xi_chain[t, ], "==") * 1)
}
psm <- psm / length(post_idx)

# --- 2. Posterior Similarity Matrix Plot ---
pdf("figures/application/realdata_posterior_similarity_matrix.pdf", width = 8, height = 7)
image(1:n, 1:n, psm, col = heat.colors(32, rev = TRUE),
      main = "Posterior Similarity Matrix (2005 - 2025)",
      xlab = "Time Index (Weeks)", ylab = "Time Index (Weeks)")
dev.off()

cp_probabilities <- numeric(n - 1)
for (t in post_idx) {
  xi_t <- xi_chain[t, ]
  cp_probabilities <- cp_probabilities + ((xi_t[-1] != xi_t[-n]) * 1)
}
cp_probabilities <- cp_probabilities / length(post_idx)

# --- 3. Marginal Change-Point Probabilities Plot ---
pdf("figures/application/realdata_change_point_probabilities.pdf", width = 10, height = 5)
plot(1:(n - 1), cp_probabilities, type = "h", lwd = 1.5, col = "darkred",
     main = "Posterior Marginal Change-Point Probabilities",
     xlab = "Time Index (Weeks)", ylab = "Probability")
dev.off()

cat("\n==================================================\n")
cat("            MCMC PERFORMANCE METRICS              \n")
cat("==================================================\n")
cat(sprintf("Effective Sample Size for alpha0: %.2f\n", effectiveSize(as.mcmc(alpha0_chain[post_idx]))))

if (sd(L_chain[post_idx]) > 0) {
  cat(sprintf("Effective Sample Size for K:      %.2f\n", effectiveSize(as.mcmc(L_chain[post_idx]))))
  cat(sprintf("Geweke Z-score for K:             %.2f\n", geweke.diag(as.mcmc(L_chain[post_idx]))$z))
} else {
  cat("Effective Sample Size for K:      0.00 (Chain settled on static regimes)\n")
  cat("Geweke Z-score for K:             NaN (Zero variance in sample)\n")
}
cat("==================================================\n")

# =============================================================================
# SECTION 6: POSTERIOR PREDICTIVE CHECKS
# =============================================================================

# -----------------------------------------------------------------------------
# 6.1 Point estimate of the partition (Dahl's method)
# -----------------------------------------------------------------------------

dahl_best_idx <- NA
dahl_best_dist <- Inf
for (t in post_idx) {
  co_t <- outer(xi_chain[t, ], xi_chain[t, ], "==") * 1
  dist_t <- sum((co_t - psm)^2)
  if (dist_t < dahl_best_dist) {
    dahl_best_dist <- dist_t
    dahl_best_idx  <- t
  }
}
xi_hat <- xi_chain[dahl_best_idx, ]

change_weeks_hat <- which(xi_hat[-1] != xi_hat[-n])
cat("\n--- DAHL POINT ESTIMATE OF THE PARTITION ---\n")
cat("Selected iteration:", dahl_best_idx, "\n")
cat("Number of clusters (K_hat):", length(unique(xi_hat)), "\n")
cat("Change point week indices:", change_weeks_hat, "\n")
cat("Corresponding approximate dates:", format(date_index[change_weeks_hat]), "\n\n")

# -----------------------------------------------------------------------------
# 6.2 Bayes factor comparison across candidate partitions
# -----------------------------------------------------------------------------

build_stats_for_block <- function(idx_range, X, mu0, n0) {
  stats_l <- init_stats(ncol(X))
  for (idx in idx_range) {
    stats_l <- add_point_stats(stats_l, X[idx, ], mu0, n0)
  }
  stats_l
}

search_map_graph <- function(stats_l, delta0, D0, pi_edge, p, n_search_iter = 300) {
  G_curr <- matrix(0, p, p)
  for (s in 1:n_search_iter) {
    move_curr    <- move_decomposable(G_curr)
    G_prop       <- move_curr$A_new
    n_moves_curr <- move_curr$n_moves
    nodes_star   <- move_curr$nodes
    type_move    <- move_curr$type
    n_moves_prop <- move_decomposable(G_prop)$n_moves
    
    if (type_move == "add") {
      cliques_prop <- mpd(G_prop)$cliques
      C.star <- as.integer(unlist(cliques_prop[sapply(lapply(cliques_prop, intersect, nodes_star), length) == 2]))
    } else {
      cliques_curr <- mpd(G_curr)$cliques
      C.star <- as.integer(unlist(cliques_curr[sapply(lapply(cliques_curr, intersect, nodes_star), length) == 2]))
    }
    
    C.u <- setdiff(C.star, nodes_star[1])
    C.v <- setdiff(C.star, nodes_star[2])
    C.0 <- setdiff(C.star, nodes_star)
    
    if (type_move == "add") {
      log_marg_ratio <- marg_S(C.star, stats_l, delta0, D0) + marg_S(C.0, stats_l, delta0, D0) -
        marg_S(C.u, stats_l, delta0, D0) - marg_S(C.v, stats_l, delta0, D0)
    } else {
      log_marg_ratio <- marg_S(C.u, stats_l, delta0, D0) + marg_S(C.v, stats_l, delta0, D0) -
        marg_S(C.star, stats_l, delta0, D0) - marg_S(C.0, stats_l, delta0, D0)
    }
    
    log_prior_ratio <- log_prior_graph(G_prop, pi_edge, p) - log_prior_graph(G_curr, pi_edge, p)
    log_alpha <- log_marg_ratio + log_prior_ratio + log(n_moves_curr) - log(n_moves_prop)
    if (!is.finite(log_alpha)) log_alpha <- -Inf
    
    if (log(runif(1)) < log_alpha) {
      G_curr <- G_prop
    }
  }
  cs_final <- mpd(G_curr)
  list(G = G_curr,
       cliques    = lapply(cs_final$cliques,    as.integer),
       separators = lapply(cs_final$separators, as.integer))
}

total_log_marginal_partition <- function(xi_vec, X, mu0, n0, delta0, D0, pi_edge, p, n_search_iter = 300) {
  labels <- unique(xi_vec)
  total <- 0
  for (lab in labels) {
    idx_range <- which(xi_vec == lab)
    stats_l <- build_stats_for_block(idx_range, X, mu0, n0)
    gsearch <- search_map_graph(stats_l, delta0, D0, pi_edge, p, n_search_iter)
    total <- total + log_marginal_G(gsearch$cliques, gsearch$separators, stats_l, delta0, D0)
  }
  total
}

crisis_2008_idx <- which(date_index >= as.Date("2008-09-01") & date_index <= as.Date("2009-03-01"))
covid_2020_idx  <- which(date_index >= as.Date("2020-02-15") & date_index <= as.Date("2020-04-30"))

week_2008 <- if (length(crisis_2008_idx) > 0) round(median(crisis_2008_idx)) else NA
week_2020 <- if (length(covid_2020_idx) > 0) round(median(covid_2020_idx)) else NA

cat("--- EDA-INFORMED BREAKPOINTS ---\n")
cat("2008 crisis midpoint week index:", week_2008, "(", format(date_index[week_2008]), ")\n")
cat("2020 COVID crash midpoint week index:", week_2020, "(", format(date_index[week_2020]), ")\n\n")

xi_eda  <- rep(1, n)
if (!is.na(week_2008)) xi_eda[week_2008:n] <- 2
if (!is.na(week_2020)) xi_eda[week_2020:n] <- 3

xi_null <- rep(1, n)

cat("--- COMPUTING TOTAL LOG-MARGINAL LIKELIHOOD FOR EACH CANDIDATE PARTITION ---\n")
cat("This step runs a short graph search per block per partition; it may take a few minutes.\n")

logmarg_hat  <- total_log_marginal_partition(xi_hat,  X, mu0, n0, delta0, D0, pi_edge, p)
logmarg_eda  <- total_log_marginal_partition(xi_eda,  X, mu0, n0, delta0, D0, pi_edge, p)
logmarg_null <- total_log_marginal_partition(xi_null, X, mu0, n0, delta0, D0, pi_edge, p)

cat("\n--- BAYES FACTOR COMPARISON (log scale) ---\n")
cat(sprintf("Sampler partition (xi_hat), K = %d: log-marginal = %.2f\n", length(unique(xi_hat)), logmarg_hat))
cat(sprintf("EDA-informed partition (xi_eda), K = %d: log-marginal = %.2f\n", length(unique(xi_eda)), logmarg_eda))
cat(sprintf("Null partition (xi_null), K = 1: log-marginal = %.2f\n", logmarg_null))
cat(sprintf("log BF(sampler vs EDA-informed) = %.2f\n", logmarg_hat - logmarg_eda))
cat(sprintf("log BF(sampler vs null)         = %.2f\n", logmarg_hat - logmarg_null))
cat(sprintf("log BF(EDA-informed vs null)    = %.2f\n", logmarg_eda - logmarg_null))
cat("Positive log BF favors the first partition; values above roughly 5 on the log scale\n")
cat("indicate strong preference for one partition over the other.\n\n")

# -----------------------------------------------------------------------------
# 6.3 Rolling realized volatility overlaid on posterior change-point probability
# -----------------------------------------------------------------------------

roll_window <- 10
realized_vol <- rowMeans(X^2)
rolling_vol <- stats::filter(realized_vol, rep(1 / roll_window, roll_window), sides = 2)

# --- 4. Volatility vs Change-Point Probabilities Plot ---
pdf("figures/application/realdata_volatility_vs_cp.pdf", width = 10, height = 8)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(date_index, rolling_vol, type = "l", col = "steelblue",
     main = "Rolling Realized Volatility (mean squared log-return, 9 sectors)",
     xlab = "Date", ylab = "Rolling volatility")
abline(v = date_index[week_2008], col = "darkorange", lty = 2)
abline(v = date_index[week_2020], col = "firebrick", lty = 2)
legend("topleft", legend = c("2008 crisis", "2020 COVID crash"),
       col = c("darkorange", "firebrick"), lty = 2, bty = "n", cex = 0.8)

plot(date_index[-1], cp_probabilities, type = "h", col = "darkred",
     main = "Posterior Marginal Change-Point Probability",
     xlab = "Date", ylab = "Probability")
abline(v = date_index[week_2008], col = "darkorange", lty = 2)
abline(v = date_index[week_2020], col = "firebrick", lty = 2)
par(mfrow = c(1, 1))
dev.off()

# -----------------------------------------------------------------------------
# 6.4 Network density comparison across the blocks of xi_hat
# -----------------------------------------------------------------------------

cat("--- NETWORK DENSITY COMPARISON ACROSS BLOCKS OF xi_hat ---\n")
labels_hat <- sort(unique(xi_hat))
for (lab in labels_hat) {
  idx_range <- which(xi_hat == lab)
  stats_l <- build_stats_for_block(idx_range, X, mu0, n0)
  gsearch <- search_map_graph(stats_l, delta0, D0, pi_edge, p)
  n_edges <- sum(gsearch$G[lower.tri(gsearch$G)])
  max_edges <- p * (p - 1) / 2
  cat(sprintf("Block %d (weeks %d-%d, n = %d): edges = %d / %d, density = %.2f\n",
              lab, min(idx_range), max(idx_range), length(idx_range),
              n_edges, max_edges, n_edges / max_edges))
}
cat("\n")

# -----------------------------------------------------------------------------
# 6.5 Simplified predictive discrepancy check around known crisis windows
# -----------------------------------------------------------------------------

discrepancy <- numeric(n)
for (lab in labels_hat) {
  idx_range <- which(xi_hat == lab)
  stats_l <- build_stats_for_block(idx_range, X, mu0, n0)
  delta_post <- delta0 + stats_l$n_l
  D_post     <- D0 + stats_l$S + stats_l$A_adj
  Sigma_pred <- D_post / (delta_post - p - 1)
  Sigma_pred_inv <- solve(Sigma_pred)
  mu_pred <- stats_l$Xbar
  
  for (idx in idx_range) {
    d <- as.numeric(X[idx, ] - mu_pred)
    discrepancy[idx] <- as.numeric(t(d) %*% Sigma_pred_inv %*% d)
  }
}

# --- 5. Predictive Discrepancy Plot ---
pdf("figures/application/realdata_predictive_discrepancy.pdf", width = 10, height = 5)
plot(date_index, discrepancy, type = "l", col = "purple",
     main = "Approximate Predictive Discrepancy by Week (assigned regime)",
     xlab = "Date", ylab = "Squared Mahalanobis-type discrepancy")
abline(v = date_index[week_2008], col = "darkorange", lty = 2)
abline(v = date_index[week_2020], col = "firebrick", lty = 2)
abline(h = qchisq(0.95, df = p), col = "gray40", lty = 3)
legend("topleft", legend = c("2008 crisis", "2020 COVID crash", "chi-sq 95% reference"),
       col = c("darkorange", "firebrick", "gray40"), lty = c(2, 2, 3), bty = "n", cex = 0.8)
dev.off()

