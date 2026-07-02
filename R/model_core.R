# =============================================================================
# MODEL CORE
# Shared primitives for the collapsed samplers (mu_l*, K_l* integrated out).
# Used by: sim2_changepoint.R, sim3_clustering.R
# log_prior_graph is also used by sim1_base_model.R
# =============================================================================

# log_prior_graph: log Bernoulli(pi_edge) prior on graph adjacency matrix A
log_prior_graph <- function(A, pi_edge, p) {
  M  <- p * (p - 1) / 2
  nE <- sum(A[lower.tri(A)])
  return(nE * log(pi_edge) + (M - nE) * log(1 - pi_edge))
}

# init_stats: initialize sufficient statistics for an empty cluster
# p is the number of covariates
init_stats <- function(p) {
  list(
    n_l   = 0,
    Xbar  = rep(0, p),
    S     = matrix(0, p, p),
    A_adj = matrix(0, p, p)
  )
}

# add_point_stats: incrementally update sufficient statistics after adding
# observation x to a cluster with current stats "stats"
add_point_stats <- function(stats, x, mu0, n0) {

  n_old <- stats$n_l
  n_new <- n_old + 1

  x <- as.numeric(x)

  if (n_old == 0) {
    Xbar_new <- x
    S_new    <- matrix(0, length(x), length(x))
  } else {
    Xbar_old <- stats$Xbar
    Xbar_new <- (n_old * Xbar_old + x) / n_new

    d_old <- x - Xbar_old

    S_new <- stats$S + (n_old / n_new) * (d_old %*% t(d_old))
  }

  dbar  <- matrix(Xbar_new - mu0, ncol = 1)
  A_adj <- (n_new * n0) / (n_new + n0) * (dbar %*% t(dbar))

  list(n_l = n_new, Xbar = Xbar_new, S = S_new, A_adj = A_adj)
}

# remove_point_stats: incrementally update sufficient statistics after
# removing observation x from a cluster with current stats "stats"
remove_point_stats <- function(stats, x, mu0, n0) {

  n_old <- stats$n_l
  n_new <- n_old - 1

  if (n_new == 0) {
    return(init_stats(length(x)))
  }

  x <- as.numeric(x)

  Xbar_old <- stats$Xbar
  Xbar_new <- (n_old * Xbar_old - x) / n_new

  d_new <- x - Xbar_new

  S_new <- stats$S - (n_new / n_old) * (d_new %*% t(d_new))

  dbar  <- matrix(Xbar_new - mu0, ncol = 1)
  A_adj <- (n_new * n0) / (n_new + n0) * (dbar %*% t(dbar))

  list(n_l = n_new, Xbar = Xbar_new, S = S_new, A_adj = A_adj)
}

# lgamma_mv: log multivariate gamma function, needed for the G-Wishart
# normalizing constant
lgamma_mv <- function(p, a) {
  (p * (p - 1) / 4) * log(pi) +
    sum(lgamma(a + (1 - seq_len(p)) / 2))
}

# marg_S: log marginal likelihood contribution of a single clique or
# separator (a set of variable indices "vars"), integrating out the
# G-Wishart precision matrix restricted to that node set
marg_S <- function(vars, stats, delta0, D0) {

  v <- length(vars)
  if (v == 0) return(0)

  n_l   <- stats$n_l
  S_l   <- stats$S
  A_adj <- stats$A_adj

  delta_post <- delta0 + n_l
  D_post     <- D0 + S_l + A_adj

  D0_v     <- D0[vars, vars, drop = FALSE]
  D_post_v <- D_post[vars, vars, drop = FALSE]

  log_norm <- function(delta, D_mat) {
    -(delta / 2) * log(det(D_mat)) +
      lgamma_mv(nrow(D_mat), delta / 2)
  }

  -(n_l * v / 2) * log(pi) +
    log_norm(delta_post, D_post_v) -
    log_norm(delta0, D0_v)
}

# log_marginal_G: log marginal likelihood p(X_l | G_l), obtained via the
# clique/separator factorization of a decomposable graph
log_marginal_G <- function(cl, se, stats, delta0, D0) {

  sum(sapply(cl, marg_S, stats = stats, delta0, D0)) -
    sum(sapply(se, marg_S, stats = stats, delta0, D0))
}
