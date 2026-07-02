# =============================================================================
# MCMC UPDATES
# Reusable update blocks shared across the collapsed samplers.
# Used by: sim2_changepoint.R, sim3_clustering.R
# Requires: model_core.R (marg_S, log_prior_graph), graph_moves.R
#           (move_decomposable), gRbase (mpd)
# =============================================================================

# update_graph_MH: runs n_graph Metropolis-Hastings steps updating the graph
# structure for a single cluster, using the simplified four-node-set
# marginal likelihood ratio (only the clique/separator touched by the
# proposed edge add/remove needs to be recomputed).
# Returns the updated graph together with its refreshed cliques/separators.
update_graph_MH <- function(G_curr, stats_l, n_graph, pi_edge, p, delta0, D0) {

  for (mh_step in 1:n_graph) {

    move_curr    <- move_decomposable(G_curr)
    G_prop       <- move_curr$A_new
    n_moves_curr <- move_curr$n_moves
    nodes_star   <- move_curr$nodes
    type_move    <- move_curr$type

    n_moves_prop <- move_decomposable(G_prop)$n_moves

    # Only four node sets are needed instead of all cliques/separators,
    # since adding/removing one edge changes the decomposition only locally
    if (type_move == "add") {

      cliques_prop <- mpd(G_prop)$cliques
      C.star <- as.integer(unlist(
        cliques_prop[sapply(lapply(cliques_prop, intersect, nodes_star), length) == 2]
      ))

    } else {

      cliques_curr <- mpd(G_curr)$cliques
      C.star <- as.integer(unlist(
        cliques_curr[sapply(lapply(cliques_curr, intersect, nodes_star), length) == 2]
      ))
    }

    C.u <- setdiff(C.star, nodes_star[1])
    C.v <- setdiff(C.star, nodes_star[2])
    C.0 <- setdiff(C.star, nodes_star)

    if (type_move == "add") {
      log_marg_ratio <- marg_S(C.star, stats_l, delta0, D0) +
        marg_S(C.0, stats_l, delta0, D0) -
        marg_S(C.u, stats_l, delta0, D0) -
        marg_S(C.v, stats_l, delta0, D0)
    } else {
      log_marg_ratio <- marg_S(C.u, stats_l, delta0, D0) +
        marg_S(C.v, stats_l, delta0, D0) -
        marg_S(C.star, stats_l, delta0, D0) -
        marg_S(C.0, stats_l, delta0, D0)
    }

    log_prior_ratio <- log_prior_graph(G_prop, pi_edge, p) -
      log_prior_graph(G_curr, pi_edge, p)

    log_alpha <- log_marg_ratio + log_prior_ratio +
      log(n_moves_curr) - log(n_moves_prop)

    if (log(runif(1)) < log_alpha) {
      G_curr <- G_prop
    }
  }

  cs_new <- mpd(G_curr)

  list(
    G          = G_curr,
    cliques    = lapply(cs_new$cliques,    as.integer),
    separators = lapply(cs_new$separators, as.integer)
  )
}

# update_alpha0_gibbs: Escobar-West auxiliary variable Gibbs step for the
# Dirichlet Process concentration parameter alpha0, given a Gamma(c_alpha,
# d_alpha) prior, the current number of active clusters K_curr, and the
# total number of observations n.
update_alpha0_gibbs <- function(alpha0, K_curr, n, c_alpha, d_alpha) {

  eta <- rbeta(1, alpha0 + 1, n)

  ratio <- (c_alpha + K_curr - 1) / (n * (-log(eta)))
  g     <- ratio / (1 + ratio)

  if (runif(1) < g) {
    rgamma(1, shape = c_alpha + K_curr,     rate = d_alpha - log(eta))
  } else {
    rgamma(1, shape = c_alpha + K_curr - 1, rate = d_alpha - log(eta))
  }
}
