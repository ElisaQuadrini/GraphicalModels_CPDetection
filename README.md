# Clustering and Change Point Detection in Bayesian Graphical Models 

R implementation for the thesis *Clustering and Change Point Detection in
Bayesian Graphical Models*.

## Structure

```
R/
  graph_moves.R        Proposal moves on the space of decomposable graphs
                        (move_decomposable).

  model_core.R          Shared collapsed-sampler primitives: sufficient
                        statistics (init_stats, add_point_stats,
                        remove_point_stats), marginal likelihood via
                        clique/separator factorization (marg_S,
                        log_marginal_G, lgamma_mv), and the graph prior
                        (log_prior_graph).

  baseline_sampling.R   Baseline graph measure and new-cluster predictive
                        (sample_baseline_graphs, normalize_weights,
                        log_pred_new). Used only by the DP mixture models.

  mcmc_updates.R        Reusable MCMC update blocks: per-cluster graph MH
                        step (update_graph_MH) and the Escobar-West
                        auxiliary variable Gibbs step for alpha0
                        (update_alpha0_gibbs).

  validation_checks.R   Shared diagnostics for the thesis Illustrations
                        chapters: graph_recovery_metrics, bayes_p_value,
                        compute_block_pip, plot_matrix, compute_ari,
                        compute_psm.

simulations/
  sim1_base_model.R              Single decomposable graph, no clustering
                                  (Chapter 1 base illustration).

  sim1_base_model_sparse.R       Sparse-graph variant of the base model
                                  illustration (chain graph 1-2-3-4-5, 4
                                  edges out of 10 possible pairs), used to
                                  assess graph recovery under lower edge
                                  density.

  sim2_clustering.R              Collapsed sampler for exchangeable DP
                                  mixtures of decomposable GGMs (Chapter 2
                                  illustration, 3 distinct clusters).

  sim2_clustering_2regimes.R     Shared-regime variant of
                                  the clustering sampler, testing partition
                                  recovery when two of the three underlying
                                  components share the same regime. Since
                                  clustering is exchangeable, cluster order
                                  is not meaningful here.

  sim3_changepoint.R             Split-and-merge sampler for ordered
                                  change-point partitions (Chapter 3
                                  illustration, 3 ordered blocks, 3
                                  distinct regimes).

  sim3_changepoint_12regimes.R   Shared-regime robustness check (adjacent
                                  case): Blocks 1 and 2 are generated from
                                  the same regime (same graph, precision
                                  matrix, and mean) and are temporally
                                  adjacent, while Block 3 is distinct. The
                                  true number of change points therefore
                                  collapses from 2 to 1. Tests whether the
                                  sampler avoids spuriously flagging a
                                  change point between two adjacent
                                  same-regime blocks.

  sim3_changepoint_2regimes.R    Shared-regime robustness check
                                  (non-adjacent case): two of the three
                                  ordered blocks share the same regime, but
                                  are separated by a block generated from a
                                  different regime (e.g. Block 1 and Block
                                  3 share a regime, with Block 2 distinct
                                  in between). Because the shared-regime
                                  blocks are not contiguous, the true
                                  number of change points remains 2, even
                                  though two of the three blocks are
                                  generated from an identical graph,
                                  precision matrix, and mean. Tests whether
                                  the sampler correctly detects both
                                  change points based on local adjacency
                                  rather than on global regime identity.
```

## Notes

- `sim1_base_model.R` keeps its own `log_marginal_full_graph()` function,
  since it computes the G-Wishart normalizing constant directly via
  `gnorm()` on the whole graph rather than through the clique/separator
  factorization used by `log_marginal_G()` in `model_core.R`. The two are
  mathematically different routes to the same quantity and are kept
  separate on purpose.
- `sim1_base_model.R` has no Dirichlet Process component, so it does not
  source `baseline_sampling.R`.
- Run each simulation script from the repository root, e.g.:
  `source("simulations/sim1_base_model.R")`.



## Citation

If you use this code in your research, please cite:

Quadrini, E. (2026). *Clustering and Change Point Detection in Bayesian Graphical Models*.
https://github.com/ElisaQuadrini/GraphicalModels_CPDetection
