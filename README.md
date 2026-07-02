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
  sim1_base_model.R      Single decomposable graph, no clustering.

  sim2_changepoint.R     Split-and-merge sampler for ordered change-point
                        partitions.

  sim3_clustering.R      Collapsed sampler for exchangeable DP mixtures of
                        decomposable GGMs.
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
