# LHS 475 b pipeline and template sensitivity (14 November 2026)

Two complementary *conditional* numerical diagnostics, separate from the released comparison figures:

1. **Shared-photon pipeline difference:** calculate the binwise difference between two reductions with variance `sigma_a² + sigma_b² - 2 rho sigma_a sigma_b`. Vary the assumed within-bin cross-pipeline correlation `rho`. The three pipelines reduce the same two transits, so independent-reduction error bars must not be assumed without justification. This calculation does not model correlations *between* wavelength bins.
2. **Offset-profiled template amplitude:** whiten a candidate shape and the observed spectrum using an explicitly supplied positive-definite covariance, project out a constant transit-depth offset, and calculate the GLS template coefficient and conditional standard error. This is a descriptive linear template calculation, not a retrieval, a Bayes factor, or an atmosphere detection. Templates selected after inspecting the data invalidate naive significance interpretation.

The functions are import-only and do not change Zenodo inputs, pipeline comparisons, or figures. Run `python -m unittest discover -s tests -p 'test_pipeline_sensitivity.py'` from the repository root. No real-data result is claimed from the synthetic regression fixtures. Full observational validation requires careful matching of released wavelength bins and justified covariance choices.

AI-assisted numerical implementation and tests are disclosed for scientific review. Biswajit Jana is the sole requested Git commit author.
