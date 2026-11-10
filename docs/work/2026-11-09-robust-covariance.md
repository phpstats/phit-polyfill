# Robust covariance, 0.9.0

Matches native phit HC0–HC3, one-way clustered HC0/HC1, and Newey–West on
formula OLS and IV. Cluster labels join shared original-row missing handling;
cluster inference uses G-1 df separately from residual df. Prediction, tidy,
and summaries use that inference convention. IV diagnostics retain iid formulas.

Shared validation has 107 sandwich cases and 134,673 numerical comparisons per
implementation; maximum absolute reference error is 8.101e-13. Cases include
weighted/no-intercept/offset/factor OLS, weighted and multi-endogenous IV, missing
instruments/clusters, G=2, day clusters, and heating NW6/NW24. Edge tests check
lossless float and 64-bit integer group labels and integral-float normalization.
Independent scientific review confirmed projected weighted IV formulas against
public APIs. No reference-package source was read.

Focused covariance/API/NaN/IV/prediction/summary tests pass (112 tests, 7,211
assertions), as do updated alias signatures, additional edge tests, and PHPCS.
The reference pin records the paired native commit. docs/vcov.md contains the
specification, docs/vcov-validation.json the actual results, and
benchmarks/vcov-results.json the robust-only benchmark measurements.
