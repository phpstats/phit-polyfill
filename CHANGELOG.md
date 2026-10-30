# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project uses [Semantic Versioning](https://semver.org/).

Each release states the phit commit it reproduces. The same commit is
pinned in `tests/fixtures/REFERENCE`, and the reference fixtures are copied
from it.

## [0.1.0] - Unreleased

Reproduces phit commit
[`f4887a4`](https://github.com/phpstats/phit/commit/f4887a4a293a3aaa129910a6006d691777ecf06f).
phit has no tagged release yet. Requires PHP 8.3 or later.

### Added

- Polyfill bootstrap that defines each function in the `Phit` namespace only
  when the phit extension does not, generated from the extension's API
  contract. Signatures, parameter names, and defaults match the extension.
- `Phit\Polyfill\Polyfill::supports()`, `implementation()`,
  `functions()`, and `nativeLoaded()`.
- Descriptive statistics: `Phit\stat_mean`, `Phit\stat_var`, `Phit\stat_sd`,
  `Phit\stat_median`, `Phit\stat_quantile`, `Phit\stat_quantile_array`, `Phit\stat_cov`,
  `Phit\stat_cor`, `Phit\stat_weighted_mean`, and `Phit\stat_weighted_var`.
- Smoothing: `Phit\stat_density` (R's `density()` with its seven kernels),
  `Phit\stat_bw_nrd0`, `Phit\stat_bw_nrd`, `Phit\stat_approx`, and
  `Phit\stat_ecdf`.
- Optimization: `Phit\optimize` and `Phit\uniroot`, R's `optimize()` and
  `uniroot()` (Brent_fmin and R_zeroin2) evaluated at the same points.
- Regression: `Phit\fit_ols`, `Phit\fit_logistic`, `Phit\fit_probit`, and
  `Phit\fit_poisson`, with weights, offsets, the automatic intercept, listwise
  deletion, and the extension's result keys, diagnostics, and exception codes.
  Results are identical to the extension's to the bit.
- Distributions: 13 families (normal, t, chi-square, F, exponential, logistic,
  log-normal, gamma, beta, binomial, Poisson, negative binomial, and negative
  binomial by mean). Each has the density or mass function, cdf, survival
  function, their logarithms, and the quantile, with `_array` variants (182
  functions).
- Formula interface: `Phit\model_matrix`, `Phit\fit_lm`, and `Phit\fit_glm`, with
  the fitted factor levels and `poly()` coefficients in a `predvars` key.
- `Phit\predict` for `Phit\fit_lm`, `Phit\fit_glm`, `Phit\fit_lmer`, and
  `Phit\fit_glmer` fits, as R's `predict()`: link or response scale, delta-method
  standard errors, confidence and prediction intervals for linear models, and
  random-effect or population predictions for mixed models.
- Seeded sampling: the 13 `*_rng` functions, with the extension's ChaCha20
  generator and inverse-CDF sampling, so the same seed gives the same draws.
- Posterior-draw kernels: `Phit\draws_moments`, `Phit\draws_covariance_map`,
  and `Phit\draws_affine`, including the fused multiply-add accumulation of
  the extension's matrix product.
- Mixed models: `Phit\fit_mixed_intercept`, `Phit\fit_mixed`, `Phit\fit_logistic_mixed`,
  `Phit\fit_lmer`, and `Phit\fit_glmer`, with exact PHP ports of the extension's
  BOBYQA (bobyqa 0.3.2), Nelder–Mead (argmin 0.11.0), and nalgebra SVD,
  symmetric eigendecomposition, and matrix products.
- A randomized differential parity suite that compares the polyfill with the
  loaded extension on seeded random inputs.
- Test suites ported from the extension and run against the R, SciPy, and
  mpmath reference fixtures. A parity suite compares the polyfill with the
  loaded extension.
- `bin/fetch-fixtures`, `bin/generate-bootstrap`, `bin/accuracy-report`, and
  `benchmarks/run.php` with recorded results in `benchmarks/results.json`.
- ACCURACY.md with the largest observed errors per function and reference.
- Content for the package site under `site/`: the manifest, Quarto docs
  rendered to HTML, two vignettes, benchmarks against the extension per model
  and distribution family, and verification files. The parity suites record
  every comparison when `PARITY_LOG` is set, and `bin/site-matrix`,
  `bin/site-data`, `bin/render-site-docs`, and `benchmarks/site.php` produce
  the files.
- CI runs the suite on PHP 8.3, 8.4, and 8.5 without OPcache, with the
  optimizer alone, with the tracing JIT, and with `opcache.jit=function` on
  Linux and macOS, and fails if PHP reports a value it cannot restore at a
  JIT entry point.

### Extension-only

`Phit\read_csv` is left undefined and needs the native extension.

### Known differences from the extension

- The listwise-deletion warning is `E_USER_WARNING`, because PHP code cannot
  raise `E_WARNING`. The text is identical.
- Scalar arguments of the wrong type follow PHP's rules for userland
  functions.
- Every polyfilled `Phit\` function throws an exception with code 1007 while
  OPcache's JIT compiles the process, on every PHP version, because the JIT
  has computed wrong results in the polyfill on PHP 8.3, 8.4, and 8.5. Use
  OPcache with `opcache.jit=disable`.

[0.1.0]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.1.0
