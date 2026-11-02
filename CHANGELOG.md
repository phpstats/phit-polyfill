# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project uses [Semantic Versioning](https://semver.org/).

Each release states the phit commit it reproduces. The same commit is
pinned in `tests/fixtures/REFERENCE`, and the reference fixtures are copied
from it.

## [Unreleased]

## [0.5.0] - 2026-11-01

Reproduces phit [v0.5.0](https://github.com/phpstats/phit/releases/tag/v0.5.0),
commit [`17729a8`](https://github.com/phpstats/phit/commit/17729a89b46c3c4234417701e8e86cc22b480607).
This release removes the deprecated names, so code that still calls them
must move to the replacements first.

### Added

- `Phit\summary()` prints a fit from any model function in the layout of R's
  `print(summary(fit))`, and `Phit\tidy()` returns its coefficient table as
  columns, with optional Wald intervals (`confInt: true`). The text is
  byte-identical to the extension's on phit's fixtures for every model
  function.
- Every model function's result has a `model` key: `linear`,
  `generalized_linear`, `mixed_gaussian`, or `mixed_logistic`. `summary()` and
  `tidy()` choose the printer from it.

### Removed

- The ten 0.3 model names, deprecated since 0.4.0, as phit 0.5.0 removes
  them: `fit_lm`, `fit_glm`, `fit_lmer`, `fit_glmer`, `fit_logistic`,
  `fit_probit`, `fit_poisson`, `fit_mixed_intercept`, `fit_mixed`, and
  `fit_logistic_mixed`. Each replacement takes the same arguments, except that
  `fit_mixed_logistic()` has no `family` or `link`.
- The matrix form of `fit_ols()`. It takes only a formula, and its `offset`,
  `intercept`, `x`, and `y` parameters are gone. Call `fit_ols_matrix()` for a
  matrix.
- `FunctionMap::DEPRECATED`, which listed the aliases.

## [0.4.2] - 2026-11-01

Reproduces phit [v0.4.4](https://github.com/phpstats/phit/releases/tag/v0.4.4),
commit [`066e2e1`](https://github.com/phpstats/phit/commit/066e2e16a4b1db644e2700451b91ddd88cd635b1).
The distribution fixtures were rerecorded with the same cases.

### Fixed

- The listwise-deletion warning names the function once, as phit 0.4.4's does:
  `Phit\fit_ols(): dropped 1 of 10 observations with missing values (listwise
  deletion)`. The prefix is still the name the caller used.

## [0.4.1] - 2026-11-01

Reproduces phit [v0.4.3](https://github.com/phpstats/phit/releases/tag/v0.4.3),
commit [`db75800`](https://github.com/phpstats/phit/commit/db7580022223bdf67daa1c7fc2e4dbd8bf568d96).
phit v0.4.1 and v0.4.2 change only the CSV reader, which the polyfill leaves
to the extension. Its reference fixtures are unchanged from v0.4.0.

### Changed

- `fit_mixed_gaussian()` and `fit_mixed_gaussian_matrix()` factor single and
  nested grouping terms block by block, as phit 0.4.3 does. The random-effect
  columns that share rows form independent blocks, one per group or outermost
  group, and each likelihood evaluation factors the blocks separately. Crossed
  terms keep the dense path. The results are identical to phit 0.4.3's. A
  random intercept on 20,000 rows in 200 groups takes about 90 ms, where the
  dense path took 9 s. Single and nested designs are no longer limited to 512
  random coefficients.
- A single or nested design whose weighted rows overflow throws code 1003,
  "weighted design exceeds floating-point range", as the dense path does.

## [0.4.0] - 2026-11-01

Reproduces phit [v0.4.0](https://github.com/phpstats/phit/releases/tag/v0.4.0),
commit [`076f019`](https://github.com/phpstats/phit/commit/076f0197287f983399ecb80c6ffee48f0975d225).
Its numeric reference fixtures are unchanged from v0.3.0, and every result is
the same as in 0.3.1.

### Changed

- The model functions have phit 0.4's names. The formula fits have the short
  names, and the matrix fits end in `_matrix`.

  | 0.3 name | 0.4 name | Input |
  |---|---|---|
  | `fit_lm` | `fit_ols` | formula |
  | `fit_glm` | `fit_generalized_linear` | formula |
  | `fit_lmer` | `fit_mixed_gaussian` | formula |
  | `fit_glmer` | `fit_mixed_logistic` | formula |
  | `fit_ols` (matrix) | `fit_ols_matrix` | matrix |
  | `fit_logistic` | `fit_logistic_matrix` | matrix |
  | `fit_probit` | `fit_probit_matrix` | matrix |
  | `fit_poisson` | `fit_poisson_matrix` | matrix |
  | `fit_mixed` | `fit_mixed_gaussian_matrix` | matrix |
  | `fit_mixed_intercept` | `fit_mixed_gaussian_intercept_matrix` | matrix |
  | `fit_logistic_mixed` | `fit_mixed_logistic_matrix` | matrix |

- `Phit\fit_ols` chooses its form from its first argument. A formula string,
  or named `formula:` and `data:`, is the formula fit. An array, or named `x:`
  and `y:`, is the matrix fit, with a deprecation. Mixed or incomplete
  arguments throw an exception with code 1001.
- `Phit\fit_mixed_logistic` takes no `family` or `link` argument.
- Messages and warnings name the 0.4 functions. The listwise-deletion
  warning's prefix is the name the caller used, as for the extension.

### Deprecated

- The ten other 0.3 names are aliases. Each raises one `E_USER_DEPRECATED`
  naming its replacement and returns exactly the replacement's result.
  `Phit\fit_glmer` still accepts `family` and `link`, and it rejects any
  family other than binomial with the logit link with code 1007.
- `Phit\fit_ols` with a matrix. Call `Phit\fit_ols_matrix` instead.
- phit 0.5 will remove the aliases, and `Phit\fit_ols` will reject a matrix.

## [0.3.1] - 2026-10-31

Reproduces phit [v0.3.0](https://github.com/phpstats/phit/releases/tag/v0.3.0)
and [v0.3.1](https://github.com/phpstats/phit/releases/tag/v0.3.1). The pin
stays at v0.3.0's commit
[`2629fe1`](https://github.com/phpstats/phit/commit/2629fe105e8aeb4a73ac96313e3ff6c05e1df7ac).
v0.3.1, commit
[`96e29e1`](https://github.com/phpstats/phit/commit/96e29e108c1527d30d959f2ec2cd0df5c52d0585),
changes only phit's packaging and documentation. Its functions, API contract,
and reference fixtures are those of v0.3.0.

### Changed

- `composer.json` requires PHP `>=8.3 <8.6`, the versions CI tests. It
  previously allowed every later 8.x release.
- The package's homepage is <https://phit.phpstats.org>.
- The documentation states that when any version of the phit extension is
  loaded, the polyfill defines nothing, and that the CSV reader needs the
  extension.

## [0.3.0] - 2026-10-31

Reproduces phit [v0.3.0](https://github.com/phpstats/phit/releases/tag/v0.3.0),
commit [`2629fe1`](https://github.com/phpstats/phit/commit/2629fe105e8aeb4a73ac96313e3ff6c05e1df7ac).

### Added

- Seeded uniform, categorical, and permutation sampling:
  `Phit\dist_uniform_rng` (R's `runif()`), `Phit\dist_categorical_rng`
  (0-based indices drawn with given weights), `Phit\dist_categorical_rows_rng`
  (one index per row of weights), and `Phit\permutation_rng` (a partial
  Fisher–Yates shuffle, as R's `sample(n, size)`). They take one uniform per
  draw from the same ChaCha20 stream as the other samplers, so the same seed
  gives the extension's draws.

## [0.2.0] - 2026-10-30

Reproduces phit [v0.2.0](https://github.com/phpstats/phit/releases/tag/v0.2.0),
commit [`ac176ff`](https://github.com/phpstats/phit/commit/ac176ff6567a27d75191685ca7e8a6876786d548).

### Extension-only

phit 0.2.0's chunked CSV reader, the `Phit\CsvReader` class, is left
undefined and needs the native extension, as `Phit\read_csv` and its new
`block_size` option do. No polyfilled function changed.

## [0.1.0] - 2026-10-30

Reproduces phit [v0.1.1](https://github.com/phpstats/phit/releases/tag/v0.1.1),
commit [`f4887a4`](https://github.com/phpstats/phit/commit/f4887a4a293a3aaa129910a6006d691777ecf06f).
Requires PHP 8.3 or later.

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
  OPcache with `opcache.jit=disable`. PHP 8.3's function JIT with hot-counter
  triggers (`opcache.jit=1235`) also crashes PHP with a segmentation fault
  when the test suite runs the implementations directly, so that one CI job
  is informational; the `Phit\` functions refuse before any of that code runs.

[Unreleased]: https://github.com/phpstats/phit-polyfill/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.5.0
[0.4.2]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.4.2
[0.4.1]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.4.1
[0.4.0]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.4.0
[0.3.1]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.3.1
[0.3.0]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.3.0
[0.2.0]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.2.0
[0.1.0]: https://github.com/phpstats/phit-polyfill/releases/tag/v0.1.0
