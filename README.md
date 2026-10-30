# phit-polyfill

A pure-PHP polyfill for the [phit](https://github.com/phpstats/phit)
native extension. It defines the extension's functions in the `Phit` namespace
only when the extension is not loaded, with the same names, parameter names,
order and defaults, result keys, warnings, and exception codes. Code written
against phit therefore runs unchanged on hosts where the extension cannot be
installed. When the extension is loaded, it is used and the polyfill does
nothing.

The package is part of the [phpstats](https://github.com/phpstats) family, alongside
phibble, phplyr, phpolars, and cmdstanphp. It follows the pattern of
`symfony/polyfill-mbstring`.

## Versions

Each polyfill release reproduces one phit commit, which is pinned in
`tests/fixtures/REFERENCE`. The reference fixtures are copied from it.

| Polyfill | phit | Notes |
|---|---|---|
| 0.1.0 | [`1bf73bd`](https://github.com/phpstats/phit/commit/1bf73bd3bdf66f825039713ac58157d0f9d6ddf9) | phit has no tagged release yet |

See [CHANGELOG.md](CHANGELOG.md) for what each release contains.

## Installation

```bash
composer require phpstats/phit-polyfill
```

PHP 8.3 or later is required. Composer loads `bootstrap.php` automatically. For
each function the polyfill implements, it declares the function in the `Phit`
namespace only if `function_exists()` reports that it is missing. Import a
function with `use function Phit\fit_ols;` or call it as `Phit\fit_ols()`.

```php
use function Phit\fit_ols;

$fit = fit_ols(
    x: [[1.2, 0.4], [2.3, 1.1], [3.1, 0.9], [4.8, 2.2], [5.0, 1.7]],
    y: [2.1, 4.2, 5.8, 8.1, 10.2],
);

$fit['coefficients'];  // intercept first, as with the extension
$fit['p_values'];
```

To check at runtime whether a function is available without the extension,
use `Phit\Polyfill\Polyfill`:

```php
use Phit\Polyfill\Polyfill;

Polyfill::supports('Phit\fit_logistic');  // true
Polyfill::supports('Phit\read_csv');      // false: needs the extension
Polyfill::nativeLoaded();                 // whether phit is loaded
```

## Coverage

| Area | Functions | Polyfill |
|---|---|---|
| Descriptive statistics | `Phit\stat_mean`, `Phit\stat_var`, `Phit\stat_sd`, `Phit\stat_median`, `Phit\stat_quantile`, `Phit\stat_quantile_array`, `Phit\stat_cov`, `Phit\stat_cor`, `Phit\stat_weighted_mean`, `Phit\stat_weighted_var` | Provided |
| Smoothing | `Phit\stat_density`, `Phit\stat_bw_nrd0`, `Phit\stat_bw_nrd`, `Phit\stat_approx`, `Phit\stat_ecdf` | Provided |
| Optimization | `Phit\optimize`, `Phit\uniroot` | Provided |
| Regression | `Phit\fit_ols`, `Phit\fit_logistic`, `Phit\fit_probit`, `Phit\fit_poisson`, with weights, offsets, the automatic intercept, and both missing-value policies | Provided |
| Formula interface | `Phit\model_matrix`, `Phit\fit_lm`, `Phit\fit_glm`, `Phit\predict` | Provided |
| Distributions | 13 families × 14 functions; see [Distributions](#distributions) | Provided |
| Mixed models | `Phit\fit_mixed_intercept`, `Phit\fit_mixed`, `Phit\fit_logistic_mixed`, `Phit\fit_lmer`, `Phit\fit_glmer` | Provided (slower than the extension; see [ACCURACY.md](ACCURACY.md#speed)) |
| Random numbers | 13 seeded `*_rng` functions, one per distribution family | Provided |
| Posterior draws | `Phit\draws_moments`, `Phit\draws_covariance_map`, `Phit\draws_affine` | Provided |
| CSV reading | `Phit\read_csv` (phibble's PHP reader is the fallback) | Extension only |

Functions marked "Extension only" stay undefined without the extension, so
`function_exists()` reports them correctly.

### Distributions

All 13 distribution families are provided, each with the full set of
operations and their `_array` variants:

| Family | Prefix | Parameters | Operations |
|---|---|---|---|
| Normal | `dist_normal_` | mean, sd | pdf, logpdf, cdf, logcdf, sf, logsf, quantile |
| Student's t | `dist_t_` | df | same |
| Chi-square | `dist_chisquare_` | df | same |
| F | `dist_f_` | df1, df2 | same |
| Exponential | `dist_exponential_` | rate | same |
| Logistic | `dist_logistic_` | location, scale | same |
| Log-normal | `dist_lognormal_` | meanlog, sdlog | same |
| Gamma | `dist_gamma_` | shape, scale | same |
| Beta | `dist_beta_` | alpha, beta | same |
| Binomial | `dist_binomial_` | size, prob | pmf, logpmf, cdf, logcdf, sf, logsf, quantile |
| Poisson | `dist_poisson_` | lambda | same |
| Negative binomial | `dist_nbinom_` | size, prob | same |
| Negative binomial by mean | `dist_nbinom_mu_` | size, mu | same |

For example, `Phit\dist_t_sf(2.5, df: 12)` and
`Phit\dist_normal_quantile_array([0.025, 0.975])`.

Each family also has a seeded sampler, such as
`Phit\dist_normal_rng(1000, seed: 42)`. The polyfill uses the extension's generator (ChaCha20 seeded as in rand_core,
with inverse-CDF sampling), so a seed gives the same draws, bit for bit, with
or without the extension.

### Formula interface

`Phit\model_matrix`, `Phit\fit_lm`, and `Phit\fit_glm` accept an R-style
formula over named columns. They port the extension's parser, term expansion, and
coding, so column names, term order, contrasts, and missing-value handling
match the extension. Where the extension differs from R, the polyfill follows
the extension (see the phit formula documentation).

```php
$data = [
    'rent' => [1879.25, 1805.03, 2022.29, 1267.85, 1509.24, 2240.10],
    'area' => [74.4, 57.4, 66.4, 49.1, 45.7, 92.4],
    'district' => ['north', 'south', 'north', 'east', 'east', 'south'],
];

$fit = Phit\fit_lm('rent ~ log(area) + district', $data);
$fit['named'];    // coefficients keyed by column name

$design = Phit\model_matrix('~ area * district', $data);
$design['columns'];

$logit = Phit\fit_glm('y ~ x + offset(log(t))', $other, family: 'binomial');

// Predictions for new data, as R's predict(): fitted factor levels and poly()
// bases are reused, and rows with a missing predictor predict null.
$new = Phit\predict($fit, ['area' => [60.0, 80.0], 'district' => ['east', 'north']],
    interval: 'prediction');
```

Supported syntax includes `+ - : * / %in% ^`, `0`/`1`/`-1`, `.`, `offset()`,
`I()`, the common transforms (`log`, `log2`, `log10`, `log1p`, `exp`, `sqrt`,
`abs`, and others), `poly()`, `factor()`, `ordered()`, and `cbind()` responses
for grouped binomial data. Families are gaussian, binomial (logit or probit),
and poisson. Formulas with random-effect bars such as `(1 | g)` are accepted
by `Phit\model_matrix`, which returns the random-effect design, and are fitted
by `Phit\fit_lmer` and `Phit\fit_glmer`. `Phit\fit_lm` and `Phit\fit_glm` reject
them with code 1007. `Phit\predict` works with all four fit functions.

## How it works

The implementations are ports of the extension's Rust code, not independent
reimplementations, so results agree to rounding in most cases.

- **Least squares.** `Phit\fit_ols` scales the design columns to unit norm and
  factors the design by Householder QR, using the same operation order as the
  extension's linear-algebra library. The singular values of the small R
  factor (computed by a Jacobi sweep) equal the design's, so the extension's
  rank rule (s ≤ ε·max(n, p)·s_max) makes the same rank decisions.
  Coefficients and covariances come from R by back-substitution. Weighted fits
  factor rows in decreasing-weight order, as the extension does.
- **Generalized linear models.** `Phit\fit_logistic` and `Phit\fit_poisson` use
  damped Newton IRLS and `Phit\fit_probit` uses Fisher scoring. Each iteration
  solves a weighted least-squares step with the same QR, halves the step
  until the likelihood does not decrease, and stops only after a small full
  step. The separation certificate, boundary handling, starting values for
  extreme offsets, and exception codes match the extension.
- **Distributions.** The extension's distribution layer is ported together
  with the parts of cdflib it calls: TOMS 708 (incomplete beta), TOMS 654
  (incomplete gamma), Cody's normal CDF, and the `dinvnr`/`dt1` quantile
  starts. Log-scale tails, saddle-point densities, and the safeguarded quantile
  search are ported as well.
- **Special functions.** The statrs routines that the extension's estimators
  call (erfc, ln Γ, the regularized incomplete beta, and the inverse error
  function) are ported from the same crate version.
- **Validation.** Inputs are validated in the extension's order, so an input
  that is invalid in several ways raises the same exception code.

Exceptions are plain `\Exception` instances whose code identifies the category,
as in the extension:

| Code | Meaning |
|---:|---|
| 1001 | Invalid input |
| 1002 | Dimension mismatch |
| 1003 | Non-finite value or overflow |
| 1004 | Rank deficiency |
| 1005 | Factorization failure |
| 1006 | Non-convergence (including separation) |
| 1007 | Unsupported input size |
| 1008 | Internal error |

## Accuracy

Every provided function passes the extension's R, SciPy, and mpmath fixtures
at the extension's own tolerances, and a parity suite that compares the
polyfill with the loaded extension. [ACCURACY.md](ACCURACY.md) lists, per
function, the largest observed error against each reference and the tested
domain.

In summary:

- Descriptive statistics agree with the extension exactly on every input
  tested, including values offset by 10¹².
- Kernel density, the bandwidths, `Phit\stat_approx`, and `Phit\stat_ecdf`
  are bit-identical to the extension and within 10⁻¹² of R 4.6.1.
- `Phit\optimize` and `Phit\uniroot` evaluate f at the same points as the
  extension and as R 4.6.1, and return the same bits.
- OLS, logistic, probit, and Poisson regression return the same bits as the
  extension on every tested input, including designs with condition numbers
  up to 10¹⁰, extreme prior weights, separation, and extreme offsets.
- Distribution functions are bit-identical to the extension on every input
  tested, and within 10⁻¹² of mpmath down to probabilities near 10⁻³⁰⁰.
- The formula functions produce identical design matrices, column names, and
  terms, and identical fits.
- Seeded samplers, posterior-draw kernels, and the five mixed-model functions
  are bit-identical to the extension in the randomized parity suite.

### Known differences

- The listwise-deletion warning uses `E_USER_WARNING` because PHP code cannot
  raise `E_WARNING`. The text is identical.
- Scalar arguments of the wrong type follow PHP's userland rules (coercion in
  non-strict files, otherwise `TypeError`). The extension raises an
  `\Exception` instead.
- Every `Phit\` function the polyfill defines throws an exception with code
  1007 while OPcache's JIT compiles the process, on every PHP version. PHP's
  JIT has computed wrong results in the polyfill without an error on PHP 8.3,
  8.4, and 8.5: wrong beta probabilities on 8.3, a mixed-model optimizer moved
  to a minimum 18% away from the correct one on 8.4, and wrong singular values
  under a function JIT with hot-counter triggers (php-src #24214). Run OPcache
  with `opcache.jit=disable`, or install the extension, which the JIT does not
  affect. `Phit\Polyfill\Polyfill::implementation()` returns the
  implementation without the check, for tests and comparisons. Under PHP
  8.3's function JIT with hot-counter triggers (`opcache.jit=1235`), the test
  suite's direct calls crash PHP with a segmentation fault, so that CI job is
  informational. See
  [docs/notes/function-jit.md](docs/notes/function-jit.md) and
  [docs/notes/php83-jit.md](docs/notes/php83-jit.md).
- Matrix products above 5×5 follow the extension's aarch64 build, which uses
  fused multiply-adds. On an x86-64 CPU without FMA the extension's products
  can differ in the last bit.

## Performance

Benchmarks compare the polyfill with the extension in the same process on an
Apple M5 Ultra with PHP 8.5.10. Each design has 10 standard normal predictors
plus an intercept. Times are the median of five runs with OPcache on and the
JIT off (`-d opcache.enable_cli=1 -d opcache.jit=disable`), the configuration
recommended above.

| Model | n | Polyfill | Extension | Ratio |
|---|---:|---:|---:|---:|
| `Phit\fit_ols` | 10,000 | 26 ms | 2.9 ms | 8.9× |
| `Phit\fit_ols` | 100,000 | 262 ms | 29 ms | 9.1× |
| `Phit\fit_logistic` | 10,000 | 160 ms | 7.0 ms | 23× |
| `Phit\fit_logistic` | 100,000 | 1.60 s | 72 ms | 22× |

The logistic fit includes the null model that every intercept model fits for
`null_deviance`.

For comparison, solving the normal equations (X'X)b = X'y in plain PHP took
87 ms at n = 100,000. That baseline returns coefficients only, without
standard errors, rank detection, or the other results. Forming X'X also
squares the condition number, so it loses accuracy on ill-conditioned designs.

The polyfill's coefficients were identical to the extension's in every case.
The recorded run is in [benchmarks/results.json](benchmarks/results.json),
with its full context: CPU, platform, PHP and JIT settings, median and minimum
times, the plain-PHP baseline for OLS, and the agreement check. Reproduce it
with:

```bash
php -d extension=/path/to/libphit.so -d opcache.enable_cli=1 -d opcache.jit=disable \
    benchmarks/run.php --sizes=10000,100000 --repeat=5 --seed=5 \
    --json=benchmarks/results.json
```

## Development

```bash
composer install
composer test        # polyfill only
composer test:ext    # with the extension loaded from .ext/libphit.dylib
```

The composer scripts run the tests with the PHP that runs Composer, so
`php85 "$(which composer)" test` picks a version. Any PHP 8.3 or later works.
`composer test:ext` needs the PHP that `.ext/libphit.dylib` was built for;
under another one the extension does not load and the parity tests skip.

- `bin/fetch-fixtures [checkout]` copies the reference fixtures from the
  phit commit pinned in `tests/fixtures/REFERENCE`. The copies are
  committed and never edited by hand.
- `bin/generate-bootstrap` regenerates `bootstrap.php` and
  `src/FunctionMap.php` from the extension's API contract
  (`tests/fixtures/php-api.json`). A function is declared only when its
  implementation exists.
- `bin/accuracy-report` prints the tables in ACCURACY.md.
- `site/` holds the content of the package site: `site.json`, the docs
  (`site/docs/*.qmd`, rendered next to each source by `bin/render-site-docs`),
  the vignettes, and the benchmark and verification files.
  `benchmarks/site.php` writes the benchmarks. `bin/site-data` writes the
  verifications from a parity log (`PARITY_LOG=… vendor/bin/pest tests/Parity`)
  and the records of `bin/site-matrix`, which runs the suite under each
  OPcache and JIT mode.
- `tests/Unit` holds the ports of the extension's suites. They call the
  polyfill directly, so they test it whether or not the extension is loaded.
  `tests/Parity` runs only with the extension and compares the two.

## License

MIT. See [LICENSE](LICENSE).
