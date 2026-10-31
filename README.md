# phit-polyfill

A pure-PHP implementation of the [phit](https://github.com/phpstats/phit)
extension's functions. When the extension is not loaded, the polyfill defines
them in the `Phit` namespace with the extension's names, parameter names,
defaults, result keys, warnings, and exception codes. Code written against
phit runs on hosts where the extension cannot be installed. When the extension
is loaded, the polyfill defines nothing.

The implementations are ports of phit's Rust code. On every input in the test
suites they return the same bits as the extension. The documentation is in
[`site/docs`](site/docs): [overview](site/docs/index.qmd),
[functions](site/docs/functions.qmd),
[compatibility and limitations](site/docs/php-versions.qmd),
[how it works](site/docs/how-it-works.qmd), and
[development](site/docs/development.qmd).

## Installation

```bash
composer require phpstats/phit-polyfill
```

The polyfill supports PHP 8.3, 8.4, and 8.5 (`composer.json` requires
`>=8.3 <8.6`). CI tests all three on Linux and macOS. No PHP extension is
required.

Composer's autoloader loads `bootstrap.php`, which defines the functions. The
polyfill's functions refuse to run while OPcache's JIT compiles the process.
They throw an exception with code 1007. Set `opcache.jit=disable` where the
polyfill runs. See [The JIT](site/docs/php-versions.qmd#the-jit).

## Usage

Import a function with `use function`, or call it by its qualified name.

```php
use function Phit\fit_ols;

$fit = fit_ols(
    x: [[1.2, 0.4], [2.3, 1.1], [3.1, 0.9], [4.8, 2.2], [5.0, 1.7]],
    y: [2.1, 4.2, 5.8, 8.1, 10.2],
);

$fit['coefficients'];
$fit['p_values'];
```

`$fit['coefficients']` is `[-0.2843, 2.5864, -1.6819]` (rounded here), with
the intercept first. The same code runs with the extension loaded.

`Phit\Polyfill\Polyfill` reports what is available:

```php
use Phit\Polyfill\Polyfill;

Polyfill::supports('Phit\fit_logistic');  // true
Polyfill::supports('Phit\read_csv');      // false: needs the extension
Polyfill::nativeLoaded();                 // whether phit is loaded
```

## With the extension

- When any version of the `phit` extension is loaded, `bootstrap.php` returns
  before it defines anything, and every `Phit\` function is the extension's.
  An extension older than the phit release the polyfill reproduces lacks the
  newer functions, and the polyfill does not define them. Upgrade the
  extension to that release to use them (see [Versions](#versions)).
- `composer.json` declares `"provide": {"ext-phit": "*"}`, so Composer treats
  a package that requires `ext-phit` as satisfied by the polyfill. The CSV
  reader, `Phit\read_csv` and `Phit\CsvReader`, needs the extension.
- `Polyfill::implementation('Phit\fit_ols')` returns the polyfill's callable
  even when the extension is loaded. It skips the JIT check. The parity tests
  and benchmarks use it to run both implementations in one process.

## Coverage

The polyfill implements 232 of phit 0.3.0's 233 functions.

| Area | Functions | Polyfill |
|---|---|---|
| Descriptive statistics | `stat_mean`, `stat_var`, `stat_sd`, `stat_median`, `stat_quantile`, `stat_quantile_array`, `stat_cov`, `stat_cor`, `stat_weighted_mean`, `stat_weighted_var` | Provided |
| Smoothing | `stat_density`, `stat_bw_nrd0`, `stat_bw_nrd`, `stat_approx`, `stat_ecdf` | Provided |
| Optimization | `optimize`, `uniroot` | Provided |
| Regression | `fit_ols`, `fit_logistic`, `fit_probit`, `fit_poisson` | Provided |
| Formula interface | `model_matrix`, `fit_lm`, `fit_glm`, `predict` | Provided |
| Mixed models | `fit_mixed_intercept`, `fit_mixed`, `fit_logistic_mixed`, `fit_lmer`, `fit_glmer` | Provided |
| Distributions | 13 families with 14 functions each | Provided |
| Random numbers | `dist_<family>_rng` for the 13 families, `dist_uniform_rng`, `dist_categorical_rng`, `dist_categorical_rows_rng`, `permutation_rng` | Provided |
| Posterior draws | `draws_moments`, `draws_covariance_map`, `draws_affine` | Provided |
| CSV reading | `read_csv` and the `CsvReader` class | Extension only |

All names are in the `Phit` namespace. The [functions page](site/docs/functions.qmd)
lists the distribution families, the formula syntax, and the mixed-model
structures. Functions and classes marked "Extension only" stay undefined
without the extension, so `function_exists()` and `class_exists()` report
them correctly. phibble's PHP CSV reader covers CSV input without the
extension.

## Versions

Each release reproduces one phit release. Its commit is pinned in
`tests/fixtures/REFERENCE`, and the reference fixtures are copied from it.

| Polyfill | phit |
|---|---|
| 0.3.0 | [v0.3.0](https://github.com/phpstats/phit/releases/tag/v0.3.0) (`2629fe1`) |
| 0.2.0 | [v0.2.0](https://github.com/phpstats/phit/releases/tag/v0.2.0) (`ac176ff`) |
| 0.1.0 | [v0.1.1](https://github.com/phpstats/phit/releases/tag/v0.1.1) (`f4887a4`) |

[CHANGELOG.md](CHANGELOG.md) lists what each release contains.

## Differences from the extension

- The listwise-deletion warning uses `E_USER_WARNING` because PHP code cannot
  raise `E_WARNING`. The text is identical.
- Scalar arguments of the wrong type follow PHP's userland rules: coercion in
  non-strict files, otherwise `TypeError`. The extension raises an
  `\Exception`.
- The polyfill's functions throw an exception with code 1007 while OPcache's
  JIT compiles the process. The JIT has computed wrong results in the
  polyfill without an error on PHP 8.3, 8.4, and 8.5. The extension runs
  under any JIT setting.
- Matrix products above 5×5 follow the extension's aarch64 build, which uses
  fused multiply-adds. On an x86-64 CPU without FMA, the extension's products
  can differ from the polyfill's in the last bit.

[ACCURACY.md](ACCURACY.md) lists, per function, the largest observed error
against the R, SciPy, and mpmath references and the tested domain.

## Performance

These timings were recorded on 2026-10-28 on an Apple M5 Ultra with PHP
8.5.10, OPcache on, and the JIT off. They compare polyfill commit `ab803fd`
with phit commit `3861d77` in the same process. Each design has 10 standard
normal predictors plus an intercept. Times are the median of five runs.

| Model | n | Polyfill | Extension | Ratio |
|---|---:|---:|---:|---:|
| `Phit\fit_ols` | 10,000 | 26 ms | 2.9 ms | 8.9× |
| `Phit\fit_ols` | 100,000 | 262 ms | 29 ms | 9.1× |
| `Phit\fit_logistic` | 10,000 | 160 ms | 7.0 ms | 23× |
| `Phit\fit_logistic` | 100,000 | 1.60 s | 72 ms | 22× |

The logistic fit includes the null model that every intercept model fits for
`null_deviance`. The coefficients were identical in every case. In the recorded mixed-model
cases the polyfill took 13 to 78 times as long as the extension. The cost of
`fit_mixed` and `fit_lmer` grows roughly with the cube of the number of random
coefficients. One fit with 120 random coefficients took 8.0 s against the
extension's 103 ms ([ACCURACY.md](ACCURACY.md#speed)).

The recorded run is in [benchmarks/results.json](benchmarks/results.json),
with its full context. Reproduce it with:

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
`php85 "$(which composer)" test` picks a version. PHP 8.3, 8.4, and 8.5 work.
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
- `site/` holds `site.json`, the docs (`site/docs/*.qmd`, rendered next to
  each source by `bin/render-site-docs`), and the examples in
  `site/vignettes`. `benchmarks/site.php` writes benchmark files. `bin/site-data` writes the
  verifications from a parity log (`PARITY_LOG=… vendor/bin/pest tests/Parity`)
  and the records of `bin/site-matrix`, which runs the suite under each
  OPcache and JIT mode.
- `tests/Unit` holds the ports of the extension's suites. They call the
  polyfill directly, so they test it whether or not the extension is loaded.
  `tests/Parity` runs only with the extension and compares the two.

## License

MIT. See [LICENSE](LICENSE).
