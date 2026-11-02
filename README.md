# phit-polyfill

A pure-PHP implementation of the [phit](https://github.com/phpstats/phit)
extension's functions, except the CSV reader (`Phit\read_csv` and
`Phit\CsvReader`), which needs the extension. When the extension is not
loaded, the polyfill defines them in the `Phit` namespace with the
extension's names, parameter names, defaults, result keys, warnings, and
exception codes. Code written against
phit runs on hosts where the extension cannot be installed. When the extension
is loaded, the polyfill defines nothing.

The implementations follow phit's Rust code operation for operation. On every
input in the test suites they return the same bits as the extension. The documentation is in
[`site/docs`](site/docs): [overview](site/docs/index.qmd),
[functions](site/docs/functions.qmd),
[compatibility and limitations](site/docs/php-versions.qmd),
[implementation](site/docs/how-it-works.qmd), and
[development](site/docs/development.qmd).

## Features

The polyfill implements 234 of phit 0.5.1's 235 functions.

| Area | Functions | Polyfill |
|---|---|---|
| Descriptive statistics | `stat_mean`, `stat_var`, `stat_sd`, `stat_median`, `stat_quantile`, `stat_quantile_array`, `stat_cov`, `stat_cor`, `stat_weighted_mean`, `stat_weighted_var` | Provided |
| Smoothing | `stat_density`, `stat_bw_nrd0`, `stat_bw_nrd`, `stat_approx`, `stat_ecdf` | Provided |
| Optimization | `optimize`, `uniroot` | Provided |
| Models from a formula | `fit_ols`, `fit_generalized_linear`, `fit_mixed_gaussian`, `fit_mixed_logistic`, `model_matrix`, `predict` | Provided |
| Models from a matrix | `fit_ols_matrix`, `fit_logistic_matrix`, `fit_probit_matrix`, `fit_poisson_matrix`, `fit_mixed_gaussian_matrix`, `fit_mixed_gaussian_intercept_matrix`, `fit_mixed_logistic_matrix` | Provided |
| Model summaries | `summary`, `tidy` | Provided |
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

This script runs in a new Composer project after
`composer require phpstats/phit-polyfill`. Save it as `fit.php` there and run
`php fit.php`. It reads mtcars from the
[PHPstats datasets repository](https://github.com/phpstats/datasets) and
prints the fit with `summary()`.

```php
<?php

require __DIR__.'/vendor/autoload.php';

use function Phit\{fit_ols, summary};

// mtcars from the PHPstats datasets repository, read into columns.
$url = 'https://raw.githubusercontent.com/phpstats/datasets/main/csv/'
    .'mtcars.csv';
$rows = array_map(
    fn (string $line) => str_getcsv($line, escape: ''),
    file($url, FILE_IGNORE_NEW_LINES),
);
$cars = [];
foreach (array_shift($rows) as $j => $name) {
    $cars[$name] = array_map(
        fn (string $v) => is_numeric($v) ? (float) $v : $v,
        array_column($rows, $j),
    );
}

$fit = fit_ols('mpg ~ wt + factor(cyl)', $cars);
echo summary($fit);
```

```text
Linear model fit by least squares
Formula: mpg ~ wt + factor(cyl)

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept)   33.9908     1.8878  18.006  < 2e-16
wt            -3.2056     0.7539  -4.252 0.000213
factor(cyl)6  -4.2556     1.3861  -3.070 0.004718
factor(cyl)8  -6.0709     1.6523  -3.674 0.000999

Residual standard error: 2.557 on 28 degrees of freedom
Multiple R-squared: 0.8374, Adjusted R-squared: 0.82
F-statistic: 48.08 on 3 and 28 DF, p-value: 3.594e-11
```

The same code runs with the extension loaded. Without the extension,
`Phit\read_csv` is undefined, so the script reads the CSV file itself.

`Phit\Polyfill\Polyfill` reports what is available:

```php
use Phit\Polyfill\Polyfill;

Polyfill::supports('Phit\fit_generalized_linear'); // true
Polyfill::supports('Phit\read_csv');               // false: extension only
Polyfill::nativeLoaded();                          // whether phit is loaded
```

## Extension

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

## Versions

Each release reproduces the phit commit pinned in `tests/fixtures/REFERENCE`,
and the reference fixtures are copied from it. phit v0.3.1 changes only
packaging and documentation, so it has the same functions as v0.3.0. phit
v0.4.1 and v0.4.2 change only the CSV reader, which the polyfill leaves to the
extension.

| Polyfill | phit |
|---|---|
| 0.5.1 | [v0.5.1](https://github.com/phpstats/phit/releases/tag/v0.5.1) (`c60fb9b`) |
| 0.5.0 | [v0.5.0](https://github.com/phpstats/phit/releases/tag/v0.5.0) (`17729a8`) |
| 0.4.2 | [v0.4.4](https://github.com/phpstats/phit/releases/tag/v0.4.4) (`066e2e1`) |
| 0.4.1 | [v0.4.3](https://github.com/phpstats/phit/releases/tag/v0.4.3) (`db75800`) |
| 0.4.0 | [v0.4.0](https://github.com/phpstats/phit/releases/tag/v0.4.0) (`076f019`) |
| 0.3.1 | [v0.3.0](https://github.com/phpstats/phit/releases/tag/v0.3.0) (`2629fe1`) and [v0.3.1](https://github.com/phpstats/phit/releases/tag/v0.3.1) (`96e29e1`) |
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
| `Phit\fit_ols_matrix` | 10,000 | 26 ms | 2.9 ms | 8.9× |
| `Phit\fit_ols_matrix` | 100,000 | 262 ms | 29 ms | 9.1× |
| `Phit\fit_logistic_matrix` | 10,000 | 160 ms | 7.0 ms | 23× |
| `Phit\fit_logistic_matrix` | 100,000 | 1.60 s | 72 ms | 22× |

The logistic fit includes the null model that every intercept model fits for
`null_deviance`. The coefficients were identical in every case. In the recorded mixed-model
cases the polyfill took 13 to 78 times as long as the extension. The cost of
`fit_mixed_gaussian_matrix` and `fit_mixed_gaussian` grows roughly with the cube of the number of random
coefficients. One fit with 120 random coefficients took 8.0 s against the
extension's 103 ms ([ACCURACY.md](ACCURACY.md#speed)).

The recorded run is in [benchmarks/results.json](benchmarks/results.json),
with its full context. Reproduce it with:

```bash
php -d extension=/path/to/libphit.so \
    -d opcache.enable_cli=1 -d opcache.jit=disable \
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
- `tests/Unit` holds the extension's test suites, rewritten for Pest. They call the
  polyfill directly, so they test it whether or not the extension is loaded.
  `tests/Parity` runs only with the extension and compares the two.

## License

Released under the MIT license. See [LICENSE](LICENSE).
