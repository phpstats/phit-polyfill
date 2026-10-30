# Accuracy

This document records how closely the polyfill reproduces two references. The
first is the R, SciPy, and mpmath values that the phit extension is
validated against. The second is the extension itself. Every number below was
measured, not estimated. A function is shipped only when it passes the
extension's own fixture tolerances and the parity suite. Anything that cannot
meet them stays undefined and is listed under [Not provided](#not-provided).

Reference: phit commit pinned in `tests/fixtures/REFERENCE`. The fixtures
are copied from that commit by `bin/fetch-fixtures` and never edited.

## How the numbers are produced

- `bin/accuracy-report` runs the polyfill on every case in the R fixtures and
  prints the tables in [Against R](#against-r). Each comparison uses the
  tolerance of the extension's own suite for that fixture, |a − b| ≤ abs + rel·|b|.
  "Tolerance used" is the largest |a − b| / (abs + rel·|b|) over all entries, so
  100% would be a failure. "Max relative error" is |a − b| / |b| over entries
  that are not rounding noise. It skips entries whose expected magnitude is
  below the absolute tolerance or below 10⁻³ of the largest value in the same
  result entry, such as the residuals of an exact fit.
- The parity suites in `tests/Parity` compare the polyfill with the loaded
  extension on the fixtures and on seeded random and adversarial inputs. Run
  them with `PARITY_REPORT=parity.json` to record the largest differences.
- Both are part of the test suite: `composer test` (polyfill only) and
  `composer test:ext` (with the extension loaded).

## Descriptive statistics

The descriptive functions are ports of the extension's algorithms, including
its two-pass and compensated summations. Against the extension they agree
exactly (largest difference 0.0) on the fixtures and on seeded random inputs.
The random inputs include values offset by 10⁹ and 10¹² (catastrophic
cancellation for a naive variance), integers, tiny values, n from 1 to 500,
and null/NaN under both missing-value policies. Exception codes are identical
for every invalid input tested.

The smoothing functions (`Phit\stat_density`, `Phit\stat_bw_nrd0`,
`Phit\stat_bw_nrd`, `Phit\stat_approx`, `Phit\stat_ecdf`) port the extension's
binning, radix-2 FFT, and interpolation step for step. They agree with the
extension to the bit on the R fixtures, on seeded random inputs (every
kernel, fixed and rule-of-thumb bandwidths, weights, grids above 512 points,
ties, signed zeros, and null, NaN, and infinite evaluation points), and on
every invalid input tested. Against R 4.6.1 (smoothing-r.json) they agree to
a relative 10⁻¹², with density heights measured against the peak.

## Regression models

`Phit\fit_ols` solves least squares by Householder QR of the column-scaled design,
then a one-sided Jacobi SVD of the small R factor. The singular values are
those of the design, so the extension's rank rule (s ≤ ε·max(n, p)·s_max)
makes the same decisions. `Phit\fit_logistic`, `Phit\fit_probit`, and
`Phit\fit_poisson` port the extension's damped Newton / Fisher-scoring loop step
for step. This includes the starting values, line search, separation
certificate, convergence rule, compensated likelihood sums, and statistics. The
statrs special functions behind p-values and likelihood constants (erfc, ln Γ,
the regularized incomplete beta, and the inverse error function) are ported
from the same crate version.

### Against R

| Function | Reference | Comparisons | Max relative error | Tolerance used | Tested domain |
|---|---|---:|---:|---:|---|
| `Phit\stat_cor` | R (descriptive-r.json) | 3 | 1.3e-16 | 0% | correlated, negative, offset |
| `Phit\stat_cov` | R (descriptive-r.json) | 6 | 2.1e-16 | 0% | correlated, negative, offset |
| `Phit\fit_logistic` | R (logistic-r.json) | 78 | 7.4e-9 | 6% | n=80, p=3, n=80, p=2, n=80, p=1 |
| `Phit\fit_logistic` | R (r-compatibility.json) | 52 | 3.5e-6 | 20% | logistic_intercept_individual, logistic_intercept_grouped, logistic_no_intercept_individual, logistic_no_intercept_grouped |
| `Phit\fit_logistic` | R (stress-r.json) | 44 | 2.5e-12 | 0% | logistic_rare, logistic_offset_-1000, logistic_offset_1000, logistic_extreme_known_predictors |
| `Phit\stat_mean` | R (descriptive-r.json) | 8 | 1.5e-16 | 0% | small, integers, offset, tiny, and 4 more |
| `Phit\stat_median` | R (descriptive-r.json) | 8 | 0.0e+0 | 0% | small, integers, offset, tiny, and 4 more |
| `Phit\fit_ols` | R (ols-r.json) | 108 | 3.5e-13 | 0% | n=8, p=2, n=8, p=3, n=5, p=2, n=10, p=3, and 1 more |
| `Phit\fit_ols` | R (r-compatibility.json) | 30 | 3.8e-12 | 0% | ols_intercept_individual, ols_no_intercept_individual |
| `Phit\fit_ols` | R (stress-r.json) | 54 | 1.6e-8 | 0% | nearly_collinear, poorly_scaled, small_outcome, large_outcome, and 2 more |
| `Phit\fit_poisson` | R (poisson-r.json) | 78 | 8.8e-11 | 0% | n=120, p=3, n=120, p=2, n=120, p=1 |
| `Phit\fit_poisson` | R (r-compatibility.json) | 26 | 8.0e-7 | 39% | poisson_intercept_individual, poisson_no_intercept_individual |
| `Phit\fit_poisson` | R (stress-r.json) | 55 | 1.7e-8 | 21% | poisson_tiny_known_means, poisson_sparse, poisson_large, poisson_offset_-1000, and 1 more |
| `Phit\fit_probit` | R (probit-r.json) | 78 | 6.3e-8 | 17% | n=120, p=3, n=120, p=2, n=120, p=1 |
| `Phit\fit_probit` | R (r-compatibility.json) | 52 | 1.8e-5 | 20% | probit_intercept_individual, probit_intercept_grouped, probit_no_intercept_individual, probit_no_intercept_grouped |
| `Phit\fit_probit` | R (stress-r.json) | 44 | 7.4e-11 | 0% | probit_rare, probit_offset_-1000, probit_offset_1000, probit_extreme_known_predictors |
| `Phit\stat_quantile_array` | R (descriptive-r.json) | 8 | 0.0e+0 | 0% | small, integers, offset, tiny, and 4 more |
| `Phit\stat_quantile` | R (descriptive-r.json) | 80 | 0.0e+0 | 0% | small, integers, offset, tiny, and 4 more |
| `Phit\stat_sd` | R (descriptive-r.json) | 8 | 3.9e-16 | 0% | small, integers, offset, tiny, and 4 more |
| `Phit\stat_var` | R (descriptive-r.json) | 16 | 5.2e-16 | 0% | small, integers, offset, tiny, and 4 more |
| `Phit\stat_weighted_mean` | R (descriptive-r.json) | 2 | 2.1e-16 | 0% | analytic, counts |
| `Phit\stat_weighted_var` | R (descriptive-r.json) | 6 | 1.9e-16 | 0% | analytic, counts |

The GLM differences from R (up to about 2·10⁻⁵ relative in the probit
r-compatibility case) come from R's own stopping rule (`glm` stops when the
relative deviance change is below 10⁻⁸). The extension shows the same
differences. Descriptive statistics are compared at a relative tolerance of
10⁻¹².

### Against the extension

The parity suite compares every result key. Each value is judged on the scale
of its quantity. Coefficients and intervals are compared relative to
|b| + standard error, t/z statistics and R² relative to |b| + 1, p-values on the
log scale, and fitted values and residuals relative to |b| + max |y|. Other keys
are compared relative to |b| + 10⁻³·max |b|.

The least-squares steps reproduce the extension's arithmetic: Householder QR
with nalgebra's operation order (including its eight-way blocked dot products
and sign handling), solves and (R'R)⁻¹ by back-substitution on R, weighted OLS
rows factored in decreasing-weight order (phit 91ee538), the R⁻¹R⁻ᵀ
product as matrixmultiply's fused kernel computes it once p exceeds 5, and the
extension's order of multiplication in every derived statistic. The two
implementations therefore take the same path through every iterative fit and
return the same bits.

| Group | Inputs | OLS | Logistic | Probit | Poisson |
|---|---|---:|---:|---:|---:|
| fixtures | all R fixtures, with and without the automatic intercept, plus the r-compatibility and stress fixtures | 0 | 0 | 0 | 0 |
| random | 40 seeded designs per model: n 12–400, 1–5 predictors, weights, offsets, missing rows | 0 | 0 | 0 | 0 |
| scaling | predictors and outcomes scaled by 10⁻⁸ to 10⁸; prior weights 10¹⁰⁰ to 10²⁰⁰ | 0 | 0 | 0 | 0 |
| collinear | condition numbers about 10⁶ to 10¹⁰ | 0 | 0 | (rejected identically) | 0 |
| separation | complete, quasi-complete, rare events (1 in 2000), grouped near-boundary proportions | | 0 | 0 | |
| offsets | binary offsets −50 to 40; per-row offsets spread over −50 to 55 (40 seeded fits each) | | 0 | 0 | 0 |
| counts | Poisson counts up to 10¹² and scaled by 10¹⁶⁰, common offsets −700 to 300, sparse counts | | | | 0 |
| errors | invalid shapes, types, domains, weights, policies, rank deficiency | 0 | same codes | same codes | 0 |

The `iterations` reported were identical in every case. Every invalid input in
the error group raises the same exception code in both, and every fit that
fails in one implementation fails with the same code in the other.

The randomized parity suite (`tests/Parity/RandomizedParityTest.php`) adds
seeded designs with 1 to 9 predictors, integer and zero entries, weights,
offsets, missing rows, and fits without an intercept, and requires every
result to be identical, including the sign of zero. Four seeds of 150 designs
per model passed.

### Inputs where the results are sensitive to rounding

With per-row offsets spread over a wide range (for example −50 to 55, so the
fitted means differ by a factor of about e¹⁰⁰), or one count larger than the
others by 10¹² or more, the IRLS working residual is dominated by a few rows.
Each step then carries rounding error close to the convergence threshold. The
polyfill and the extension still agree, because they perform the same
arithmetic. Their common answer can differ from R's by more than usual, and in
some such cases both report non-convergence (1006).

### Weighted least squares with a dominating weight

Before phit 91ee538 the extension's OLS was wrong when one prior weight was
far larger than the rest. For y = (1, 2, 2.5, 4, 4.5, 6) on x = 1, …, 6 with
weights (1, 2, 3, 4, 5, 10¹⁸), it returned (319.6, −52.27) instead of the
exact (−0.4, 1.0667). The extension now orders rows by weight and solves R by
back-substitution, and the polyfill does the same. Both return
(−0.4000000000000019, 1.0666666666666669), and both raise rank deficiency
(1004) at a weight of 10³⁰.

## Distributions

All 13 families in the extension are provided: normal, t, chi-square, F,
exponential, logistic, log-normal, gamma, beta, binomial, Poisson, negative
binomial, and negative binomial by mean. Each family has 14 functions. These
are the density or mass function and its logarithm, cdf, logcdf, sf, logsf,
and quantile, each with an `_array` variant. That makes 182 functions in all.

The implementation ports the extension's distribution layer and the parts of
the cdflib crate it calls. The cdflib parts are TOMS 708 (incomplete beta),
TOMS 654 (incomplete gamma), Cody's normal CDF, and the `dinvnr`/`dt1`
quantile starts. The extension's log-scale tails, Loader's saddle-point
densities, and its safeguarded Newton quantile search are ported as well.

### Against the references

The fixture comparisons use the extension's own tolerances. These are
1e-12 relative for probabilities and densities (with a 1e-300 floor),
1e-12·max(1, |value|) for log-scale values, and 1e-15 + 1e-12·|value| for
quantiles. The table gives the largest error observed against each reference,
measured the same way. R and SciPy are measured only at points where that
reference itself agrees with mpmath, and the bracketed count is the number of
points excluded because it did not.

| Family | mpmath | R | SciPy | Tested domain |
|---|---:|---:|---:|---|
| normal | 2.9e-14 | 1.8e-14 | 1.1e-13 | sd 1e-8 to 1e3, \|x\| to 1e13 |
| t | 1.8e-13 | 2.6e-13 [4] | 6.6e-13 [40] | df 0.5 to 1e6, \|x\| to 1e10 |
| chi-square | 8.3e-13 | 8.3e-13 [3] | 8.9e-13 [64] | df 0.5 to 1e6, x to 1e7 |
| F | 1.6e-13 | 5.2e-13 [15] | 3.4e-13 [46] | df1 0.5 to 100, df2 0.5 to 1e6, x to 1e20 |
| exponential | 1.5e-14 | 2.1e-16 | 2.1e-16 | rate 1e-5 to 1e5 |
| logistic | 1.5e-15 | 2.0e-16 | 2.0e-16 | scale 1e-5 to 2 |
| log-normal | 9.1e-14 | 2.1e-14 [1] | 1.1e-13 | sdlog 0.5 to 3 |
| gamma | 7.5e-13 | 7.5e-13 | 5.6e-13 [40] | shape 0.01 to 1e4, scale 0.01 to 1e5 |
| beta | 1.5e-13 | 5.7e-14 | 2.6e-13 [15] | shapes 0.1 to 1e3 |
| binomial | 3.3e-14 | 5.7e-14 [3] | 9.9e-13 [31] | size up to 1e6 |
| Poisson | 3.7e-13 | 3.7e-13 [6] | 3.7e-13 [57] | lambda 0 to 1e6 |
| negative binomial | 6.8e-14 | 2.5e-14 [4] | 4.5e-14 [35] | size 0.5 to 1e5 |
| negative binomial (mu) | 5.1e-14 | 3.8e-13 [10] | 9.5e-13 [70] | mu 0 to 1e5, size 0.5 to 1e6 |

Probabilities are tested down to about 1e-300, and the log-scale functions
beyond the double range. `tests/Unit/DistributionsStressTest.php` adds
comparisons with closed forms and summed series in the extreme tails. These
cover the normal asymptotics to z = 1e150, Cauchy and t with df = 2 to 1e300,
chi-square with df 1 and 2, the incomplete gamma for integer shapes up to 1e6,
beta with a unit shape, F with df 2, and the support edges of the binomial and
negative binomial. The suite also checks densities at extreme parameters and
quantile/cdf round trips for p from 1e-300 to 1 − 2⁻⁵².

### Against the extension

Every family is bit-identical to the extension: the largest difference is 0 on
about 8,600 fixture points and on seeded sweeps over tails, ±INF evaluation
points, p in {0, 1, 5e-324, 1e-300, 1 − 2⁻⁵²}, and extreme parameters.
Exception codes are the same for every invalid input tested.

### Limits shared with the extension

For gamma, chi-square, and Poisson, log tails below 1e-280 lose accuracy at
shape or lambda of about 1e6 and above. The relative error is 2e-12 at 1e6,
7e-12 at 1e7, and 5e-10 at 1e8, from cancellation in the extension's
formula. The polyfill reproduces the extension's values exactly there, and a
stress test records the bound.

### Extreme parameters

cdflib's incomplete beta and gamma routines never return for some extreme
arguments and panic in the deep asymptotic regime of the incomplete gamma.
Since phit bd98bf7 the extension returns a failure up front for those
regimes, and the polyfill applies the same conditions. They are an overflowing
a + b, a beta parameter above 1e300 with x or 1 − x below 1e-300, and a gamma
shape or point above 1e300. A panic inside either routine also becomes a
failure. The functions report every such failure as exception 1006, and a
quantile search stops with 1006 as soon as a tail evaluation fails.
`tests/Parity/DistributionsExtremeParityTest.php` compares both
implementations on a grid of parameters from 1e-300 to `PHP_FLOAT_MAX` for F,
beta, gamma, chi-square, t, Poisson, and the negative binomial by mean. The
values or exception codes agree on every point.

For `Phit\dist_t_quantile` at very small df, the quantile search used to step from
a very large starting point by ±2 and could take seconds per call. Both the
extension (since 0aafc26) and the polyfill now search that run of points by
exponential probing and bisection. The visited points are the same, so the
results equal those of the original walk.

## Formula interface

`Phit\model_matrix`, `Phit\fit_lm`, and `Phit\fit_glm` port the extension's formula
parser, term expansion, evaluation, and coding. The fits reuse the polyfill's
OLS and GLM code described above.

### Against R

The 62 cases in `formula-r.json` cover 33 model matrices, 14 `lm` fits, and 15
`glm` fits. The tolerances are those of the extension's suite: 1e-12 for
design matrices (1e-9 + 1e-7·|b| for the orthogonal `poly()` case) and
1e-9 + 1e-7·|b| for fits. Column names, terms, and retained rows must match
exactly.

| Function | Comparisons | Max relative error | Tolerance used | Datasets |
|---|---:|---:|---:|---|
| `Phit\model_matrix` | 266 | 6.7e-16 | 0% | housing, schools, cbpp, synthetic, small |
| `Phit\fit_lm` | 238 | 2.7e-12 | 0% | housing, schools, synthetic |
| `Phit\fit_glm` | 240 | 5.9e-10 | 1% | learning, cbpp, clinic, synthetic |

### Against the extension

`tests/Parity/FormulaParityTest.php` compares the two on every fixture case.
It also covers about 170 syntax forms, about 120 malformed formulas, factor
level and collation edge cases, both missing-value policies, about 28 weight
shapes, 20 family and link combinations, about 50 data shapes, and 30 random
designs. Design matrices, column names, terms, `assign`, row accounting, key
order, and warning text are identical, and every invalid input raises the same
exception code with the same message. The fits are identical to the bit in
every group.

### Prediction

`Phit\predict` ports the extension's prediction design (the fitted factor
levels and `poly()` recurrence coefficients from `predvars`) and its
prediction code. The extension's `predict.php` passes unchanged against the
polyfill: 25 cases against R's `predict()` for lm, glm, lmer, and glmer fits,
with 877 comparisons at the extension's tolerances. Against the extension,
441 calls over the fixture cases (each with seven option sets on new and
training data, plus malformed fits and invalid options) and the randomized
parity suite's 1,390 predictions per seed are identical, including errors.

## Random numbers

Each `*_rng` call seeds a fresh ChaCha20 stream with rand_core's
`seed_from_u64`, draws open-interval uniforms as rand 0.9's `Open01` does, and
maps them through the family's quantile function, as the extension does. The
keystream comes from OpenSSL's `chacha20` cipher when it is available and from
a pure-PHP ChaCha20 otherwise. The two give identical streams.

Against the extension, every draw is identical, as are the exception codes and
messages. This holds for the 5 fixture seeds checked during development and
for the randomized parity suite (all 13 families, random parameters, sizes,
and seeds, including invalid ones). Against the R and NumPy reference vectors
in `rng-python-r.json`, continuous draws agree within 1e-12 relative error and
counts agree exactly.

## Posterior draws

`Phit\draws_moments`, `Phit\draws_covariance_map`, and `Phit\draws_affine`
follow the extension's operation order: Neumaier sums, nalgebra's Cholesky
factorization, triangular solve, and partially pivoted LU. For
`Phit\draws_affine` with more than five rows and columns, nalgebra multiplies
through matrixmultiply, whose aarch64 kernel accumulates each product with a
fused multiply-add. PHP has no `fma()`, so the polyfill emulates a correctly
rounded one (Boldo and Melquiond's round-to-odd algorithm). The emulation
matched exact rational arithmetic on 200,000 random cases.

The results are identical to the extension's in the randomized parity suite.
Without the fused accumulation, roughly 40 to 90 percent of the
entries in a typical affine transform would differ in the last bit. The
numerical-rank check in `Phit\draws_affine` compares singular values with a
tolerance. The polyfill computes them with a one-sided Jacobi method, not
nalgebra's bidiagonal SVD, so a matrix whose smallest singular value lies
within rounding of the threshold could be classified differently. Against the
NumPy references in `draws-python.json`, all cases agree within 1e-11.

## Mixed models

`Phit\fit_mixed_intercept`, `Phit\fit_mixed`, `Phit\fit_logistic_mixed`, `Phit\fit_lmer`, and
`Phit\fit_glmer` port the extension's stats-core code (mixed.rs, lmm.rs,
logistic_mixed.rs, and the formula front end) together with the numerical
libraries it relies on, each reproduced to the bit and checked against the
compiled Rust library on randomized problems:

- BOBYQA from the `bobyqa` 0.3.2 crate (2,224 problems, every evaluated point
  identical).
- Nelder–Mead and the executor loop from argmin 0.11.0 (3,276 problems).
- nalgebra 0.35.0's SVD, symmetric eigendecomposition (for the Gauss–Hermite
  nodes), Householder QR, Cholesky, triangular solves, norms, and matrix
  products, including matrixmultiply's fused kernel (about 30,000
  decompositions and products).

### Against the extension

Every result key is identical to the extension's, bit for bit, including
theta, the variance components, conditional modes and variances, evaluation
counts, and the diagnostics. This holds for the extension's R fixtures and for
the randomized parity suite, which fits seeded designs with 2 to 8 groups,
string and integer group labels, missing rows, random intercepts and slopes
(correlated and uncorrelated), ML and REML, both optimizers, Laplace and
adaptive quadrature, and the formula interfaces. Four seeds of 50 designs
each passed.

### Against R

The extension's own integration scripts (`lmm.php`, `lmm-expanded.php`,
`logistic-mixed.php`, and `formula-mixed.php`) pass unchanged against the
polyfill, apart from the warning level noted below. They compare with lme4
and nlme at the extension's tolerances. The largest differences from lme4 on
the main quantities:

| Function | Reference | Coefficients | Std. errors | Variance components | Log-likelihood |
|---|---|---:|---:|---:|---:|
| `Phit\fit_mixed_intercept` | lmm-r.json (6 fits) | 5.4e-9 | 2.8e-8 | 7.8e-8 | 2.4e-16 |
| `Phit\fit_mixed_intercept` | lmm-public-r.json (4 public datasets) | 1.7e-9 | 4.5e-8 | 1.2e-7 | 4.0e-16 |
| `Phit\fit_logistic_mixed` | logistic-mixed-r.json (28 fits) | 1.7e-6 | 2.2e-5 | 1.2e-6 | 3.1e-9 |

The binomial standard errors differ from lme4 by up to 2·10⁻⁵ relative because
both programs differentiate the approximate likelihood numerically, with
different step rules. `lmm-expanded.php` covers 15 general models against four
lme4 optimizers, and `formula-mixed.php` makes 782 comparisons over 7 `lmer`
and 3 `glmer` fits.

### Speed

The mixed models are the slowest part of the polyfill. Measured on an Apple M5
Ultra with PHP 8.5.10, in the same process as the extension:

| Fit | Extension | Polyfill |
|---|---:|---:|
| `Phit\fit_mixed_intercept`, sleepstudy (n = 180) | 0.2 ms | 2.6 ms |
| `Phit\fit_logistic_mixed`, n = 185, Laplace | 1.3 ms | 26 ms |
| `Phit\fit_logistic_mixed`, n = 185, 9 quadrature points | 2.1 ms | 41 ms |
| `Phit\fit_mixed`, random slopes, q = 24 | 0.7 ms | 37 ms |
| `Phit\fit_mixed`, nested, q = 84 | 7.8 ms | 604 ms |
| `Phit\fit_mixed`, four nested levels, q = 120 | 103 ms | 8.0 s |

Here q is the number of random coefficients. `Phit\fit_mixed` and `Phit\fit_lmer`
repeat a dense QR of a matrix with about n + q rows and q + p columns at every
optimizer step, so their cost grows roughly with q³. Like every polyfilled
function, they throw an exception with code 1007 while a JIT compiles the
process; PHP's JIT has moved their optimizer to a wrong minimum without an
error. See
[docs/notes/function-jit.md](docs/notes/function-jit.md).

## Not provided

`Phit\read_csv` needs the native extension and is left undefined; phibble's
PHP reader is the fallback. `Phit\Polyfill\Polyfill::supports()` reports
which functions are available.

## Known differences from the extension

- The listwise-deletion warning is raised with `E_USER_WARNING`, because PHP
  code cannot raise `E_WARNING`. The message text is identical, including the
  `Phit\fit_ols(): ` prefix.
- A scalar argument of the wrong type follows PHP's rules for a userland
  function. In files without `strict_types`, numeric strings and booleans are
  coerced, and other values raise `TypeError`. The extension raises an
  `\Exception` with the message "Invalid value given for argument" instead.
- The parity suites ran on aarch64 (Apple silicon), where matrixmultiply
  accumulates with fused multiply-adds, as the polyfill does. On an x86-64 CPU
  without FMA the extension's products of matrices larger than 5×5 can differ
  from the polyfill's in the last bit.
