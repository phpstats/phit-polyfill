<?php

declare(strict_types=1);

// Polyfill vs native extension for Phit\fit_ols_matrix/logistic/probit/poisson: same keys,
// values within tolerance, same exception codes and the same warning text.
// Runs only when phit is loaded. Set PARITY_REPORT=path to write the
// largest observed differences per model as JSON (used for ACCURACY.md).

const MODEL_FUNCTIONS = ['Phit\fit_ols_matrix', 'Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix', 'Phit\fit_poisson_matrix'];

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

/**
 * Largest scaled difference for one result key. Each value is judged against
 * the natural scale of its quantity, so rounding noise in a near-zero number
 * is not mistaken for disagreement:
 * - coefficients and conf_int: |a - b| / (|b| + standard error);
 * - t/z statistics and R-squared: |a - b| / (|b| + 1);
 * - p_values: |ln a - ln b| / max(1, |ln b|) (relative error on the log scale);
 * - fitted and residuals: |a - b| / (|b| + max |y|), y = fitted + residual;
 * - everything else: |a - b| / (|b| + 1e-3 max |b|).
 */
function modelDifference(array $polyfill, array $native, string $key): float
{
    $flatten = function (mixed $value) use (&$flatten): array {
        if (! is_array($value)) {
            return [$value];
        }

        return array_merge(...array_map($flatten, array_values($value)));
    };
    $a = $flatten($polyfill[$key]);
    $b = $flatten($native[$key]);
    if (count($a) !== count($b)) {
        return INF;
    }
    $largest = 0.0;
    foreach ($b as $value) {
        if (is_float($value) && is_finite($value)) {
            $largest = max($largest, abs($value));
        }
    }
    $se = $native['std_errors'] ?? [];
    $floor = function (int $i) use ($key, $native, $se, $largest): float {
        return match ($key) {
            'coefficients' => $se[$i],
            'conf_int' => $se[intdiv($i, 2)],
            't_statistics', 'z_statistics', 'r_squared', 'adj_r_squared' => 1.0,
            'fitted', 'residuals' => max(array_map(fn ($f, $r) => abs($f + $r), $native['fitted'], $native['residuals'])),
            default => 1e-3 * $largest,
        };
    };
    $max = 0.0;
    foreach ($b as $i => $value) {
        if (! is_int($value) && ! is_float($value)) {
            $max = max($max, $a[$i] === $value ? 0.0 : INF);

            continue;
        }
        $x = (float) $a[$i];
        $y = (float) $value;
        if (is_nan($x) || is_nan($y) || is_infinite($x) || is_infinite($y)) {
            $max = max($max, ($x === $y || (is_nan($x) && is_nan($y))) ? 0.0 : INF);

            continue;
        }
        if ($key === 'p_values' && $x > 0.0 && $y > 0.0) {
            $max = max($max, abs(log($x) - log($y)) / max(1.0, abs(log($y))));

            continue;
        }
        $max = max($max, abs($x - $y) / (abs($y) + $floor($i) + PHP_FLOAT_MIN));
    }

    return $max;
}

function recordParity(string $fn, string $group, float $difference): void
{
    $GLOBALS['model_parity'][$fn][$group] = max($GLOBALS['model_parity'][$fn][$group] ?? 0.0, $difference);
}

/**
 * Compare polyfill and native on one input. Values must agree within $tolerance
 * (relative, see modelDifference); errors must carry the same code.
 */
function assertModelParity(string $fn, array $arguments, float $tolerance, string $group): void
{
    if ($fn !== 'Phit\fit_ols_matrix') {
        // Iterative fits stop once the full step is below 1e-8 relative (1e-11 for
        // probit), so two implementations agree to about that, not to rounding;
        // slow boundary fits and extreme tails amplify it a little.
        $tolerance = max($tolerance, 1e-7);
    }
    $call = function (callable $f) use ($arguments): array {
        $warnings = [];
        set_error_handler(function (int $severity, string $message) use (&$warnings): bool {
            $warnings[] = $message;

            return true;
        });
        try {
            return ['result' => $f(...$arguments), 'warnings' => $warnings];
        } catch (\Exception $e) {
            return ['code' => $e->getCode(), 'message' => $e->getMessage(), 'warnings' => $warnings];
        } finally {
            restore_error_handler();
        }
    };
    $polyfill = $call(\Phit\Polyfill\Polyfill::implementation($fn));
    $native = $call($fn);
    parityLog('models', $fn, $polyfill, $native, isset($native['code']) || isset($polyfill['code']));
    if (isset($native['code']) || isset($polyfill['code'])) {
        expect($polyfill['code'] ?? 'result')->toBe($native['code'] ?? 'result', "$fn [$group]: ".($native['message'] ?? $polyfill['message'] ?? ''));
        recordParity($fn, 'errors', 0.0);

        return;
    }
    $a = $polyfill['result'];
    $b = $native['result'];
    expect(array_keys($a))->toBe(array_keys($b))
        ->and($polyfill['warnings'])->toBe($native['warnings']);
    foreach (array_keys($b) as $key) {
        if ($key === 'iterations') {
            // A converged fit may stop one iteration apart when the last full
            // step sits on the convergence threshold. In the stress groups the
            // likelihood can be flat to rounding (e.g. offsets near 300), so the
            // stopping iteration is noise; only the estimates must agree there.
            $GLOBALS['model_parity'][$fn]['iteration_difference'] = max($GLOBALS['model_parity'][$fn]['iteration_difference'] ?? 0, abs($a[$key] - $b[$key]));
            if (in_array($group, ['fixtures', 'random', 'scaling', 'errors'], true)) {
                expect(abs($a[$key] - $b[$key]))->toBeLessThanOrEqual(1, "$fn [$group] iterations");
            }

            continue;
        }
        if ($key === 'diagnostics') {
            $da = $a[$key];
            $db = $b[$key];
            unset($da['iterations'], $db['iterations']);
            expect($da)->toBe($db);

            continue;
        }
        if (in_array($key, ['n_input', 'n_observations', 'n_dropped', 'retained_rows', 'dropped_rows', 'df_residual', 'converged'], true)) {
            expect($a[$key])->toBe($b[$key]);

            continue;
        }
        $difference = modelDifference($a, $b, $key);
        recordParity($fn, $group, $difference);
        expect($difference)->toBeLessThanOrEqual($tolerance, "$fn [$group] $key differs by $difference");
    }
}

/** Seeded uniform in [0, 1). */
function uniform(): float
{
    return mt_rand() / (mt_getrandmax() + 1.0);
}

function gaussian(): float
{
    return sqrt(-2.0 * log(1.0 - uniform())) * cos(2.0 * M_PI * uniform());
}

/** Design with $p predictors (no intercept column) and optional collinearity. */
function randomDesign(int $n, int $p, float $collinearity = 0.0, float $scale = 1.0): array
{
    $x = [];
    for ($i = 0; $i < $n; $i++) {
        $row = [];
        for ($j = 0; $j < $p; $j++) {
            $row[] = gaussian() * $scale;
        }
        if ($collinearity > 0.0 && $p > 1) {
            // Last column = first column + tiny noise: condition number ~ 1/collinearity.
            $row[$p - 1] = $row[0] + $collinearity * gaussian() * $scale;
        }
        $x[] = $row;
    }

    return $x;
}

function linearPredictor(array $x, array $beta): array
{
    return array_map(function (array $row) use ($beta) {
        $eta = $beta[0];
        foreach ($row as $j => $v) {
            $eta += $v * $beta[$j + 1];
        }

        return $eta;
    }, $x);
}

function randomOutcome(string $fn, array $eta): array
{
    return array_map(fn (float $e) => match ($fn) {
        'Phit\fit_ols_matrix' => $e + gaussian(),
        'Phit\fit_logistic_matrix' => (float) (uniform() < 1.0 / (1.0 + exp(-$e))),
        'Phit\fit_probit_matrix' => (float) (gaussian() < $e),
        'Phit\fit_poisson_matrix' => (float) poissonDraw(exp($e)),
    }, $eta);
}

function poissonDraw(float $mean): int
{
    if ($mean > 50) {
        return max(0, (int) round($mean + sqrt($mean) * gaussian()));
    }
    $l = exp(-$mean);
    $k = 0;
    $p = 1.0;
    do {
        $k++;
        $p *= uniform();
    } while ($p > $l);

    return $k - 1;
}

afterAll(function () {
    $path = getenv('PARITY_REPORT');
    if ($path && isset($GLOBALS['model_parity'])) {
        $existing = is_file($path) ? json_decode(file_get_contents($path), true) : [];
        file_put_contents($path, json_encode(array_replace_recursive($existing, $GLOBALS['model_parity']), JSON_PRETTY_PRINT));
    }
});

it('agrees with the extension on the reference fixtures', function () {
    $files = ['Phit\fit_ols_matrix' => 'ols-r.json', 'Phit\fit_logistic_matrix' => 'logistic-r.json', 'Phit\fit_probit_matrix' => 'probit-r.json', 'Phit\fit_poisson_matrix' => 'poisson-r.json'];
    foreach ($files as $fn => $file) {
        foreach (loadFixture($file)['cases'] as $case) {
            // A perfect fit's t statistics are rounding noise divided by rounding noise.
            if ($case['name'] === 'perfect') {
                continue;
            }
            assertModelParity($fn, [$case['x'], $case['y'], 'drop', null, null, false], 1e-10, 'fixtures');
            $x = array_map(fn ($row) => array_slice($row, 1), $case['x']);
            assertModelParity($fn, [$x, $case['y']], 1e-10, 'fixtures');
        }
    }
    foreach ([...loadFixture('r-compatibility.json')['cases'], ...loadFixture('stress-r.json')['cases']] as $case) {
        assertModelParity('Phit\fit_'.$case['model'].'_matrix', [$case['x'], $case['y'], 'drop', $case['weights'], $case['offset'], false], 1e-8, 'fixtures');
    }
});

it('agrees with the extension on random designs', function (string $fn) {
    mt_srand(crc32($fn));
    for ($trial = 0; $trial < 40; $trial++) {
        $n = [12, 30, 100, 400][$trial % 4];
        $p = 1 + $trial % 5;
        $x = randomDesign($n, $p);
        $beta = array_map(fn () => gaussian() * 0.5, range(0, $p));
        $y = randomOutcome($fn, linearPredictor($x, $beta));
        $weights = $trial % 3 === 0 ? array_map(fn () => (float) mt_rand(1, 5), $y) : null;
        $offset = $trial % 4 === 1 ? array_map(fn () => gaussian() * 0.3, $y) : null;
        if ($trial % 5 === 2) {
            $x[3][0] = null;
            $y[7] = NAN;
        }
        assertModelParity($fn, [$x, $y, 'drop', $weights, $offset], 1e-9, 'random');
    }
})->with(MODEL_FUNCTIONS);

it('agrees on near-collinear designs (condition numbers 1e6 to 1e10)', function (string $fn) {
    mt_srand(7 + crc32($fn));
    foreach ([1e-6, 1e-7, 1e-8, 1e-9, 1e-10] as $collinearity) {
        $x = randomDesign(200, 3, $collinearity);
        $y = randomOutcome($fn, linearPredictor($x, [0.2, 0.5, -0.3, 0.4]));
        // Two backward-stable solvers agree to about cond * eps on ill-conditioned coefficients.
        assertModelParity($fn, [$x, $y], max(1e-9, 1e3 * PHP_FLOAT_EPSILON / $collinearity), 'collinear');
    }
})->with(MODEL_FUNCTIONS);

it('agrees on badly scaled predictors and outcomes', function (string $fn) {
    mt_srand(11 + crc32($fn));
    foreach ([1e-8, 1e-3, 1e3, 1e8] as $scale) {
        $x = randomDesign(150, 3, 0.0, $scale);
        $beta = [0.1, 0.4 / $scale, -0.3 / $scale, 0.2 / $scale];
        $y = randomOutcome($fn, linearPredictor($x, $beta));
        if ($fn === 'Phit\fit_ols_matrix') {
            $y = array_map(fn ($v) => $v * $scale, $y);
        }
        assertModelParity($fn, [$x, $y], 1e-9, 'scaling');
    }
})->with(MODEL_FUNCTIONS);

it('agrees on separation and near-separation', function (string $fn) {
    mt_srand(13 + crc32($fn));
    $x = array_map(fn ($i) => [(float) $i], range(-10, 9));
    // Complete separation, quasi-complete separation, and one overlapping pair.
    $separated = array_map(fn ($row) => $row[0] >= 0 ? 1.0 : 0.0, $x);
    $quasi = $separated;
    $x2 = $x;
    $x2[] = [0.0];
    $quasi[] = 0.0;
    $overlap = $separated;
    $overlap[9] = 1.0;
    $overlap[10] = 0.0;
    assertModelParity($fn, [$x, $separated], 1e-8, 'separation');
    assertModelParity($fn, [$x2, $quasi], 1e-8, 'separation');
    assertModelParity($fn, [$x, $overlap], 1e-8, 'separation');
    // Rare events: one success in 2000 trials, and grouped near-boundary proportions.
    $rare = array_fill(0, 2000, 0.0);
    $rare[1234] = 1.0;
    assertModelParity($fn, [array_map(fn () => [gaussian()], $rare), $rare], 1e-8, 'separation');
    $grouped = array_map(fn ($i) => [$i / 10], range(0, 19));
    $trials = array_fill(0, 20, 1000.0);
    $y = array_map(fn ($i) => $i < 10 ? 1 / 1000 : 999 / 1000, range(0, 19));
    assertModelParity($fn, [$grouped, $y, 'drop', $trials], 1e-8, 'separation');
})->with(['Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix']);

it('agrees on huge counts, extreme offsets and sparse counts in Poisson', function () {
    mt_srand(17);
    $x = randomDesign(120, 2);
    foreach ([1e3, 1e6, 1e9, 1e12] as $level) {
        $y = array_map(fn ($row) => (float) round($level * exp(0.1 * $row[0] - 0.05 * $row[1]) + sqrt($level) * gaussian()), $x);
        assertModelParity('Phit\fit_poisson_matrix', [$x, $y], 1e-8, 'counts');
    }
    foreach ([-700.0, -100.0, 30.0, 300.0] as $shift) {
        $offset = array_map(fn () => $shift + gaussian(), $x);
        $y = array_map(fn ($row) => (float) poissonDraw(2.0 * exp(0.3 * $row[0])), $x);
        assertModelParity('Phit\fit_poisson_matrix', [$x, $y, 'drop', null, $offset], 1e-8, 'counts');
    }
    $sparse = array_fill(0, 500, 0.0);
    $sparse[17] = 1.0;
    $sparse[402] = 2.0;
    assertModelParity('Phit\fit_poisson_matrix', [array_map(fn () => [gaussian()], $sparse), $sparse], 1e-8, 'counts');
});

it('agrees on extreme offsets for binary models', function (string $fn) {
    mt_srand(19 + crc32($fn));
    $x = randomDesign(150, 2);
    foreach ([-50.0, -12.0, 9.0, 40.0] as $shift) {
        $offset = array_map(fn () => $shift + gaussian(), $x);
        $y = array_map(fn ($i) => (float) ($i % 2), array_keys($x));
        assertModelParity($fn, [$x, $y, 'drop', null, $offset], 1e-8, 'offsets');
    }
})->with(['Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix']);

it('raises the same exception codes on invalid input', function (string $fn) {
    $binary = in_array($fn, ['Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix'], true);
    $cases = [
        [[], []], [[[]], [1]], [[[1], [1, 2]], [1, 2]], [[[1], [2]], [1]],
        [[[1], [2], [3]], ['a', 1, 0]], [[[1], [2], [3]], [true, 1, 0]], [[[1], [2], [3]], [1, INF, 0]],
        [[[1], ['x'], [3]], [1, 0, 1]], [[[1], [INF], [3]], [1, 0, 1]], [[[1], 2, [3]], [1, 0, 1]],
        [[[1, 1], [1, 1], [1, 1], [1, 1]], [1, 0, 1, 0]],
        [[[1], [2], [3], [4]], [0.5, 0, 1, 0]], [[[1], [2], [3], [4]], [-1, 0, 1, 0]],
        [[[1], [2], [3], [4]], [2, 0, 1, 0]], [[[1], [2], [3], [4]], [0, 0, 0, 0]],
        [[[null], [null], [null]], [1, 0, 1]],
        [[[1], [2], [3], [4]], [1, 0, 1, 0], 'raise', [1, 1, 1]],
        [[[1], [2], [3], [4]], [1, 0, 1, 0], 'drop', [1, 0, 1, 1]],
        [[[1], [2], [3], [4]], [1, 0, 1, 0], 'drop', [1, 1.5, 1, 1]],
        [[[1], [2], [3], [4]], [1, 0, 1, 0], 'drop', [1, 2e6, 1, 1]],
        [[[1], [2], [3], [4]], [1, 0, 1, 0], 'drop', null, [0, INF, 0, 0]],
        [[[1], [2], [3], [4]], [1, 0, 1, 0], 'drop', null, [0, '1', 0, 0]],
        [[[1], [2], [3], [4]], [1, 0, 1, 0], 'nope'],
        [[[1], [null], [3], [4]], [1, 0, 1, 0], 'raise'],
        [array_fill(0, 4, []), [1, 0, 1, 0], 'drop', null, null, false],
        [[[1], [2]], [1, 0]],
        [[[1e200], [2e200], [3e200]], [1, 0, 1]],
    ];
    if ($binary) {
        $cases[] = [[[1], [2], [3], [4]], [0.5, 0, 1, 0], 'drop', [3, 1, 1, 1]];
        $cases[] = [[[-2], [-1], [1], [2]], [0, 0, 1, 1]];
    }
    foreach ($cases as $arguments) {
        assertModelParity($fn, $arguments, 1e-9, 'errors');
    }
})->with(MODEL_FUNCTIONS);

it('agrees with the extension on huge weights and counts', function (string $fn) {
    // Review finding: column norms above about 1e77 overflowed the Jacobi rotation test.
    $x = [[1.0], [2.0], [3.0], [4.5], [5.0], [6.0]];
    $y = $fn === 'Phit\fit_poisson_matrix' ? [1, 2, 1, 4, 3, 6] : [0, 1, 0, 1, 1, 0];
    foreach ([1e100, 1e160, 1e200] as $w) {
        assertModelParity($fn, [$x, $y, 'drop', array_fill(0, 6, $w)], 1e-9, 'scaling');
        assertModelParity($fn, [$x, $y, 'drop', [$w, 1, 1, 1, 1, 1]], 1e-9, 'scaling');
    }
    if ($fn === 'Phit\fit_poisson_matrix') {
        foreach ([1e100, 1e160] as $scale) {
            // Likelihoods near 1e160 are flat to rounding at the optimum, so the
            // stopping iteration is noise; the 'counts' group checks estimates only.
            assertModelParity($fn, [$x, array_map(fn ($v) => $v * $scale, $y)], 1e-9, 'counts');
        }
    }
})->with(['Phit\fit_ols_matrix', 'Phit\fit_poisson_matrix']);

it('agrees on widely spread per-row offsets', function (string $fn) {
    // Review finding: with offsets spread over about -50 to 55 the IRLS steps are
    // dominated by rounding, so agreement needs the extension's exact QR
    // arithmetic (phit 91ee538: Householder QR, triangular solves).
    mt_srand(23 + crc32($fn));
    for ($trial = 0; $trial < 40; $trial++) {
        $n = mt_rand(3, 35);
        $x = array_fill(0, $n, []);
        $offset = array_map(fn () => -50.0 + 105.0 * uniform(), range(1, $n));
        $y = $fn === 'Phit\fit_poisson_matrix'
            ? array_map(fn () => (float) mt_rand(0, 3), range(1, $n))
            : array_map(fn () => (float) mt_rand(0, 1), range(1, $n));
        $weights = $trial % 3 === 0 ? array_map(fn () => (float) mt_rand(1, 4), range(1, $n)) : null;
        assertModelParity($fn, [$x, $y, 'drop', $weights, $offset], 1e-9, 'offsets');
    }
})->with(['Phit\fit_logistic_matrix', 'Phit\fit_poisson_matrix']);
