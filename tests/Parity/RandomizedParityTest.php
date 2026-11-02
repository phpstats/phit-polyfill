<?php

declare(strict_types=1);

// Randomized differential tests: seeded random inputs (including invalid and
// non-finite ones) go to the polyfill and to the loaded extension, and the
// results must be identical bit for bit, or both must throw the same exception
// class, code, and message. Skipped without the extension.
//
// PARITY_SEED picks the generator seed (default 20261016) and
// PARITY_ITERATIONS the rounds per function group (default 60).

use Phit\Polyfill\Polyfill;

if (! function_exists('randomParityOutcome')) {
    function randomParityOutcome(callable $call): array
    {
        try {
            return ['value', $call()];
        } catch (\Throwable $e) {
            return ['error', get_class($e), $e->getCode(), $e->getMessage()];
        }
    }

    /** Path and values of the first difference between two outcomes, or null. */
    function randomParityFirstDifference(mixed $a, mixed $b, string $path = ''): ?string
    {
        if (is_array($a) && is_array($b)) {
            if (array_keys($a) !== array_keys($b)) {
                return "$path keys ".json_encode(array_keys($a)).' vs '.json_encode(array_keys($b));
            }
            foreach ($a as $k => $v) {
                $d = randomParityFirstDifference($v, $b[$k], "$path/$k");
                if ($d !== null) {
                    return $d;
                }
            }

            return null;
        }

        return serialize($a) === serialize($b) ? null : sprintf('%s: %s vs %s', $path, var_export($a, true), var_export($b, true));
    }

    function randomParityCompare(string $fn, array $arguments): void
    {
        $native = randomParityOutcome(fn () => $fn(...$arguments));
        $polyfill = randomParityOutcome(fn () => Polyfill::implementation($fn)(...$arguments));
        parityLog('randomized', $fn, $polyfill, $native, $native[0] === 'error' || $polyfill[0] === 'error');
        if (serialize($native) !== serialize($polyfill)) {
            $show = fn (mixed $v) => substr(json_encode($v, JSON_PARTIAL_OUTPUT_ON_ERROR | JSON_PRESERVE_ZERO_FRACTION), 0, 400);
            throw new \RuntimeException("$fn differs at ".randomParityFirstDifference($native, $polyfill)."\nfor ".$show($arguments));
        }
    }

    function randomUniform(float $lo, float $hi): float
    {
        return $lo + ($hi - $lo) * mt_rand() / mt_getrandmax();
    }

    function randomParityIterations(): int
    {
        return (int) (getenv('PARITY_ITERATIONS') ?: 60);
    }

    function randomParitySeed(string $group): void
    {
        mt_srand(crc32($group) ^ (int) (getenv('PARITY_SEED') ?: 20261016));
    }

    /** A covariance matrix A A' + D with a random positive diagonal. */
    function randomCovariance(int $p): array
    {
        $a = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $a[$i][$j] = randomUniform(-1, 1);
            }
        }
        $c = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $s = 0.0;
                for ($k = 0; $k < $p; $k++) {
                    $s += $a[$i][$k] * $a[$j][$k];
                }
                $c[$i][$j] = $s + ($i === $j ? randomUniform(0, 0.5) : 0.0);
            }
        }

        return $c;
    }
}

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

it('samples exactly like the extension', function () {
    randomParitySeed('rng');
    $families = [
        'normal' => fn () => [randomUniform(-5, 5), randomUniform(0.01, 10)],
        't' => fn () => [mt_rand(0, 3) > 0 ? randomUniform(0.1, 50) : randomUniform(1e-3, 1e6)],
        'chisquare' => fn () => [randomUniform(0.05, 100)],
        'f' => fn () => [randomUniform(0.1, 50), randomUniform(0.1, 50)],
        'exponential' => fn () => [randomUniform(0.01, 10)],
        'logistic' => fn () => [randomUniform(-5, 5), randomUniform(0.01, 10)],
        'lognormal' => fn () => [randomUniform(-3, 3), randomUniform(0.01, 3)],
        'gamma' => fn () => [randomUniform(0.01, 50), randomUniform(0.01, 10)],
        'beta' => fn () => [randomUniform(0.01, 20), randomUniform(0.01, 20)],
        'binomial' => fn () => [mt_rand(0, 500), randomUniform(0, 1)],
        'poisson' => fn () => [randomUniform(0, 300)],
        'nbinom' => fn () => [randomUniform(0.01, 50), randomUniform(0.001, 1)],
        'nbinom_mu' => fn () => [randomUniform(0.01, 50), randomUniform(0, 100)],
    ];
    $calls = 0;
    for ($round = 0; $round < randomParityIterations(); $round++) {
        foreach ($families as $family => $parameters) {
            $params = $parameters();
            if (mt_rand(0, 20) === 0) {
                $k = array_rand($params);
                // Integer parameters stay integers: the extension rejects floats
                // for them where PHP would coerce (a documented difference).
                $params[$k] = is_int($params[$k])
                    ? [-1, 0, PHP_INT_MAX][mt_rand(0, 2)]
                    : [-1.0, 0.0, NAN, INF, 1e308][mt_rand(0, 4)];
            }
            $n = [0, 1, 7, 64, 333][mt_rand(0, 4)];
            if (mt_rand(0, 40) === 0) {
                $n = [-1, 1000001][mt_rand(0, 1)];
            }
            $seed = mt_rand(0, 3) > 0 ? mt_rand(0, PHP_INT_MAX) : [0, -1, 1][mt_rand(0, 2)];
            randomParityCompare("Phit\\dist_{$family}_rng", [$n, $seed, ...$params]);
            $calls++;
        }
    }
    expect($calls)->toBeGreaterThan(0);
});

it('computes draw kernels exactly like the extension', function () {
    randomParitySeed('draws');
    for ($round = 0; $round < randomParityIterations() * 3; $round++) {
        $n = mt_rand(1, 12);
        $p = mt_rand(1, 9);
        $scale = 10 ** mt_rand(-3, 6);
        $offset = mt_rand(0, 3) > 0 ? 0 : 10 ** mt_rand(2, 9);
        $draws = [];
        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $draws[$i][$j] = $offset + $scale * randomUniform(-1, 1);
            }
        }
        $weights = mt_rand(0, 2) > 0 ? null : array_map(fn () => mt_rand(0, 4) > 0 ? randomUniform(0, 5) : 0.0, range(1, $n));
        randomParityCompare('Phit\draws_moments', [$draws, $weights, ['population', 'unbiased'][mt_rand(0, 1)]]);

        randomParityCompare('Phit\draws_covariance_map', [randomCovariance($p), randomCovariance(mt_rand(0, 8) > 0 ? $p : $p + 1)]);

        $matrix = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $matrix[$i][$j] = randomUniform(-2, 2) * 10 ** mt_rand(-2, 2);
            }
        }
        if ($p > 1 && mt_rand(0, 10) === 0) {
            $matrix[1] = $matrix[0];
        }
        $shift = array_map(fn () => randomUniform(-10, 10), range(1, $p));
        // More than five rows and columns takes the extension's fused kernel.
        $rows = mt_rand(0, 1) === 1 ? mt_rand(6, 40) : mt_rand(1, 5);
        $input = [];
        for ($i = 0; $i < $rows; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $input[$i][$j] = randomUniform(-3, 3);
            }
        }
        randomParityCompare('Phit\draws_affine', [$input, $matrix, $shift]);
    }
    expect(true)->toBeTrue();
});

if (! function_exists('randomDesign')) {
    /** @return array{0: list<list<float|int|null>>, 1: list<float>, 2: list<int>, 3: list<int>, 4: list<string|int|null>} */
    function randomDesign(int $n, int $p, int $groups): array
    {
        $x = $y = $binary = $counts = $labels = [];
        $effects = array_map(fn () => randomUniform(-1, 1), range(0, $groups - 1));
        $scale = 10 ** mt_rand(-2, 3);
        for ($i = 0; $i < $n; $i++) {
            $row = [];
            for ($j = 0; $j < $p; $j++) {
                $row[] = mt_rand(0, 9) === 0 ? mt_rand(-3, 3) : $scale * randomUniform(-1, 1) * ($j + 1);
            }
            $g = mt_rand(0, $groups - 1);
            $signal = 0.0;
            foreach ($row as $j => $v) {
                $signal += $v / ($scale * ($j + 1)) * (($j % 3) - 1) * 0.7;
            }
            $signal = max(-6.0, min(6.0, $signal + $effects[$g]));
            $x[] = $row;
            $y[] = 2.0 + $signal + randomUniform(-1, 1);
            $binary[] = randomUniform(0, 1) < 1 / (1 + exp(-$signal)) ? 1 : 0;
            $counts[] = (int) floor(exp(0.3 * $signal) * randomUniform(0, 4));
            $labels[] = mt_rand(0, 1) === 1 ? "g$g" : $g;
        }
        if (mt_rand(0, 4) === 0) {
            $x[mt_rand(0, $n - 1)][0] = null; // listwise deletion
        }

        return [$x, $y, $binary, $counts, $labels];
    }
}

it('fits regressions exactly like the extension', function () {
    randomParitySeed('models');
    for ($round = 0; $round < randomParityIterations(); $round++) {
        $n = mt_rand(8, 120);
        $p = mt_rand(1, 9);
        [$x, $y, $binary, $counts] = randomDesign($n, $p, 4);
        $weights = mt_rand(0, 3) === 0 ? array_map(fn () => (float) mt_rand(1, 4), range(1, $n)) : null;
        $offset = mt_rand(0, 3) === 0 ? array_map(fn () => randomUniform(-0.5, 0.5), range(1, $n)) : null;
        $intercept = mt_rand(0, 5) > 0;
        set_error_handler(fn () => true);
        try {
            randomParityCompare('Phit\fit_ols_matrix', [$x, $y, 'drop', $weights, $offset, $intercept]);
            randomParityCompare('Phit\fit_logistic_matrix', [$x, $binary, 'drop', $weights, $offset, $intercept]);
            randomParityCompare('Phit\fit_probit_matrix', [$x, $binary, 'drop', $weights, $offset, $intercept]);
            randomParityCompare('Phit\fit_poisson_matrix', [$x, $counts, 'drop', $weights, $offset, $intercept]);
            $values = array_map(fn () => randomUniform(-5, 5) * 10 ** mt_rand(-3, 6), range(1, $n));
            randomParityCompare('Phit\stat_mean', [$values]);
            randomParityCompare('Phit\stat_var', [$values]);
            randomParityCompare('Phit\stat_quantile', [$values, randomUniform(0, 1)]);
            randomParityCompare('Phit\stat_cor', [$values, $y]);
        } finally {
            restore_error_handler();
        }
    }
    expect(true)->toBeTrue();
});

it('fits mixed models exactly like the extension', function () {
    skipMixedFitsUnderJit();
    randomParitySeed('mixed');
    for ($round = 0; $round < max(1, intdiv(randomParityIterations(), 3)); $round++) {
        $groups = mt_rand(2, 8);
        $n = mt_rand(3 * $groups, 12 * $groups);
        $p = mt_rand(1, 3);
        [$x, $y, $binary, , $labels] = randomDesign($n, $p, $groups);
        set_error_handler(fn () => true);
        try {
            randomParityCompare('Phit\fit_mixed_gaussian_intercept_matrix', [$x, $y, $labels]);
            $intercepts = array_fill(0, $n, [1.0]);
            $slopes = array_map(fn ($row) => [1.0, $row[0] ?? 0.0], $x);
            randomParityCompare('Phit\fit_mixed_gaussian_matrix', [$x, $y, [['z' => $intercepts, 'groups' => $labels]], 'reml' => mt_rand(0, 1) === 1]);
            randomParityCompare('Phit\fit_mixed_gaussian_matrix', [$x, $y, [['z' => $slopes, 'groups' => $labels, 'correlated' => mt_rand(0, 1) === 1]],
                'optimizer' => mt_rand(0, 2) === 0 ? 'nelder_mead' : 'bobyqa']);
            randomParityCompare('Phit\fit_mixed_logistic_matrix', [$x, $binary, $labels, 'quadraturePoints' => [1, 1, 3, 9][mt_rand(0, 3)]]);
            $data = ['y' => $y, 'b' => $binary, 'g' => array_map('strval', $labels)];
            foreach ($x as $i => $row) {
                foreach ($row as $j => $v) {
                    $data["x$j"][$i] = $v;
                }
            }
            $rhs = implode(' + ', array_map(fn ($j) => "x$j", range(0, $p - 1)));
            randomParityCompare('Phit\fit_mixed_gaussian', ["y ~ $rhs + (1 | g)", $data]);
            randomParityCompare('Phit\fit_mixed_logistic', ["b ~ $rhs + (1 | g)", $data]);
        } finally {
            restore_error_handler();
        }
    }
    expect(true)->toBeTrue();
});

it('fits nested and crossed Gaussian mixed models exactly like the extension', function () {
    // Nested terms split into a block per outer group (phit 0.4.3); crossed
    // terms usually connect every random column and keep the dense summary.
    skipMixedFitsUnderJit();
    randomParitySeed('mixed-blocks');
    for ($round = 0; $round < max(1, intdiv(randomParityIterations(), 6)); $round++) {
        $outer = mt_rand(2, 6);
        $inner = mt_rand(2, 4);
        $n = mt_rand(3 * $outer * $inner, 6 * $outer * $inner);
        $p = mt_rand(1, 2);
        [$x, $y, , , $labels] = randomDesign($n, $p, $outer);
        $nested = array_map(fn ($label) => "$label/".mt_rand(0, $inner - 1), $labels);
        $crossed = array_map(fn () => mt_rand(0, $inner), $labels);
        $intercepts = array_fill(0, $n, [1.0]);
        $slopes = array_map(fn ($row) => [1.0, $row[0] ?? 0.0], $x);
        $options = ['reml' => mt_rand(0, 1) === 1, 'optimizer' => mt_rand(0, 2) === 0 ? 'nelder_mead' : 'bobyqa'];
        set_error_handler(fn () => true);
        try {
            randomParityCompare('Phit\fit_mixed_gaussian_matrix', [$x, $y, [['z' => $intercepts, 'groups' => $labels], ['z' => $intercepts, 'groups' => $nested]], ...$options]);
            randomParityCompare('Phit\fit_mixed_gaussian_matrix', [$x, $y, [['z' => $slopes, 'groups' => $labels, 'correlated' => mt_rand(0, 1) === 1], ['z' => $intercepts, 'groups' => $nested]], ...$options]);
            randomParityCompare('Phit\fit_mixed_gaussian_matrix', [$x, $y, [['z' => $intercepts, 'groups' => $labels], ['z' => $intercepts, 'groups' => $crossed]], ...$options]);
        } finally {
            restore_error_handler();
        }
    }
    expect(true)->toBeTrue();
});

it('predicts exactly like the extension', function () {
    randomParitySeed('predict');
    for ($round = 0; $round < randomParityIterations(); $round++) {
        $n = mt_rand(12, 60);
        $make = function (int $rows, bool $missing) {
            $data = ['x' => [], 'z' => [], 'g' => [], 'f' => [], 'y' => [], 'b' => []];
            for ($i = 0; $i < $rows; $i++) {
                $x = randomUniform(-3, 3);
                $g = ['a', 'b', 'c'][mt_rand(0, 2)];
                $data['x'][] = $missing && mt_rand(0, 9) === 0 ? null : $x;
                $data['z'][] = mt_rand(0, 1) === 1;
                $data['g'][] = $g;
                $data['f'][] = 'u'.mt_rand(1, 4);
                $eta = 0.5 * $x + ($g === 'b' ? 0.7 : 0.0);
                $data['y'][] = 1.0 + $eta + randomUniform(-1, 1);
                $data['b'][] = randomUniform(0, 1) < 1 / (1 + exp(-$eta)) ? 1 : 0;
            }

            return $data;
        };
        $train = $make($n, false);
        $new = $make(mt_rand(1, 15), true);
        $rhs = ['x + g', 'poly(x, 2) + g + z', 'x * g + offset(0.1 * x)', 'poly(x, 3) + I(x > 0)', 'log(abs(x) + 1) + z'][mt_rand(0, 4)];
        $fits = [];
        set_error_handler(fn () => true);
        try {
            foreach ([
                fn () => \Phit\fit_ols("y ~ $rhs", $train),
                fn () => \Phit\fit_generalized_linear("b ~ $rhs", $train, family: 'binomial'),
                fn () => \Phit\fit_generalized_linear("b ~ $rhs", $train, family: 'binomial', link: 'probit'),
                fn () => \Phit\fit_mixed_gaussian('y ~ x + (1 | f)', $train),
                fn () => \Phit\fit_mixed_logistic('b ~ x + (1 | f)', $train),
            ] as $fit) {
                try {
                    $fits[] = $fit();
                } catch (\Throwable) {
                    // An unidentifiable random design; the fit errors match elsewhere.
                }
            }
        } finally {
            restore_error_handler();
        }
        foreach ($fits as $fit) {
            foreach ([[], ['type' => 'response', 'seFit' => true], ['interval' => 'prediction', 'level' => 0.9], ['random' => false], ['allowNewLevels' => true]] as $options) {
                randomParityCompare('Phit\predict', [$fit, $new, ...$options]);
            }
        }
    }
    expect(true)->toBeTrue();
});

it('prints summaries and tidies fits exactly like the extension', function () {
    // Synthetic fits of every kind, with numbers chosen to hit R's formatting
    // edges (ties, widening at a power of ten, tiny p-values, zeros, and
    // non-finite values) and keys removed or given the wrong type.
    randomParitySeed('summary');
    $number = function (): mixed {
        return match (mt_rand(0, 14)) {
            0 => 0.0,
            1 => -0.0,
            2 => NAN,
            3 => INF,
            4 => null,
            5 => [9999.6, 999.96, 99999.6, 99996.0, 0.00005, 0.125, 2.5, 1e-300, 1e300][mt_rand(0, 8)],
            6 => mt_rand(-50, 50),
            7 => randomUniform(-1, 1) * 10 ** mt_rand(-320, 300),
            default => randomUniform(-1, 1) * 10 ** mt_rand(-12, 6),
        };
    };
    $pvalue = fn (): mixed => match (mt_rand(0, 6)) {
        0 => 10 ** -mt_rand(16, 300) * randomUniform(0, 1),
        1 => NAN,
        2 => 10 ** -mt_rand(3, 6) * randomUniform(0.5, 1),
        default => randomUniform(0, 1),
    };
    // range(1, 0) counts down, so build lists with a loop.
    $list = function (int $n, \Closure $f): array {
        $out = [];
        for ($i = 0; $i < $n; $i++) {
            $out[] = $f($i);
        }

        return $out;
    };
    $maybe = fn (mixed $v): mixed => mt_rand(0, 25) === 0 ? ['x', true, '1', null][mt_rand(0, 3)] : $v;
    for ($round = 0; $round < randomParityIterations() * 4; $round++) {
        $p = mt_rand(0, 6);
        $fit = [
            'coefficients' => $list($p, fn () => $number()),
            'std_errors' => $list($p, fn () => mt_rand(0, 3) === 0 ? $number() : abs(randomUniform(0, 1) * 10 ** mt_rand(-8, 4))),
            (mt_rand(0, 1) ? 't_statistics' : 'z_statistics') => $list($p, fn () => $number()),
            'p_values' => $list($p, fn () => $pvalue()),
            'df_residual' => mt_rand(0, 9) === 0 ? (float) mt_rand(0, 40) : mt_rand(0, 400),
            'n_observations' => mt_rand(1, 500),
        ];
        if (mt_rand(0, 2) > 0) {
            $fit['columns'] = $list($p, fn ($i) => ['(Intercept)', 'wt', 'hp', 'factor(cyl)6', 'x:z', 'é'][$i % 6].($i >= 6 ? $i : ''));
        }
        foreach (['formula' => 'y ~ x + (1 | g)', 'family' => 'binomial', 'link' => 'logit'] as $key => $value) {
            if (mt_rand(0, 2) > 0) {
                $fit[$key] = $value;
            }
        }
        if (mt_rand(0, 3) > 0) {
            $fit['intercept'] = mt_rand(0, 1) === 1;
        }
        $kind = mt_rand(0, 3);
        $fit['model'] = mt_rand(0, 30) === 0 ? ['quantile', null, 3][mt_rand(0, 2)]
            : ['linear', 'generalized_linear', $fit['family'] ?? null ? 'mixed_logistic' : 'mixed_gaussian'][min($kind, 2)];
        if ($kind === 0) {
            $fit += ['sigma' => $number(), 'r_squared' => randomUniform(0, 1), 'adj_r_squared' => randomUniform(-0.1, 1)];
        } elseif ($kind === 1) {
            $fit += ['null_deviance' => abs(randomUniform(0, 1)) * 10 ** mt_rand(-2, 6), 'deviance' => abs(randomUniform(0, 1)) * 10 ** mt_rand(-2, 6),
                'aic' => randomUniform(0, 1) * 10 ** mt_rand(0, 6), 'residual_variance' => $number()];
        } elseif ($kind === 2) {
            $fit += ['log_likelihood' => randomUniform(-1, 0) * 10 ** mt_rand(0, 5), 'method' => mt_rand(0, 1) ? 'REML' : 'ML',
                'residual_variance' => abs($number() ?? 1.0), 'n_groups' => mt_rand(1, 30)];
            foreach (['aic', 'bic'] as $key) {
                if (mt_rand(0, 2) > 0) {
                    $fit[$key] = randomUniform(0, 1) * 10 ** mt_rand(0, 5);
                }
            }
            if (mt_rand(0, 2) === 0 || $fit['model'] === 'mixed_logistic') {
                $fit['integration'] = mt_rand(0, 1) ? 'laplace' : 'agq';
                $fit['quadrature_points'] = mt_rand(1, 25);
            }
            $fit[mt_rand(0, 1) ? 'singular' : 'boundary'] = mt_rand(0, 1) === 1;
            if (mt_rand(0, 2) === 0) {
                $fit['random_intercept_variance'] = abs(randomUniform(0, 2));
                if (mt_rand(0, 3) > 0) {
                    $groups = $fit['n_groups'];
                    $fit['random_effects'] = $list($groups, fn () => randomUniform(-2, 2));
                    $fit['random_effects_cond_var'] = $list($groups + (mt_rand(0, 12) === 0 ? 1 : 0), fn () => randomUniform(0, 1));
                    $fit['group_labels'] = $list($groups + (mt_rand(0, 12) === 0 ? -1 : 0), fn ($i) => mt_rand(0, 1) ? $i : "s$i");
                }
            } else {
                $fit['random'] = [];
                for ($t = mt_rand(1, 3); $t > 0; $t--) {
                    $q = mt_rand(1, 3);
                    $covariance = randomCovariance($q);
                    if (mt_rand(0, 15) === 0) {
                        // A covariance smaller than its columns, or ragged.
                        $covariance = mt_rand(0, 1) ? [[1.0]] : [...array_slice($covariance, 0, -1), [1.0]];
                    }
                    $term = ['covariance' => $covariance, 'group_labels' => $list(mt_rand(1, 20), fn ($i) => "g$i")];
                    if (mt_rand(0, 3) > 0) {
                        $term['group'] = ['Subject', 'herd', 'g'][mt_rand(0, 2)];
                    }
                    if (mt_rand(0, 2) > 0) {
                        $term['columns'] = array_slice(['(Intercept)', 'Days', 'x'], 0, $q);
                    }
                    if (mt_rand(0, 2) === 0) {
                        $term['correlated'] = mt_rand(0, 1) === 1;
                    }
                    if (mt_rand(0, 4) > 0) {
                        // Conditional modes, sometimes the wrong size, with string, int, or float labels.
                        $groups = count($term['group_labels']) + (mt_rand(0, 12) === 0 ? 1 : 0);
                        $term['effects'] = $list($groups, fn () => $list($q + (mt_rand(0, 20) === 0 ? 1 : 0), fn () => randomUniform(-2, 2)));
                        $term['effects_cond_var'] = $list($groups + (mt_rand(0, 12) === 0 ? -1 : 0), fn () => randomCovariance($q));
                        $term['group_labels'] = $list(count($term['group_labels']), fn ($i) => [
                            "g$i", $i, $i + 0.5, randomUniform(-1, 1) * 10 ** mt_rand(-30, 30), (float) $i,
                        ][mt_rand(0, 4)]);
                    }
                    $fit['random'][] = $term;
                }
            }
        }
        $key = array_rand($fit);
        $fit[$key] = $maybe($fit[$key]);
        randomParityCompare('Phit\summary', [$fit]);
        randomParityCompare('Phit\tidy', [$fit, mt_rand(0, 1) === 1, [0.95, 0.5, 0.99, 1.0, 0.0, NAN][mt_rand(0, 5)]]);
        randomParityCompare('Phit\tidy', [$fit, mt_rand(0, 1) === 1, [0.95, 0.5, 1.0][mt_rand(0, 2)], ['ran_vals', 'ran_vals', 'fixed', 'random'][mt_rand(0, 3)]]);
        if (mt_rand(0, 9) === 0) {
            unset($fit['model']);
            randomParityCompare('Phit\summary', [$fit]);
        }
    }

    // Real fits from random designs.
    randomParitySeed('summary-fits');
    for ($round = 0; $round < randomParityIterations(); $round++) {
        $n = mt_rand(8, 80);
        [$x, $y, $binary, $counts, $labels] = randomDesign($n, mt_rand(1, 4), 4);
        $intercept = mt_rand(0, 5) > 0;
        set_error_handler(fn () => true);
        try {
            foreach ([
                fn () => Phit\fit_ols_matrix($x, $y, intercept: $intercept),
                fn () => Phit\fit_logistic_matrix($x, $binary, intercept: $intercept),
                fn () => Phit\fit_poisson_matrix($x, $counts, intercept: $intercept),
                fn () => Phit\fit_probit_matrix($x, $binary, intercept: $intercept),
                fn () => Phit\fit_mixed_gaussian_intercept_matrix($x, $y, $labels),
                fn () => Phit\fit_mixed_gaussian_matrix($x, $y, [['z' => array_map(fn ($row) => [1.0, $row[0]], $x), 'groups' => $labels]]),
                fn () => Phit\fit_mixed_logistic_matrix($x, $binary, $labels),
            ] as $make) {
                try {
                    $fit = $make();
                } catch (\Exception) {
                    continue;
                }
                randomParityCompare('Phit\summary', [$fit]);
                randomParityCompare('Phit\tidy', [$fit, true, randomUniform(0.5, 0.999)]);
                randomParityCompare('Phit\tidy', [$fit, true, randomUniform(0.5, 0.999), 'ran_vals']);
            }
        } finally {
            restore_error_handler();
        }
    }
    expect(true)->toBeTrue();
});
