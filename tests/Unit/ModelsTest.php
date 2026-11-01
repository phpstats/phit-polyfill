<?php

declare(strict_types=1);

// Ports of the extension's tests/ols.php, logistic.php, probit.php, poisson.php,
// r-compatibility.php, stability.php and intercept.php, run against the polyfill.

const MODEL_FIXTURES = [
    'Phit\fit_ols_matrix' => 'ols-r.json',
    'Phit\fit_logistic_matrix' => 'logistic-r.json',
    'Phit\fit_probit_matrix' => 'probit-r.json',
    'Phit\fit_poisson_matrix' => 'poisson-r.json',
];

/** Call a polyfilled model, collecting the warnings it raises. */
function fitQuietly(string $fn, mixed ...$arguments): array
{
    warnings(fn () => pf($fn, ...$arguments), $result);

    return $result;
}

function expectCode(callable $f, int $code): void
{
    expect(thrownCode($f))->toBe($code);
}

function expectThrows(callable $f): void
{
    expect(thrownCode($f))->not->toBeNull();
}

/** Select rows of a list by original position. */
function withoutRows(array $values, array $drop): array
{
    return array_values(array_filter($values, fn ($i) => ! in_array($i, $drop, true), ARRAY_FILTER_USE_KEY));
}

describe('R reference fits', function () {
    foreach (MODEL_FIXTURES as $fn => $file) {
        foreach (loadFixture($file)['cases'] as $case) {
            it("$fn matches R on {$case['name']}", function () use ($fn, $case) {
                $result = pf($fn, x: $case['x'], y: $case['y'], intercept: false);
                foreach ($case['expected'] as $key => $expected) {
                    // Perfect and all-zero OLS fits have undefined inference; checked separately.
                    if ($fn === 'Phit\fit_ols_matrix' && in_array($case['name'], ['perfect', 'zero'], true)
                        && in_array($key, ['t_statistics', 'p_values', 'log_likelihood', 'r_squared', 'adj_r_squared', 'conf_int'], true)) {
                        continue;
                    }
                    expect($result[$key])->toBeClose($expected, 1e-9, 1e-7);
                }
                if ($fn !== 'Phit\fit_ols_matrix') {
                    expect($result['converged'])->toBeTrue()
                        ->and($result['iterations'])->toBeGreaterThan(0)->toBeLessThanOrEqual(100);
                }
            });
        }
    }

    it('gives explicit undefined inference for an all-zero outcome', function () {
        $case = array_values(array_filter(loadFixture('ols-r.json')['cases'], fn ($c) => $c['name'] === 'zero'))[0];
        $result = pf('Phit\fit_ols_matrix', $case['x'], $case['y'], intercept: false);
        expect(is_nan($result['r_squared']))->toBeTrue()
            ->and(is_nan($result['t_statistics'][0]))->toBeTrue()
            ->and($result['log_likelihood'])->toBe(INF)
            ->and($result['conf_int'])->toBe([[0.0, 0.0], [0.0, 0.0]]);
    });

    foreach (loadFixture('r-compatibility.json')['cases'] as $case) {
        it("Phit\\fit_{$case['model']}_matrix matches base R with weights and offsets on {$case['name']}", function () use ($case) {
            $result = pf('Phit\fit_'.$case['model'].'_matrix', x: $case['x'], y: $case['y'], weights: $case['weights'], offset: $case['offset'], intercept: false);
            foreach ($case['expected'] as $key => $expected) {
                expect($result[$key])->toBeClose($expected, 1e-7, 1e-6);
            }
        });
    }

    foreach (loadFixture('stress-r.json')['cases'] as $case) {
        it("Phit\\fit_{$case['model']}_matrix is stable on {$case['name']}", function () use ($case) {
            $result = pf('Phit\fit_'.$case['model'].'_matrix', $case['x'], $case['y'], weights: $case['weights'], offset: $case['offset'], intercept: false);
            foreach ($case['expected'] as $key => $expected) {
                $actualKey = $key === 'statistics' ? ($case['model'] === 'ols' ? 't_statistics' : 'z_statistics') : $key;
                $factor = match ($key) {
                    'covariance' => $case['response_scale'] ** 2,
                    'coefficients', 'std_errors', 'conf_int', 'fitted', 'residuals' => $case['response_scale'],
                    default => 1,
                };
                expect($result[$actualKey])->toBeClose($expected, $case['tolerance']['absolute'] * $factor, $case['tolerance']['relative']);
            }
            expect($result['diagnostics']['converged'])->toBeTrue()
                ->and($result['diagnostics']['warnings'])->toBe([]);
        });
    }
});

describe('validation', function () {
    it('rejects invalid OLS inputs', function (array $x, array $y) {
        expectThrows(fn () => pf('Phit\fit_ols_matrix', $x, $y, intercept: false));
    })->with([
        [[], []], [[[]], [1]], [[[1], [1, 2]], [1, 2]],
        [[[1], [2]], [1]], [[[1, 2], [1, 3]], [1, 2]],
        [[[1, 1], [1, 1], [1, 1]], [1, 2, 3]],
        [[[NAN], [1]], [1, 2]], [[[1], [2]], [INF, 2]],
        [[[true], [2]], [1, 2]], [[['1'], [2]], [1, 2]],
        [[1, 2], [1, 2]], [[[1], [2]], [null, 2]],
    ]);

    it('rejects invalid binary inputs and separation', function (string $fn, array $x, array $y) {
        expectThrows(fn () => pf($fn, x: $x, y: $y, intercept: false));
    })->with(['Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix'])->with([
        [[], []], [[[]], [0]], [[[1], [1, 2], [1]], [0, 1, 0]],
        [[[1], [1], [1]], [0, 0.5, 1]], [[[1], [1], [1]], [0, true, 1]],
        [[[1], [1], [1]], [0, '1', 1]], [[[1], [1], [1]], [0, INF, 1]],
        [[[1], [1], [1]], [0, 0, 0]], [[[1, 1], [1, 1], [1, 1]], [0, 1, 0]],
        [[[1, -2], [1, -1], [1, 1], [1, 2]], [0, 0, 1, 1]],
        [[[1, -2], [1, -1], [1, 0], [1, 0], [1, 1], [1, 2]], [0, 0, 0, 1, 1, 1]],
        [[[null], [1], [1], [1]], [2, 0, 1, 0]],
        [[[null, INF], [1, 1], [1, 2], [1, 3]], [0, 0, 1, 0]],
        [[[1], [null], [1]], [0, 1, 0]],
    ]);

    it('rejects invalid Poisson inputs and boundary fits', function (array $x, array $y) {
        expectThrows(fn () => pf('Phit\fit_poisson_matrix', x: $x, y: $y, intercept: false));
    })->with([
        [[], []], [[[]], [0]], [[[1], [1, 2], [1]], [0, 1, 0]],
        [[[1], [1], [1]], [0, -1, 1]], [[[1], [1], [1]], [0, true, 1]],
        [[[1], [1], [1]], [0, '1', 1]], [[[1], [1], [1]], [0, INF, 1]],
        [[[1], [1], [1]], [0, 0, 0]], [[[1, 1], [1, 1], [1, 1]], [0, 1, 0]],
        [[[1, 0], [1, 0], [1, 0], [1, 1], [1, 1], [1, 1]], [1, 2, 3, 0, 0, 0]],
        [[[null], [1], [1], [1]], [0.5, 0, 1, 0]],
        [[[null], [1], [1], [1]], [-2, 0, 1, 0]],
        [[[null, INF], [1, 1], [1, 2], [1, 3]], [0, 0, 1, 0]],
        [[[1], [null], [1]], [0, 1, 0]],
    ]);

    it('validates every row before deletion can hide it', function (array $x, array $y, string $policy) {
        expectThrows(fn () => pf('Phit\fit_ols_matrix', x: $x, y: $y, missing: $policy, intercept: false));
    })->with([
        [[[1, 1], [1, null], [1, 2], [1, 3], [1, 4], [1, NAN], [1, 5]], [2.1, 99, 4.2, null, 8.1, 99, 10.2], 'raise'],
        [[[1, 1], [1, null], [1, 2], [1, 3], [1, 4], [1, NAN], [1, 5]], [2.1, 99, 4.2, null, 8.1, 99, 10.2], 'unknown'],
        [[[1, null], [1, 2]], [null, null], 'drop'],
        [[[1, null], [1, 2], [1, 3]], [1, 2, 3], 'drop'],
        [[[1, null], [1, 1], [1, 1], [1, 1]], [1, 2, 3, 4], 'drop'],
        [[[null, 'bad'], [1, 2], [1, 3], [1, 4]], [1, 2, 3, 4], 'drop'],
        [[[null, INF], [1, 2], [1, 3], [1, 4]], [1, 2, 3, 4], 'drop'],
        [[[1, null], [1], [1, 3], [1, 4]], [1, 2, 3, 4], 'drop'],
    ]);

    it('uses the extension exception categories', function (callable $call, int $code) {
        expectCode($call, $code);
    })->with([
        'overflowing weighted outcome' => [fn () => pf('Phit\fit_ols_matrix', [[1, 1], [1, 2], [1, 3], [1, 4]], [1e155, 2e155, 3e155, 4e155], weights: [1, 2, 3, 4], intercept: false), 1003],
        'string outcome' => [fn () => pf('Phit\fit_ols_matrix', [[1], [1], [1]], ['oops', 2, 3], intercept: false), 1001],
        'short outcome' => [fn () => pf('Phit\fit_ols_matrix', [[1], [1], [1]], [1], intercept: false), 1002],
        'infinite outcome' => [fn () => pf('Phit\fit_ols_matrix', [[1], [1], [1]], [1, INF, 3], intercept: false), 1003],
        'duplicate columns' => [fn () => pf('Phit\fit_ols_matrix', [[1, 1], [1, 1], [1, 1]], [1, 2, 3], intercept: false), 1004],
        'separation' => [fn () => pf('Phit\fit_logistic_matrix', [[1, -2], [1, -1], [1, 1], [1, 2]], [0, 0, 1, 1], intercept: false), 1006],
        'too many trials' => [fn () => pf('Phit\fit_logistic_matrix', array_fill(0, 4, [1]), [0, 1, 0, 1], weights: array_fill(0, 4, 1000001), intercept: false), 1007],
        'underflowing variance' => [fn () => pf('Phit\fit_ols_matrix', array_fill(0, 4, [1]), [1e-170, 2e-170, 3e-170, 4e-170], intercept: false), 1003],
        'logistic y below 0' => [fn () => pf('Phit\fit_logistic_matrix', [[null], [1], [1], [1]], [-0.1, 1, 0, 1], intercept: false), 1001],
        'logistic y above 1' => [fn () => pf('Phit\fit_logistic_matrix', [[null], [1], [1], [1]], [1.1, 1, 0, 1], intercept: false), 1001],
        'probit y below 0' => [fn () => pf('Phit\fit_probit_matrix', [[null], [1], [1], [1]], [-0.1, 1, 0, 1], intercept: false), 1001],
        'probit y above 1' => [fn () => pf('Phit\fit_probit_matrix', [[null], [1], [1], [1]], [1.1, 1, 0, 1], intercept: false), 1001],
        'negative count' => [fn () => pf('Phit\fit_poisson_matrix', [[null], [1], [1], [1]], [-1, 1, 2, 3], intercept: false), 1001],
        'fractional count' => [fn () => pf('Phit\fit_poisson_matrix', [[null], [1], [1], [1]], [0.5, 1, 2, 3], intercept: false), 1001],
        'all rows dropped' => [fn () => pf('Phit\fit_ols_matrix', [[null], [null]], [1, 2]), 1003],
        'missing policy' => [fn () => pf('Phit\fit_logistic_matrix', [[1]], [1], missing: 'bogus'), 1001],
    ]);

    it('canonicalizes binary counts and requires both classes', function (string $fn) {
        $x = array_fill(0, 8, [1]);
        $y = [0, 1, 0, 1, 0, 1, 0, 1];
        $near = array_map(fn ($v) => $v === 0 ? 1e-10 : 1 - 1e-10, $y);
        expect(pf($fn, $x, $y, intercept: false)['coefficients'])->toBe(pf($fn, $x, $near, intercept: false)['coefficients']);
        expectCode(fn () => pf($fn, $x, array_fill(0, 8, 1e-10), intercept: false), 1001);
    })->with(['Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix']);

    it('rejects a non-array x', function () {
        expect(fn () => pf('Phit\fit_ols_matrix', 'bad', [], intercept: false))->toThrow(TypeError::class);
    });
});

describe('listwise deletion', function () {
    it('drops incomplete rows, keeps original positions and warns once', function (string $fn) {
        $file = MODEL_FIXTURES[$fn];
        if ($fn === 'Phit\fit_ols_matrix') {
            $x = [[1, 1], [1, null], [1, 2], [1, 3], [1, 4], [1, NAN], [1, 5]];
            $y = [2.1, 99, 4.2, null, 8.1, 99, 10.2];
            $drop = [1, 3, 5];
        } else {
            $case = loadFixture($file)['cases'][0];
            $x = $case['x'];
            $y = $case['y'];
            $x[1][1] = null;
            $y[4] = NAN;
            $drop = [1, 4];
        }
        $n = count($y);
        $warnings = warnings(fn () => pf($fn, $x, $y, intercept: false), $result);
        $complete = pf($fn, withoutRows($x, $drop), withoutRows($y, $drop), intercept: false);
        expect($result['coefficients'])->toBeClose($complete['coefficients'])
            ->and($result['fitted'])->toBeClose($complete['fitted'])
            ->and($result['dropped_rows'])->toBe($drop)
            ->and($result['retained_rows'])->toBe(array_values(array_diff(range(0, $n - 1), $drop)))
            ->and($result['n_input'])->toBe($n)
            ->and($result['n_observations'])->toBe($n - count($drop))
            ->and($result['n_dropped'])->toBe(count($drop))
            ->and($result['df_residual'])->toBe($n - count($drop) - count($x[0]))
            ->and($warnings)->toHaveCount(1)
            ->and($warnings[0][0])->toBe(E_USER_WARNING)
            ->and($warnings[0][1])->toBe("$fn(): dropped ".count($drop)." of $n observations with missing values (listwise deletion)")
            ->and($result['diagnostics']['warnings'])->toBe(['listwise deletion removed '.count($drop)." of $n original rows"])
            ->and($complete['n_dropped'])->toBe(0)
            ->and($complete['dropped_rows'])->toBe([]);

        foreach (['raise', 'unknown'] as $policy) {
            expectThrows(fn () => pf($fn, $x, $y, $policy, intercept: false));
        }
        // An application can promote the warning to an exception.
        set_error_handler(function (int $severity, string $message): never {
            throw new ErrorException($message, 0, $severity);
        });
        try {
            expect(fn () => pf($fn, $x, $y, intercept: false))->toThrow(ErrorException::class);
        } finally {
            restore_error_handler();
        }
    })->with(array_keys(MODEL_FIXTURES));

    it('deletes rows with missing weights or offsets', function (string $model) {
        $c = array_values(array_filter(loadFixture('r-compatibility.json')['cases'], fn ($c) => $c['model'] === $model))[0];
        $c['weights'][3] = null;
        $c['offset'][5] = NAN;
        $fn = 'Phit\fit_'.$model.'_matrix';
        $warnings = warnings(fn () => pf($fn, $c['x'], $c['y'], weights: $c['weights'], offset: $c['offset'], intercept: false), $result);
        expect($result['dropped_rows'])->toBe([3, 5])->and($warnings)->toHaveCount(1);
        $complete = pf($fn, withoutRows($c['x'], [3, 5]), withoutRows($c['y'], [3, 5]), weights: withoutRows($c['weights'], [3, 5]), offset: withoutRows($c['offset'], [3, 5]), intercept: false);
        expect($result['coefficients'])->toBeClose($complete['coefficients'], 1e-7, 1e-6);
        $base = ['x' => $c['x'], 'y' => $c['y'], 'weights' => $c['weights'], 'offset' => $c['offset'], 'intercept' => false];
        foreach ([['weights' => [0]], ['weights' => array_fill(0, count($c['y']), -1)], ['offset' => [0]], ['weights' => ['1']], ['missing' => 'raise']] as $options) {
            expectThrows(fn () => pf($fn, ...array_replace($base, $options)));
        }
    })->with(['ols', 'logistic', 'probit', 'poisson']);
});

describe('result contract', function () {
    it('returns exactly the extension result keys', function (string $fn) {
        $case = loadFixture(MODEL_FIXTURES[$fn])['cases'][0];
        $result = pf($fn, $case['x'], $case['y'], intercept: false);
        $keys = array_keys($result);
        sort($keys);
        expect($keys)->toBe(loadFixture('php-api.json')[$fn]['result_keys'])
            ->and($result['diagnostics']['converged'])->toBeTrue()
            ->and($result['diagnostics']['warnings'])->toBeArray();
        if (array_key_exists('iterations', $result)) {
            expect($result['diagnostics']['iterations'])->toBe($result['iterations']);
        }
    })->with(array_keys(MODEL_FIXTURES));

    it('reports the extension solver names', function () {
        $x = [[1], [2], [3], [4], [5], [6]];
        expect(pf('Phit\fit_ols_matrix', $x, [1, 3, 2, 5, 4, 6])['diagnostics']['solver'])->toBe('scaled_svd')
            ->and(pf('Phit\fit_logistic_matrix', $x, [0, 1, 0, 1, 1, 0])['diagnostics']['solver'])->toBe('newton_irls')
            ->and(pf('Phit\fit_poisson_matrix', $x, [0, 1, 0, 1, 1, 0])['diagnostics']['solver'])->toBe('newton_irls')
            ->and(pf('Phit\fit_probit_matrix', $x, [0, 1, 0, 1, 1, 0])['diagnostics']['solver'])->toBe('fisher_scoring');
    });
});

describe('automatic intercept', function () {
    it('matches a literal ones column exactly and does not mutate inputs', function (string $fn) {
        $case = loadFixture(MODEL_FIXTURES[$fn])['cases'][0];
        $literal = pf($fn, $case['x'], $case['y'], intercept: false);
        $x = array_map(fn ($row) => array_slice($row, 1), $case['x']);
        $copy = $x;
        expect(pf($fn, $x, $case['y']))->toEqual($literal)
            ->and($x)->toBe($copy);
    })->with(array_keys(MODEL_FIXTURES));

    it('fits intercept-only models from empty rows', function () {
        $x = array_fill(0, 6, []);
        $r = pf('Phit\fit_ols_matrix', $x, [2, 3, 7, 4, 8, 6]);
        expect($r['coefficients'])->toHaveCount(1)
            ->and(abs($r['coefficients'][0] - 5))->toBeLessThan(1e-12)
            ->and($r['df_residual'])->toBe(5);
        foreach (['Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix', 'Phit\fit_poisson_matrix'] as $fn) {
            $y = $fn === 'Phit\fit_poisson_matrix' ? [1, 3, 2, 5, 4, 2] : [0, 1, 0, 0, 1, 0];
            expect(pf($fn, $x, $y))->toEqual(pf($fn, array_fill(0, 6, [1]), $y, intercept: false));
        }
    });

    it('keeps rank errors for a duplicated intercept and rejects empty designs', function () {
        expectCode(fn () => pf('Phit\fit_ols_matrix', array_fill(0, 6, [1]), [2, 3, 7, 4, 8, 6]), 1004);
        foreach ([[[], [], true], [[[1], []], [1, 2], true], [array_fill(0, 6, []), [1, 2, 3, 4, 5, 6], false]] as [$x, $y, $intercept]) {
            expect(thrownCode(fn () => pf('Phit\fit_ols_matrix', $x, $y, intercept: $intercept)))->toBeIn([1001, 1002]);
        }
    });

    it('uses uncentered R-squared without an intercept', function () {
        $y = [5, 6, 9, 15, 20];
        $r = pf('Phit\fit_ols_matrix', [[2], [3], [5], [7], [11]], $y, intercept: false);
        $tss = array_sum(array_map(fn ($v) => $v * $v, $y));
        expect($r['coefficients'])->toHaveCount(1)
            ->and($r['df_residual'])->toBe(4)
            ->and(abs($r['r_squared'] - (1 - $r['rss'] / $tss)))->toBeLessThan(1e-12);
    });

    it('reproduces the closed-form Poisson exposure fit', function (array $weights) {
        $exposure = [1, 2, 4, 8, 16, 3, 6, 12];
        $counts = [1, 4, 3, 12, 23, 2, 11, 14];
        $offset = array_map(log(...), $exposure);
        $events = $time = 0;
        foreach ($counts as $i => $count) {
            $events += $weights[$i] * $count;
            $time += $weights[$i] * $exposure[$i];
        }
        $rate = $events / $time;
        $r = pf('Phit\fit_poisson_matrix', array_fill(0, 8, []), $counts, weights: $weights, offset: $offset);
        expect(abs($r['coefficients'][0] - log($rate)))->toBeLessThan(1e-10)
            ->and(abs($r['std_errors'][0] - 1 / sqrt($events)))->toBeLessThan(1e-10);
        $ll = 0;
        foreach ($counts as $i => $count) {
            $mu = $rate * $exposure[$i];
            $logFactorial = array_sum(array_map(log(...), range(1, $count)));
            $ll += $weights[$i] * ($count * log($mu) - $mu - $logFactorial);
            expect(abs($r['fitted'][$i] - $mu))->toBeLessThan(1e-9);
        }
        expect(abs($r['log_likelihood'] - $ll))->toBeLessThan(1e-9);
        $scaled = pf('Phit\fit_poisson_matrix', array_fill(0, 8, []), $counts, weights: $weights, offset: array_map(fn ($o) => $o + log(2), $offset));
        expect(abs($scaled['coefficients'][0] - $r['coefficients'][0] + log(2)))->toBeLessThan(1e-10);
        foreach ($r['fitted'] as $i => $mu) {
            expect(abs($scaled['fitted'][$i] - $mu))->toBeLessThan(1e-9);
        }
    })->with([[array_fill(0, 8, 1)], [[1, 2, 1, 3, 2, 1, 2, 1]]]);

    it('matches the R weighted Poisson fixture with the default intercept', function () {
        foreach (loadFixture('r-compatibility.json')['cases'] as $case) {
            if ($case['model'] !== 'poisson' || $case['x'][0][0] !== 1) {
                continue;
            }
            $r = pf('Phit\fit_poisson_matrix', array_map(fn ($row) => array_slice($row, 1), $case['x']), $case['y'], weights: $case['weights'], offset: $case['offset']);
            foreach ($case['expected'] as $key => $value) {
                expect($r[$key])->toBeClose($value, 1e-7, 1e-6);
            }
        }
    });

    it('drops a missing offset and enforces missing=raise', function () {
        $offset = array_map(log(...), [1, 2, 4, 8, 16, 3, 6, 12]);
        $offset[2] = null;
        $counts = [1, 4, 3, 12, 23, 2, 11, 14];
        $warnings = warnings(fn () => pf('Phit\fit_poisson_matrix', array_fill(0, 8, []), $counts, offset: $offset), $r);
        expect($r['dropped_rows'])->toBe([2])->and($r['df_residual'])->toBe(6)->and($warnings)->toHaveCount(1);
        expectCode(fn () => pf('Phit\fit_poisson_matrix', array_fill(0, 8, []), $counts, missing: 'raise', offset: $offset), 1001);
    });
});

describe('rounding near the optimum', function () {
    // Plain summation of the likelihood accumulates systematic rounding in
    // intercept-only fits, which made the exact Newton step look like a decrease
    // (extension issue fixed in 4ba33e4 with compensated sums).
    it('fits intercept-only binary and count models on large samples', function (string $fn) {
        mt_srand(3);
        foreach ([1000, 10000, 20000] as $n) {
            for ($trial = 0; $trial < 5; $trial++) {
                $q = mt_rand(5, 95) / 100;
                $y = [];
                for ($i = 0; $i < $n; $i++) {
                    $y[] = mt_rand() / mt_getrandmax() < $q ? 1 : 0;
                }
                $result = pf($fn, array_fill(0, $n, []), $y);
                $mean = array_sum($y) / $n;
                // The intercept-only MLE reproduces the sample mean on the response scale.
                expect(abs($result['fitted'][0] - $mean))->toBeLessThan(1e-9);
            }
        }
    })->with(['Phit\fit_logistic_matrix', 'Phit\fit_probit_matrix', 'Phit\fit_poisson_matrix']);
});

describe('extreme magnitudes', function () {
    // Review finding: the Jacobi SVD's rotation test multiplied two squared
    // column norms, which overflowed for heavily weighted IRLS steps and stopped
    // every rotation, giving wrong fits or spurious 1006 errors.
    $x = [[1.0], [2.0], [3.0], [4.5], [5.0], [6.0]];
    $y = [1, 2, 1, 4, 3, 6];

    it('gives the unweighted Poisson fit for equal huge weights', function () use ($x, $y) {
        $plain = pf('Phit\fit_poisson_matrix', $x, $y);
        $weighted = pf('Phit\fit_poisson_matrix', $x, $y, weights: array_fill(0, 6, 1e160));
        expect($weighted['coefficients'])->toBeClose($plain['coefficients'], 1e-9, 1e-9);
    });

    it('shifts only the Poisson intercept when counts are scaled by 1e160', function () use ($x, $y) {
        $plain = pf('Phit\fit_poisson_matrix', $x, $y);
        $scaled = pf('Phit\fit_poisson_matrix', $x, array_map(fn ($v) => $v * 1e160, $y));
        // The log-likelihood is about 1e162 here, so the fit converges only to about
        // 1e-8 relative (as in the extension); compare at the extension's R tolerance.
        expect($scaled['coefficients'][1])->toBeClose($plain['coefficients'][1], 1e-9, 1e-7)
            ->and($scaled['coefficients'][0])->toBeClose($plain['coefficients'][0] + log(1e160), 1e-9, 1e-7);
    });

    it('rejects a single dominating weight as the extension does', function () use ($x, $y) {
        expect(thrownCode(fn () => pf('Phit\fit_poisson_matrix', $x, $y, weights: [1e160, 1, 1, 1, 1, 1])))->toBe(1006);
    });
});
