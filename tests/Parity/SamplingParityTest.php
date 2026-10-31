<?php

declare(strict_types=1);

// Compares the uniform, categorical and permutation samplers with the loaded
// extension on the same inputs: identical draws, or the same exception.

function samplingCompare(string $fn, array $args): void
{
    $run = function (callable $call) use ($fn, $args): array {
        try {
            return ['ok', $call($fn, ...$args)];
        } catch (\Throwable $e) {
            return ['error', get_class($e), $e->getCode(), $e->getMessage()];
        }
    };
    $native = $run(native(...));
    $polyfill = $run(pf(...));
    parityLog('sampling', $fn, $polyfill, $native, $native[0] === 'error' || $polyfill[0] === 'error');
    if (serialize($native) !== serialize($polyfill)) {
        throw new \RuntimeException("$fn differs for ".json_encode($args)."\nnative ".var_export($native, true)."\npolyfill ".var_export($polyfill, true));
    }
    expect(true)->toBeTrue();
}

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

it('agrees on seeded random inputs', function () {
    mt_srand(20261031);
    for ($case = 0; $case < 200; $case++) {
        $seed = mt_rand(0, PHP_INT_MAX);
        $n = mt_rand(0, 300);
        $k = mt_rand(1, 8);
        $prob = [];
        for ($i = 0; $i < $k; $i++) {
            $prob[] = match (mt_rand(0, 3)) {
                0 => 0,
                1 => mt_rand(1, 9),
                default => mt_rand() / mt_getrandmax() * 10 ** mt_rand(-300, 300),
            };
        }
        $min = (mt_rand() / mt_getrandmax() - .5) * 10 ** mt_rand(-5, 300);
        samplingCompare('Phit\dist_uniform_rng', [$n, $seed, $min, $min + abs($min) * mt_rand(0, 3)]);
        samplingCompare('Phit\dist_categorical_rng', [$n, $seed, $prob]);
        $rows = [];
        for ($i = 0; $i < $n; $i++) {
            $rows[] = array_map(fn ($w) => mt_rand(0, 9) ? $w * mt_rand(1, 3) : 0, $prob);
        }
        samplingCompare('Phit\dist_categorical_rows_rng', [$seed, $rows]);
        samplingCompare('Phit\permutation_rng', [$n, $seed, mt_rand(0, 1) ? null : mt_rand(0, $n + 1)]);
    }
});

it('agrees on edge cases and errors', function () {
    foreach ([
        ['Phit\dist_uniform_rng', [5, 1, -PHP_FLOAT_MAX, PHP_FLOAT_MAX]], ['Phit\dist_uniform_rng', [5, 1, 4, 4]],
        ['Phit\dist_uniform_rng', [5, 1, 1, 0]], ['Phit\dist_uniform_rng', [1000001, 1, NAN]], ['Phit\dist_uniform_rng', [1000001, 1]],
        ['Phit\dist_uniform_rng', [-1, -1, NAN]], ['Phit\dist_uniform_rng', [5, PHP_INT_MAX, -0.0, 0.0]],
        ['Phit\dist_categorical_rng', [5, 1, []]], ['Phit\dist_categorical_rng', [5, 1, [1 => 1]]],
        ['Phit\dist_categorical_rng', [5, 1, ['a' => 1]]], ['Phit\dist_categorical_rng', [5, 1, ['1']]],
        ['Phit\dist_categorical_rng', [5, 1, [null]]], ['Phit\dist_categorical_rng', [5, 1, [true]]],
        ['Phit\dist_categorical_rng', [5, 1, [INF, 'x']]], ['Phit\dist_categorical_rng', [5, 1, [-1, INF]]],
        ['Phit\dist_categorical_rng', [5, 1, [NAN]]], ['Phit\dist_categorical_rng', [5, 1, [-0.0, 0]]],
        ['Phit\dist_categorical_rng', [5, 1, [PHP_FLOAT_MAX, PHP_FLOAT_MAX]]], ['Phit\dist_categorical_rng', [1000001, 1, [1]]],
        ['Phit\dist_categorical_rng', [1000001, 1, ['x']]], ['Phit\dist_categorical_rng', [1000001, 1, [-1]]],
        ['Phit\dist_categorical_rng', [50, 3, [5e-324, 0, 0]]], ['Phit\dist_categorical_rng', [50, 3, [1e308, 1e308, 0]]],
        ['Phit\dist_categorical_rows_rng', [1, [[1], 1]]], ['Phit\dist_categorical_rows_rng', [1, [1 => [1]]]],
        ['Phit\dist_categorical_rows_rng', [1, [[1], [1 => 1]]]], ['Phit\dist_categorical_rows_rng', [1, [[-1], ['x']]]],
        ['Phit\dist_categorical_rows_rng', [1, [[1], []]]], ['Phit\dist_categorical_rows_rng', [-1, 1 => [1]]],
        ['Phit\dist_categorical_rows_rng', [1, [[1], [0, 0]]]], ['Phit\dist_categorical_rows_rng', [1, []]],
        ['Phit\permutation_rng', [3, 1, 4]], ['Phit\permutation_rng', [3, 1, -1]], ['Phit\permutation_rng', [-1, 1, -1]],
        ['Phit\permutation_rng', [1000001, 1, -1]], ['Phit\permutation_rng', [1000001, 1, 2000000]], ['Phit\permutation_rng', [1, 0]],
    ] as [$fn, $args]) {
        samplingCompare($fn, $args);
    }
});
