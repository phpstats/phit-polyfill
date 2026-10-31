<?php

declare(strict_types=1);

// Port of the "Uniform, categorical and permutation" section of the
// extension's tests/rng.php: the documented maps applied to the fixture's
// independent Python ChaCha20 uniforms, then validation in the extension's order.

$fixture = loadFixture('rng-python-r.json');

it('maps the reference uniforms', function () use ($fixture) {
    $uniform = [];
    foreach ($fixture['cases'] as $c) {
        if ($c['distribution'] === 'uniform') {
            $uniform[$c['seed']] = $c['expected'];
        }
    }
    expect($uniform)->not->toBe([]);
    foreach ($uniform as $seed => $u) {
        $u = array_slice($u, 0, 96);
        expect(pf('Phit\dist_uniform_rng', 96, $seed))->toBe($u, "runif seed $seed")
            ->and(pf('Phit\dist_uniform_rng', 96, $seed, min: -2, max: 3))->toBe(array_map(fn ($u) => -2 + 5 * $u, $u));
        $prob = [.2, 0, .5, .3];
        $draws = pf('Phit\dist_categorical_rng', 96, $seed, $prob);
        expect($draws)->toBe(array_map(fn ($u) => $u < .2 ? 0 : ($u < .7 ? 2 : 3), $u), "categorical seed $seed")
            ->and(pf('Phit\dist_categorical_rng', 96, $seed, [2, 0, 5, 3]))->toBe($draws)
            ->and(pf('Phit\dist_categorical_rows_rng', $seed, array_fill(0, 96, $prob)))->toBe($draws);
        $order = range(0, 49);
        for ($i = 0; $i < 50; $i++) {
            $j = $i + min((int) ($u[$i] * (50 - $i)), 49 - $i);
            [$order[$i], $order[$j]] = [$order[$j], $order[$i]];
        }
        expect(pf('Phit\permutation_rng', 50, $seed))->toBe($order, "permutation seed $seed")
            ->and(pf('Phit\permutation_rng', 50, $seed, size: 7))->toBe(array_slice($order, 0, 7));
    }
});

it('handles degenerate, empty and zero-weight cases', function () {
    expect(pf('Phit\dist_uniform_rng', 3, 1, min: 4, max: 4))->toBe([4.0, 4.0, 4.0])
        ->and(pf('Phit\dist_categorical_rows_rng', 9, [[1, 0], [0, 1], [0, 0, 1]]))->toBe([0, 1, 2])
        ->and(pf('Phit\dist_categorical_rows_rng', 9, []))->toBe([])
        ->and(pf('Phit\permutation_rng', 0, 9))->toBe([])
        ->and(pf('Phit\permutation_rng', 5, 9, size: 0))->toBe([])
        ->and(pf('Phit\dist_categorical_rng', 1000, 7, [0, 1, 0]))->toBe(array_fill(0, 1000, 1));
    $sorted = pf('Phit\permutation_rng', 1000, 5);
    sort($sorted);
    expect($sorted)->toBe(range(0, 999));
    $counts = array_count_values(pf('Phit\dist_categorical_rng', 50000, 175, [.1, .6, .3]));
    expect(abs($counts[0] / 50000 - .1))->toBeLessThan(.01)
        ->and(abs($counts[1] / 50000 - .6))->toBeLessThan(.01);
});

it('validates in the extension order', function () {
    $cases = [
        1001 => [
            fn () => pf('Phit\dist_uniform_rng', 0, 1, min: 1, max: 0), fn () => pf('Phit\dist_categorical_rng', 0, 1, []),
            fn () => pf('Phit\dist_categorical_rng', 0, 1, [-1, 2]), fn () => pf('Phit\dist_categorical_rng', 0, 1, [0, 0]),
            fn () => pf('Phit\dist_categorical_rng', 0, 1, ['a' => 1]), fn () => pf('Phit\dist_categorical_rng', 0, 1, ['1']),
            fn () => pf('Phit\dist_categorical_rows_rng', 1, [[1], 1]), fn () => pf('Phit\dist_categorical_rows_rng', 1, [[1], []]),
            fn () => pf('Phit\dist_categorical_rows_rng', -1, [[1]]), fn () => pf('Phit\permutation_rng', 3, 1, size: 4),
            fn () => pf('Phit\permutation_rng', 3, 1, size: -1), fn () => pf('Phit\permutation_rng', -1, 1),
            // An element error precedes the draw-count limit.
            fn () => pf('Phit\dist_categorical_rng', 1000001, 1, [1, 'x']),
        ],
        1003 => [
            fn () => pf('Phit\dist_uniform_rng', 0, 1, min: NAN), fn () => pf('Phit\dist_uniform_rng', 1, 1, min: -PHP_FLOAT_MAX, max: PHP_FLOAT_MAX),
            fn () => pf('Phit\dist_categorical_rng', 0, 1, [INF]), fn () => pf('Phit\dist_categorical_rng', 0, 1, [PHP_FLOAT_MAX, PHP_FLOAT_MAX]),
            // Parameters are checked before the draw-count limit.
            fn () => pf('Phit\dist_uniform_rng', 1000001, 1, max: INF),
        ],
        1007 => [
            fn () => pf('Phit\dist_uniform_rng', 1000001, 1), fn () => pf('Phit\dist_categorical_rng', 1000001, 1, [1]),
            fn () => pf('Phit\permutation_rng', 1000001, 1, size: 0), fn () => pf('Phit\dist_categorical_rng', 1000001, 1, [-1]),
            fn () => pf('Phit\permutation_rng', 1000001, 1, size: 2000000),
        ],
    ];
    foreach ($cases as $code => $calls) {
        foreach ($calls as $i => $call) {
            expect(thrownCode($call))->toBe($code, "case $i for $code");
        }
    }
    $message = function (callable $call): string {
        try {
            $call();
        } catch (\Throwable $e) {
            return $e->getMessage();
        }

        return '';
    };
    expect($message(fn () => pf('Phit\dist_categorical_rows_rng', 1, [[1], [1, -1]])))->toBe('prob row 1 probabilities must be nonnegative')
        ->and($message(fn () => pf('Phit\dist_categorical_rows_rng', 1, [[1], [1 => 1]])))->toBe('prob row 1 must be a list')
        ->and($message(fn () => pf('Phit\dist_categorical_rng', 1, 1, [0.0, 0])))->toBe('prob probabilities must have a positive sum');
});
