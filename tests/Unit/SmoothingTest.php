<?php

declare(strict_types=1);

// Port of the extension's tests/smoothing.php: stat_density, stat_bw_nrd0,
// stat_bw_nrd, stat_approx and stat_ecdf against R 4.6.1 at relative 1e-12,
// with density heights compared on the scale of their peak.

/** @param list<float> $actual */
function smoothClose(array $actual, array $expected, string $label, float $scale = 0.0): void
{
    expect(array_is_list($actual))->toBeTrue("$label shape")->and($actual)->toHaveCount(count($expected));
    foreach ($expected as $i => $e) {
        if ($e === null) {
            expect(is_nan($actual[$i]))->toBeTrue("$label/$i is not NAN");

            continue;
        }
        expect(abs($actual[$i] - $e))->toBeLessThanOrEqual(1e-12 * max(abs($e), $scale), "$label/$i");
    }
}

$fixture = loadFixture('smoothing-r.json');

it('matches R for density', function () use ($fixture) {
    foreach ($fixture['density'] as $c) {
        $d = pf('Phit\stat_density', $c['x'], ...$c['args']);
        $e = $c['expected'];
        expect(array_keys($d))->toBe(['x', 'y', 'bw', 'n'])
            ->and($d['n'])->toBe($e['n']);
        smoothClose([$d['bw']], [$e['bw']], "{$c['name']} bw");
        smoothClose($d['x'], (array) $e['x'], "{$c['name']} x");
        smoothClose($d['y'], (array) $e['y'], "{$c['name']} y", max((array) $e['y']));
    }
});

it('matches R for the bandwidths', function () use ($fixture) {
    foreach ($fixture['bandwidth'] as $c) {
        smoothClose([pf('Phit\stat_bw_nrd0', $c['x'])], [$c['nrd0']], "{$c['name']} nrd0");
        smoothClose([pf('Phit\stat_bw_nrd', $c['x'])], [$c['nrd']], "{$c['name']} nrd");
    }
});

it('matches R for approx and ecdf', function () use ($fixture) {
    foreach ($fixture['approx'] as $c) {
        smoothClose(pf('Phit\stat_approx', $c['x'], $c['y'], $c['xout'], ...$c['args']), $c['expected'], $c['name']);
    }
    foreach ($fixture['ecdf'] as $c) {
        expect(pf('Phit\stat_ecdf', $c['x'], $c['q']))->toBe(array_map('floatval', $c['expected']));
    }
});

it('reads null and NAN evaluation points as NAN and infinities as outside the data', function () {
    $out = pf('Phit\stat_approx', [1, 2], [10, 20], [null, NAN, -INF, INF], rule: 2);
    expect(is_nan($out[0]) && is_nan($out[1]))->toBeTrue()
        ->and(array_slice($out, 2))->toBe([10.0, 20.0]);
    $out = pf('Phit\stat_ecdf', [1, 2], [null, -INF, INF]);
    expect(is_nan($out[0]))->toBeTrue()
        ->and(array_slice($out, 1))->toBe([0.0, 1.0]);
});

it('rejects missing data by default and removes it with missing: drop', function () {
    expect(thrownCode(fn () => pf('Phit\stat_density', [1, null, 3])))->toBe(1001)
        ->and(pf('Phit\stat_density', [1.0, null, 3.0, NAN, 4.0], missing: 'drop'))->toBe(pf('Phit\stat_density', [1.0, 3.0, 4.0]))
        ->and(pf('Phit\stat_approx', [1, null, 3], [1, 5, 3], [2], missing: 'drop'))->toBe([2.0])
        ->and(pf('Phit\stat_ecdf', [1, NAN, 3], [2], missing: 'drop'))->toBe([0.5]);
});

it('renormalizes weights that summed to one after dropping, as R does', function () {
    $w = pf('Phit\stat_density', [1.0, null, 3.0], bw: 1.0, weights: [0.25, 0.25, 0.5], missing: 'drop');

    expect($w['y'])->toEqualWithDelta(pf('Phit\stat_density', [1.0, 3.0], bw: 1.0, weights: [1 / 3, 2 / 3])['y'], 1e-15);
});

it('rejects invalid input with the extension codes', function () {
    foreach ([
        'density infinite' => [fn () => pf('Phit\stat_density', [1, INF]), 1003],
        'automatic bandwidth needs two values' => [fn () => pf('Phit\stat_density', [1]), 1001],
        'unsupported bandwidth rule' => [fn () => pf('Phit\stat_density', [1, 2], bw: 'ucv'), 1001],
        'negative bandwidth' => [fn () => pf('Phit\stat_density', [1, 2], bw: -1.0), 1001],
        'non-finite adjust' => [fn () => pf('Phit\stat_density', [1, 2], adjust: NAN), 1003],
        'unknown kernel' => [fn () => pf('Phit\stat_density', [1, 2], kernel: 'box'), 1001],
        'n zero' => [fn () => pf('Phit\stat_density', [1, 2], n: 0), 1001],
        'from above to' => [fn () => pf('Phit\stat_density', [1, 2], from: 3.0, to: 1.0), 1001],
        'weights length' => [fn () => pf('Phit\stat_density', [1, 2], weights: [1.0]), 1002],
        'negative weights' => [fn () => pf('Phit\stat_density', [1, 2], weights: [-1.0, 2.0]), 1001],
        'rule' => [fn () => pf('Phit\stat_approx', [1, 2], [1, 2], [1], rule: 3), 1001],
        'method' => [fn () => pf('Phit\stat_approx', [1, 2], [1, 2], [1], method: 'spline'), 1001],
        'ties' => [fn () => pf('Phit\stat_approx', [1, 2], [1, 2], [1], ties: 'ordered'), 1001],
        'f' => [fn () => pf('Phit\stat_approx', [1, 2], [1, 2], [1], method: 'constant', f: 2.0), 1001],
        'linear needs two distinct x' => [fn () => pf('Phit\stat_approx', [1, 1], [1, 2], [1]), 1001],
        'approx lengths' => [fn () => pf('Phit\stat_approx', [1, 2], [1], [1]), 1002],
        'xout type' => [fn () => pf('Phit\stat_approx', [1, 2], [1, 2], ['a']), 1001],
        'ecdf empty' => [fn () => pf('Phit\stat_ecdf', [], [1]), 1001],
    ] as $label => [$call, $code]) {
        expect(thrownCode($call))->toBe($code, $label);
    }
});

it('falls back as bw.nrd0 does when the spread is zero', function () {
    expect(pf('Phit\stat_bw_nrd0', [3, 3]))->toBe(0.9 * 3.0 * 2.0 ** -0.2)
        ->and(pf('Phit\stat_bw_nrd0', [0, 0]))->toBe(0.9 * 2.0 ** -0.2);
});

it('steps a constant approx with weight f on the right-hand value', function () {
    expect(pf('Phit\stat_approx', [1, 2], [5, 7], [0, 1.5, 9], method: 'constant', rule: 2))->toBe([5.0, 5.0, 7.0])
        ->and(pf('Phit\stat_approx', [1, 2], [5, 7], [1.5], method: 'constant', f: 0.25))->toBe([5.5]);
});
