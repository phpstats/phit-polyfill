<?php

declare(strict_types=1);

// Compares the smoothing functions with the loaded extension; skipped without
// it. Results and exceptions must be identical to the bit. The extension's
// FFT calls sin and cos separately since phit 3dc8bb2; before that, Apple's
// combined sincos gave different twiddles on macOS.

/** @return array{0: string, 1: mixed} */
function smoothingOutcome(callable $call): array
{
    try {
        return ['ok', $call()];
    } catch (\Throwable $e) {
        return ['error', [get_class($e), $e->getCode(), $e->getMessage()]];
    }
}

function smoothingCompare(string $fn, array $args): void
{
    $native = smoothingOutcome(fn () => native($fn, ...$args));
    $polyfill = smoothingOutcome(fn () => pf($fn, ...$args));
    parityLog('smoothing', $fn, $polyfill, $native, $native[0] === 'error' || $polyfill[0] === 'error');
    $label = fn () => $fn.' '.substr((string) json_encode($args, JSON_PARTIAL_OUTPUT_ON_ERROR), 0, 300);
    if (serialize($native) !== serialize($polyfill)) {
        throw new \RuntimeException('mismatch: '.$label()."\nnative ".substr(var_export($native, true), 0, 400)."\npolyfill ".substr(var_export($polyfill, true), 0, 400));
    }
    expect(true)->toBeTrue();
}

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

it('agrees on the R fixtures', function () {
    $fixture = loadFixture('smoothing-r.json');
    foreach ($fixture['density'] as $c) {
        smoothingCompare('Phit\stat_density', [$c['x'], ...$c['args']]);
    }
    foreach ($fixture['bandwidth'] as $c) {
        smoothingCompare('Phit\stat_bw_nrd0', [$c['x']]);
        smoothingCompare('Phit\stat_bw_nrd', [$c['x']]);
    }
    foreach ($fixture['approx'] as $c) {
        smoothingCompare('Phit\stat_approx', [$c['x'], $c['y'], $c['xout'], ...$c['args']]);
    }
    foreach ($fixture['ecdf'] as $c) {
        smoothingCompare('Phit\stat_ecdf', [$c['x'], $c['q']]);
    }
});

it('agrees on seeded random inputs', function () {
    mt_srand(crc32('smoothing') ^ (int) (getenv('PARITY_SEED') ?: 20261030));
    $u = fn (float $lo, float $hi): float => $lo + ($hi - $lo) * mt_rand() / mt_getrandmax();
    $kernels = ['gaussian', 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine', 'optcosine'];
    for ($it = 0, $iterations = (int) (getenv('PARITY_ITERATIONS') ?: 60) * 4; $it < $iterations; $it++) {
        $m = mt_rand(1, 60);
        $x = [];
        for ($i = 0; $i < $m; $i++) {
            $x[] = mt_rand(0, 9) === 0 ? (float) mt_rand(-3, 3) : $u(-50, 50) * (mt_rand(0, 1) === 1 ? 1.0 : 1e-6);
        }
        if (mt_rand(0, 5) === 0) {
            $x[mt_rand(0, $m - 1)] = null;
        }
        $w = mt_rand(0, 2) === 0 ? array_map(fn () => $u(0, 2), $x) : null;
        if ($w !== null && mt_rand(0, 1) === 1) {
            $total = array_sum($w);
            $w = array_map(fn ($v) => $v / $total, $w);
        }
        smoothingCompare('Phit\stat_density', [
            $x, ['nrd0', 'nrd', 'NRD0', $u(0.01, 5), -1.0, 0, 'x'][mt_rand(0, 6)], $u(0.2, 3), $kernels[mt_rand(0, 6)],
            [1, 2, 3, 50, 512, 513, 700, 1024, 2000][mt_rand(0, 8)],
            mt_rand(0, 3) === 0 ? $u(-80, 0) : null, mt_rand(0, 3) === 0 ? $u(0, 80) : null, $u(0, 4), $w,
            mt_rand(0, 2) === 0 ? 'raise' : 'drop',
        ]);
        smoothingCompare('Phit\stat_bw_nrd0', [$x, 'drop']);
        smoothingCompare('Phit\stat_bw_nrd', [$x, 'drop']);
        // Rounded x for ties, signed zeros in x and y, and every kind of point.
        $y = array_map(fn () => mt_rand(0, 4) > 0 ? $u(-10, 10) : mt_rand(-1, 1) * 0.0, $x);
        $xr = array_map(fn ($v) => $v === null ? null : (mt_rand(0, 2) > 0 ? round($v) : $v) * (mt_rand(0, 9) > 0 ? 1.0 : -0.0), $x);
        $points = [];
        for ($i = 0; $i < 20; $i++) {
            $points[] = [null, NAN, INF, -INF, $u(-60, 60), (float) mt_rand(-50, 50), -0.0][mt_rand(0, 6)];
        }
        smoothingCompare('Phit\stat_approx', [$xr, $y, $points, ['linear', 'constant'][mt_rand(0, 1)], mt_rand(1, 2), ['mean', 'min', 'max'][mt_rand(0, 2)], [0.0, 1.0, $u(0, 1)][mt_rand(0, 2)], 'drop']);
        smoothingCompare('Phit\stat_ecdf', [$x, $points, 'drop']);
    }
});

it('rejects the same invalid input with the same exception', function () {
    foreach ([
        ['Phit\stat_density', [[1, INF]]], ['Phit\stat_density', [[1]]], ['Phit\stat_density', [[1, 2], 'ucv']],
        ['Phit\stat_density', [[1, 2], -1.0]], ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'box']],
        ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'gaussian', 0]], ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'gaussian', 512, 3.0, 1.0]],
        ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'gaussian', 512, null, null, 3.0, [1.0]]],
        ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'gaussian', 512, null, null, 3.0, [-1.0, 2.0]]],
        ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'gaussian', 512, null, null, 3.0, [1.0, INF]]],
        ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'gaussian', 512, null, null, 3.0, [1.0, null]]],
        ['Phit\stat_density', [[1, 2], 'nrd0', NAN]], ['Phit\stat_density', [[1, 2], 'nrd0', 1.0, 'gaussian', 512, NAN]],
        ['Phit\stat_density', [[], 1.0]], ['Phit\stat_density', [[1, 'a']]], ['Phit\stat_density', [[1, null], 'ucv']],
        ['Phit\stat_approx', [[1, 2], [1, 2], [1], 'linear', 3]], ['Phit\stat_approx', [[1, 2], [1, 2], [1], 'spline']],
        ['Phit\stat_approx', [[1, 2], [1, 2], [1], 'linear', 1, 'ordered']], ['Phit\stat_approx', [[1, 1], [1, 2], [1]]],
        ['Phit\stat_approx', [[1, 2], [1], [1]]], ['Phit\stat_approx', [[1, 2], [1, 2], ['a']]],
        ['Phit\stat_approx', [[1, 2], [1, 2], [1], 'constant', 1, 'mean', 2.0]], ['Phit\stat_approx', [[1, 2], [1, 2], [1], 'constant', 1, 'mean', NAN]],
        ['Phit\stat_approx', [[], [], [1], 'constant']], ['Phit\stat_approx', [[1], [5], [0, 1, 2], 'constant', 2]],
        ['Phit\stat_approx', [[1, 2], [1, 2], [1], 'linear', 1, 'mean', 0.0, 'bogus']],
        ['Phit\stat_approx', [[1e308, -1e308, 1e308], [1e308, 1, 1e308], [0]]],
        ['Phit\stat_ecdf', [[], [1]]], ['Phit\stat_ecdf', [[1], ['x']]], ['Phit\stat_ecdf', [[1, INF], [1]]],
        ['Phit\stat_bw_nrd0', [[0, 0]]], ['Phit\stat_bw_nrd0', [[3, 3]]], ['Phit\stat_bw_nrd', [[3, 3]]], ['Phit\stat_bw_nrd0', [[1]]],
    ] as [$fn, $args]) {
        smoothingCompare($fn, $args);
    }
});
