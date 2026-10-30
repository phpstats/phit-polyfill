<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Fma;

// Port of the extension's tests/draws.php against the NumPy reference cases
// in draws-python.json (1e-11 absolute plus relative).

function drawClose(mixed $actual, mixed $expected, string $name): void
{
    if (is_array($expected)) {
        expect(is_array($actual) && array_keys($actual) === array_keys($expected))->toBeTrue("$name shape/keys");
        foreach ($expected as $key => $v) {
            drawClose($actual[$key], $v, "$name.$key");
        }

        return;
    }
    expect(is_finite((float) $actual) && abs($actual - $expected) <= 1e-11 + 1e-11 * abs($expected))
        ->toBeTrue("$name: $actual != $expected");
}

function drawTranspose(array $a): array
{
    $out = [];
    for ($j = 0; $j < count($a[0]); $j++) {
        $out[] = array_column($a, $j);
    }

    return $out;
}

function drawMultiply(array $a, array $b): array
{
    $out = [];
    foreach ($a as $i => $row) {
        foreach (drawTranspose($b) as $j => $col) {
            $out[$i][$j] = array_sum(array_map(fn ($x, $y) => $x * $y, $row, $col));
        }
    }

    return $out;
}

$fixture = loadFixture('draws-python.json');

it('matches the reference draw cases', function () use ($fixture) {
    foreach ($fixture['cases'] as $case) {
        $name = $case['name'];
        switch ($case['operation']) {
            case 'moments':
                $result = pf('Phit\draws_moments', draws: $case['draws'], weights: $case['weights'], normalization: $case['normalization']);
                drawClose($result, $case['expected'], $name);
                break;
            case 'covariance_map':
                $result = pf('Phit\draws_covariance_map', sourceCovariance: $case['source_covariance'], targetCovariance: $case['target_covariance']);
                drawClose($result, $case['expected'], $name);
                $a = $result['matrix'];
                $cov = drawMultiply(drawMultiply($a, $case['source_covariance']), drawTranspose($a));
                drawClose($cov, $case['target_covariance'], "$name identity");
                foreach ($cov as $i => $row) {
                    foreach ($row as $j => $v) {
                        $ref = $case['target_covariance'][$i][$j];
                        if ($ref != 0) {
                            expect(abs($v / $ref - 1))->toBeLessThan(1e-12, "$name relative identity");
                        }
                    }
                }
                $reverse = pf('Phit\draws_covariance_map', $case['target_covariance'], $case['source_covariance']);
                drawClose($result['log_abs_det'] + $reverse['log_abs_det'], 0, "$name reverse logdet");
                break;
            case 'affine':
                $result = pf('Phit\draws_affine', draws: $case['draws'], matrix: $case['matrix'], shift: $case['shift']);
                drawClose($result, $case['expected'], $name);
                foreach ($result['draws'] as $i => $row) {
                    foreach ($row as $j => $v) {
                        $ref = $case['expected']['draws'][$i][$j];
                        if ($ref != 0) {
                            expect(abs($v / $ref - 1))->toBeLessThan(1e-12, "$name relative output");
                        }
                    }
                }
                $back = pf('Phit\draws_affine', $result['draws'], $case['inverse'], $case['inverse_shift']);
                drawClose($back['draws'], $case['draws'], "$name round trip");
                drawClose($result['log_abs_det'] + $back['log_abs_det'], 0, "$name round trip logdet");
                break;
            default:
                throw new RuntimeException('Unknown fixture operation');
        }
    }
});

it('returns the result keys of the API contract', function () {
    $x = [[1, 2], [3, 4]];
    $id = [[1, 0], [0, 1]];
    $schema = loadFixture('php-api.json');
    foreach ([
        'Phit\draws_moments' => pf('Phit\draws_moments', $x),
        'Phit\draws_covariance_map' => pf('Phit\draws_covariance_map', $id, $id),
        'Phit\draws_affine' => pf('Phit\draws_affine', $x, $id, [0, 0]),
    ] as $name => $result) {
        $keys = array_keys($result);
        sort($keys);
        expect($keys)->toBe($schema[$name]['result_keys']);
    }
    drawClose(pf('Phit\draws_moments', $x), pf('Phit\draws_moments', $x, [1, 1], 'population'), 'defaults');
    $subnormal = pf('Phit\draws_affine', [[1], [2]], [[1e-310]], [0]);
    expect($subnormal['draws'])->toBe([[1e-310], [2e-310]]);
    drawClose($subnormal['log_abs_det'], log(1e-310), 'subnormal log determinant');
    drawClose(pf('Phit\draws_affine', ['second' => [9 => 3, 2 => 4], 'first' => [1, 2]], $id, [0, 0])['draws'],
        [[3, 4], [1, 2]], 'draw order and literal columns');
});

it('applies the strict input contract', function () {
    $x = [[1, 2], [3, 4]];
    $id = [[1, 0], [0, 1]];
    $code = fn (callable $f) => thrownCode($f);
    foreach ([null, '1', true, []] as $bad) {
        expect($code(fn () => pf('Phit\draws_moments', [[1, 2], [3, $bad]], [1, 0])))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_moments', $x, [1, $bad])))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_covariance_map', [[1, $bad], [0, 1]], $id)))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_covariance_map', $id, [[1, 0], [$bad, 1]])))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_affine', [[1, $bad]], $id, [0, 0])))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_affine', $x, [[1, 0], [$bad, 1]], [0, 0])))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_affine', $x, $id, [0, $bad])))->toBe(1001);
    }
    foreach ([NAN, INF, -INF] as $bad) {
        expect($code(fn () => pf('Phit\draws_moments', [[1, 2], [3, $bad]], [1, 0])))->toBe(1003)
            ->and($code(fn () => pf('Phit\draws_moments', $x, [1, $bad])))->toBe(1003)
            ->and($code(fn () => pf('Phit\draws_covariance_map', [[1, 0], [0, $bad]], $id)))->toBe(1003)
            ->and($code(fn () => pf('Phit\draws_covariance_map', $id, [[1, 0], [0, $bad]])))->toBe(1003)
            ->and($code(fn () => pf('Phit\draws_affine', [[1, $bad]], $id, [0, 0])))->toBe(1003)
            ->and($code(fn () => pf('Phit\draws_affine', $x, [[1, 0], [0, $bad]], [0, 0])))->toBe(1003)
            ->and($code(fn () => pf('Phit\draws_affine', $x, $id, [0, $bad])))->toBe(1003);
    }
    foreach ([[], [[]]] as $bad) {
        expect($code(fn () => pf('Phit\draws_moments', $bad)))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_covariance_map', $bad, $id)))->toBe(1001)
            ->and($code(fn () => pf('Phit\draws_affine', $bad, $id, [0, 0])))->toBe(1001);
    }
    $cases = [
        [fn () => pf('Phit\draws_moments', [1, 2]), 1001],
        [fn () => pf('Phit\draws_moments', [[1], [2, 3]]), 1002],
        [fn () => pf('Phit\draws_moments', $x, [1]), 1002],
        [fn () => pf('Phit\draws_moments', $x, [0, 0]), 1001],
        [fn () => pf('Phit\draws_moments', $x, [-1, 2]), 1001],
        [fn () => pf('Phit\draws_moments', $x, [1, 0], 'unbiased'), 1001],
        [fn () => pf('Phit\draws_moments', $x, [1, 1e-30], 'unbiased'), 1001],
        [fn () => pf('Phit\draws_moments', [[1]], normalization: 'unbiased'), 1001],
        [fn () => pf('Phit\draws_moments', $x, normalization: 'sample'), 1001],
        [fn () => pf('Phit\draws_covariance_map', [[1]], $id), 1002],
        [fn () => pf('Phit\draws_covariance_map', [[1, 0]], $id), 1002],
        [fn () => pf('Phit\draws_covariance_map', [[1, 0], [0]], $id), 1002],
        [fn () => pf('Phit\draws_covariance_map', [[1, 0.2], [0.3, 1]], $id), 1001],
        [fn () => pf('Phit\draws_covariance_map', [[1, 0.2 + 1e-12], [0.2, 1]], $id), 1001],
        [fn () => pf('Phit\draws_affine', $x, [[1, 2], [2, 4]], [0, 0]), 1004],
        [fn () => pf('Phit\draws_affine', $x, [[1, 1], [1, 1 + 1e-15]], [0, 0]), 1004],
        [fn () => pf('Phit\draws_affine', $x, [[1, 2]], [0, 0]), 1002],
        [fn () => pf('Phit\draws_affine', $x, [[1, 0], [0]], [0, 0]), 1002],
        [fn () => pf('Phit\draws_affine', [[1], [2, 3]], $id, [0, 0]), 1002],
        [fn () => pf('Phit\draws_affine', $x, $id, [0]), 1002],
        [fn () => pf('Phit\draws_affine', [[PHP_FLOAT_MAX]], [[2]], [0]), 1003],
        [fn () => pf('Phit\draws_moments', [[-1e200], [1e200]]), 1003],
        [fn () => pf('Phit\draws_moments', [[1], [null]]), 1001],
    ];
    foreach ($cases as $i => [$f, $expected]) {
        expect($code($f))->toBe($expected, "case $i");
    }
    foreach ([[[1, 1], [1, 1]], [[1, 2], [2, 1]], [[-1, 0], [0, 1]], [[0, 0], [0, 1]], [[1, 1 - 1e-15], [1 - 1e-15, 1]]] as $singular) {
        expect($code(fn () => pf('Phit\draws_covariance_map', $singular, $id)))->toBe(1005)
            ->and($code(fn () => pf('Phit\draws_covariance_map', $id, $singular)))->toBe(1005);
    }
    expect(pf('Phit\draws_covariance_map', [[1, 0.2 + 1e-15], [0.2, 1]], $id))->toBeArray();
    drawClose(pf('Phit\draws_moments', [[1, 2], [3, 4]], [1, 0]), ['mean' => [1, 2], 'covariance' => [[0, 0], [0, 0]]], 'one effective draw');
});

it('emulates a correctly rounded fused multiply-add', function () {
    // a*b + c where the separately rounded product loses the answer entirely.
    $a = 1.0 + 2 ** -30;
    $b = 1.0 - 2 ** -30;
    expect(Fma::fma($a, $b, -1.0))->toBe(-(2 ** -60))
        ->and($a * $b - 1.0)->toBe(0.0)
        ->and(Fma::fma(0.1, 10.0, -1.0))->toBe(5.551115123125783e-17)
        ->and(Fma::fma(2.0, 3.0, 4.0))->toBe(10.0)
        ->and(Fma::fma(1e308, 10.0, -INF))->toBe(-INF);
});
