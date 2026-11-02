<?php

declare(strict_types=1);

use Phit\Polyfill\Descriptive;
use Phit\Polyfill\Optimize;
use Phit\Polyfill\Smoothing;

// The tables of phit's specs/cleanroom/data (copied to tests/fixtures/cleanroom)
// for Brent's methods (brent.md) and for interpolation, the normal density and
// the type-7 quantile (smoothing.md), checked with the spec's tolerances. Each
// check reports its largest difference on STDERR.

if (! function_exists('numericCell')) {
    /** A %.17g cell of the R tables; "NA" and "NaN" are NaN. */
    function numericCell(string|int|float $cell): float
    {
        if (! is_string($cell)) {
            return (float) $cell;
        }

        return match ($cell) {
            'Inf' => INF,
            '-Inf' => -INF,
            'NaN', 'NA' => NAN,
            default => (float) $cell,
        };
    }

    /** The bits of a double, so that comparisons are exact and see the sign of zero. */
    function numericBits(float $x): string
    {
        return bin2hex(pack('E', $x));
    }

    /** |got - expected| / |expected| (absolute at zero); 0 when both are the same non-finite value. */
    function numericDifference(float $got, float $expected): float
    {
        if (is_nan($expected) || is_infinite($expected)) {
            return ($got === $expected || (is_nan($got) && is_nan($expected))) ? 0.0 : INF;
        }
        $diff = abs($got - $expected);

        return $expected == 0.0 ? $diff : $diff / abs($expected);
    }

    /** A double for a failure message (interpolating NaN into a string warns). */
    function numericText(float $x): string
    {
        return var_export($x, true);
    }

    function numericReport(string $table, int $rows, float $worst): void
    {
        fwrite(STDERR, sprintf("%s: %d rows, largest relative difference %.3e\n", $table, $rows, $worst));
    }

    /** A private static method of a polyfill class, as a closure. */
    function numericPrivate(string $class, string $method): Closure
    {
        return (new ReflectionMethod($class, $method))->getClosure();
    }

    /** The functions of brent.md. */
    function brentFunction(string $name): Closure
    {
        return match ($name) {
            'parabola' => fn (float $x): float => ($x - 1 / 3) * ($x - 1 / 3),
            'sine' => fn (float $x): float => sin($x),
            'kink' => fn (float $x): float => abs($x - 0.7),
            'flat' => fn (float $x): float => 1.0,
            'linear' => fn (float $x): float => $x,
            'decay' => fn (float $x): float => exp(-$x),
            'double_well' => fn (float $x): float => ($x * $x - 1) * ($x * $x - 1) + 0.3 * $x,
            'wall' => fn (float $x): float => $x < 0.45 ? INF : ($x - 0.5) * ($x - 0.5),
            'far' => fn (float $x): float => ($x - 1e6) * ($x - 1e6),
            'hump' => fn (float $x): float => $x * exp(-$x),
            'sqrt2' => fn (float $x): float => $x * $x - 2,
            'fixed_point' => fn (float $x): float => cos($x) - $x,
            'exp5' => fn (float $x): float => exp($x) - 5,
            'at_one' => fn (float $x): float => $x - 1,
            'at_three' => fn (float $x): float => $x - 3,
            'step' => fn (float $x): float => $x < 0.3 ? -1.0 : 1.0,
            'cubic' => fn (float $x): float => $x * $x * $x - 2 * $x - 5,
            'log' => fn (float $x): float => $x == 0.0 ? -INF : log($x),
            'cliff' => fn (float $x): float => $x < 0.5 ? -INF : $x - 0.6,
        };
    }

    /** Runs $search with the replacement warnings of Optimize silenced, and returns its result. */
    function brentQuietly(Closure $search): mixed
    {
        set_error_handler(fn (int $no, string $message): bool => str_contains($message, 'replaced')
            || str_contains($message, 'not converged'), E_USER_WARNING);
        try {
            return $search();
        } finally {
            restore_error_handler();
        }
    }
}

it('matches R optimize() bit for bit', function () {
    $cases = loadFixture('cleanroom/optimize.json')['optimize'];
    expect($cases)->toHaveCount(12);
    $worst = 0.0;
    foreach ($cases as $case) {
        $f = brentFunction($case['f']);
        $calls = 0;
        $counted = function (float $x) use ($f, &$calls): float {
            $calls++;

            return $f($x);
        };
        $maximum = $case['maximum'] ?? false;
        $result = brentQuietly(fn () => Optimize::optimize($counted, (float) $case['lower'], (float) $case['upper'], isset($case['tol']) ? (float) $case['tol'] : null, $maximum));
        $x = $result[$maximum ? 'maximum' : 'minimum'];
        $label = json_encode($case);
        expect(numericBits($x))->toBe(numericBits((float) $case['expected']['x']), $label.' got '.numericText($x));
        // R's count includes the final call for the objective, as does the extension's.
        expect($calls)->toBe($case['expected']['evaluations'], $label);
        expect($result['objective'])->toBe((float) $case['expected']['objective'], $label);
        $worst = max($worst, numericDifference($x, (float) $case['expected']['x']));
    }
    numericReport('optimize', count($cases), $worst);
});

it('matches R uniroot() bit for bit', function () {
    $cases = loadFixture('cleanroom/optimize.json')['uniroot'];
    expect($cases)->toHaveCount(11);
    $worst = 0.0;
    foreach ($cases as $case) {
        $f = brentFunction($case['f']);
        $calls = 0;
        $counted = function (float $x) use ($f, &$calls): float {
            $calls++;

            return $f($x);
        };
        $expected = $case['expected'];
        $result = brentQuietly(fn () => Optimize::uniroot(
            $counted,
            (float) $case['lower'],
            (float) $case['upper'],
            isset($case['tol']) ? (float) $case['tol'] : null,
            $case['maxiter'] ?? 1000,
        ));
        $label = json_encode($case);
        expect(numericBits($result['root']))->toBe(numericBits((float) $expected['root']), $label.' got '.numericText($result['root']));
        expect(numericBits($result['estim_prec']))->toBe(numericBits((float) $expected['estim_prec']), $label.' got '.numericText($result['estim_prec']));
        // Without convergence the binding reports maxiter, as R does.
        expect($result['iter'])->toBe($expected['iter'], $label);
        // Both endpoints, the search, and the root, as R counts them.
        expect($calls)->toBe($expected['evaluations'], $label);
        expect($result['f_root'])->toBe((float) $expected['f_root'], $label);
        $worst = max($worst, numericDifference($result['root'], (float) $expected['root']), numericDifference($result['estim_prec'], (float) $expected['estim_prec']));
    }
    numericReport('uniroot', count($cases), $worst);
});

it('matches R approx() to 1e-15, and exactly at the data points', function () {
    $interpolate = numericPrivate(Smoothing::class, 'interpolate');
    $cases = loadFixture('cleanroom/approx.json');
    expect($cases)->toHaveCount(80);
    $worst = 0.0;
    $points = 0;
    foreach ($cases as $case) {
        $x = array_map(numericCell(...), $case['x']);
        $y = array_map(numericCell(...), $case['y']);
        $constant = $case['method'] === 'constant';
        [$low, $high] = $case['rule'] === 2 ? [$y[0], $y[count($y) - 1]] : [NAN, NAN];
        foreach ($case['xout'] as $i => $cell) {
            $v = numericCell($cell);
            $expected = numericCell($case['out'][$i]);
            $got = $interpolate($x, $y, $v, $constant, (float) $case['f'], $low, $high);
            $label = json_encode($case)." at $cell got ".numericText($got);
            $diff = numericDifference($got, $expected);
            if (in_array($v, $x, true)) {
                expect(numericBits($got))->toBe(numericBits($expected), $label);
            }
            expect($diff)->toBeLessThanOrEqual(1e-15, $label);
            $worst = max($worst, $diff);
            $points++;
        }
    }
    numericReport('approx', $points, $worst);
});

it('matches R dnorm() to 1e-14', function () {
    $table = loadFixture('cleanroom/dnorm.json');
    $worst = 0.0;
    $subnormalUlps = 0.0;
    foreach ($table['z'] as $i => $cell) {
        $got = Smoothing::standardNormalDensity(numericCell($cell));
        $expected = numericCell($table['out'][$i]);
        if ($expected > 0.0 && $expected < PHP_FLOAT_MIN) {
            // Subnormal: the spec's relative bound covers normal results only, and the
            // density may underflow as the true value does; it must not exceed R's.
            $ulps = abs($got - $expected) / 2.0 ** -1074;
            expect($ulps)->toBeLessThanOrEqual(1.0, "dnorm($cell) got ".numericText($got))
                ->and($got)->toBeLessThanOrEqual($expected, "dnorm($cell) got ".numericText($got));
            $subnormalUlps = max($subnormalUlps, $ulps);

            continue;
        }
        $diff = numericDifference($got, $expected);
        expect($diff)->toBeLessThanOrEqual(1e-14, "dnorm($cell) got ".numericText($got));
        $worst = max($worst, $diff);
    }
    numericReport('dnorm', count($table['z']), $worst);
    fwrite(STDERR, sprintf("dnorm: largest difference on subnormal results %g of the smallest subnormal\n", $subnormalUlps));

    $scaled = $table['scaled'];
    $worst = 0.0;
    foreach ($scaled['x'] as $i => $cell) {
        $sigma = numericCell($scaled['sigma'][$i]);
        $got = Smoothing::standardNormalDensity(numericCell($cell) / $sigma) / $sigma;
        $diff = numericDifference($got, numericCell($scaled['out'][$i]));
        expect($diff)->toBeLessThanOrEqual(1e-14, "dnorm($cell, 0, ".numericText($sigma).") got ".numericText($got));
        $worst = max($worst, $diff);
    }
    numericReport('dnorm scaled', count($scaled['x']), $worst);
});

it('keeps the normal density symmetric, and zero beyond the tails', function () {
    foreach ([0.5, 1.0, 5.25, 12.0, 37.9, 38.5] as $z) {
        expect(numericBits(Smoothing::standardNormalDensity(-$z)))->toBe(numericBits(Smoothing::standardNormalDensity($z)));
    }
    expect(Smoothing::standardNormalDensity(INF))->toBe(0.0)
        ->and(Smoothing::standardNormalDensity(-INF))->toBe(0.0)
        ->and(Smoothing::standardNormalDensity(40.5))->toBe(0.0)
        ->and(is_nan(Smoothing::standardNormalDensity(NAN)))->toBeTrue();
});

it('matches R quantile(type = 7) to 1e-15', function () {
    $quantile = numericPrivate(Descriptive::class, 'quantileType7');
    $cases = loadFixture('cleanroom/quantile7.json');
    expect($cases)->toHaveCount(150);
    $worst = 0.0;
    $points = 0;
    foreach ($cases as $case) {
        $sorted = array_map(numericCell(...), $case['x']);
        usort($sorted, Descriptive::totalOrder(...));
        $public = Descriptive::quantileArray(array_map(numericCell(...), $case['x']), array_map(numericCell(...), $case['p']));
        foreach ($case['p'] as $i => $cell) {
            $got = $quantile($sorted, numericCell($cell));
            $expected = numericCell($case['out'][$i]);
            $label = json_encode($case)." at p = $cell got ".numericText($got);
            $diff = numericDifference($got, $expected);
            expect($diff)->toBeLessThanOrEqual(1e-15, $label);
            expect(numericBits($public[$i]))->toBe(numericBits($got), $label);
            $worst = max($worst, $diff);
            $points++;
        }
    }
    numericReport('quantile7', $points, $worst);
});
