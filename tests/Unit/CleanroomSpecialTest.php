<?php

declare(strict_types=1);

use Phit\Polyfill\Distribution\Saddle;
use Phit\Polyfill\Internal\Statrs;

// The tables of phit's specs/cleanroom/data/special.json (copied to
// tests/fixtures/cleanroom), checked against the saddle-point and incomplete
// gamma and beta routines with the tolerances of specs/cleanroom/special.md.
// Each row lists the inputs and then the expected value.

if (! function_exists('cleanroomCell')) {
    function cleanroomCell(mixed $cell): float
    {
        if (is_string($cell)) {
            return match (strtolower($cell)) {
                'inf' => INF,
                '-inf' => -INF,
                'nan' => NAN,
                default => (float) $cell,
            };
        }

        return (float) $cell;
    }

    /**
     * Checks $routine on every row of a table and returns the largest relative
     * difference. A scaled check measures it against max(1, |value|).
     */
    function cleanroomCheck(string $table, float $tolerance, bool $scaled, Closure $routine): float
    {
        $rows = loadFixture('cleanroom/special.json')[$table];
        expect($rows)->not->toBeEmpty();
        $worst = 0.0;
        foreach ($rows as $row) {
            $cells = array_map(cleanroomCell(...), $row);
            $expected = array_pop($cells);
            $got = $routine(...$cells);
            $label = $table.' '.json_encode($row);
            if (is_infinite($expected) || is_nan($expected)) {
                expect($got === $expected || (is_nan($got) && is_nan($expected)))->toBeTrue($label);

                continue;
            }
            expect(is_finite($got))->toBeTrue($label.' got '.$got);
            $scale = $scaled ? max(abs($expected), 1.0) : abs($expected);
            $diff = $scale == 0.0 ? abs($got - $expected) : abs($got - $expected) / $scale;
            expect($diff)->toBeLessThanOrEqual($tolerance, $label.' got '.$got);
            $worst = max($worst, $diff);
        }
        fwrite(STDERR, sprintf("%s: %d rows, largest relative difference %.3e\n", $table, count($rows), $worst));

        return $worst;
    }
}

it('matches the stirling_error table', function () {
    cleanroomCheck('stirling_error', 1e-14, false, Saddle::stirlingError(...));
});

it('matches the deviance_term table', function () {
    cleanroomCheck('deviance_term', 1e-13, false, Saddle::devianceTerm(...));
});

it('matches the ln_poisson_raw table', function () {
    cleanroomCheck('ln_poisson_raw', 1e-13, true, Saddle::lnPoissonRaw(...));
});

it('matches the ln_binomial_raw table', function () {
    cleanroomCheck('ln_binomial_raw', 1e-13, true, Saddle::lnBinomialRaw(...));
});

it('matches the ln_gamma_lower table', function () {
    cleanroomCheck('ln_gamma_lower', 1e-13, true, Saddle::lnGammaLowerSeries(...));
});

it('matches the ln_gamma_upper table', function () {
    cleanroomCheck('ln_gamma_upper', 1e-13, true, Saddle::lnGammaUpperFraction(...));
});

it('matches the ln_beta_lower table', function () {
    cleanroomCheck(
        'ln_beta_lower',
        1e-13,
        true,
        fn (float $a, float $b, float $x): float => Saddle::lnBetaLowerFraction($a, $b, $x, log($x), log1p(-$x)),
    );
});

it('matches the ln_beta_lower_tiny table, where only ln x carries an underflowed x', function () {
    cleanroomCheck('ln_beta_lower_tiny', 1e-13, true, Saddle::lnBetaLowerFraction(...));
});

it('follows ln x when x is subnormal', function () {
    $lnX = -720.0;
    $x = exp($lnX);
    expect($x > 0.0 && $x < PHP_FLOAT_MIN)->toBeTrue();
    [$a, $b] = [2.5, 3.0];
    $lnABeta = Statrs::lnGamma($a + 1.0) + Statrs::lnGamma($b) - Statrs::lnGamma($a + $b);
    $want = $a * $lnX - $lnABeta;
    $got = Saddle::lnBetaLowerFraction($a, $b, $x, $lnX, log1p(-$x));
    expect(abs($got - $want))->toBeLessThanOrEqual(1e-13 * abs($want));
});

it('keeps the deviance term finite and right near the top of the double range', function () {
    // x + np and 2x overflow here; the series once took 0 · ∞ = NaN for its first term and never settled.
    // Expected values are x ln(x / np) + np - x at 60 digits (mpmath) for the exact doubles.
    $cases = [
        [1e308, 9.9e307, 5.033585350144054e303],
        [PHP_FLOAT_MAX, 1.7e308, 2.7547670212519258e305],
        [4.4e307, 6.6e305, 1.4144702342671678e308], // x ln(x / np) alone overflows, but D does not
        [5e-324, 1e308, 1e308], // a quarter of x underflows to zero
        [1e-300, 1e308, 1e308],
    ];
    foreach ($cases as [$x, $np, $want]) {
        $got = Saddle::devianceTerm($x, $np);
        expect(abs($got - $want) / $want)->toBeLessThanOrEqual(1e-13, "D($x, $np) = $got");
    }
    expect(Saddle::devianceTerm(PHP_FLOAT_MAX, PHP_FLOAT_MAX))->toBe(0.0)
        ->and(Saddle::devianceTerm(1e308, 1e-300))->toBe(INF);
});

it('scales the deviance term exactly by powers of two', function () {
    foreach ([[1.9, 1.8], [1.0, 1.05], [1.9, 0.2], [0.3, 1.9]] as [$x, $np]) {
        $small = Saddle::devianceTerm($x, $np);
        foreach ([1000, 1021, 1022, 1023] as $k) {
            $c = 2.0 ** $k;
            expect(Saddle::devianceTerm($x * $c, $np * $c))->toBe($small * $c, "D($x · 2^$k, $np · 2^$k)");
        }
    }
});

it('gives the beta log density with a shape parameter near the largest double', function () {
    // ln f(x) = (a - 1) ln x + (b - 1) ln(1 - x) - ln B(a, b), at 400 digits (mpmath).
    $logpdf = \Phit\Polyfill\Polyfill::implementation('Phit\dist_beta_logpdf');
    foreach ([[0.01, -1.0050335853501442e306], [1.1868334525134252e-307, 707.0426934330717]] as [$x, $want]) {
        $got = $logpdf($x, 12.2, 1e308);
        expect(abs(($got - $want) / $want))->toBeLessThanOrEqual(1e-13, "logpdf($x) = $got");
    }
});
