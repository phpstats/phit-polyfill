<?php

declare(strict_types=1);

// Ported from phit's tests/summary.php and the tests in stats-core's summary.rs:
// Phit\summary() and Phit\tidy() against tests/fixtures/summary/summary.json,
// the summary text byte for byte, the tidy columns exactly, and each printed
// coefficient row against R's printCoefmat() for the same model.

use Phit\Polyfill\Summary\Format;

// The same printed form, allowing one unit in the last printed place: phit and R stop their iterations at slightly
// different points, so a sixth significant digit can differ. The width and the number of decimals must match.
function summarySameToken(string $a, string $b): bool
{
    if ($a === $b) {
        return true;
    }
    $fixed = '/^-?\d+\.(\d+)$/';
    if (strlen($a) !== strlen($b) || ! preg_match($fixed, $a, $ma) || ! preg_match($fixed, $b, $mb)) {
        return false;
    }

    return strlen($ma[1]) === strlen($mb[1]) && abs((float) $a - (float) $b) <= 1.000001 * 10 ** -strlen($ma[1]);
}

// The same tidy columns: terms exactly, numbers to 1e-5 relative. The fixtures come from one platform, and the
// iterative fits can differ elsewhere: by about 1e-14 for the GLMs and up to 7e-7 for the logistic mixed fits,
// where BOBYQA takes another path on Linux. The printed summary must still match byte for byte.
function summarySameTidy(array $ours, array $expected): bool
{
    if (array_keys($ours) !== array_keys($expected)) {
        return false;
    }
    foreach ($expected as $column => $values) {
        if (in_array($column, ['term', 'group', 'level'], true)) {
            if ($ours[$column] !== $values) {
                return false;
            }

            continue;
        }
        foreach ($values as $i => $v) {
            $a = $ours[$column][$i];
            if (is_nan($v) ? ! is_nan($a) : abs($a - $v) > 1e-5 * max(abs($v), 1e-300)) {
                return false;
            }
        }
    }

    return true;
}

/** The fit of a fixture case, from the polyfill. */
function summaryCaseFit(array $case): array
{
    $fixture = loadFixture('summary/summary.json');
    $args = $case['args'];
    if (isset($args['data'])) {
        $args['data'] = $fixture['data'][$args['data']];
    }

    return pf('Phit\\'.$case['function'], ...$args);
}

describe('Phit\summary and Phit\tidy on the fixture models', function () {
    beforeEach(fn () => skipMixedFitsUnderJit());

    it('prints each model as the extension does', function (array $case) {
        $fit = summaryCaseFit($case);
        $text = pf('Phit\summary', $fit);
        expect($text)->toBe($case['summary']);
        expect(summarySameTidy(pf('Phit\tidy', $fit, confInt: true), $case['tidy']))->toBeTrue("{$case['name']}: tidy differs");

        // Random effects: the same rows as the fixture, and lme4's ranef(condVar = TRUE) to 1e-4 (scaled by
        // max(1, |x|)), which absorbs lme4's optimizer tolerance. Matrix fits name their groups and columns generically.
        if (isset($case['ran_vals'])) {
            $values = pf('Phit\tidy', $fit, effects: 'ran_vals');
            expect(summarySameTidy($values, $case['ran_vals']))->toBeTrue("{$case['name']}: ran_vals differ");
            $r = $case['r_ran_vals'];
            expect($values['level'])->toBe($r['level']);
            if (! str_contains($case['function'], 'matrix')) {
                expect($values['group'])->toBe($r['group'])
                    ->and($values['term'])->toBe($r['term']);
            }
            foreach (['estimate', 'std_error'] as $column) {
                foreach ($r[$column] as $i => $v) {
                    expect(abs($values[$column][$i] - $v))->toBeLessThanOrEqual(1e-4 * max(1.0, abs($v)), "{$case['name']}: ran_vals $column $i");
                }
            }
        }

        // Every number R prints in a coefficient row prints the same here.
        $lines = explode("\n", $text);
        $start = array_search(true, array_map(fn ($l) => in_array($l, ['Coefficients:', 'Fixed effects:'], true), $lines), true);
        foreach ($case['r_rows'] as $i => $rRow) {
            $ours = preg_split('/\s{1,}(?=[-<\d])/', trim($lines[$start + 2 + $i]));
            $theirs = preg_split('/\s{1,}(?=[-<\d])/', trim($rRow));
            $numbers = min(count($ours), count($theirs)) - 1;
            foreach (array_map(null, array_slice($ours, 1, $numbers), array_slice($theirs, 1, $numbers)) as [$a, $b]) {
                expect(summarySameToken($a, $b))->toBeTrue("{$case['name']}: row $i '{$lines[$start + 2 + $i]}' vs R '$rRow'");
            }
        }
    })->with(fn () => array_map(fn (array $case) => [$case], array_column(loadFixture('summary/summary.json')['cases'], null, 'name')));

    it('has a case for every model function', function () {
        $covered = array_unique(array_column(loadFixture('summary/summary.json')['cases'], 'function'));
        foreach (Phit\Polyfill\Polyfill::functions() as $function) {
            if (str_starts_with($function, 'Phit\\fit_')) {
                expect($covered)->toContain(substr($function, 5));
            }
        }
    });

    it('gives tidy columns without intervals by default', function () {
        $fit = summaryCaseFit(loadFixture('summary/summary.json')['cases'][0]);
        expect(array_keys(pf('Phit\tidy', $fit)))->toBe(['term', 'estimate', 'std_error', 'statistic', 'p_value']);
    });
});

describe('Phit\summary and Phit\tidy errors', function () {
    it('rejects arrays that are not fits', function (array $bad, string $message) {
        expect(fn () => pf('Phit\summary', $bad))->toThrow(Exception::class, "fit must be a result of a Phit model function: $message");
    })->with([
        'empty' => [[], 'missing model'],
        'unknown model' => [['model' => 'quantile', 'coefficients' => [1.0]], "unknown model 'quantile'"],
        'coefficients only' => [['model' => 'linear', 'coefficients' => [1.0]], 'missing z_statistics'],
        'strings' => [['model' => 'linear', 'coefficients' => ['a'], 't_statistics' => [1.0]], 'coefficients must be numbers'],
        'lengths' => [['model' => 'linear', 'coefficients' => [1.0], 't_statistics' => [1.0], 'std_errors' => [1.0, 2.0], 'p_values' => [0.5]],
            'coefficient lists differ in length'],
        'random term' => [['model' => 'mixed_gaussian', 'log_likelihood' => -1.0, 'df_residual' => 3, 'n_observations' => 4, 'random' => [1]],
            'random terms must be arrays'],
        'logistic integration' => [['model' => 'mixed_logistic'], 'missing integration'],
    ]);

    it('rejects random terms whose covariance does not fit their columns', function () {
        skipMixedFitsUnderJit();
        $fit = pf('Phit\fit_mixed_gaussian', 'Reaction ~ Days + (Days | Subject)', loadFixture('summary/summary.json')['data']['sleepstudy']);
        $short = $fit;
        $short['random'][0]['covariance'] = [[1.0]];
        $ragged = $fit;
        $ragged['random'][0]['covariance'][1] = [1.0];

        expect(fn () => pf('Phit\summary', $short))->toThrow(Exception::class, "a random term's columns and covariance differ in size")
            ->and(fn () => pf('Phit\summary', $ragged))->toThrow(Exception::class, "a random term's covariance must be square");
    });

    it('rejects other effects, and ran_vals without random effects', function () {
        $ols = pf('Phit\fit_ols', 'mpg ~ wt', loadFixture('summary/summary.json')['data']['mtcars']);

        expect(fn () => pf('Phit\tidy', $ols, effects: 'random'))->toThrow(Exception::class, "effects must be 'fixed' or 'ran_vals'")
            ->and(fn () => pf('Phit\tidy', $ols, effects: 'ran_vals'))->toThrow(Exception::class, "effects 'ran_vals' needs a mixed model fit");
    });

    it('gives ran_vals normal intervals', function () {
        skipMixedFitsUnderJit();
        $fit = pf('Phit\fit_mixed_gaussian', 'Reaction ~ Days + (Days | Subject)', loadFixture('summary/summary.json')['data']['sleepstudy']);
        $values = pf('Phit\tidy', $fit, confInt: true, effects: 'ran_vals');

        expect(abs($values['conf_high'][0] - $values['estimate'][0] - 1.959963984540054 * $values['std_error'][0]))->toBeLessThan(1e-9);
    });

    it('needs the random effects, their variances, and the labels to agree in size', function () {
        skipMixedFitsUnderJit();
        $fit = pf('Phit\fit_mixed_gaussian', 'Reaction ~ Days + (Days | Subject)', loadFixture('summary/summary.json')['data']['sleepstudy']);
        $missing = $fit;
        unset($missing['random'][0]['effects']);
        $variances = $fit;
        array_pop($variances['random'][0]['effects_cond_var']);
        $labels = $fit;
        array_pop($labels['random'][0]['group_labels']);
        $prefix = 'fit must be a result of a Phit model function: ';

        expect(pf('Phit\summary', $missing))->toBe(pf('Phit\summary', $fit))
            ->and(fn () => pf('Phit\tidy', $missing, effects: 'ran_vals'))->toThrow(Exception::class, $prefix.'missing random effects')
            ->and(fn () => pf('Phit\summary', $variances))->toThrow(Exception::class, $prefix.'random effects and their variances differ in size')
            ->and(fn () => pf('Phit\summary', $labels))->toThrow(Exception::class, $prefix.'random effects and group labels differ in size');
    });

    it('throws code 1001', function () {
        expect(fn () => pf('Phit\summary', []))->toThrow(function (Exception $e) {
            expect($e->getCode())->toBe(1001);
        });
    });

    it('checks the level before the fit', function (float $level) {
        expect(fn () => pf('Phit\tidy', [], true, $level))->toThrow(Exception::class, 'level must be in (0, 1)');
    })->with([0.0, 1.0, -0.5, NAN]);
});

describe('summary number formatting', function () {
    it('formats reals as R does', function (array $xs, int $digits, array $expected) {
        expect(Format::formatReal($xs, $digits))->toBe($expected);
    })->with([
        [[33.9908, -3.2056, 1.8878, 0.7539], 4, ['33.9908', '-3.2056', ' 1.8878', ' 0.7539']],
        [[1378.2, 960.5], 4, ['1378.2', ' 960.5']],
        [[2.220446049250313e-16], 1, ['2e-16']],
        [[0.0, -0.0], 4, ['0', '0']],
        [[1.0, NAN, INF], 4, ['  1', 'NaN', 'Inf']],
        [[123456789.0], 4, ['123456789']],
        [[0.000123456], 4, ['0.0001235']],
        [[1e-10, 0.5], 4, ['1e-10', '5e-01']],
        [[-0.00001234, 12.0], 4, ['-1.234e-05', ' 1.200e+01']],
        [[0.1 + 0.2], 4, ['0.3']],
        [[100.0, 1e-8], 4, ['1e+02', '1e-08']],
        [[9999.6], 4, ['10000']],
        [[999.96], 4, ['1000']],
        [[99999.6], 4, ['1e+05']],
        [[99999.6, 1.5], 4, ['1.0e+05', '1.5e+00']],
        [[99996.0], 4, ['99996']],
        [[9999.6, 1.0], 4, ['10000', '    1']],
    ]);

    it('formats p-values as format.pval does', function (array $ps, int $digits, array $expected) {
        expect(Format::formatPval($ps, $digits))->toBe($expected);
    })->with([
        [[1e-20, 0.000213, 0.004718, 0.000999], 3, ['< 2e-16', '0.000213', '0.004718', '0.000999']],
        [[3.594e-11], 4, ['3.594e-11']],
        [[1e-20], 3, ['<2e-16']],
        [[0.5, 1e-5, 1e-20], 3, ['0.5', '1e-05', '<2e-16']],
        [[0.04, 0.3], 3, ['0.04', '0.30']],
        [[NAN], 3, ['NA']],
    ]);

    it('formats %g as C does', function (float $x, string $expected) {
        expect(Format::formatG($x, 4))->toBe($expected);
    })->with([
        [0.83743, '0.8374'],
        [0.82, '0.82'],
        [48.08, '48.08'],
        [1.5e-7, '1.5e-07'],
        [123456.0, '1.235e+05'],
        [0.0001234, '0.0001234'],
        [100.0, '100'],
    ]);

    it('rounds beyond sprintf\'s 53 decimals as the exact decimal would', function (float $x, int $decimals, float $expected) {
        expect(Format::round($x, $decimals))->toBe($expected);
    })->with([
        [1.2345, 300, 1.2345],
        [1.23456789e-60, 62, 1.23e-60],
        [-1.25e-60, 61, -1.3e-60],
        [6e-61, 60, 1e-60],
        [4e-61, 60, 0.0],
        [1e-62, 60, 0.0],
    ]);
});
