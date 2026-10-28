<?php

declare(strict_types=1);

// Port of the extension's tests/logistic-mixed.php: Phit\fit_logistic_mixed against
// 28 lme4 glmer fits (logistic-mixed-r.json) with Laplace and adaptive
// Gauss-Hermite quadrature, plus union deletion, invalid inputs and typed group
// labels. The polyfill raises E_USER_WARNING where the extension raises
// E_WARNING; the text is the same.

if (! function_exists('mixedCompare')) {
    /**
     * Compare as the extension's mixed-model scripts do: arrays by count and key,
     * booleans exactly, numbers finite and within abs + rel * |expected|.
     */
    function mixedCompare(mixed $actual, mixed $expected, string $label, float $abs = 1e-7, float $rel = 1e-6): void
    {
        if (is_array($expected)) {
            expect(is_array($actual) && count($actual) === count($expected))->toBeTrue("$label shape");
            foreach ($expected as $key => $value) {
                expect(array_key_exists($key, $actual))->toBeTrue("$label/$key missing");
                mixedCompare($actual[$key], $value, "$label/$key", $abs, $rel);
            }

            return;
        }
        if (is_bool($expected)) {
            expect($actual)->toBe($expected, "$label bool");

            return;
        }
        $numeric = is_int($actual) || is_float($actual);
        $shown = $numeric ? (string) $actual : get_debug_type($actual);
        expect($numeric && is_finite((float) $actual) && abs($actual - $expected) <= $abs + $rel * abs($expected))
            ->toBeTrue("$label: $shown != $expected");
    }
}

/** The grouped_offset_1 case with string school labels, as the script's deletion checks use it. */
function logisticMixedOffsetCase(): array
{
    $c = array_values(array_filter(
        loadFixture('logistic-mixed-r.json')['cases'],
        fn ($c) => $c['name'] === 'grouped_offset_1',
    ))[0];
    $c['groups'] = array_map(fn ($g) => 'school-'.$g, $c['groups']);
    $c['x'][2][1] = null;
    $c['y'][4] = NAN;
    $c['groups'][6] = null;
    $c['weights'][8] = null;
    $c['offset'][10] = NAN;

    return $c;
}

describe('Phit\fit_logistic_mixed against lme4', function () {
    foreach (loadFixture('logistic-mixed-r.json')['cases'] as $case) {
        it("matches glmer on {$case['name']}", function () use ($case) {
            skipMixedFitsUnderJit();
            $r = pf('Phit\fit_logistic_mixed', $case['x'], $case['y'], $case['groups'],
                weights: $case['weights'], offset: $case['offset'],
                quadraturePoints: $case['quadrature_points'], intercept: false);
            foreach ($case['expected'] as $key => $expected) {
                [$abs, $rel] = in_array($key, ['covariance', 'std_errors', 'z_statistics', 'conf_int', 'p_values'], true)
                    ? [2e-5, 2e-3]
                    : (in_array($key, ['log_likelihood', 'neg2_log_likelihood', 'aic', 'bic'], true) ? [1e-6, 1e-8] : [2e-5, 2e-4]);
                mixedCompare($r[$key], $expected, $case['name'].'/'.$key, $abs, $rel);
            }
            expect($r['n_parameters'])->toBe(count($case['x'][0]) + 1)
                ->and($r['df_residual'])->toBe(count($case['y']) - $r['n_parameters'])
                ->and($r['group_labels'])->toBe(array_values(array_unique($case['groups'])))
                ->and($r['method'])->toBe('ML')
                ->and($r['family'])->toBe('binomial')
                ->and($r['link'])->toBe('logit')
                ->and($r['converged'])->toBeTrue()
                ->and($r['integration'])->toBe($case['quadrature_points'] === 1 ? 'laplace' : 'adaptive_gauss_hermite');
        });
    }
});

describe('Phit\fit_logistic_mixed inputs', function () {
    it('deletes the union of missing rows with one warning and keeps alignment', function () {
        skipMixedFitsUnderJit();
        $c = logisticMixedOffsetCase();
        $warnings = warnings(fn () => pf('Phit\fit_logistic_mixed', $c['x'], $c['y'], $c['groups'],
            weights: $c['weights'], offset: $c['offset'], quadraturePoints: 9, intercept: false), $r);
        expect($r['dropped_rows'])->toBe([2, 4, 6, 8, 10])
            ->and($warnings)->toHaveCount(1)
            ->and($warnings[0][0])->toBe(E_USER_WARNING);

        $select = fn ($v) => array_values(array_filter($v, fn ($i) => ! in_array($i, [2, 4, 6, 8, 10], true), ARRAY_FILTER_USE_KEY));
        $complete = pf('Phit\fit_logistic_mixed', $select($c['x']), $select($c['y']), $select($c['groups']),
            weights: $select($c['weights']), offset: $select($c['offset']), quadraturePoints: 9, intercept: false);
        mixedCompare($r['fitted'], $complete['fitted'], 'retained alignment', 1e-10, 1e-10);
    });

    it('rejects invalid inputs', function () {
        $c = logisticMixedOffsetCase();
        $n = count($c['y']);
        foreach ([
            ['missing' => 'raise'], ['missing' => 'bad'], ['quadraturePoints' => 0], ['quadraturePoints' => 62],
            ['maxEvaluations' => 4], ['tolerance' => NAN],
            ['groups' => array_fill(0, $n, 1)], ['groups' => array_fill(0, $n, 1.0)], ['weights' => [1]], ['offset' => [0]],
            ['weights' => array_fill(0, $n, 0.5)], ['weights' => array_fill(0, $n, 1000001)], ['y' => array_fill(0, $n, 1)],
        ] as $options) {
            $arguments = array_replace(['x' => $c['x'], 'y' => $c['y'], 'groups' => $c['groups'], 'weights' => $c['weights'], 'offset' => $c['offset']], $options);
            expect(thrownCode(fn () => pf('Phit\fit_logistic_mixed', ...$arguments, intercept: false)))
                ->not->toBeNull('invalid logistic mixed inputs: '.implode(',', array_keys($options)));
        }
    });

    it('validates values even on a row that would be dropped', function () {
        $c = logisticMixedOffsetCase();
        $bad = $c['x'];
        $bad[2][0] = INF;
        expect(thrownCode(fn () => pf('Phit\fit_logistic_mixed', $bad, $c['y'], $c['groups'],
            weights: $c['weights'], offset: $c['offset'], intercept: false)))->not->toBeNull();
    });

    it('keeps typed IDs distinct and reports an exact boundary', function () {
        skipMixedFitsUnderJit();
        $r = pf('Phit\fit_logistic_mixed', array_fill(0, 8, [1]), [0, 1, 0, 1, 0, 1, 0, 1], [1, 1, '1', '1', 'A', 'A', 'B', 'B'], intercept: false);
        expect($r['group_labels'])->toBe([1, '1', 'A', 'B'])
            ->and($r['boundary'])->toBeTruthy();
    });
});
