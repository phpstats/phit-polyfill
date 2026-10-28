<?php

declare(strict_types=1);

// Port of the extension's tests/lmm-expanded.php: Phit\fit_mixed with random
// slopes, crossed, nested, weighted and boundary models against lme4's four
// optimizers (lmm-expanded-r.json), each fitted with BOBYQA and with the native
// Nelder-Mead simplex, plus union deletion, invalid inputs and typed group
// labels. The comparison helper follows tests/r-compatibility.php. The polyfill
// raises E_USER_WARNING where the extension raises E_WARNING; the text is the same.
//
// The whole fixture takes several minutes in PHP without JIT, so only the
// models in LMM_EXPANDED_ALWAYS run by default; set POLYFILL_SLOW=1 for the rest.

const LMM_EXPANDED_ALWAYS = ['random_slopes_ml', 'crossed', 'sleepstudy_reml'];

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

/** Reference terms: (x || g) keeps only conditional-variance diagonals, as lme4 reports them. */
function lmmExpandedTerms(array $random, array $input): array
{
    return array_map(function ($t, $in) {
        $out = ['covariance' => $t['covariance'], 'effects' => $t['effects']];
        if ($in['correlated'] ?? true) {
            $out['effects_cond_var'] = $t['effects_cond_var'];
        } else {
            $out['effects_cond_var_diagonal'] = array_map(fn ($b) => array_map(fn ($j) => $b[$j][$j], array_keys($b)), $t['effects_cond_var']);
        }

        return $out;
    }, $random, $input);
}

/** The crossed case with prefixed labels and one missing value in each of four inputs. */
function lmmExpandedMissingCase(): array
{
    $c = loadFixture('lmm-expanded-r.json')['cases'][3];
    $c['terms'][0]['groups'] = array_map(fn ($g) => 'school-'.$g, $c['terms'][0]['groups']);
    $c['terms'][1]['groups'] = array_map(fn ($g) => 'teacher-'.$g, $c['terms'][1]['groups']);
    $c['terms'][0]['z'][1][1] = null;
    $c['terms'][1]['groups'][3] = null;
    $c['weights'][5] = NAN;
    $c['offset'][7] = null;

    return $c;
}

describe('Phit\fit_mixed against lme4 optimizers', function () {
    foreach (loadFixture('lmm-expanded-r.json')['cases'] as $c) {
        $slow = ! in_array($c['name'], LMM_EXPANDED_ALWAYS, true);
        it("matches lme4 on {$c['name']} with BOBYQA and Nelder-Mead", function () use ($c) {
            skipMixedFitsUnderJit();
            $r = pf('Phit\fit_mixed', x: $c['x'], y: $c['y'], random: $c['terms'], reml: $c['reml'],
                weights: $c['weights'], offset: $c['offset'], intercept: false);
            $r['terms'] = lmmExpandedTerms($r['random'], $c['terms']);
            foreach ($c['references'] as $optimizer => $ref) {
                foreach ($ref['result'] as $k => $v) {
                    [$abs, $rel] = in_array($k, ['log_likelihood', 'aic', 'bic'], true)
                        ? [5e-6, 1e-9]
                        : ($optimizer === 'bobyqa' ? [1e-5, 1e-5] : [1e-3, 3e-4]);
                    mixedCompare($r[$k], $v, $c['name'].'/'.$optimizer.'/'.$k, $abs, $rel);
                }
            }

            $nm = pf('Phit\fit_mixed', $c['x'], $c['y'], $c['terms'], reml: $c['reml'],
                weights: $c['weights'], offset: $c['offset'], optimizer: 'nelder_mead', intercept: false);
            $nm['terms'] = lmmExpandedTerms($nm['random'], $c['terms']);
            foreach ($c['references']['Nelder_Mead']['result'] as $k => $v) {
                [$abs, $rel] = in_array($k, ['log_likelihood', 'aic', 'bic'], true) ? [5e-6, 1e-9] : [1e-3, 3e-4];
                mixedCompare($nm[$k], $v, $c['name'].'/native_Nelder_Mead/'.$k, $abs, $rel);
            }

            expect($nm['optimizer'])->toBe('nelder_mead')
                ->and($nm['converged'])->toBeTrue()
                ->and($r['method'])->toBe($c['reml'] ? 'REML' : 'ML')
                ->and($r['optimizer'])->toBe('bobyqa')
                ->and($r['converged'])->toBeTrue();
            foreach ($c['terms'] as $i => $term) {
                expect($r['random'][$i]['group_labels'])->toBe(array_values(array_unique($term['groups'])), "term $i group order");
            }
        })->skip($slow && getenv('POLYFILL_SLOW') !== '1', 'slow lme4 model; set POLYFILL_SLOW=1');
    }
});

describe('Phit\fit_mixed inputs', function () {
    it('deletes the union of missing rows across all inputs with one warning', function () {
        skipMixedFitsUnderJit();
        $c = lmmExpandedMissingCase();
        $warnings = warnings(fn () => pf('Phit\fit_mixed', $c['x'], $c['y'], $c['terms'],
            weights: $c['weights'], offset: $c['offset'], intercept: false), $r);
        expect($r['dropped_rows'])->toBe([1, 3, 5, 7])
            ->and($warnings)->toHaveCount(1)
            ->and($warnings[0][0])->toBe(E_USER_WARNING);

        $select = fn ($v) => array_values(array_filter($v, fn ($i) => ! in_array($i, [1, 3, 5, 7], true), ARRAY_FILTER_USE_KEY));
        $terms = array_map(fn ($t) => ['z' => $select($t['z']), 'groups' => $select($t['groups']), 'correlated' => $t['correlated']], $c['terms']);
        $complete = pf('Phit\fit_mixed', $select($c['x']), $select($c['y']), $terms,
            weights: $select($c['weights']), offset: $select($c['offset']), intercept: false);
        mixedCompare($r['fitted'], $complete['fitted'], 'all retained rows aligned');
    });

    it('rejects invalid inputs', function () {
        $c = lmmExpandedMissingCase();
        $n = count($c['y']);
        foreach ([
            ['random' => []],
            ['random' => [['z' => [[1]], 'groups' => [1]]]],
            ['random' => [['z' => array_fill(0, $n, [1]), 'groups' => array_fill(0, $n, true)]]],
            ['random' => [['z' => array_fill(0, $n, [1]), 'groups' => array_fill(0, $n, 1), 'correlated' => 'yes']]],
            ['optimizer' => 'nloptwrap'], ['objectiveTolerance' => NAN], ['missing' => 'raise'], ['tolerance' => NAN],
            ['start' => [-1]], ['maxEvaluations' => -1], ['maxEvaluations' => 4],
        ] as $i => $options) {
            $arguments = array_replace(['x' => $c['x'], 'y' => $c['y'], 'random' => $c['terms'], 'weights' => $c['weights'], 'offset' => $c['offset']], $options);
            expect(thrownCode(fn () => pf('Phit\fit_mixed', ...$arguments, intercept: false)))
                ->not->toBeNull("invalid general mixed inputs case $i");
        }
    });

    it('keeps integer and string labels distinct', function () {
        skipMixedFitsUnderJit();
        $typed = pf('Phit\fit_mixed', array_fill(0, 8, [1]), [-3, -2, -1, 1, 2, 3, 4, 5],
            [['z' => array_fill(0, 8, [1]), 'groups' => [1, 1, '1', '1', 'A', 'A', 'B', 'B']]], intercept: false);
        expect($typed['random'][0]['group_labels'])->toBe([1, '1', 'A', 'B']);
    });
});
