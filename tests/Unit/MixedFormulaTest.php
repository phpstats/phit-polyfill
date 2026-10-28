<?php

declare(strict_types=1);

// Port of the extension's tests/formula-mixed.php: Phit\fit_lmer and Phit\fit_glmer
// formulas against lme4 (formula-lme4-r.json, from
// tests/generate-formula-lme4-fixtures.R), with tolerances that follow
// lmm-expanded-r (lmer) and logistic-mixed-r (glmer). Also checks that formula
// fits equal the matrix interfaces they wrap, the error codes, and the deletion
// warning. The polyfill raises E_USER_WARNING where the extension raises
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

/**
 * lme4 ranef(condVar = TRUE) postVar, keyed by lme4 level, against a term's
 * per-group blocks. $blocks maps a group label to its k x k block; (x || g)
 * records only the diagonal because lme4 splits it into separate one-column terms.
 */
function mixedCondVarCompare(array $blocks, array $labels, array $reference, string $label, float $abs, float $rel): void
{
    $index = array_flip(array_map('strval', $labels));
    expect(count($index))->toBe(count($reference['levels']), "$label levels");
    foreach ($reference['levels'] as $i => $level) {
        $block = $blocks[$index[(string) $level] ?? -1] ?? null;
        expect($block)->toBeArray("$label/$level missing");
        if (isset($reference['cond_var'])) {
            mixedCompare($block, $reference['cond_var'][$i], "$label/$level", $abs, $rel);
        } else {
            mixedCompare(array_map(fn ($j) => $block[$j][$j], array_keys($block)), $reference['cond_var_diagonal'][$i], "$label/$level", $abs, $rel);
        }
    }
}

describe('Phit\fit_lmer and Phit\fit_glmer against lme4', function () {
    $fixture = loadFixture('formula-lme4-r.json');
    foreach ($fixture['cases'] as $case) {
        it("matches {$case['function']} on {$case['name']}", function () use ($fixture, $case) {
            skipMixedFitsUnderJit();
            $data = $fixture['datasets'][$case['dataset']];
            $name = $case['name'];
            warnings(fn () => $case['function'] === 'lmer'
                ? pf('Phit\fit_lmer', $case['formula'], $data, reml: $case['reml'])
                : pf('Phit\fit_glmer', $case['formula'], $data, quadraturePoints: $case['quadrature_points']), $result);
            $expected = $case['expected'];
            expect($result['columns'])->toBe($expected['columns'], "$name: columns")
                ->and($result['n_observations'])->toBe($expected['n_observations'], "$name: observations")
                ->and(array_column($result['random'], 'group'))->toBe(array_column($case['bars'], 'group'), "$name: random groups")
                ->and(array_map(fn ($t) => (array) $t['columns'], $result['random']))
                ->toBe(array_map(fn ($b) => (array) $b['columns'], $case['bars']), "$name: random columns");
            foreach ($expected as $key => $value) {
                if (in_array($key, ['columns', 'n_observations'], true)) {
                    continue;
                }
                if ($key === 'random_cond_var') {
                    [$abs, $rel] = $case['function'] === 'lmer' ? [1e-5, 1e-5] : [2e-5, 2e-4];
                    foreach ($value as $t => $reference) {
                        [$blocks, $labels] = $case['function'] === 'lmer'
                            ? [$result['random'][$t]['effects_cond_var'], $result['random'][$t]['group_labels']]
                            : [array_map(fn ($v) => [[$v]], $result['random_effects_cond_var']), $result['group_labels']];
                        mixedCondVarCompare($blocks, $labels, $reference, "$name/random_cond_var/$t", $abs, $rel);
                    }

                    continue;
                }
                $likelihood = in_array($key, ['log_likelihood', 'aic', 'bic'], true);
                if ($case['function'] === 'lmer') {
                    [$abs, $rel] = $likelihood ? [5e-6, 1e-9] : [1e-5, 1e-5];
                    $actual = $key === 'random_covariances' ? array_column($result['random'], 'covariance') : $result[$key];
                } else {
                    [$abs, $rel] = $likelihood ? [1e-6, 1e-8] : (in_array($key, ['std_errors', 'covariance'], true) ? [2e-5, 2e-3] : [2e-5, 2e-4]);
                    $actual = $result[$key];
                }
                mixedCompare($actual, $value, "$name/$key", $abs, $rel);
            }
            expect(array_keys($result['named']))->toBe($expected['columns'], "$name: named");
        });
    }
});

describe('mixed formula interfaces', function () {
    it('fits Phit\fit_lmer exactly as the Phit\fit_mixed matrix interface it wraps', function () {
        skipMixedFitsUnderJit();
        $schools = loadFixture('formula-lme4-r.json')['datasets']['schools'];
        $lmer = pf('Phit\fit_lmer', 'score ~ study_hours + (1 | school)', $schools, reml: false);
        $matrix = pf('Phit\fit_mixed', array_map(fn ($h) => [$h], $schools['study_hours']), $schools['score'],
            [['z' => array_fill(0, count($schools['score']), [1]), 'groups' => $schools['school']]]);
        mixedCompare($lmer['coefficients'], $matrix['coefficients'], 'lmer equals Phit\fit_mixed', 1e-12, 1e-12);
        expect($lmer['method'])->toBe('ML')
            ->and($lmer['random'][0]['label'])->toBe('1 | school')
            ->and($lmer['random'][0]['group_labels'])->toBe($matrix['random'][0]['group_labels']);
    });

    it('reports the same Phit\fit_mixed_intercept conditional variances as lme4', function () {
        $fixture = loadFixture('formula-lme4-r.json');
        $schools = $fixture['datasets']['schools'];
        $ml = array_values(array_filter($fixture['cases'], fn ($c) => $c['name'] === 'lmer_schools_intercept_ml'))[0];
        $lmm = pf('Phit\fit_mixed_intercept', array_map(null, $schools['study_hours'], $schools['prior_score']), $schools['score'], $schools['school']);
        mixedCondVarCompare(array_map(fn ($v) => [[$v]], $lmm['random_effects_cond_var']), $lmm['group_labels'],
            $ml['expected']['random_cond_var'][0], 'Phit\fit_mixed_intercept/random_effects_cond_var', 1e-5, 1e-5);
    });

    it('reports Phit\fit_glmer metadata', function () {
        skipMixedFitsUnderJit();
        $panel = loadFixture('formula-lme4-r.json')['datasets']['panel'];
        warnings(fn () => pf('Phit\fit_glmer', 'yb ~ x + (1 | subject)', $panel), $glmer);
        expect($glmer['family'])->toBe('binomial')
            ->and($glmer['random'][0]['group'])->toBe('subject')
            ->and($glmer['n_dropped'])->toBe(2);
    });

    it('throws the documented codes and messages', function () {
        $fixture = loadFixture('formula-lme4-r.json');
        $schools = $fixture['datasets']['schools'];
        $panel = $fixture['datasets']['panel'];
        foreach ([
            ['no random term', fn () => pf('Phit\fit_lmer', 'score ~ study_hours', $schools), 1001, 'random-effect term'],
            ['glmer slopes', fn () => pf('Phit\fit_glmer', 'yb ~ x + (x | subject)', $panel), 1007, 'exactly one random intercept'],
            ['glmer two terms', fn () => pf('Phit\fit_glmer', 'yb ~ x + (1 | subject) + (1 | item)', $panel), 1007, 'exactly one random intercept'],
            ['glmer family', fn () => pf('Phit\fit_glmer', 'yb ~ x + (1 | subject)', $panel, 'poisson'), 1007, 'binomial'],
            ['glmer link', fn () => pf('Phit\fit_glmer', 'yb ~ x + (1 | subject)', $panel, link: 'probit'), 1007, 'logit'],
            ['lmer optimizer', fn () => pf('Phit\fit_lmer', 'y ~ x + (1 | subject)', $panel, optimizer: 'nloptwrap'), 1007, 'optimizer'],
            ['lmer raise', fn () => pf('Phit\fit_lmer', 'y ~ x + (1 | subject)', $panel, missing: 'raise'), 1001, "missing='raise'"],
        ] as [$label, $call, $code, $fragment]) {
            $thrown = null;
            try {
                warnings($call);
            } catch (\Throwable $e) {
                $thrown = $e;
            }
            expect($thrown)->not->toBeNull("$label: no exception")
                ->and($thrown->getCode())->toBe($code, "$label: ".$thrown->getMessage())
                ->and($thrown->getMessage())->toContain($fragment);
        }
    });

    it('warns once when Phit\fit_lmer deletes rows', function () {
        skipMixedFitsUnderJit();
        $panel = loadFixture('formula-lme4-r.json')['datasets']['panel'];
        $warnings = warnings(fn () => pf('Phit\fit_lmer', 'y ~ x + (1 | subject)', $panel), $r);
        expect($warnings)->toHaveCount(1)
            ->and($warnings[0][0])->toBe(E_USER_WARNING)
            ->and($warnings[0][1])->toContain('Phit\fit_lmer dropped 5 of 384')
            ->and($r['dropped_rows'])->toBe([4, 39, 76, 210, 299]);
    });
});
