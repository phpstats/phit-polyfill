<?php

declare(strict_types=1);

// Port of the extension's tests/lmm.php: Phit\fit_mixed_gaussian_intercept_matrix against lme4 (lmm-r.json) and
// against lme4 and nlme on published datasets (lmm-public-r.json), plus invalid
// inputs, listwise deletion and typed group labels. The polyfill raises
// E_USER_WARNING where the extension raises E_WARNING; the text is the same.

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

describe('Phit\fit_mixed_gaussian_intercept_matrix against lme4', function () {
    foreach (loadFixture('lmm-r.json')['cases'] as $case) {
        it("matches lme4 on {$case['name']}", function () use ($case) {
            $r = pf('Phit\fit_mixed_gaussian_intercept_matrix', $case['x'], $case['y'], $case['groups'], intercept: false);
            foreach ($case['reference']['result'] as $key => $expected) {
                mixedCompare($r[$key], $expected, $case['name'].'/'.$key);
            }
            expect($r['method'])->toBe('ML')
                ->and($r['converged'])->toBeTrue()
                ->and($r['evaluations'])->toBeGreaterThanOrEqual(66)
                ->and($r['group_labels'])->toBe(array_values(array_unique($case['groups'])));
        });
    }
});

describe('Phit\fit_mixed_gaussian_intercept_matrix against lme4 and nlme on published data', function () {
    foreach (loadFixture('lmm-public-r.json')['cases'] as $case) {
        it("matches both backends on {$case['name']}", function () use ($case) {
            $groups = array_map(fn ($g) => $case['group_labels'][$g], $case['groups']);
            $r = pf('Phit\fit_mixed_gaussian_intercept_matrix', $case['x'], $case['y'], $groups, intercept: false);
            foreach (['reference', 'nlme_reference'] as $backend) {
                foreach ($case[$backend]['result'] as $key => $expected) {
                    $abs = $case['tolerance']['quantity_absolute'][$key] ?? $case['tolerance']['absolute'];
                    $rel = $case['tolerance']['quantity_relative'][$key] ?? $case['tolerance']['relative'];
                    mixedCompare($r[$key], $expected, $case['name'].'/'.$backend.'/'.$key, $abs, $rel);
                }
            }
            expect($r['group_labels'])->toBe($case['group_labels'])
                ->and($r['n_observations'])->toBe(count($case['y']))
                ->and($r['n_dropped'])->toBe(0);
        });
    }
});

describe('Phit\fit_mixed_gaussian_intercept_matrix inputs', function () {
    it('rejects invalid inputs and group boundaries', function () {
        foreach ([
            [[], [], []],
            [[[1], [1], [1]], [1, 2, 3], [1, 2]],
            [[[1], [1], [1]], [1, 2, 3], [1, 1, 1]],
            [[[1], [1], [1]], [1, 2, 3], [1, 2, 3]],
            [[[1], [1], [1], [1]], [1, 2, 3, 4], [1, 1, true, 2]],
            [[[null], [1], [1], [1]], [1, 2, 3, 4], [1.5, 1, 2, 2]],
            [[[1], [1], [1], [1]], [1, 2, 3, 4], [1, 1, 2.0, 2]],
            [[[1], [1], [1], [1]], [2, 2, 2, 2], [1, 1, 2, 2]],
            [[[1], [1], [1], [1]], [1, 1, 2, 2], [1, 1, 2, 2]],
        ] as $i => [$x, $y, $g]) {
            expect(thrownCode(fn () => pf('Phit\fit_mixed_gaussian_intercept_matrix', $x, $y, $g, intercept: false)))
                ->not->toBeNull("invalid/group boundary case $i must fail");
        }
    });

    it('deletes rows listwise across x, y and groups with one warning', function () {
        $c = loadFixture('lmm-r.json')['cases'][0];
        [$x, $y, $g] = [$c['x'], $c['y'], $c['groups']];
        $x[1][1] = null;
        $y[4] = NAN;
        $g[8] = null;
        $warnings = warnings(fn () => pf('Phit\fit_mixed_gaussian_intercept_matrix', $x, $y, $g, intercept: false), $r);
        expect($warnings)->toHaveCount(1)
            ->and($warnings[0][0])->toBe(E_USER_WARNING)
            ->and($warnings[0][1])->toContain('Phit\fit_mixed_gaussian_intercept_matrix')
            ->and($r['dropped_rows'])->toBe([1, 4, 8])
            ->and($r['n_dropped'])->toBe(3);

        $xx = $yy = $gg = [];
        foreach ($x as $i => $row) {
            if (! in_array($i, [1, 4, 8], true)) {
                $xx[] = $row;
                $yy[] = $y[$i];
                $gg[] = $g[$i];
            }
        }
        $complete = pf('Phit\fit_mixed_gaussian_intercept_matrix', $xx, $yy, $gg, intercept: false);
        mixedCompare($r['fitted'], $complete['fitted'], 'retained ordering');

        foreach (['raise', 'unknown'] as $policy) {
            expect(thrownCode(fn () => pf('Phit\fit_mixed_gaussian_intercept_matrix', $x, $y, $g, $policy, intercept: false)))
                ->not->toBeNull("missing policy $policy");
        }
    });

    it('keeps integer and string labels distinct even when their text matches', function () {
        $c = loadFixture('lmm-r.json')['cases'][5];
        $ids = [1, 1, 1, '1', '1', '1', 'A', 'A', 'A', 'B', 'B', 'B'];
        $r = pf('Phit\fit_mixed_gaussian_intercept_matrix', $c['x'], $c['y'], $ids, intercept: false);
        expect($r['group_labels'])->toBe([1, '1', 'A', 'B']);
    });
});
