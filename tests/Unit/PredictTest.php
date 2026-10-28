<?php

declare(strict_types=1);

// Port of the extension's tests/predict.php: Phit\predict versus R's predict()
// (predict-r.json). Fixed-model tolerances follow r-compatibility; GLM standard
// errors allow the documented final-versus-penultimate IRLS weight difference;
// mixed tolerances follow formula-mixed.

function predictCompare(mixed $actual, mixed $expected, string $label, float $abs, float $rel): void
{
    if (is_array($expected)) {
        expect(is_array($actual) && count($actual) === count($expected))->toBeTrue("$label shape");
        foreach ($expected as $i => $value) {
            predictCompare($actual[$i] ?? null, $value, "$label/$i", $abs, $rel);
        }

        return;
    }
    if ($expected === null) {
        expect($actual)->toBeNull("$label: expected a missing prediction");

        return;
    }
    expect(is_float($actual) || is_int($actual))->toBeTrue("$label: missing prediction");
    expect(abs($actual - $expected) <= $abs + $rel * abs($expected))->toBeTrue("$label: $actual != $expected");
}

function predictFit(array $case, array $data): array
{
    $args = $case['args'];
    $result = null;
    warnings(function () use ($case, $data, $args, &$result) {
        $result = match ($case['function']) {
            'lm' => pf('Phit\fit_lm', $case['formula'], $data, weights: $args['weights'] ?? null),
            'glm' => pf('Phit\fit_glm', $case['formula'], $data, family: $args['family'], link: $args['link']),
            'lmer' => pf('Phit\fit_lmer', $case['formula'], $data, reml: $args['reml'] ?? true),
            'glmer' => pf('Phit\fit_glmer', $case['formula'], $data),
        };
    });

    return $result;
}

function predictMessage(callable $f): ?string
{
    try {
        $f();
    } catch (\Exception $e) {
        return $e->getMessage();
    }

    return null;
}

$fixture = loadFixture('predict-r.json');

foreach ($fixture['cases'] as $case) {
    it("matches R's predict() for {$case['name']}", function () use ($case, $fixture) {
        if (in_array($case['function'], ['lmer', 'glmer'], true)) {
            skipMixedFitsUnderJit();
        }
        $fit = predictFit($case, $fixture['datasets'][$case['dataset']]);
        $options = [];
        foreach ($case['options'] as $option => $value) {
            $options[['se_fit' => 'seFit', 'allow_new_levels' => 'allowNewLevels'][$option] ?? $option] = $value;
        }
        $result = pf('Phit\predict', $fit, $fixture['newdata'][$case['dataset']], ...$options);
        foreach ($case['expected'] as $key => $expected) {
            [$abs, $rel] = match (true) {
                $case['function'] === 'lmer' => [1e-5, 1e-5],
                $case['function'] === 'glmer' => [2e-5, 2e-4],
                // R's penultimate-weight vcov; near separation it differs by about 1e-4.
                $case['function'] === 'glm' && $key === 'se_fit' => [1e-6, 1e-3],
                default => [1e-7, 1e-6],
            };
            $actual = $key === 'se_fit_final' ? $result['se_fit'] : ($result[$key] ?? null);
            expect($actual !== null || $expected === null)->toBeTrue("{$case['name']}: missing $key");
            predictCompare($actual, $expected, "{$case['name']}/$key", $abs, $rel);
        }
        $missing = array_keys(array_filter($case['expected']['fit'], fn ($v) => $v === null));
        expect($result['missing_rows'])->toBe($missing);
    });
}

it('reproduces fitted values on the training data', function () use ($fixture) {
    $cars = $fixture['datasets']['cars'];
    $fit = pf('Phit\fit_lm', 'mpg ~ poly(hp, 3) + factor(cyl) + wt:am', $cars);
    predictCompare(pf('Phit\predict', $fit, $cars)['fit'], $fit['fitted'], 'lm training data', 1e-12, 1e-12);
    $fit = pf('Phit\fit_glm', 'am ~ wt + hp', $cars, family: 'binomial');
    predictCompare(pf('Phit\predict', $fit, $cars, type: 'response')['fit'], $fit['fitted'], 'glm training data', 1e-12, 1e-10);
    skipMixedFitsUnderJit();
    $schools = $fixture['datasets']['schools'];
    $fit = pf('Phit\fit_lmer', 'score ~ study_hours + (1 | school)', $schools);
    predictCompare(pf('Phit\predict', $fit, $schools)['fit'], $fit['fitted'], 'lmer training data', 1e-10, 1e-10);
    predictCompare(pf('Phit\predict', $fit, $schools, random: false)['fit'], $fit['fixed_fitted'], 'lmer population', 1e-10, 1e-10);
});

it("raises R's errors", function () use ($fixture) {
    skipMixedFitsUnderJit();
    $cars = $fixture['datasets']['cars'];
    $schools = $fixture['datasets']['schools'];
    $fit = pf('Phit\fit_lm', 'mpg ~ wt + gear_label', $cars);
    $lmer = pf('Phit\fit_lmer', 'score ~ study_hours + (1 | school)', $schools);
    $glm = pf('Phit\fit_glm', 'am ~ wt', $cars, family: 'binomial');
    $cases = [
        'factor gear_label has new level six' => fn () => pf('Phit\predict', $fit, ['wt' => [3.0], 'gear_label' => ['six']]),
        'was a factor' => fn () => pf('Phit\predict', $fit, ['wt' => [3.0], 'gear_label' => [4]]),
        'was numeric' => fn () => pf('Phit\predict', $fit, ['wt' => ['heavy'], 'gear_label' => ['four']]),
        "unknown column 'wt'" => fn () => pf('Phit\predict', $fit, ['gear_label' => ['four']]),
        'type must be' => fn () => pf('Phit\predict', $fit, ['wt' => [3.0], 'gear_label' => ['four']], type: 'probability'),
        'new levels detected in newdata: school=nowhere' => fn () => pf('Phit\predict', $lmer, ['study_hours' => [3.0], 'school' => ['nowhere']]),
        'not available for mixed' => fn () => pf('Phit\predict', $lmer, ['study_hours' => [3.0], 'school' => ['nowhere']], seFit: true),
        'intervals are available for Phit\fit_lm' => fn () => pf('Phit\predict', $glm, $cars, interval: 'confidence'),
        'fit must be a result' => fn () => pf('Phit\predict', ['formula' => 'y ~ x'], $cars),
    ];
    foreach ($cases as $needle => $f) {
        expect((string) predictMessage($f))->toContain($needle);
    }
    // A response column is optional; population predictions need no grouping column.
    $p = pf('Phit\predict', $lmer, ['study_hours' => [1.0, 2.0]], random: false);
    expect($p['fit'])->toHaveCount(2)->and($p['missing_rows'])->toBe([]);
});

it('reports the fitted factor levels and poly() coefficients', function () {
    $data = [
        'y' => [1, 2, 3, 4, 5, 6.5, 7, 9, 8, 10],
        'x' => [1, 2, 3, 5, 8, 9, 4, 6, 7, 2],
        'a' => ['p', 'q', 'p', 'r', 'q', 'p', 'r', 'r', 'q', 'p'],
    ];
    $fit = pf('Phit\fit_lm', 'y ~ poly(x, 2) + a + I(x > 3)', $data);
    expect($fit['predvars']['levels'])->toBe(['a' => ['p', 'q', 'r'], 'I(x > 3)' => ['FALSE', 'TRUE']])
        ->and(array_keys($fit['predvars']['poly']))->toBe(['poly(x, 2)'])
        ->and($fit['predvars']['poly']['poly(x, 2)']['norm2'][0])->toBe(1.0)
        ->and($fit['predvars']['poly']['poly(x, 2)']['alpha'])->toHaveCount(2);
    // Out-of-range values use the fitted basis instead of refitting it.
    $p = pf('Phit\predict', $fit, ['x' => [20.0], 'a' => ['q']]);
    expect($p['fit'])->toHaveCount(1)->and($p['n_predicted'])->toBe(1);
});
