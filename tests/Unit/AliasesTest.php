<?php

declare(strict_types=1);

// phit 0.4's transition (the extension's tests/aliases.php): every 0.3 model
// name is a deprecated alias that returns exactly what its replacement returns,
// with one E_USER_DEPRECATED per call, and Phit\fit_ols dispatches on its
// first argument without changing any call's meaning.

use Phit\Polyfill\Polyfill;

// The datasets and cases are in tests/Pest.php, shared with the parity suite.

describe('deprecated aliases', function () {
    foreach (aliasCases() as [$old, $new, $positional, $named]) {
        it("Phit\\$old returns exactly Phit\\$new's result, with one deprecation", function () use ($old, $new, $positional, $named) {
            [$expected, $none] = deprecations(fn () => pf("Phit\\$new", ...$positional, ...$named));
            expect($none)->toBe([]);
            [$actual, $messages] = deprecations(fn () => pf("Phit\\$old", ...$positional, ...$named));
            expect($actual)->toBe($expected)
                ->and($messages)->toBe(["Phit\\$old(): Deprecated since phit 0.4.0; use Phit\\$new() instead"]);
        });
    }

    it('lists every alias with its replacement', function () {
        $expected = [];
        foreach (aliasCases() as [$old, $new]) {
            $expected["Phit\\$old"] = "Phit\\$new";
        }
        ksort($expected);
        $listed = Phit\Polyfill\FunctionMap::DEPRECATED;
        ksort($listed);

        expect($listed)->toBe($expected);
    });

    it('keeps fit_glmer\'s family and link, which fit_mixed_logistic drops', function () {
        ['herds' => $herds] = aliasData();
        [$glmer] = deprecations(fn () => pf('Phit\fit_glmer', 'cbind(incidence, healthy) ~ period + (1 | herd)', $herds, 'binomial', 'logit'));
        $parameters = array_map(fn ($p) => $p->getName(), (new ReflectionMethod(...Polyfill::implementation('Phit\fit_mixed_logistic')))->getParameters());

        expect($glmer)->toBe(pf('Phit\fit_mixed_logistic', 'cbind(incidence, healthy) ~ period + (1 | herd)', $herds))
            ->and($parameters)->not->toContain('family')
            ->and($parameters)->not->toContain('link');
    });

    it('raises the deprecation before fit_glmer rejects a family', function () {
        ['herds' => $herds] = aliasData();
        $thrown = null;
        [, $messages] = deprecations(function () use ($herds, &$thrown) {
            try {
                pf('Phit\fit_glmer', 'cbind(incidence, healthy) ~ period + (1 | herd)', $herds, 'poisson');
            } catch (Exception $e) {
                $thrown = $e;
            }
        });

        expect($thrown?->getCode())->toBe(1007)
            ->and($thrown?->getMessage())->toBe("Phit\\fit_glmer supports family 'binomial' with the logit link only (got 'poisson'/'')")
            ->and($messages)->toBe(['Phit\fit_glmer(): Deprecated since phit 0.4.0; use Phit\fit_mixed_logistic() instead']);
    });

    it('prefixes the listwise-deletion warning with the name the caller used', function () {
        $x = [[1.0], [2.0], [3.0], [null], [5.0], [6.0]];
        $y = [1.0, 2.1, 2.9, 4.2, 5.1, 6.3];
        $warning = fn (callable $call) => array_values(array_filter(warnings($call), fn ($w) => $w[0] === E_USER_WARNING))[0][1];

        expect($warning(fn () => pf('Phit\fit_logistic', $x, [0, 1, 0, 1, 1, 0])))
            ->toBe('Phit\fit_logistic(): dropped 1 of 6 observations with missing values (listwise deletion)')
            ->and($warning(fn () => pf('Phit\fit_ols', $x, $y)))
            ->toBe('Phit\fit_ols(): dropped 1 of 6 observations with missing values (listwise deletion)')
            ->and($warning(fn () => pf('Phit\fit_ols_matrix', $x, $y)))
            ->toBe('Phit\fit_ols_matrix(): dropped 1 of 6 observations with missing values (listwise deletion)');
    });

    it('raises the deprecation from the bootstrap\'s functions', function () {
        ['x' => $x, 'pass' => $pass] = aliasData();
        [$result, $messages] = deprecations(fn () => Phit\fit_probit($x, $pass));

        expect($result)->toBe(pf('Phit\fit_probit_matrix', $x, $pass))
            ->and($messages)->toBe(['Phit\fit_probit(): Deprecated since phit 0.4.0; use Phit\fit_probit_matrix() instead']);
    })->skip(fn () => extension_loaded('phit') || Phit\Polyfill\Internal\Jit::compiling(), 'the bootstrap defines nothing when the extension is loaded, and its functions refuse to run under the JIT');
});

describe('fit_ols dispatch', function () {
    it('fits a formula, positional or named, with no deprecation', function () {
        ['schools' => $schools] = aliasData();
        $formula = pf('Phit\fit_ols', 'score ~ study_hours + prior_score', $schools);
        foreach ([
            fn () => pf('Phit\fit_ols', 'score ~ study_hours + prior_score', $schools),
            fn () => pf('Phit\fit_ols', formula: 'score ~ study_hours + prior_score', data: $schools),
            fn () => pf('Phit\fit_ols', data: $schools, formula: 'score ~ study_hours + prior_score'),
        ] as $call) {
            [$result, $messages] = deprecations($call);
            expect($result)->toBe($formula)->and($messages)->toBe([]);
        }
        expect($formula['columns'])->toBe(['(Intercept)', 'study_hours', 'prior_score'])
            ->and(deprecations(fn () => pf('Phit\fit_lm', 'score ~ study_hours + prior_score', $schools))[0])->toBe($formula);
    });

    it('fits a matrix, positional or named, with one deprecation', function () {
        ['x' => $x, 'y' => $y, 'weights' => $weights] = aliasData();
        $matrix = pf('Phit\fit_ols_matrix', $x, $y, weights: $weights);
        foreach ([
            fn () => pf('Phit\fit_ols', $x, $y, weights: $weights),
            fn () => pf('Phit\fit_ols', x: $x, y: $y, weights: $weights),
            fn () => pf('Phit\fit_ols', $x, y: $y, weights: $weights),
            fn () => pf('Phit\fit_ols', $x, $y, 'drop', $weights, null, true),
        ] as $call) {
            [$result, $messages] = deprecations($call);
            expect($result)->toBe($matrix)
                ->and($messages)->toBe(['Phit\fit_ols(): Passing a matrix is deprecated since phit 0.4.0; use Phit\fit_ols_matrix() instead']);
        }
        $offset = array_fill(0, count($y), 1.0);
        [$result] = deprecations(fn () => pf('Phit\fit_ols', $x, $y, intercept: false, offset: $offset));
        expect($result)->toBe(pf('Phit\fit_ols_matrix', $x, $y, intercept: false, offset: $offset));
    });

    it('gives the matrix and formula fits the same numbers', function () {
        ['schools' => $schools, 'x' => $x, 'y' => $y] = aliasData();

        expect(pf('Phit\fit_ols_matrix', $x, $y)['coefficients'])
            ->toBe(pf('Phit\fit_ols', 'score ~ study_hours + prior_score', $schools)['coefficients']);
    });

    foreach (fitOlsErrorCases() as $label => [$arguments, $fragment]) {
        it("rejects $label with code 1001 and no deprecation", function () use ($arguments, $fragment) {
            $thrown = null;
            [, $messages] = deprecations(function () use ($arguments, &$thrown) {
                try {
                    pf('Phit\fit_ols', ...$arguments);
                } catch (Exception $e) {
                    $thrown = $e;
                }
            });

            expect($thrown?->getCode())->toBe(1001)
                ->and($thrown?->getMessage())->toContain('Phit\fit_ols')
                ->and($thrown?->getMessage())->toContain($fragment)
                ->and($messages)->toBe([]);
        });
    }
});
