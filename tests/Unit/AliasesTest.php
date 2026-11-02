<?php

declare(strict_types=1);

// The extension's tests/aliases.php: the 0.3 model names are gone, and
// Phit\fit_ols takes only a formula (both removed in phit 0.5.0).

use Phit\Polyfill\Polyfill;

// The datasets are in tests/Pest.php, shared with the parity suite.

describe('removed aliases', function () {
    it('neither implements nor defines the 0.3 names', function () {
        foreach (removedAliases() as $old) {
            expect(Polyfill::supports("Phit\\$old"))->toBeFalse("Phit\\$old is polyfilled")
                ->and(function_exists("Phit\\$old"))->toBeFalse("Phit\\$old is defined");
        }
    });

    it('gives fit_mixed_logistic no family or link', function () {
        $parameters = array_map(fn ($p) => $p->getName(), (new ReflectionMethod(...Polyfill::implementation('Phit\fit_mixed_logistic')))->getParameters());

        expect($parameters)->not->toContain('family')
            ->and($parameters)->not->toContain('link');
    });
});

describe('fit_ols', function () {
    it('takes a formula, positional or named', function () {
        ['schools' => $schools] = aliasData();
        $formula = pf('Phit\fit_ols', 'score ~ study_hours + prior_score', $schools);

        expect(pf('Phit\fit_ols', data: $schools, formula: 'score ~ study_hours + prior_score'))->toBe($formula)
            ->and($formula['columns'])->toBe(['(Intercept)', 'study_hours', 'prior_score']);
    });

    it('has the extension\'s parameters', function () {
        $parameters = array_map(fn ($p) => $p->getName(), (new ReflectionMethod(...Polyfill::implementation('Phit\fit_ols')))->getParameters());

        expect($parameters)->toBe(['formula', 'data', 'missing', 'weights']);
    });

    it('gives the same numbers as fit_ols_matrix on the formula\'s model matrix', function () {
        ['schools' => $schools, 'x' => $x, 'y' => $y] = aliasData();

        expect(pf('Phit\fit_ols_matrix', $x, $y)['coefficients'])
            ->toBe(pf('Phit\fit_ols', 'score ~ study_hours + prior_score', $schools)['coefficients']);
    });

    it('rejects a matrix', function () {
        ['x' => $x, 'y' => $y] = aliasData();

        // The formula is a string parameter, and x and y are not parameters.
        expect(fn () => pf('Phit\fit_ols', $x, $y))->toThrow(TypeError::class);
        expect(fn () => pf('Phit\fit_ols', x: $x, y: $y))->toThrow(Error::class);
    });

    it('prefixes the listwise-deletion warning with the function', function () {
        $x = [[1.0], [2.0], [3.0], [null], [5.0], [6.0]];
        $warning = fn (callable $call) => array_values(array_filter(warnings($call), fn ($w) => $w[0] === E_USER_WARNING))[0][1];

        expect($warning(fn () => pf('Phit\fit_logistic_matrix', $x, [0, 1, 0, 1, 1, 0])))
            ->toBe('Phit\fit_logistic_matrix(): dropped 1 of 6 observations with missing values (listwise deletion)')
            ->and($warning(fn () => pf('Phit\fit_ols', 'y ~ x', ['y' => [1.0, 2.0, 3.1, null, 5.2, 5.9], 'x' => [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])))
            ->toBe('Phit\fit_ols(): dropped 1 of 6 observations with missing values (listwise deletion)');
    });
});
