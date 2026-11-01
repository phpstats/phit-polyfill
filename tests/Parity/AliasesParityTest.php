<?php

declare(strict_types=1);

// Compares phit 0.4's deprecated aliases and Phit\fit_ols's dispatch with the
// loaded extension: the same result or exception, and the same deprecations
// and warnings, in order.

/** Runs $fn through $call and records its result or exception and every notice it raised. */
function aliasRun(callable $call, string $fn, array $arguments): array
{
    $notices = [];
    set_error_handler(function (int $level, string $message) use (&$notices): bool {
        // The polyfill raises E_USER_WARNING where the extension raises E_WARNING.
        $notices[] = [$level === E_USER_WARNING ? E_WARNING : $level, $message];

        return true;
    });
    try {
        $outcome = ['ok', $call($fn, ...$arguments)];
    } catch (\Throwable $e) {
        $outcome = ['error', get_class($e), $e->getCode(), $e->getMessage()];
    } finally {
        restore_error_handler();
    }

    return [$outcome, $notices];
}

function aliasCompare(string $fn, array $arguments): void
{
    $native = aliasRun(native(...), $fn, $arguments);
    $polyfill = aliasRun(pf(...), $fn, $arguments);
    parityLog('aliases', $fn, $polyfill, $native, $native[0][0] === 'error');
    if (serialize($native) !== serialize($polyfill)) {
        throw new \RuntimeException("$fn differs\nnative ".var_export([$native[0][0], $native[0][3] ?? null, $native[1]], true)
            ."\npolyfill ".var_export([$polyfill[0][0], $polyfill[0][3] ?? null, $polyfill[1]], true));
    }
    expect(true)->toBeTrue();
}

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

describe('aliases parity', function () {
    foreach (aliasCases() as [$old, $new, $positional, $named]) {
        it("matches Phit\\$old and Phit\\$new", function () use ($old, $new, $positional, $named) {
            aliasCompare("Phit\\$old", [...$positional, ...$named]);
            aliasCompare("Phit\\$new", [...$positional, ...$named]);
        });
    }

    it('matches fit_glmer\'s family and link checks', function () {
        ['herds' => $herds] = aliasData();
        foreach ([['binomial', 'logit'], ['poisson', ''], ['binomial', 'probit']] as [$family, $link]) {
            aliasCompare('Phit\fit_glmer', ['cbind(incidence, healthy) ~ period + (1 | herd)', $herds, $family, $link]);
        }
    });

    it('matches fit_ols on formulas and matrices, positional and named', function () {
        ['schools' => $schools, 'x' => $x, 'y' => $y, 'weights' => $weights] = aliasData();
        foreach ([
            ['score ~ study_hours + prior_score', $schools],
            ['formula' => 'score ~ study_hours + prior_score', 'data' => $schools, 'weights' => $weights],
            [$x, $y, 'weights' => $weights],
            ['x' => $x, 'y' => $y],
            [$x, 'y' => $y, 'intercept' => false, 'offset' => array_fill(0, count($y), 1.0)],
        ] as $arguments) {
            aliasCompare('Phit\fit_ols', $arguments);
        }
    });

    it('matches fit_ols\'s errors', function () {
        foreach (fitOlsErrorCases() as [$arguments]) {
            aliasCompare('Phit\fit_ols', $arguments);
        }
    });

    it('matches the listwise-deletion warnings through aliases and fit_ols', function () {
        $x = [[1.0], [2.0], [3.0], [null], [5.0], [6.0]];
        $y = [1.0, 2.1, 2.9, 4.2, 5.1, 6.3];
        $data = ['y' => [1.0, 2.0, 3.1, null, 5.2, 5.9], 'x' => [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]];
        aliasCompare('Phit\fit_ols', [$x, $y]);
        aliasCompare('Phit\fit_ols', ['y ~ x', $data]);
        aliasCompare('Phit\fit_lm', ['y ~ x', $data]);
        aliasCompare('Phit\fit_glm', ['y ~ x', $data]);
        aliasCompare('Phit\fit_ols_matrix', [$x, $y]);
        aliasCompare('Phit\fit_logistic', [$x, [0, 1, 0, 1, 1, 0]]);
        aliasCompare('Phit\fit_poisson', [$x, [1, 2, 3, 4, 5, 6]]);
    });
});
