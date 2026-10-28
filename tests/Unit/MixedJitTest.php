<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Jit;

// The mixed-model optimizers refuse to run while a JIT compiles the process:
// PHP 8.4.26's tracing JIT moved BOBYQA to a wrong minimum without an error
// (docs/notes/function-jit.md). These tests run only in the CI jobs that
// enable a JIT; everywhere else the mixed-model suites fit the same models.

$fixture = loadFixture('formula-lme4-r.json');

it('refuses every optimizer-backed mixed model under the JIT', function () use ($fixture) {
    if (! Jit::compiling()) {
        test()->markTestSkipped('no JIT is compiling this process');
    }
    $schools = $fixture['datasets']['schools'];
    $panel = $fixture['datasets']['panel'];
    $x = array_map(fn ($h) => [$h], $schools['study_hours']);
    $intercepts = [['z' => array_fill(0, count($schools['score']), [1]), 'groups' => $schools['school']]];
    foreach ([
        'fit_mixed bobyqa' => fn () => pf('Phit\fit_mixed', $x, $schools['score'], $intercepts),
        'fit_mixed nelder_mead' => fn () => pf('Phit\fit_mixed', $x, $schools['score'], $intercepts, optimizer: 'nelder_mead'),
        'fit_lmer bobyqa' => fn () => pf('Phit\fit_lmer', 'score ~ study_hours + (1 | school)', $schools),
        'fit_lmer nelder_mead' => fn () => pf('Phit\fit_lmer', 'score ~ study_hours + (1 | school)', $schools, optimizer: 'nelder_mead'),
        'fit_logistic_mixed' => fn () => pf('Phit\fit_logistic_mixed', array_map(fn ($v) => [$v], $panel['x']), $panel['yb'], $panel['subject']),
        'fit_glmer' => fn () => pf('Phit\fit_glmer', 'yb ~ x + (1 | subject)', $panel),
    ] as $label => $call) {
        $thrown = null;
        try {
            warnings($call);
        } catch (\Throwable $e) {
            $thrown = $e;
        }
        expect($thrown)->not->toBeNull("$label: no exception")
            ->and($thrown->getCode())->toBe(1007, "$label: ".$thrown->getMessage())
            ->and($thrown->getMessage())->toContain('opcache.jit=disable');
    }
});

it('still fits Phit\fit_mixed_intercept under the JIT, which uses no optimizer', function () use ($fixture) {
    if (! Jit::compiling()) {
        test()->markTestSkipped('no JIT is compiling this process');
    }
    $schools = $fixture['datasets']['schools'];
    $r = pf('Phit\fit_mixed_intercept', array_map(fn ($h) => [$h], $schools['study_hours']), $schools['score'], $schools['school']);

    expect($r['converged'])->toBeTrue();
});

it('validates input before refusing', function () use ($fixture) {
    $schools = $fixture['datasets']['schools'];
    $x = array_map(fn ($h) => [$h], $schools['study_hours']);
    $intercepts = [['z' => array_fill(0, count($schools['score']), [1]), 'groups' => $schools['school']]];

    expect(thrownCode(fn () => pf('Phit\fit_mixed', $x, $schools['score'], $intercepts, tolerance: NAN)))->toBe(1003);
});
