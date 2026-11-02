<?php

declare(strict_types=1);

// Port of the extension's tests/optimize.php: Phit\optimize and Phit\uniroot
// against R's optimize() and uniroot(), bit for bit, with the same number of
// function evaluations.

require_once __DIR__.'/../Support/optimize.php';

/** $call with $f counted: the result, the number of calls and the warnings. */
function optimizeCounted(callable $f, callable $call): array
{
    $calls = 0;
    $warnings = [];
    set_error_handler(function (int $no, string $message) use (&$warnings) {
        $warnings[] = $message;

        return true;
    });
    try {
        $result = $call(function (float $x) use ($f, &$calls) {
            $calls++;

            return $f($x);
        });
    } finally {
        restore_error_handler();
    }

    return [$result, $calls, $warnings];
}

$fixture = loadFixture('optimize-r.json');

it('matches R for optimize bit for bit, with the same evaluations', function () use ($fixture) {
    $functions = optimizeFunctions();
    foreach ($fixture['optimize'] as $i => $c) {
        $label = "optimize/$i {$c['f']}";
        $maximum = $c['maximum'] ?? false;
        [$r, $calls, $warnings] = optimizeCounted($functions[$c['f']], fn ($f) => pf('Phit\optimize', $f, $c['lower'], $c['upper'], $c['tol'] ?? null, $maximum));
        $key = $maximum ? 'maximum' : 'minimum';
        expect(array_keys($r))->toBe([$key, 'objective'], $label)
            ->and($r[$key])->toBe((float) $c['expected']['x'], $label)
            ->and($r['objective'])->toBe((float) $c['expected']['objective'], $label)
            ->and($calls)->toBe($c['expected']['evaluations'], $label)
            ->and($warnings)->toHaveCount($c['f'] === 'wall' ? 1 : 0);
    }
});

it('matches R for uniroot bit for bit, with the same evaluations', function () use ($fixture) {
    $functions = optimizeFunctions();
    foreach ($fixture['uniroot'] as $i => $c) {
        $label = "uniroot/$i {$c['f']}";
        [$r, $calls, $warnings] = optimizeCounted($functions[$c['f']], fn ($f) => pf('Phit\uniroot', $f, $c['lower'], $c['upper'], $c['tol'] ?? null, $c['maxiter'] ?? 1000));
        $e = $c['expected'];
        expect(array_keys($r))->toBe(['root', 'f_root', 'iter', 'estim_prec'], $label)
            ->and($r['root'])->toBe((float) $e['root'], $label)
            ->and($r['f_root'])->toBe((float) $e['f_root'], $label)
            ->and($r['iter'])->toBe((int) $e['iter'], $label)
            ->and($r['estim_prec'])->toBe((float) $e['estim_prec'], $label)
            ->and($calls)->toBe($e['evaluations'], $label)
            ->and(array_filter($warnings, fn ($w) => str_contains($w, 'not converged')) === [])->toBe($e['converged'], $label);
    }
});

it('takes any callable and lets its exceptions through', function () {
    expect(pf('Phit\optimize', 'abs', -1, 2)['minimum'])->toBe(pf('Phit\optimize', fn ($x) => abs($x), -1, 2)['minimum'])
        ->and(fn () => pf('Phit\optimize', function () {
            throw new DomainException('stop');
        }, 0, 1))->toThrow(DomainException::class, 'stop')
        ->and(fn () => pf('Phit\uniroot', function (float $x) {
            if ($x > 0.5) {
                throw new DomainException('stop');
            }

            return $x - 0.75;
        }, 0, 1))->toThrow(DomainException::class, 'stop');
});

it('reads int results as numbers and rejects anything else', function () {
    expect(pf('Phit\uniroot', fn ($x) => $x < 2 ? -1 : 1, 0, 4)['f_root'])->toBe(-1.0)
        ->and(thrownCode(fn () => pf('Phit\optimize', fn ($x) => '1', 0, 1)))->toBe(1001)
        ->and(thrownCode(fn () => pf('Phit\optimize', fn ($x) => null, 0, 1)))->toBe(1001)
        ->and(thrownCode(fn () => pf('Phit\optimize', fn ($x) => true, 0, 1)))->toBe(1001);
});

it('warns once that nonfinite values were replaced, with the extension text', function () {
    [, , $warnings] = optimizeCounted(fn ($x) => $x < 0.5 ? NAN : $x, fn ($f) => pf('Phit\optimize', $f, 0, 1));

    expect($warnings)->toBe(['Phit\optimize(): nonfinite function values were replaced by the largest finite value']);
});

it('rejects invalid input with code 1001', function () {
    $sq = optimizeFunctions()['sqrt2'];
    foreach ([
        'optimize empty interval' => fn () => pf('Phit\optimize', $sq, 1, 1),
        'optimize infinite upper' => fn () => pf('Phit\optimize', $sq, 0, INF),
        'optimize zero tol' => fn () => pf('Phit\optimize', $sq, 0, 1, tol: 0.0),
        'optimize NAN tol' => fn () => pf('Phit\optimize', $sq, 0, 1, tol: NAN),
        'uniroot reversed' => fn () => pf('Phit\uniroot', $sq, 2, 0),
        'uniroot same signs' => fn () => pf('Phit\uniroot', $sq, 2, 3),
        'uniroot NAN at lower' => fn () => pf('Phit\uniroot', fn ($x) => $x < 0.5 ? NAN : 1.0, 0, 1),
        'uniroot maxiter' => fn () => pf('Phit\uniroot', $sq, 0, 2, maxiter: 0),
        'uniroot negative tol' => fn () => pf('Phit\uniroot', $sq, 0, 2, tol: -1.0),
    ] as $label => $call) {
        expect(thrownCode($call))->toBe(1001, $label);
    }
});
