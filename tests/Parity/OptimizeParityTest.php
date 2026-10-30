<?php

declare(strict_types=1);

// Compares Phit\optimize and Phit\uniroot with the loaded extension; skipped
// without it. The two must evaluate f at the same points in the same order
// and return the same bits, warnings (whose level differs: the polyfill can
// only raise E_USER_WARNING) and exceptions.

/** Every point f was evaluated at, the outcome and the warning texts. */
function optimizeTrace(callable $f, callable $call): array
{
    $points = [];
    $warnings = [];
    set_error_handler(function (int $no, string $message) use (&$warnings) {
        $warnings[] = $message;

        return true;
    });
    try {
        $outcome = ['ok', $call(function (float $x) use ($f, &$points) {
            $points[] = $x;

            return $f($x);
        })];
    } catch (\Throwable $e) {
        $outcome = ['error', get_class($e), $e->getCode(), $e->getMessage()];
    } finally {
        restore_error_handler();
    }

    return [$outcome, $points, $warnings];
}

function optimizeCompare(string $fn, callable $f, array $args): void
{
    $native = optimizeTrace($f, fn ($g) => native($fn, $g, ...$args));
    $polyfill = optimizeTrace($f, fn ($g) => pf($fn, $g, ...$args));
    parityLog('optimize', $fn, $polyfill[0], $native[0], $native[0][0] === 'error' || $polyfill[0][0] === 'error');
    if (serialize($native) !== serialize($polyfill)) {
        throw new \RuntimeException("$fn differs for ".json_encode($args)."\nnative ".var_export($native, true)."\npolyfill ".var_export($polyfill, true));
    }
    expect(true)->toBeTrue();
}

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

it('agrees on the R fixtures', function () {
    $fixture = loadFixture('optimize-r.json');
    $functions = optimizeFunctions();
    foreach ($fixture['optimize'] as $c) {
        optimizeCompare('Phit\optimize', $functions[$c['f']], [$c['lower'], $c['upper'], $c['tol'] ?? null, $c['maximum'] ?? false]);
    }
    foreach ($fixture['uniroot'] as $c) {
        optimizeCompare('Phit\uniroot', $functions[$c['f']], [$c['lower'], $c['upper'], $c['tol'] ?? null, $c['maxiter'] ?? 1000]);
    }
});

it('agrees on seeded random functions and intervals', function () {
    mt_srand(crc32('optimize') ^ (int) (getenv('PARITY_SEED') ?: 20261030));
    $u = fn (float $lo, float $hi): float => $lo + ($hi - $lo) * mt_rand() / mt_getrandmax();
    for ($it = 0, $iterations = (int) (getenv('PARITY_ITERATIONS') ?: 60) * 4; $it < $iterations; $it++) {
        // A random cubic with a sine term, sometimes nonfinite or integer-valued
        // on part of the interval.
        [$a, $b, $c, $d, $s] = [$u(-3, 3), $u(-3, 3), $u(-3, 3), $u(-1, 1), $u(0, 2)];
        $cut = $u(-5, 5);
        $bad = [null, INF, -INF, NAN][mt_rand(0, 3)];
        $integer = mt_rand(0, 9) === 0;
        $f = function (float $x) use ($a, $b, $c, $d, $s, $cut, $bad, $integer) {
            if ($bad !== null && $x < $cut) {
                return $bad;
            }
            $y = (($a * $x + $b) * $x + $c) * $x + $d + $s * sin(3 * $x);

            return $integer ? (int) round($y) : $y;
        };
        $lower = $u(-5, 4);
        $upper = mt_rand(0, 19) === 0 ? $lower : $u($lower, 5);
        $tol = [null, null, 1e-10, $u(1e-6, 0.1), 0.0][mt_rand(0, 4)];
        optimizeCompare('Phit\optimize', $f, [$lower, $upper, $tol, mt_rand(0, 1) === 1]);
        optimizeCompare('Phit\uniroot', $f, [$lower, $upper, $tol, [1, 2, 5, 1000][mt_rand(0, 3)]]);
    }
});

it('rejects the same invalid input and lets the same exceptions through', function () {
    $sq = fn (float $x) => $x * $x - 2;
    foreach ([
        ['Phit\optimize', $sq, [1, 1]], ['Phit\optimize', $sq, [0, INF]], ['Phit\optimize', $sq, [NAN, 1]],
        ['Phit\optimize', $sq, [0, 1, 0.0]], ['Phit\optimize', $sq, [0, 1, INF]],
        ['Phit\uniroot', $sq, [2, 0]], ['Phit\uniroot', $sq, [2, 3]], ['Phit\uniroot', $sq, [0, 2, null, 0]],
        ['Phit\uniroot', $sq, [0, 2, -1.0]], ['Phit\uniroot', fn ($x) => $x < 0.5 ? NAN : 1.0, [0, 1]],
        ['Phit\uniroot', fn ($x) => $x < 0.5 ? -INF : INF, [0, 1]], ['Phit\uniroot', fn ($x) => 0, [0, 1]],
        ['Phit\optimize', fn ($x) => '1', [0, 1]], ['Phit\optimize', fn ($x) => null, [0, 1]], ['Phit\optimize', fn ($x) => true, [0, 1]],
        ['Phit\uniroot', fn ($x) => [1], [0, 1]], ['Phit\uniroot', fn ($x) => $x < 0.5 ? 'bad' : 1.0, [0, 1]],
        ['Phit\uniroot', fn ($x) => $x < 0.5 ? 'bad' : throw new DomainException('upper'), [0, 1]],
        ['Phit\uniroot', fn ($x) => $x < 0.5 ? throw new DomainException('lower') : 'bad', [0, 1]],
        ['Phit\uniroot', fn ($x) => $x < 0.5 ? -1.0 : 'bad', [0, 1]],
        ['Phit\optimize', function () {
            throw new DomainException('stop', 7);
        }, [0, 1]],
        ['Phit\uniroot', function (float $x) {
            if ($x > 0.5) {
                throw new DomainException('stop', 7);
            }

            return $x - 0.75;
        }, [0, 1]],
    ] as [$fn, $f, $args]) {
        optimizeCompare($fn, $f, $args);
    }
});
