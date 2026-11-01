<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Jit;
use Phit\Polyfill\Polyfill;

/*
 * The ported suites call the polyfill implementation directly through pf(),
 * so they test the polyfill whether or not the extension is loaded. The
 * Parity suite compares those implementations with the loaded extension.
 */

/** Decoded JSON fixture from tests/fixtures. */
function loadFixture(string $name): array
{
    static $cache = [];

    return $cache[$name] ??= json_decode(
        file_get_contents(__DIR__.'/fixtures/'.$name),
        true,
        flags: JSON_THROW_ON_ERROR,
    );
}

/** The functions of the API contract, without its classes. */
function apiFunctions(): array
{
    return array_filter(loadFixture('php-api.json'), fn (array $spec) => ! isset($spec['methods']));
}

/** The classes of the API contract, which are all extension-only. */
function apiClasses(): array
{
    return array_filter(loadFixture('php-api.json'), fn (array $spec) => isset($spec['methods']));
}

/** Call the polyfill implementation of a Phit\ function (named arguments allowed). */
function pf(string $fn, mixed ...$arguments): mixed
{
    return Polyfill::implementation($fn)(...$arguments);
}

/** Call the native extension's function; only valid when it is loaded. */
function native(string $fn, mixed ...$arguments): mixed
{
    return $fn(...$arguments);
}

/** |actual - expected| <= abs + rel * |expected|, elementwise, with matching NaN/INF. */
function close(mixed $actual, mixed $expected, float $abs = 1e-9, float $rel = 1e-7, string $label = ''): bool
{
    if (is_array($expected)) {
        if (! is_array($actual) || count($actual) !== count($expected)) {
            return false;
        }
        foreach ($expected as $key => $value) {
            if (! array_key_exists($key, $actual) || ! close($actual[$key], $value, $abs, $rel)) {
                return false;
            }
        }

        return true;
    }
    if (! is_int($actual) && ! is_float($actual)) {
        return $actual === $expected;
    }
    if (is_string($expected)) {
        // JSON fixtures may spell non-finite values as strings.
        $expected = match (strtolower($expected)) {
            'nan' => NAN, 'inf', 'infinity' => INF, '-inf', '-infinity' => -INF,
            default => (float) $expected,
        };
    }
    if ($expected === null) {
        return is_nan((float) $actual);
    }
    $actual = (float) $actual;
    $expected = (float) $expected;
    if (is_nan($expected) || is_nan($actual)) {
        return is_nan($expected) && is_nan($actual);
    }
    if (is_infinite($expected) || is_infinite($actual)) {
        return $actual === $expected;
    }

    return abs($actual - $expected) <= $abs + $rel * abs($expected);
}

/**
 * Record one polyfill-versus-extension comparison when PARITY_LOG names a file.
 * Each suite|function pair keeps the number of calls, numeric values and
 * exceptions compared, the largest absolute and relative differences, and
 * whether every result was identical bit for bit. bin/site-data reads the file.
 */
function parityLog(string $suite, string $fn, mixed $polyfill, mixed $native, bool $threw = false): void
{
    static $log = null;
    $path = getenv('PARITY_LOG');
    if (! $path) {
        return;
    }
    if ($log === null) {
        $log = [];
        register_shutdown_function(function () use (&$log, $path): void {
            // JSON has no infinity; a non-finite difference is written as "inf".
            array_walk_recursive($log, function (&$v) {
                $v = is_float($v) && ! is_finite($v) ? 'inf' : $v;
            });
            file_put_contents($path, json_encode(array_values($log), JSON_PRETTY_PRINT | JSON_PRESERVE_ZERO_FRACTION)."\n");
        });
    }
    $entry = &$log["$suite|$fn"];
    $entry ??= ['suite' => $suite, 'function' => $fn, 'calls' => 0, 'values' => 0, 'exceptions' => 0,
        'max_abs_error' => 0.0, 'max_rel_error' => 0.0, 'identical' => true];
    $entry['calls']++;
    if ($threw) {
        $entry['exceptions']++;
        $entry['identical'] = $entry['identical'] && serialize($polyfill) === serialize($native);

        return;
    }
    $entry['identical'] = $entry['identical'] && serialize($polyfill) === serialize($native);
    $leaves = function (mixed $value): array {
        $out = [];
        is_array($value) ? array_walk_recursive($value, function ($v) use (&$out) { $out[] = $v; }) : $out[] = $value;

        return $out;
    };
    $a = $leaves($polyfill);
    $b = $leaves($native);
    foreach ($b as $i => $y) {
        if (! is_float($y) && ! is_int($y)) {
            continue;
        }
        $entry['values']++;
        $x = $a[$i] ?? NAN;
        if ((is_float($x) && is_nan($x) && is_nan((float) $y)) || $x === $y) {
            continue;
        }
        $difference = abs((float) $x - (float) $y);
        $entry['max_abs_error'] = max($entry['max_abs_error'], is_nan($difference) ? INF : $difference);
        if ((float) $y != 0.0) {
            $entry['max_rel_error'] = max($entry['max_rel_error'], is_nan($difference) ? INF : $difference / abs((float) $y));
        }
    }
}

/** Collect PHP warnings raised while $callback runs. */
function warnings(callable $callback, mixed &$result = null): array
{
    $warnings = [];
    set_error_handler(function (int $severity, string $message) use (&$warnings): bool {
        $warnings[] = [$severity, $message];

        return true;
    });
    try {
        $result = $callback();
    } finally {
        restore_error_handler();
    }

    return $warnings;
}

/**
 * Skip a test that fits Phit\fit_mixed_gaussian_matrix, Phit\fit_mixed_gaussian, Phit\fit_mixed_logistic_matrix
 * or Phit\fit_mixed_logistic while a JIT compiles the process. The tests call the
 * implementations, which do not refuse, and PHP 8.4's tracing JIT has moved
 * their optimizer to a wrong minimum. JitTest checks the Phit\ functions'
 * refusal.
 */
function skipMixedFitsUnderJit(): void
{
    if (Jit::compiling()) {
        test()->markTestSkipped('mixed models do not run under the JIT (docs/notes/function-jit.md)');
    }
}

/** The exception code thrown by $callback, or null. */
function thrownCode(callable $callback): ?int
{
    try {
        warnings($callback);
    } catch (\Throwable $e) {
        return $e->getCode();
    }

    return null;
}

expect()->extend('toBeClose', function (mixed $expected, float $abs = 1e-9, float $rel = 1e-7) {
    $ok = close($this->value, $expected, $abs, $rel);
    expect($ok)->toBeTrue(sprintf(
        "Values differ beyond abs %.1e + rel %.1e:\nactual:   %s\nexpected: %s",
        $abs,
        $rel,
        json_encode($this->value, JSON_PARTIAL_OUTPUT_ON_ERROR | JSON_PRESERVE_ZERO_FRACTION),
        json_encode($expected, JSON_PARTIAL_OUTPUT_ON_ERROR | JSON_PRESERVE_ZERO_FRACTION),
    ));

    return $this;
});

// phit 0.4's deprecated aliases and Phit\fit_ols dispatch (tests/Unit/AliasesTest.php
// and tests/Parity/AliasesParityTest.php), from the extension's tests/aliases.php.

/** The cases: [old name, new name, positional arguments, named arguments]. */
function aliasCases(): array
{
    $d = aliasData();

    return [
        ['fit_lm', 'fit_ols', ['score ~ study_hours + prior_score', $d['schools']], ['weights' => $d['weights']]],
        ['fit_glm', 'fit_generalized_linear', ['pass ~ study_hours + prior_score', $d['schools']], ['family' => 'binomial', 'link' => 'probit']],
        ['fit_lmer', 'fit_mixed_gaussian', ['score ~ study_hours + (1 | school)', $d['schools']], ['reml' => false]],
        ['fit_glmer', 'fit_mixed_logistic', ['cbind(incidence, healthy) ~ period + (1 | herd)', $d['herds']], ['quadraturePoints' => 5]],
        ['fit_logistic', 'fit_logistic_matrix', [$d['x'], $d['pass']], ['offset' => array_fill(0, count($d['pass']), 0.1)]],
        ['fit_probit', 'fit_probit_matrix', [$d['x'], $d['pass']], []],
        ['fit_poisson', 'fit_poisson_matrix', [$d['x'], $d['count']], ['intercept' => true]],
        ['fit_mixed_intercept', 'fit_mixed_gaussian_intercept_matrix', [$d['x'], $d['y'], $d['groups']], []],
        ['fit_mixed', 'fit_mixed_gaussian_matrix', [$d['x'], $d['y'], $d['random']], ['reml' => true]],
        ['fit_logistic_mixed', 'fit_mixed_logistic_matrix', [$d['hx'], $d['hy'], $d['herds']['herd']], ['weights' => $d['herds']['size']]],
    ];
}

/** The fit_ols calls that must fail with code 1001: label => [arguments, message fragment]. */
function fitOlsErrorCases(): array
{
    ['schools' => $schools, 'x' => $x, 'y' => $y] = aliasData();

    return [
        'formula with x' => [['score ~ study_hours', $schools, 'x' => $x], 'a formula and data, or a matrix x and y, not both'],
        'array with x' => [[$x, 'x' => $x, 'y' => $y], 'a formula and data, or a matrix x and y, not both'],
        'data with y' => [['x' => $x, 'data' => $y, 'y' => $y], 'a formula and data, or a matrix x and y, not both'],
        'formula with offset' => [['score ~ study_hours', $schools, 'offset' => $y], 'no offset or intercept arguments'],
        'formula without intercept' => [['score ~ study_hours', $schools, 'intercept' => false], 'no offset or intercept arguments'],
        'formula without data' => [['score ~ study_hours'], 'with a formula needs data'],
        'nothing' => [[], 'needs a formula and data'],
        'x without y' => [[$x], 'needs a formula and data'],
        'number' => [[5, $y], 'needs a formula string, or a matrix x as an array'],
        'matrix with string y' => [[$x, 'score'], 'needs y as an array'],
        'matrix with column-name weights' => [[$x, $y, 'weights' => 'w'], 'needs weights as an array'],
    ];
}

/** The datasets of the extension's alias tests, read as Phit\read_csv reads them. */
function aliasData(): array
{
    static $data = null;
    if ($data !== null) {
        return $data;
    }
    $read = function (string $file, array $types): array {
        $rows = array_map(fn (string $line) => str_getcsv($line, escape: ''), file(__DIR__.'/fixtures/data/'.$file, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES));
        $names = array_shift($rows);
        $columns = [];
        foreach ($names as $j => $name) {
            $columns[$name] = array_map(fn (array $row) => match ($types[$j]) {
                'f' => (float) $row[$j],
                'i' => (int) $row[$j],
                'c' => $row[$j],
            }, $rows);
        }

        return $columns;
    };
    $schools = $read('schools.csv', ['c', 'f', 'f', 'f']);
    $herds = $read('cbpp.csv', ['c', 'i', 'i', 'c']);
    $herds['healthy'] = array_map(fn ($n, $k) => $n - $k, $herds['size'], $herds['incidence']);
    $y = $schools['score'];
    $pass = array_map(fn ($s) => $s >= 70.0 ? 1 : 0, $y);
    $schools['pass'] = $pass;

    return $data = [
        'schools' => $schools,
        'herds' => $herds,
        'x' => array_map(fn ($h, $p) => [$h, $p], $schools['study_hours'], $schools['prior_score']),
        'y' => $y,
        'pass' => $pass,
        'count' => array_map(fn ($s) => (int) floor($s / 10), $y),
        'groups' => $schools['school'],
        'weights' => array_map(fn ($i) => 1.0 + ($i % 3) / 2, array_keys($y)),
        'random' => [['z' => array_fill(0, count($y), [1.0]), 'groups' => $schools['school']]],
        'hx' => array_map(fn ($p) => [$p === '2' ? 1.0 : 0.0, $p === '3' ? 1.0 : 0.0, $p === '4' ? 1.0 : 0.0], $herds['period']),
        'hy' => array_map(fn ($k, $n) => $k / $n, $herds['incidence'], $herds['size']),
    ];
}

/** Runs $call and returns its result and the E_USER_DEPRECATED messages it raised. */
function deprecations(callable $call): array
{
    $messages = [];
    set_error_handler(function (int $level, string $message) use (&$messages): bool {
        if ($level === E_USER_DEPRECATED) {
            $messages[] = $message;

            return true;
        }

        return false;
    });
    try {
        $result = $call();
    } finally {
        restore_error_handler();
    }

    return [$result, $messages];
}
