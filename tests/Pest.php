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
 * Skip a test that fits Phit\fit_mixed, Phit\fit_lmer, Phit\fit_logistic_mixed
 * or Phit\fit_glmer while a JIT compiles the process. The tests call the
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
