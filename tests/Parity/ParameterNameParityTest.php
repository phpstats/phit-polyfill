<?php

declare(strict_types=1);

use Phit\Polyfill\FunctionMap;
use Phit\Polyfill\Polyfill;

/*
 * phit and phit-polyfill must expose identical parameter names, so a caller's
 * named arguments work against either implementation.
 *
 * tests/fixtures/php-api.json is the contract. phit generates it from the
 * built extension with scripts/php-stubs.php, and bin/generate-bootstrap
 * builds bootstrap.php from this copy. These tests read the names out of it
 * and check both implementations, so renaming a parameter in one package
 * without the other fails the suite:
 *
 * - the polyfill's backing methods are checked on every run;
 * - the native extension is checked whenever it is loaded, which is what
 *   catches phit drifting away from this package's copy of the contract.
 *   Run `composer test:ext` to include it.
 */

/** @return array<string, list<string>> Contract parameter names, per function, in order. */
function contractParameterNames(): array
{
    static $names = null;
    if ($names === null) {
        $names = [];
        foreach (apiFunctions() as $fn => $spec) {
            $names[$fn] = array_map(fn (array $parameter) => $parameter['name'], $spec['parameters']);
        }
    }

    return $names;
}

/** A value of the declared type. The call only has to bind its names, not succeed. */
function contractArgumentValue(string $type): mixed
{
    return match (ltrim($type, '?')) {
        'int' => 1,
        'float' => 1.0,
        'bool' => false,
        'string' => '',
        'array' => [],
        default => null,
    };
}

/** Every contract parameter as a named argument, optional ones at their declared default. */
function contractNamedArguments(array $parameters): array
{
    $arguments = [];
    foreach ($parameters as $parameter) {
        $arguments[$parameter['name']] = array_key_exists('default', $parameter)
            ? $parameter['default']
            : contractArgumentValue($parameter['type']);
    }

    return $arguments;
}

/**
 * Call $callable with named arguments and report only a name mismatch.
 * An invalid-value failure means the names bound, which is all this checks.
 */
function unknownNamedParameter(callable $callable, array $arguments): ?string
{
    set_error_handler(static fn () => true);
    try {
        $callable(...$arguments);
    } catch (Error $e) {
        if (str_contains($e->getMessage(), 'Unknown named parameter')) {
            return $e->getMessage();
        }
    } catch (Throwable) {
        // A statistical or input error: the parameter names were accepted.
    } finally {
        restore_error_handler();
    }

    return null;
}

it('declares the contract parameter names on every polyfilled implementation', function () {
    $problems = [];
    foreach (contractParameterNames() as $fn => $expected) {
        if (! Polyfill::supports($fn)) {
            continue;
        }
        [$class, $method] = FunctionMap::MAP[$fn];
        $actual = array_map(
            fn (ReflectionParameter $p) => $p->getName(),
            (new ReflectionMethod($class, $method))->getParameters(),
        );
        foreach ($expected as $i => $name) {
            if (($actual[$i] ?? null) !== $name) {
                $problems[] = "$fn parameter $i: contract \$$name, {$class}::{$method}() \$".($actual[$i] ?? 'missing');
            }
        }
    }

    expect($problems)->toBe([]);
});

it('binds every contract parameter name on the polyfill', function () {
    $problems = [];
    foreach (apiFunctions() as $fn => $spec) {
        if (! Polyfill::supports($fn)) {
            continue;
        }
        $message = unknownNamedParameter(Polyfill::implementation($fn), contractNamedArguments($spec['parameters']));
        if ($message !== null) {
            $problems[] = "$fn: $message";
        }
    }

    expect($problems)->toBe([]);
});

it('declares the contract parameter names on the native extension', function () {
    if (! extension_loaded('phit')) {
        test()->markTestSkipped('phit extension not loaded');
    }
    $problems = [];
    foreach (contractParameterNames() as $fn => $expected) {
        if (! function_exists($fn)) {
            $problems[] = "$fn: the extension does not define it";

            continue;
        }
        $actual = array_map(
            fn (ReflectionParameter $p) => $p->getName(),
            (new ReflectionFunction($fn))->getParameters(),
        );
        if ($actual !== $expected) {
            $problems[] = "$fn: contract (".implode(', ', $expected).'), extension ('.implode(', ', $actual).')';
        }
    }

    expect($problems)->toBe([]);
});

it('binds every contract parameter name on the native extension', function () {
    if (! extension_loaded('phit')) {
        test()->markTestSkipped('phit extension not loaded');
    }
    $problems = [];
    foreach (apiFunctions() as $fn => $spec) {
        $message = unknownNamedParameter($fn, contractNamedArguments($spec['parameters']));
        if ($message !== null) {
            $problems[] = "$fn: $message";
        }
    }

    expect($problems)->toBe([]);
});

it('agrees on the parameter names the two implementations share', function () {
    if (! extension_loaded('phit')) {
        test()->markTestSkipped('phit extension not loaded');
    }
    $problems = [];
    foreach (contractParameterNames() as $fn => $contract) {
        if (! Polyfill::supports($fn) || ! function_exists($fn)) {
            continue;
        }
        [$class, $method] = FunctionMap::MAP[$fn];
        $native = array_map(
            fn (ReflectionParameter $p) => $p->getName(),
            (new ReflectionFunction($fn))->getParameters(),
        );
        $polyfill = array_map(
            fn (ReflectionParameter $p) => $p->getName(),
            (new ReflectionMethod($class, $method))->getParameters(),
        );
        foreach ($contract as $i => $name) {
            if (($native[$i] ?? null) !== ($polyfill[$i] ?? null)) {
                $problems[] = "$fn parameter $i (\$$name): extension \$".($native[$i] ?? 'missing')
                    .", polyfill \$".($polyfill[$i] ?? 'missing');
            }
        }
    }

    expect($problems)->toBe([]);
});
