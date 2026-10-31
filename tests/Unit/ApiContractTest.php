<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Jit;
use Phit\Polyfill\Polyfill;

// The Phit functions the polyfill defines must have the extension's exact
// signatures (tests/fixtures/php-api.json), so code runs unchanged on either.

$api = apiFunctions();

it('declares every polyfilled function with the extension signature', function () use ($api) {
    foreach (Polyfill::functions() as $fn) {
        expect($api)->toHaveKey($fn);
        $spec = $api[$fn];
        $reflection = extension_loaded('phit')
            ? new ReflectionMethod(...Polyfill::implementation($fn))
            : new ReflectionFunction($fn);
        $parameters = $reflection->getParameters();
        expect($parameters)->toHaveCount(count($spec['parameters']), $fn);
        foreach ($spec['parameters'] as $i => $expected) {
            $actual = $parameters[$i];
            expect($actual->getName())->toBe($expected['name'], "$fn parameter $i")
                ->and((string) $actual->getType())->toBe($expected['type'], "$fn \${$expected['name']} type")
                ->and($actual->isOptional())->toBe($expected['optional'], "$fn \${$expected['name']} optional");
            if ($expected['optional']) {
                expect($actual->getDefaultValue())->toEqual($expected['default'], "$fn \${$expected['name']} default");
            }
        }
        expect((string) $reflection->getReturnType())->toBe($spec['return_type'], "$fn return type");
    }
});

it('defines the Phit functions only when the extension is absent', function () {
    foreach (Polyfill::functions() as $fn) {
        expect(function_exists($fn))->toBeTrue($fn);
        $function = new ReflectionFunction($fn);
        expect($function->isInternal())->toBe(extension_loaded('phit'), $fn);
    }
});

it('leaves extension-only functions undefined', function () use ($api) {
    $extensionOnly = array_diff(array_keys($api), Polyfill::functions());
    foreach ($extensionOnly as $fn) {
        expect(Polyfill::supports($fn))->toBeFalse($fn);
        if (! extension_loaded('phit')) {
            expect(function_exists($fn))->toBeFalse($fn);
        }
    }
    foreach (['Phit\read_csv'] as $fn) {
        expect(Polyfill::supports($fn))->toBeFalse($fn);
    }
});

it('leaves the contract classes to the extension', function () {
    $classes = apiClasses();
    expect($classes)->toHaveKey('Phit\CsvReader');
    foreach ($classes as $class => $spec) {
        expect(Polyfill::supports($class))->toBeFalse($class);
        if (! extension_loaded('phit')) {
            expect(class_exists($class))->toBeFalse($class);

            continue;
        }
        foreach ($spec['methods'] as $name => $method) {
            $reflection = new ReflectionMethod($class, $name);
            $parameters = array_map(fn (ReflectionParameter $p) => $p->getName(), $reflection->getParameters());
            expect($parameters)->toBe(array_column($method['parameters'], 'name'), "$class::$name");
        }
    }
});

it('defines no global stats_* aliases', function () {
    $global = array_filter(get_defined_functions()['user'], fn ($fn) => str_starts_with($fn, 'stats_'));
    expect($global)->toBe([])
        ->and(array_filter(Polyfill::functions(), fn ($fn) => ! str_starts_with($fn, 'Phit\\')))->toBe([]);
});

it('routes the Phit functions to the polyfill implementation', function () {
    if (extension_loaded('phit')) {
        $this->markTestSkipped('the extension provides the Phit functions');
    }
    if (Jit::compiling()) {
        $this->markTestSkipped('the Phit functions refuse to run under the JIT; JitTest checks that');
    }
    $x = [[1.0], [2.0], [3.0], [4.0], [5.0]];
    $y = [2.1, 4.2, 5.8, 8.1, 10.2];
    expect(\Phit\fit_ols($x, $y))->toBe(Polyfill::implementation('Phit\fit_ols')($x, $y))
        ->and(\Phit\stat_mean([1, 2, 3]))->toBe(2.0);
});
