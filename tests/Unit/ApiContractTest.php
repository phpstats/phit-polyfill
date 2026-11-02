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
            ? new ReflectionMethod(...Phit\Polyfill\FunctionMap::MAP[$fn])
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
        // The bootstrap wraps a model function's array in a fit object and tidy()'s columns in a table.
        $wrapped = str_starts_with($fn, 'Phit\\fit_') || $fn === 'Phit\\tidy';
        $return = extension_loaded('phit') && $wrapped ? 'array' : $spec['return_type'];
        expect((string) $reflection->getReturnType())->toBe($return, "$fn return type");
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

it('leaves Phit\\CsvReader to the extension', function () {
    $classes = apiClasses();
    expect($classes)->toHaveKey('Phit\CsvReader')
        ->and(Polyfill::supports('Phit\CsvReader'))->toBeFalse();
    if (! extension_loaded('phit')) {
        expect(class_exists('Phit\CsvReader'))->toBeFalse();
    }
});

it('defines the fit and table classes with the extension\'s shape', function () {
    $classes = array_filter(apiClasses(), fn (string $class) => $class !== 'Phit\CsvReader', ARRAY_FILTER_USE_KEY);
    expect(array_keys($classes))->toEqualCanonicalizing(['Phit\Fit', 'Phit\LinearFit', 'Phit\GeneralizedLinearFit',
        'Phit\MixedGaussianFit', 'Phit\MixedLogisticFit', 'Phit\Table']);
    foreach ($classes as $class => $spec) {
        $reflection = new ReflectionClass($class);
        expect($reflection->isInternal())->toBe(extension_loaded('phit'), $class);
        if (isset($spec['interface'])) {
            expect($reflection->isInterface())->toBeTrue($class)
                ->and($reflection->getInterfaceNames())->toEqualCanonicalizing($spec['extends'], $class)
                ->and($reflection->getMethods())->toHaveCount(count(array_unique(array_merge(
                    ...array_map(fn ($i) => array_map(fn ($m) => $m->name, (new ReflectionClass($i))->getMethods()), $spec['extends'])
                ))), "$class declares no methods of its own");

            continue;
        }
        expect($reflection->isFinal())->toBeTrue($class);
        $methods = array_map(fn (ReflectionMethod $m) => $m->name, $reflection->getMethods(ReflectionMethod::IS_PUBLIC));
        expect($methods)->toEqualCanonicalizing(array_keys($spec['methods']), $class);
        foreach ($spec['methods'] as $name => $method) {
            $m = new ReflectionMethod($class, $name);
            $parameters = $m->getParameters();
            expect($parameters)->toHaveCount(count($method['parameters']), "$class::$name");
            foreach ($method['parameters'] as $i => $expected) {
                expect($parameters[$i]->getName())->toBe($expected['name'], "$class::$name parameter $i")
                    ->and((string) $parameters[$i]->getType())->toBe($expected['type'], "$class::$name \${$expected['name']} type")
                    ->and($parameters[$i]->isOptional())->toBe($expected['optional'], "$class::$name \${$expected['name']} optional");
                if ($expected['optional']) {
                    expect($parameters[$i]->getDefaultValue())->toEqual($expected['default'], "$class::$name \${$expected['name']} default");
                }
            }
            expect((string) $m->getReturnType())->toBe($method['return_type'], "$class::$name return type");
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
    expect(\Phit\fit_ols_matrix($x, $y)->toArray())->toBe(Polyfill::implementation('Phit\fit_ols_matrix')($x, $y))
        ->and(\Phit\stat_mean([1, 2, 3]))->toBe(2.0);
});
