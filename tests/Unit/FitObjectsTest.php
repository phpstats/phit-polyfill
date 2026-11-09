<?php

declare(strict_types=1);

// The extension's tests/fit-objects.php: every model function returns the
// class its model key selects, reads like its array, is read-only, prints its
// summary, and agrees with the free functions; tidy() returns a Phit\Table.
// These call the Phit functions, so they test the polyfill's classes without
// the extension and check the same expectations against it when it is loaded.

use Phit\Polyfill\Internal\Jit;

/** Whether $call throws an instance of $class, with $message when given. */
function fitObjectThrows(callable $call, string $class, ?string $message = null): bool
{
    try {
        $call();
    } catch (Throwable $e) {
        return $e instanceof $class && ($message === null || $e->getMessage() === $message);
    }

    return false;
}

beforeEach(function () {
    if (Jit::compiling()) {
        $this->markTestSkipped('the Phit functions refuse to run under the JIT; JitTest checks that');
    }
});

describe('fit objects', function () {
    it('behave as their arrays for every model function', function (array $case) {
        if (str_contains($case['function'], 'mixed')) {
            skipMixedFitsUnderJit();
        }
        $fixture = loadFixture('summary/summary.json');
        $classes = ['iv' => 'Phit\IvFit', 'linear' => 'Phit\LinearFit', 'generalized_linear' => 'Phit\GeneralizedLinearFit',
            'mixed_gaussian' => 'Phit\MixedGaussianFit', 'mixed_logistic' => 'Phit\MixedLogisticFit'];
        $readOnly = 'A Phit fit is read-only; call toArray() for an array copy';
        $args = $case['args'];
        if (isset($args['data'])) {
            $args['data'] = $fixture['data'][$args['data']];
        }
        $fit = ('Phit\\'.$case['function'])(...$args);
        $array = $fit->toArray();

        expect(get_class($fit))->toBe($classes[$array['model']])
            ->and((new ReflectionClass($fit))->isFinal())->toBeTrue();
        foreach (['Phit\Fit', ArrayAccess::class, Iterator::class, Countable::class, JsonSerializable::class, Stringable::class] as $interface) {
            expect($fit)->toBeInstanceOf($interface);
        }

        // Printing and the methods agree with the free functions, which also take the array.
        expect((string) $fit)->toBe($case['summary'])
            ->and($fit->summary())->toBe($case['summary'])
            ->and(Phit\summary($array))->toBe(Phit\summary($fit))
            ->and((string) $fit->tidy(confInt: true))->toBe((string) Phit\tidy($fit, confInt: true))
            ->and(Phit\tidy($array, confInt: true)->toArray())->toBe($fit->tidy(confInt: true)->toArray());

        // Reading works as on the array: keys, nested values, isset, count, iteration, and JSON.
        foreach ($array as $key => $value) {
            expect($fit[$key])->toBe($value)
                ->and(isset($fit[$key]))->toBe(isset($array[$key]));
        }
        expect(isset($fit['no_such_key']))->toBeFalse()
            ->and($fit['no_such_key'] ?? 'default')->toBe('default')
            ->and(count($fit))->toBe(count($array))
            ->and(iterator_to_array($fit))->toBe($array);
        $keys = [];
        foreach ($fit as $key => $value) {
            $keys[] = $key;
        }
        expect($keys)->toBe(array_map('strval', array_keys($array)))
            ->and(json_encode($fit))->toBe(json_encode($array))
            ->and($fit->toMimeBundle())->toBe(['text/plain' => $case['summary'], 'application/json' => $array])
            ->and((clone $fit)->toArray())->toBe($array);

        // Writes are refused.
        expect(fitObjectThrows(function () use ($fit) {
            $fit['coefficients'] = [];
        }, LogicException::class, $readOnly))->toBeTrue()
            ->and(fitObjectThrows(function () use ($fit) {
                unset($fit['coefficients']);
            }, LogicException::class, $readOnly))->toBeTrue();

        // Formula fits predict through the method as through the function.
        if (! str_contains($case['function'], 'matrix')) {
            expect($fit->predict($args['data']))->toBe(Phit\predict($fit, $args['data']))
                ->and(Phit\predict($array, $args['data']))->toBe(Phit\predict($fit, $args['data']));
        }
    })->with(fn () => array_map(fn (array $case) => [$case], array_column(loadFixture('summary/summary.json')['cases'], null, 'name')));

    it('warns on a missing key as an array does', function () {
        $fit = Phit\fit_ols('mpg ~ wt', loadFixture('summary/summary.json')['data']['mtcars']);
        $warnings = [];
        set_error_handler(function (int $level, string $message) use (&$warnings): bool {
            // The polyfill raises E_USER_WARNING where the extension raises E_WARNING.
            $warnings[] = [$level === E_USER_WARNING ? E_WARNING : $level, $message];

            return true;
        });
        try {
            $missing = $fit['no_such_key'];
            $fit[1.5];
            $fit[7];
        } finally {
            restore_error_handler();
        }

        expect($missing)->toBeNull()
            ->and($warnings)->toBe([
                [E_WARNING, 'Phit\LinearFit::offsetGet(): Undefined array key "no_such_key"'],
                [E_WARNING, 'Phit\LinearFit::offsetGet(): Undefined array key ""'],
                [E_WARNING, 'Phit\LinearFit::offsetGet(): Undefined array key "7"'],
            ]);
    });

    it('cannot be constructed', function () {
        foreach (['Phit\LinearFit', 'Phit\GeneralizedLinearFit', 'Phit\MixedGaussianFit', 'Phit\MixedLogisticFit', 'Phit\Table'] as $class) {
            expect(fitObjectThrows(fn () => new $class(), Exception::class, 'You cannot instantiate this class from PHP.'))->toBeTrue($class);
        }
    });

    it('rejects other values in the free functions', function () {
        $message = 'fit must be a result of a Phit model function: not a fit object or array';

        expect(fitObjectThrows(fn () => Phit\summary(new stdClass()), Exception::class, $message))->toBeTrue()
            ->and(fitObjectThrows(fn () => Phit\tidy('fit'), Exception::class, $message))->toBeTrue()
            ->and(fitObjectThrows(fn () => Phit\predict(5, []), Exception::class,
                'fit must be a result of Phit\fit_ols, Phit\fit_generalized_linear, Phit\fit_mixed_gaussian or Phit\fit_mixed_logistic: not a fit object or array'))->toBeTrue();
    });
});

describe('table objects', function () {
    it('hold read-only columns', function () {
        $fit = Phit\fit_ols('mpg ~ wt', loadFixture('summary/summary.json')['data']['mtcars']);
        $table = $fit->tidy(confInt: true);
        $columns = $table->toArray();

        expect($table)->toBeInstanceOf('Phit\Table')
            ->and((new ReflectionClass($table))->isFinal())->toBeTrue()
            ->and($table->toColumns())->toBe($columns)
            ->and(array_keys($columns))->toBe(['term', 'estimate', 'std_error', 'statistic', 'p_value', 'conf_low', 'conf_high'])
            ->and(count($table))->toBe(2)
            ->and($table['term'])->toBe(['(Intercept)', 'wt'])
            ->and(iterator_to_array($table))->toBe($columns)
            ->and(json_encode($table))->toBe(json_encode($columns))
            ->and($table->toMimeBundle())->toBe(['text/plain' => (string) $table, 'application/json' => $columns])
            ->and(fitObjectThrows(function () use ($table) {
                $table['term'] = [];
            }, LogicException::class, 'A Phit table is read-only; call toArray() for an array copy'))->toBeTrue();
    });
});
