<?php

declare(strict_types=1);

namespace Phit\Polyfill\Objects;

/**
 * Builds the fit and table objects the Phit functions return, and reads the
 * array behind a Fit|array argument.
 *
 * @internal
 */
final class Wrap
{
    /** The fit classes, by the result's `model` key. */
    public const array FIT_CLASSES = [
        'linear' => 'Phit\LinearFit',
        'generalized_linear' => 'Phit\GeneralizedLinearFit',
        'mixed_gaussian' => 'Phit\MixedGaussianFit',
        'mixed_logistic' => 'Phit\MixedLogisticFit',
    ];

    /** The fit object for a model function's result, by its `model` key. */
    public static function fit(array $data): object
    {
        $class = self::FIT_CLASSES[$data['model'] ?? null] ?? throw new \LogicException('model functions set a known model key');

        return self::holding($class, $data);
    }

    public static function table(array $columns): object
    {
        return self::holding('Phit\Table', $columns);
    }

    /** The result array behind a fit argument: a Phit fit object, or an array in the same shape. */
    public static function data(mixed $fit): ?array
    {
        if (is_array($fit)) {
            return $fit;
        }
        foreach (self::FIT_CLASSES as $class) {
            if ($fit instanceof $class) {
                return $fit->toArray();
            }
        }

        return null;
    }

    private static function holding(string $class, array $data): object
    {
        // The classes have no usable constructor, as the extension's cannot be built from PHP.
        $object = (new \ReflectionClass($class))->newInstanceWithoutConstructor();
        (function (array $data): void {
            $this->data = $data;
        })->call($object, $data);

        return $object;
    }
}
