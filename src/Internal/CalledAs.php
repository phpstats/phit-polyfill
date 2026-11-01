<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

/**
 * The name a model function was called by, when it differs from the function
 * that does the work: a deprecated alias such as Phit\fit_lm, or Phit\fit_ols
 * with a matrix. PHP prefixes the extension's warnings with that name, so the
 * listwise-deletion warning uses it.
 */
final class CalledAs
{
    private static ?string $name = null;

    /**
     * @template T
     *
     * @param  callable(): T  $fit
     * @return T
     */
    public static function run(string $name, callable $fit): mixed
    {
        $previous = self::$name;
        self::$name = $name;
        try {
            return $fit();
        } finally {
            self::$name = $previous;
        }
    }

    public static function name(): ?string
    {
        return self::$name;
    }
}
