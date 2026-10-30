<?php

declare(strict_types=1);

namespace Phit\Polyfill;

/**
 * Which phit functions this package implements, and where.
 */
final class Polyfill
{
    /**
     * Whether the polyfill implements $fn. Functions it does not implement
     * (Phit\read_csv) need the native extension.
     */
    public static function supports(string $fn): bool
    {
        return isset(FunctionMap::MAP[$fn]);
    }

    /**
     * The implementation behind a polyfilled function. It is used whether or
     * not the extension is loaded, so the two can be compared side by side.
     *
     * @return callable(mixed ...): mixed
     */
    public static function implementation(string $fn): callable
    {
        $target = FunctionMap::MAP[$fn] ?? throw new \InvalidArgumentException("$fn is not provided by the polyfill");

        return $target;
    }

    /** @return list<string> */
    public static function functions(): array
    {
        return array_keys(FunctionMap::MAP);
    }

    /** Whether the native phit extension is loaded. */
    public static function nativeLoaded(): bool
    {
        return \extension_loaded('phit');
    }
}
