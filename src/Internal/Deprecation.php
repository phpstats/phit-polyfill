<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

/**
 * The E_USER_DEPRECATED notices of phit 0.4's renamed functions. The extension
 * raises them with php_error, which prefixes the function name, so the
 * polyfill writes the prefix itself.
 */
final class Deprecation
{
    /** A deprecated alias, for example Phit\fit_lm for Phit\fit_ols. */
    public static function alias(string $alias, string $replacement): void
    {
        trigger_error("$alias(): Deprecated since phit 0.4.0; use $replacement() instead", E_USER_DEPRECATED);
    }

    /** Phit\fit_ols called with a matrix x and y. */
    public static function olsMatrix(): void
    {
        trigger_error('Phit\fit_ols(): Passing a matrix is deprecated since phit 0.4.0; use Phit\fit_ols_matrix() instead', E_USER_DEPRECATED);
    }
}
