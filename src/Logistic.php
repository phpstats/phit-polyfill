<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Internal\Family\LogisticFamily;

/**
 * Binary or grouped-binomial logistic regression by damped Newton IRLS (Phit\fit_logistic_matrix).
 */
final class Logistic
{
    public static function fit(
        array $x,
        array $y,
        string $missing = 'drop',
        ?array $weights = null,
        ?array $offset = null,
        bool $intercept = true,
    ): array {
        return (new LogisticFamily())->run($x, $y, $missing, $weights, $offset, $intercept, 'binary', 'Phit\fit_logistic_matrix');
    }
}
