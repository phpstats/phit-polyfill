<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Internal\Family\ProbitFamily;

/**
 * Binary or grouped-binomial probit regression by Fisher scoring (Phit\fit_probit_matrix).
 */
final class Probit
{
    public static function fit(
        array $x,
        array $y,
        string $missing = 'drop',
        ?array $weights = null,
        ?array $offset = null,
        bool $intercept = true,
    ): array {
        return (new ProbitFamily())->run($x, $y, $missing, $weights, $offset, $intercept, 'binary', 'Phit\fit_probit_matrix');
    }
}
