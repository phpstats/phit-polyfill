<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Internal\Family\PoissonFamily;

/**
 * Poisson count regression with exposure offsets by damped Newton IRLS (Phit\fit_poisson_matrix).
 */
final class Poisson
{
    public static function fit(
        array $x,
        array $y,
        string $missing = 'drop',
        ?array $weights = null,
        ?array $offset = null,
        bool $intercept = true,
    ): array {
        return (new PoissonFamily())->run($x, $y, $missing, $weights, $offset, $intercept, 'count', 'Phit\fit_poisson_matrix');
    }
}
