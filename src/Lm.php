<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Frontend;

/**
 * Linear model from a formula: the Phit\fit_ols_matrix core plus coefficient names (Phit\fit_ols).
 */
final class Lm
{
    public static function fit(string $formula, mixed $data, string $missing = 'drop', mixed $weights = null): array
    {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($missing, 'missing');
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        Frontend::rejectRandom($design, 'Phit\fit_ols');
        $y = Frontend::numericResponse($design, 'Phit\fit_ols');
        $offset = $design->offset ?? array_fill(0, count($y), 0.0);
        $output = Ols::fitWithData($design->x, $y, $w, $offset);
        $output = Frontend::addNames($output, $design);

        return Frontend::finish($output, $design, 'Phit\fit_ols', 'scaled_svd');
    }
}
