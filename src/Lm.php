<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Frontend;

/**
 * Phit\fit_ols: a linear model from a formula, the Phit\fit_ols_matrix core plus coefficient names.
 */
final class Lm
{
    public static function fit(string $formula, mixed $data, string $missing = 'drop', mixed $weights = null, mixed $vcov = 'iid'): array
    {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($missing, 'missing');
        $spec=Vcov::spec($vcov,$data);
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true, additional:$spec['cluster']===null?[]:[$spec['cluster']]);
        Frontend::rejectRandom($design, 'Phit\fit_ols');
        $y = Frontend::numericResponse($design, 'Phit\fit_ols');
        $offset = $design->offset ?? array_fill(0, count($y), 0.0);
        $output = Ols::fitWithData($design->x, $y, $w, $offset);
        $output=Vcov::apply($output,$design,null,$w,$spec);
        $output = Frontend::addNames($output, $design);

        return Frontend::finish($output, $design, 'Phit\fit_ols', 'scaled_svd');
    }
}
