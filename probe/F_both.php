<?php

declare(strict_types=1);

namespace Probe;

final class F_both
{
    private const float EPS = PHP_FLOAT_EPSILON;

    private const int MAX_ITER = 1_000_000;

    // D: B and C together.
    public static function both(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a;
        $halfEps = self::EPS * 0.5; $cap = self::MAX_ITER;
        for ($i = 0; $i < $cap; $i++) {
            $denom += 1.0;
            $ratio = $x / $denom;
            $term *= $ratio;
            $y = $term - $comp;
            $t = $sum + $y;
            $grown = $t - $sum;
            $comp = $grown - $y;
            $sum = $t;
            if ($term <= $sum * $halfEps) { return $sum; }
        }
        return NAN;
    }
}
