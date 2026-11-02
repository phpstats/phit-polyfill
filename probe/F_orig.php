<?php

declare(strict_types=1);

namespace Probe;

final class F_orig
{
    private const float EPS = PHP_FLOAT_EPSILON;

    private const int MAX_ITER = 1_000_000;

    public static function orig(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a; $converged = false;
        for ($i = 0; $i < self::MAX_ITER; $i++) {
            $denom += 1.0;
            $term *= $x / $denom;
            $y = $term - $comp;
            $t = $sum + $y;
            $comp = ($t - $sum) - $y;
            $sum = $t;
            if ($term <= $sum * self::EPS * 0.5) { $converged = true; break; }
        }
        return $converged ? $sum : NAN;
    }
}
