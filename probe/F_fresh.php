<?php

declare(strict_types=1);

namespace Probe;

final class F_fresh
{
    private const float EPS = PHP_FLOAT_EPSILON;

    private const int MAX_ITER = 1_000_000;

    // B: fresh temporaries for each step of the compensation.
    public static function fresh(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a; $converged = false;
        for ($i = 0; $i < self::MAX_ITER; $i++) {
            $denom += 1.0;
            $ratio = $x / $denom;
            $term *= $ratio;
            $y = $term - $comp;
            $t = $sum + $y;
            $grown = $t - $sum;
            $comp = $grown - $y;
            $sum = $t;
            if ($term <= $sum * self::EPS * 0.5) { $converged = true; break; }
        }
        return $converged ? $sum : NAN;
    }
}
