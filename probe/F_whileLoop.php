<?php

declare(strict_types=1);

namespace Probe;

final class F_whileLoop
{
    private const float EPS = PHP_FLOAT_EPSILON;

    private const int MAX_ITER = 1_000_000;

    // E: a while loop over the state, with no counter variable reused as a float.
    public static function whileLoop(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a; $left = self::MAX_ITER;
        while ($left-- > 0) {
            $denom += 1.0;
            $term *= $x / $denom;
            $y = $term - $comp;
            $t = $sum + $y;
            $comp = ($t - $sum) - $y;
            $sum = $t;
            if ($term <= $sum * self::EPS * 0.5) { return $sum; }
        }
        return NAN;
    }
}
