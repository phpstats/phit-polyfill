<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/**
 * Follows stats-core's distributions/quantile.rs: safeguarded Newton on the
 * log scale for continuous quantiles and a bracketing search for discrete ones.
 *
 * @internal
 */
final class Quantile
{
    /**
     * Solve F(x) = p on the support (lower, upper), starting near $start.
     *
     * @param  \Closure(float): float  $logCdf
     * @param  \Closure(float): float  $logSf
     * @param  \Closure(float): float  $logPdf
     */
    public static function continuous(
        float $p,
        float $start,
        float $lower,
        float $upper,
        \Closure $logCdf,
        \Closure $logSf,
        \Closure $logPdf,
    ): float {
        if ($p == 0.0) {
            return $lower;
        }
        if ($p == 1.0) {
            return $upper;
        }
        $useLower = $p <= 0.5;
        $target = $useLower ? log($p) : log(1.0 - $p);
        $g = $useLower
            ? static fn (float $x): float => $logCdf($x) - $target
            : static fn (float $x): float => $target - $logSf($x);
        $tail = $useLower ? $logCdf : $logSf;
        $slope = static fn (float $x): float => exp($logPdf($x) - $tail($x));

        if (is_finite($start) && $start > $lower && $start < $upper) {
            $x = $start;
        } elseif (is_finite($lower) && is_finite($upper)) {
            $x = 0.5 * ($lower + $upper);
        } elseif (is_finite($lower)) {
            $x = $lower + 1.0;
        } elseif (is_finite($upper)) {
            $x = $upper - 1.0;
        } else {
            $x = 0.0;
        }
        $gx = $g($x);
        if ($gx == 0.0) {
            return $x;
        }
        // A failed tail evaluation (NaN) must surface as a failure, not steer the search.
        if (is_nan($gx)) {
            return NAN;
        }
        $a = $lower;
        $b = $upper;
        if ($gx < 0.0) {
            $a = $x;
            if (is_infinite($upper)) {
                // The upstream walk stops at the first point where g is NaN or
                // reaches the target; a NaN stop is a failure.
                // NaN is tested explicitly: opcache's optimizer may rewrite !(a < b) as a >= b.
                $hit = self::bracket($x, 1.0, static function (float $v) use ($g): bool {
                    $gv = $g($v);

                    return is_nan($gv) || $gv >= 0.0;
                });
                if ($hit === null) {
                    return $upper;
                }
                if (is_nan($g($hit[0]))) {
                    return NAN;
                }
                [$b, $a] = $hit;
            }
        } else {
            $b = $x;
            if (is_infinite($lower)) {
                $hit = self::bracket($x, -1.0, static function (float $v) use ($g): bool {
                    $gv = $g($v);

                    return is_nan($gv) || $gv <= 0.0;
                });
                if ($hit === null) {
                    return $lower;
                }
                if (is_nan($g($hit[0]))) {
                    return NAN;
                }
                [$a, $b] = $hit;
            } elseif ($lower == 0.0 && $g(PHP_FLOAT_MIN) > 0.0) {
                return 0.0;
            }
        }

        for ($i = 0; $i < 400; $i++) {
            $step = fdiv($gx, $slope($x));
            $next = $x - $step;
            // Written without a negation: opcache folds !(x > a) into x <= a, which is wrong for NaN.
            if (is_finite($next) && $next > $a && $next < $b) {
                // Keep the Newton step.
            } else {
                $next = self::midpoint($a, $b);
            }
            if ($next == $x || abs($next - $x) <= 4.0 * PHP_FLOAT_EPSILON * abs($next)) {
                return $next;
            }
            $x = $next;
            $gx = $g($x);
            if ($gx == 0.0) {
                return $x;
            }
            if (is_nan($gx)) {
                return NAN;
            }
            if ($gx < 0.0) {
                $a = $x;
            } else {
                $b = $x;
            }
            $absA = abs($a);
            $absB = abs($b);
            if ($b - $a <= 4.0 * PHP_FLOAT_EPSILON * ($absA > $absB ? $absA : $absB)) {
                return 0.5 * ($a + $b);
            }
        }

        return NAN;
    }

    /**
     * Upstream's bracket expansion: step outward from $x with expand() until
     * $stop holds, returning [the first point where it holds, the point before
     * it] (the start counts as "before"), or null once expand() reaches ±∞.
     *
     * While the walk is additive (steps of ±2 between |x| = 2 and 2^54, where
     * x ± 2k is exact) the first stopping step is located by exponential and
     * binary search instead of one step at a time: the extension can take
     * ~10^8 steps there (Student t with tiny df). The points are the same grid,
     * so the bracket is identical whenever $stop is monotone along it.
     *
     * @param  \Closure(float): bool  $stop
     * @return array{0: float, 1: float}|null
     */
    private static function bracket(float $x, float $direction, \Closure $stop): ?array
    {
        $limit = 18_014_398_509_481_984.0; // 2^54
        $next = $x;
        $previous = $x;
        while (true) {
            $distance = $direction < 0.0 ? $next - 2.0 : -2.0 - $next;
            if ($distance >= 0.0 && abs($next) < $limit) {
                // Steps k = 1..K keep the preceding point at |v| >= 2, so each is exact.
                $steps = floor($distance / 2.0) + 1.0;
                $origin = $next;
                $at = static fn (float $k): float => $origin + $direction * (2.0 * $k);
                $miss = 0.0;
                $probe = 1.0;
                while ($probe <= $steps && ! $stop($at($probe))) {
                    $miss = $probe;
                    $probe *= 2.0;
                }
                if ($probe > $steps) {
                    $probe = ($miss < $steps && $stop($at($steps))) ? $steps : INF;
                }
                if ($probe !== INF) {
                    while ($probe - $miss > 1.0) {
                        $mid = floor($miss + ($probe - $miss) / 2.0);
                        if ($stop($at($mid))) {
                            $probe = $mid;
                        } else {
                            $miss = $mid;
                        }
                    }

                    return [$at($probe), $probe == 1.0 ? $previous : $at($probe - 1.0)];
                }
                $next = $at($steps);
                $previous = $next;

                continue;
            }
            $step = self::expand($next, $direction);
            if ($step == $direction * INF) {
                return null;
            }
            if ($step == $next) {
                // No progress (|x| >= 2^54 in the additive phase): upstream never
                // terminates here; report a failed evaluation instead.
                return [NAN, NAN];
            }
            if ($stop($step)) {
                return [$step, $next];
            }
            $next = $step;
            $previous = $step;
        }
    }

    private static function expand(float $x, float $direction): float
    {
        if ($x * $direction <= 2.0) {
            return $x + 2.0 * $direction;
        }
        if (abs($x) == PHP_FLOAT_MAX) {
            return $direction * INF;
        }
        $next = $x * abs($x);

        return is_finite($next) ? $next : $direction * PHP_FLOAT_MAX;
    }

    private static function midpoint(float $a, float $b): float
    {
        $aFloor = $a > PHP_FLOAT_MIN ? $a : PHP_FLOAT_MIN;
        if ($a >= 0.0 && $b > 4.0 * $aFloor) {
            return sqrt($aFloor) * sqrt($b);
        }
        $bCeil = $b < -PHP_FLOAT_MIN ? $b : -PHP_FLOAT_MIN;
        if ($b <= 0.0 && $a < 4.0 * $bCeil) {
            $absB = abs($b);

            return -(sqrt(abs($a)) * sqrt($absB > PHP_FLOAT_MIN ? $absB : PHP_FLOAT_MIN));
        }

        return $a + 0.5 * ($b - $a);
    }

    /**
     * Smallest integer k in [0, max] with F(k) ≥ p, searching from $start.
     * p is reduced by 64 machine epsilons (relative), as in R.
     */
    public static function discrete(
        float $p,
        float $start,
        float $max,
        \Closure $cdf,
        \Closure $sf,
        \Closure $logCdf,
    ): float {
        if ($p == 0.0) {
            return 0.0;
        }
        if ($p == 1.0) {
            return $max;
        }
        $fuzz = 64.0 * PHP_FLOAT_EPSILON;
        // A NaN tail compares as "not reached"; remember it and fail at the end.
        $failed = false;
        if ($p > 0.5) {
            $bound = (1.0 - $p) * (1.0 + $fuzz);
            $reached = static function (float $k) use ($sf, $bound, &$failed): bool {
                $v = $sf($k);
                $failed = $failed || is_nan($v);

                return $v <= $bound;
            };
        } elseif ($p < Special::TINY) {
            $bound = log($p) + log1p(-$fuzz);
            $reached = static function (float $k) use ($logCdf, $bound, &$failed): bool {
                $v = $logCdf($k);
                $failed = $failed || is_nan($v);

                return $v >= $bound;
            };
        } else {
            $bound = $p * (1.0 - $fuzz);
            $reached = static function (float $k) use ($cdf, $bound, &$failed): bool {
                $v = $cdf($k);
                $failed = $failed || is_nan($v);

                return $v >= $bound;
            };
        }
        $limit = is_finite($max) ? $max : 9.007_199_254_740_992e15;
        if (is_finite($start)) {
            $k = floor($start);
            $k = $k < 0.0 ? 0.0 : ($k > $limit ? $limit : $k);
        } else {
            $k = 0.0;
        }

        if ($reached($k)) {
            $hi = $k;
            $step = 1.0;
            while (true) {
                $candidate = $hi - $step;
                if ($candidate < 0.0) {
                    $lo = -1.0;
                    break;
                }
                if ($reached($candidate)) {
                    $hi = $candidate;
                    $step *= 2.0;
                } else {
                    $lo = $candidate;
                    break;
                }
            }
        } else {
            $lo = $k;
            $step = 1.0;
            while (true) {
                $candidate = $lo + $step;
                if ($candidate > $limit) {
                    $candidate = $limit;
                }
                if ($reached($candidate)) {
                    $hi = $candidate;
                    break;
                }
                if ($candidate >= $limit) {
                    if ($failed) {
                        return NAN;
                    }

                    return is_finite($max) ? $max : INF;
                }
                $lo = $candidate;
                $step *= 2.0;
            }
        }
        while ($hi - $lo > 1.0) {
            $k = floor($lo + 0.5 * ($hi - $lo));
            if ($reached($k)) {
                $hi = $k;
            } else {
                $lo = $k;
            }
        }

        return $failed ? NAN : $hi;
    }

    /** Cornish–Fisher starting value for the discrete quantile search. */
    public static function cornishFisher(float $p, float $mean, float $sd, float $skewness): float
    {
        $z = Normal::stdQuantile($p);
        // f64::clamp(-40, 40)
        if ($z < -40.0) {
            $z = -40.0;
        }
        if ($z > 40.0) {
            $z = 40.0;
        }

        return $mean + $sd * ($z + $skewness * ($z * $z - 1.0) / 6.0);
    }

    /** A nonnegative finite integer (k.fract() == 0). */
    public static function isCount(float $k): bool
    {
        return $k >= 0.0 && is_finite($k) && $k == floor($k);
    }
}
