<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

use Phit\Polyfill\Internal\Statrs;

/**
 * Loader's saddle-point densities and the incomplete gamma and beta series and
 * continued fractions, mirroring phit's stats-core distributions/saddle.rs.
 *
 * Written from the clean-room specification specs/cleanroom/special.md and the
 * clean-room Rust it describes; sources: Loader (2000), Abramowitz and Stegun,
 * Lentz (1976), and Thompson and Barnett (1986).
 *
 * @internal
 */
final class Saddle
{
    /** ln(2π) */
    private const float LN_2PI = 1.8378770664093455;

    private const float EPS = PHP_FLOAT_EPSILON;

    /** Iteration cap for the series and continued fractions. */
    private const int MAX_ITER = 1_000_000;

    /** Lentz's replacement for a zero denominator. */
    private const float LENTZ_TINY = 1e-300;

    /** δ(k/2) for k = 1..=30 (n = 0.5, 1.0, ..., 15.0), from mpmath at 60 digits. */
    private const array DELTA_HALF = [
        0.15342640972002736,   // 0.5
        0.08106146679532726,   // 1.0
        0.05481412105191765,   // 1.5
        0.0413406959554093,    // 2.0
        0.03316287351993629,   // 2.5
        0.02767792568499834,   // 3.0
        0.023746163656297496,  // 3.5
        0.020790672103765093,  // 4.0
        0.018488450532673187,  // 4.5
        0.016644691189821193,  // 5.0
        0.015134973221917378,  // 5.5
        0.013876128823070748,  // 6.0
        0.012810465242920227,  // 6.5
        0.01189670994589177,   // 7.0
        0.011104559758206917,  // 7.5
        0.010411265261972096,  // 8.0
        0.009799416126158804,  // 8.5
        0.009255462182712733,  // 9.0
        0.008768700134139386,  // 9.5
        0.00833056343336287,   // 10.0
        0.00793411456431402,   // 10.5
        0.007573675487951841,  // 11.0
        0.007244554301320383,  // 11.5
        0.00694284010720953,   // 12.0
        0.006665247032707682,  // 12.5
        0.006408994188004207,  // 13.0
        0.006171712263039458,  // 13.5
        0.0059513701127588475, // 14.0
        0.0057462165130101155, // 14.5
        0.005554733551962801,  // 15.0
    ];

    /** Coefficients B_{2k} / (2k (2k - 1)) of the asymptotic series δ(n) ~ Σ c_k n^{1-2k}, k = 1..=9. */
    private const array STIRLING_COEFS = [
        1.0 / 12.0,
        -1.0 / 360.0,
        1.0 / 1260.0,
        -1.0 / 1680.0,
        1.0 / 1188.0,
        -691.0 / 360_360.0,
        1.0 / 156.0,
        -3617.0 / 122_400.0,
        43_867.0 / 244_188.0,
    ];

    /** Asymptotic series for δ(n), n >= 15. */
    private static function stirlingAsymptotic(float $n): float
    {
        $r = 1.0 / $n;
        $r2 = $r * $r;
        $acc = 0.0;
        for ($i = count(self::STIRLING_COEFS) - 1; $i >= 0; $i--) {
            $acc = $acc * $r2 + self::STIRLING_COEFS[$i];
        }

        return $acc * $r;
    }

    /**
     * δ(n) - δ(n + 1) = (n + 1/2) ln(1 + 1/n) - 1 = atanh(t)/t - 1 with t = 1/(2n + 1),
     * = t²/3 + t⁴/5 + t⁶/7 + ..., a sum of positive terms with no cancellation.
     */
    private static function stirlingStep(float $n): float
    {
        $t = 1.0 / (2.0 * $n + 1.0);
        if ($t > 0.5) {
            // Only for n < 1/2: atanh(t)/t >= 1.098, so the subtraction loses at most one digit.
            $atanh = 0.5 * log1p(fdiv(2.0 * $t, 1.0 - $t));

            return $atanh / $t - 1.0;
        }
        $t2 = $t * $t;
        $pow = $t2;
        $sum = 0.0;
        $k = 3.0;
        while (true) {
            $next = $sum + $pow / $k;
            if ($next == $sum) {
                break;
            }
            $sum = $next;
            $pow *= $t2;
            $k += 2.0;
        }

        return $sum;
    }

    /** δ(n) = ln Γ(n + 1) - (n + 1/2) ln n + n - ln √(2π), for n > 0. */
    public static function stirlingError(float $n): float
    {
        if (is_nan($n) || $n < 0.0) {
            return NAN;
        }
        if ($n == 0.0) {
            return INF;
        }
        if ($n >= 15.0) {
            return self::stirlingAsymptotic($n);
        }
        $twice = 2.0 * $n;
        if ($twice == floor($twice)) {
            return self::DELTA_HALF[(int) $twice - 1];
        }
        // Climb with the step relation until n + k >= 15, then use the asymptotic series.
        $m = $n;
        $steps = 0.0;
        while ($m < 15.0) {
            $steps += self::stirlingStep($m);
            $m += 1.0;
        }

        return self::stirlingAsymptotic($m) + $steps;
    }

    /**
     * δ(t) + ½ ln t, regrouped for t < 1 as ln Γ(t + 1) - t ln t + t - ln √(2π)
     * so that the two cancelling terms do not appear.
     */
    private static function stirlingPlusHalfLog(float $t): float
    {
        if ($t < 1.0) {
            return Statrs::lnGamma($t + 1.0) - $t * log($t) + $t - 0.5 * self::LN_2PI;
        }

        return self::stirlingError($t) + 0.5 * log($t);
    }

    /** x ln(x / np) + np - x for x >= 0 and np > 0, accurate when x ≈ np. */
    public static function devianceTerm(float $x, float $np): float
    {
        if ($x == 0.0) {
            return $np;
        }
        $diff = $x - $np;
        $total = $x + $np;
        if (abs($diff) < 0.1 * $total) {
            // With v = (x - np)/(x + np): D = (x - np) v + 2x (v³/3 + v⁵/5 + ...).
            $v = $diff / $total;
            $v2 = $v * $v;
            $pow = 2.0 * $x * $v;
            $sum = $diff * $v;
            $k = 3.0;
            while (true) {
                $pow *= $v2;
                $next = $sum + $pow / $k;
                if ($next == $sum) {
                    break;
                }
                $sum = $next;
                $k += 2.0;
            }

            return $sum;
        }
        $ratio = fdiv($x, $np);
        $logRatio = is_finite($ratio) && $ratio > 0.0 ? log($ratio) : log($x) - log($np);

        return $x * $logRatio + $np - $x;
    }

    /** ln(λ^x e^-λ / Γ(x + 1)) for real x >= 0 and λ >= 0, by the saddle-point form. */
    public static function lnPoissonRaw(float $x, float $lambda): float
    {
        if ($lambda == 0.0) {
            return $x == 0.0 ? 0.0 : -INF;
        }
        if ($x == 0.0) {
            return -$lambda;
        }

        // -δ(x) - D(x, λ) - ½ ln(2πx), with δ(x) grouped with its ½ ln x.
        return -self::stirlingPlusHalfLog($x) - self::devianceTerm($x, $lambda) - 0.5 * self::LN_2PI;
    }

    /** n ln q, computed from whichever of p and q = 1 - p is the more accurate argument. */
    private static function nLnComplement(float $n, float $p, float $q): float
    {
        if ($p < 0.5) {
            return $n * log1p(-$p);
        }

        return $n * log($q);
    }

    /** Shared body: `x` successes and `m` failures out of `n` (n = x + m). */
    private static function lnBinomialCore(float $x, float $m, float $n, float $p, float $q): float
    {
        if ($p == 0.0) {
            return $x == 0.0 ? 0.0 : -INF;
        }
        if ($q == 0.0) {
            return $m == 0.0 ? 0.0 : -INF;
        }
        if ($x == 0.0) {
            return self::nLnComplement($n, $p, $q);
        }
        if ($m == 0.0) {
            return self::nLnComplement($n, $q, $p);
        }
        // δ(n) - δ(x) - δ(m) - D(x, np) - D(m, nq) - ½ ln(2π x m / n), each δ(t) grouped with its ½ ln t.
        $stirling = self::stirlingPlusHalfLog($n) - self::stirlingPlusHalfLog($x) - self::stirlingPlusHalfLog($m);
        $deviance = self::devianceTerm($x, $n * $p) + self::devianceTerm($m, $n * $q);

        return $stirling - $deviance - 0.5 * self::LN_2PI;
    }

    /** ln(Γ(n + 1) / (Γ(x + 1) Γ(n - x + 1)) p^x q^(n-x)) for real 0 <= x <= n, with q = 1 - p given accurately. */
    public static function lnBinomialRaw(float $x, float $n, float $p, float $q): float
    {
        return self::lnBinomialCore($x, $n - $x, $n, $p, $q);
    }

    /** As lnBinomialRaw with n = x + m, for callers that know m = n - x exactly. */
    public static function lnBinomialRawSplit(float $x, float $m, float $p, float $q): float
    {
        return self::lnBinomialCore($x, $m, $x + $m, $p, $q);
    }

    /** ln(x^a e^-x / Γ(a + 1)): the Poisson saddle-point form for a >= 1, the direct sum otherwise. */
    private static function lnGammaPrefactor(float $a, float $x): float
    {
        if ($a >= 1.0) {
            return self::lnPoissonRaw($a, $x);
        }

        return $a * log($x) - $x - Statrs::lnGamma($a + 1.0);
    }

    /**
     * ln P(a, x), the regularized lower incomplete gamma function, by its power
     * series (A&S 6.5.29); used for 0 < x < a + 1.
     */
    public static function lnGammaLowerSeries(float $a, float $x): float
    {
        if ($x == 0.0) {
            return -INF;
        }
        $term = 1.0;
        $sum = 1.0;
        $comp = 0.0; // compensated (Kahan) summation
        $denom = $a;
        $converged = false;
        for ($i = 0; $i < self::MAX_ITER; $i++) {
            $denom += 1.0;
            $term *= $x / $denom;
            $y = $term - $comp;
            $t = $sum + $y;
            $comp = ($t - $sum) - $y;
            $sum = $t;
            if ($term <= $sum * self::EPS * 0.5) {
                $converged = true;
                break;
            }
        }
        if ($converged === false) {
            return NAN;
        }

        return self::lnGammaPrefactor($a, $x) + log($sum);
    }

    /**
     * Evaluates b0 + a1/(b1 + a2/(b2 + ...)) with the modified Lentz method.
     * `$terms(j)` returns [a_j, b_j] for j >= 1. Returns NaN at the iteration cap.
     *
     * @param  callable(int): array{0: float, 1: float}  $terms
     */
    private static function lentz(float $b0, callable $terms): float
    {
        $f = $b0 == 0.0 ? self::LENTZ_TINY : $b0;
        $c = $f;
        $d = 0.0;
        for ($j = 1; $j <= self::MAX_ITER; $j++) {
            [$aj, $bj] = $terms($j);
            $d = $bj + $aj * $d;
            if ($d == 0.0) {
                $d = self::LENTZ_TINY;
            }
            $c = $bj + fdiv($aj, $c);
            if ($c == 0.0) {
                $c = self::LENTZ_TINY;
            }
            $d = fdiv(1.0, $d);
            $delta = $c * $d;
            $f *= $delta;
            if (abs($delta - 1.0) <= self::EPS) {
                return $f;
            }
        }

        return NAN;
    }

    /**
     * ln Q(a, x), the regularized upper incomplete gamma function, by the
     * contracted continued fraction of A&S 6.5.31 (modified Lentz); used for x > a + 1.
     */
    public static function lnGammaUpperFraction(float $a, float $x): float
    {
        $b1 = $x + 1.0 - $a;
        $denom = self::lentz($b1, static function (int $j) use ($a, $b1): array {
            $k = (float) $j;

            return [-$k * ($k - $a), $b1 + 2.0 * $k];
        });
        if (is_nan($denom)) {
            return NAN;
        }

        return log($a) + self::lnGammaPrefactor($a, $x) - log($denom);
    }

    /**
     * ln I_x(a, b), the regularized incomplete beta function, by its continued
     * fraction (A&S 26.5.8, modified Lentz); used for 0 < x < (a + 1)/(a + b + 2).
     * The caller passes ln x and ln(1 - x), computed without cancellation.
     */
    public static function lnBetaLowerFraction(float $a, float $b, float $x, float $lnX, float $lnY): float
    {
        $denom = self::lentz(1.0, static function (int $j) use ($a, $b, $x): array {
            if ($j % 2 === 1) {
                $m = (float) intdiv($j - 1, 2);
                $d = fdiv(-($a + $m) * ($a + $b + $m) * $x, ($a + 2.0 * $m) * ($a + 2.0 * $m + 1.0));
            } else {
                $m = (float) intdiv($j, 2);
                $d = fdiv($m * ($b - $m) * $x, ($a + 2.0 * $m - 1.0) * ($a + 2.0 * $m));
            }

            return [$d, 1.0];
        });
        if (is_nan($denom)) {
            return NAN;
        }
        if ($x < PHP_FLOAT_MIN) {
            // x is 0 or subnormal and no longer carries the value that ln x does: use the direct log form.
            $prefactor = $a * $lnX + $b * $lnY - self::lnABeta($a, $b);
        } else {
            $y = exp($lnY);
            $prefactor = self::lnBinomialRawSplit($a, $b, $x, $y) + log($b / ($a + $b));
        }

        return $prefactor - log($denom);
    }

    /** ln(a B(a, b)) = ln Γ(a + 1) + ln Γ(b) - ln Γ(a + b). */
    private static function lnABeta(float $a, float $b): float
    {
        return Statrs::lnGamma($a + 1.0) + Statrs::lnGamma($b) - Statrs::lnGamma($a + $b);
    }
}
