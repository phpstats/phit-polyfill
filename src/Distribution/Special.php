<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/**
 * Follows stats-core's distributions/special.rs: log-scale tails beyond the
 * double range, Loader's saddle-point densities and logistic-scale helpers.
 *
 * @internal
 */
final class Special
{
    /** Probabilities below this are recomputed on the log scale. */
    public const float TINY = 1e-280;

    public const float LN_SQRT_2PI = 0.918_938_533_204_672_8;

    private const int MAX_TERMS = 100_000;

    private const float FPMIN = 1e-300;

    /** Stirling-series error at n = k/2, k = 1..30 (mpmath, 50 digits). */
    private const array STIRLING_HALVES = [
        1.534_264_097_200_273_5e-1,
        8.106_146_679_532_726e-2,
        5.481_412_105_191_765_4e-2,
        4.134_069_595_540_929_4e-2,
        3.316_287_351_993_628_7e-2,
        2.767_792_568_499_834e-2,
        2.374_616_365_629_749_6e-2,
        2.079_067_210_376_509_3e-2,
        1.848_845_053_267_318_5e-2,
        1.664_469_118_982_119_2e-2,
        1.513_497_322_191_737_9e-2,
        1.387_612_882_307_074_8e-2,
        1.281_046_524_292_022_7e-2,
        1.189_670_994_589_177e-2,
        1.110_455_975_820_691_7e-2,
        1.041_126_526_197_209_6e-2,
        9.799_416_126_158_803e-3,
        9.255_462_182_712_733e-3,
        8.768_700_134_139_385e-3,
        8.330_563_433_362_871e-3,
        7.934_114_564_314_02e-3,
        7.573_675_487_951_841e-3,
        7.244_554_301_320_383e-3,
        6.942_840_107_209_53e-3,
        6.665_247_032_707_682e-3,
        6.408_994_188_004_207e-3,
        6.171_712_263_039_458e-3,
        5.951_370_112_758_848e-3,
        5.746_216_513_010_116e-3,
        5.554_733_551_962_801e-3,
    ];

    /**
     * ln(lower) from independently computed complementary probabilities;
     * $extreme supplies ln(lower) once lower is below TINY.
     */
    public static function logTail(float $lower, float $upper, \Closure $extreme): float
    {
        if ($lower > 0.5) {
            return log1p(-$upper);
        }
        if ($lower >= self::TINY) {
            return log($lower);
        }

        return $extreme();
    }

    /** Replace a probability below TINY by exp of its log-scale value. */
    public static function refine(float $value, \Closure $lnValue): float
    {
        if ($value >= self::TINY) {
            return $value;
        }
        $ln = $lnValue();

        return is_nan($ln) ? $value : exp($ln);
    }

    /** ln Φ(z) for z < -8 from the continued fraction for the Mills ratio. */
    public static function normalLogLowerTail(float $z): float
    {
        $x = -$z;
        if (is_infinite($x)) {
            return -INF;
        }
        $denominator = $x;
        for ($k = 100; $k >= 1; $k--) {
            $denominator = $x + $k / $denominator;
        }

        return -0.5 * $x * $x - self::LN_SQRT_2PI - log($denominator);
    }

    /** ln P(a, x) by the power series; valid for x < a + 1. */
    public static function lnGammaLowerSeries(float $a, float $x): float
    {
        if ($x == 0.0) {
            return -INF;
        }
        if ($x >= $a + 1.0) {
            return NAN;
        }
        $term = 1.0 / $a;
        $sum = $term;
        for ($n = 1; $n < self::MAX_TERMS; $n++) {
            $term *= $x / ($a + $n);
            $sum += $term;
            if ($term < $sum * PHP_FLOAT_EPSILON) {
                return $a * log($x) - $x - Cdflib::gammaLog($a) + log($sum);
            }
        }

        return NAN;
    }

    /** ln Q(a, x) by Lentz's continued fraction; valid for x > a + 1. */
    public static function lnGammaUpperFraction(float $a, float $x): float
    {
        if (is_infinite($x)) {
            return -INF;
        }
        if ($x <= $a + 1.0) {
            return NAN;
        }
        $b = $x + 1.0 - $a;
        $c = 1.0 / self::FPMIN;
        $d = 1.0 / $b;
        $h = $d;
        for ($i = 1; $i < self::MAX_TERMS; $i++) {
            $fi = (float) $i;
            $an = -$fi * ($fi - $a);
            $b += 2.0;
            $d = $an * $d + $b;
            if (abs($d) < self::FPMIN) {
                $d = self::FPMIN;
            }
            $c = $b + $an / $c;
            if (abs($c) < self::FPMIN) {
                $c = self::FPMIN;
            }
            $d = 1.0 / $d;
            $delta = $d * $c;
            $h *= $delta;
            if (abs($delta - 1.0) < PHP_FLOAT_EPSILON) {
                return -$x + $a * log($x) - Cdflib::gammaLog($a) + log($h);
            }
        }

        return NAN;
    }

    /**
     * ln I_x(a, b) by the continued fraction; valid for x < (a + 1)/(a + b + 2).
     * Callers pass ln x and ln(1 - x) computed without cancellation.
     */
    public static function lnBetaLowerFraction(float $a, float $b, float $x, float $lnX, float $lnY): float
    {
        if ($lnX == -INF) {
            return -INF;
        }
        if ($x >= ($a + 1.0) / ($a + $b + 2.0)) {
            return NAN;
        }
        $qab = $a + $b;
        $qap = $a + 1.0;
        $qam = $a - 1.0;
        $c = 1.0;
        $d = 1.0 - $qab * $x / $qap;
        if (abs($d) < self::FPMIN) {
            $d = self::FPMIN;
        }
        $d = 1.0 / $d;
        $h = $d;
        for ($i = 1; $i < self::MAX_TERMS; $i++) {
            $m = (float) $i;
            $m2 = 2.0 * $m;
            $aa = $m * ($b - $m) * $x / (($qam + $m2) * ($a + $m2));
            $d = 1.0 + $aa * $d;
            if (abs($d) < self::FPMIN) {
                $d = self::FPMIN;
            }
            $c = 1.0 + $aa / $c;
            if (abs($c) < self::FPMIN) {
                $c = self::FPMIN;
            }
            $d = 1.0 / $d;
            $h *= $d * $c;
            $aa = -($a + $m) * ($qab + $m) * $x / (($a + $m2) * ($qap + $m2));
            $d = 1.0 + $aa * $d;
            if (abs($d) < self::FPMIN) {
                $d = self::FPMIN;
            }
            $c = 1.0 + $aa / $c;
            if (abs($c) < self::FPMIN) {
                $c = self::FPMIN;
            }
            $d = 1.0 / $d;
            $delta = $d * $c;
            $h *= $delta;
            if (abs($delta - 1.0) < PHP_FLOAT_EPSILON) {
                return $a * $lnX + $b * $lnY - Cdflib::betaLog($a, $b) - log($a) + log($h);
            }
        }

        return NAN;
    }

    /** Stirling-series error δ(n) = ln Γ(n + 1) - (n + 1/2) ln n + n - ln √(2π). */
    public static function stirlingError(float $n): float
    {
        if ($n <= 15.0) {
            $twice = $n + $n;
            if ($twice == floor($twice) && $twice >= 1.0) {
                return self::STIRLING_HALVES[(int) $twice - 1];
            }

            return Cdflib::gammaLog($n + 1.0) - ($n + 0.5) * log($n) + $n - self::LN_SQRT_2PI;
        }
        $S0 = 1.0 / 12.0;
        $S1 = 1.0 / 360.0;
        $S2 = 1.0 / 1260.0;
        $S3 = 1.0 / 1680.0;
        $S4 = 1.0 / 1188.0;
        $nn = $n * $n;
        if ($n > 500.0) {
            return ($S0 - $S1 / $nn) / $n;
        }
        if ($n > 80.0) {
            return ($S0 - ($S1 - $S2 / $nn) / $nn) / $n;
        }
        if ($n > 35.0) {
            return ($S0 - ($S1 - ($S2 - $S3 / $nn) / $nn) / $nn) / $n;
        }

        return ($S0 - ($S1 - ($S2 - ($S3 - $S4 / $nn) / $nn) / $nn) / $nn) / $n;
    }

    /** Deviance term x ln(x / np) + np - x without cancellation when x ≈ np. */
    public static function devianceTerm(float $x, float $np): float
    {
        if (abs($x - $np) < 0.1 * ($x + $np)) {
            $v = ($x - $np) / ($x + $np);
            $s = ($x - $np) * $v;
            if (abs($s) < PHP_FLOAT_MIN) {
                return $s;
            }
            $ej = 2.0 * $x * $v;
            $v *= $v;
            for ($j = 1; $j < 1000; $j++) {
                $ej *= $v;
                $next = $s + $ej / (float) (2 * $j + 1);
                if ($next == $s) {
                    return $next;
                }
                $s = $next;
            }
        }

        return $x * log(fdiv($x, $np)) + $np - $x;
    }

    /** ln of λ^x e^-λ / Γ(x + 1) for real x ≥ 0 (Loader's saddle-point form). */
    public static function lnPoissonRaw(float $x, float $lambda): float
    {
        if ($lambda == 0.0) {
            return $x == 0.0 ? 0.0 : -INF;
        }
        if (! is_finite($lambda) || $x < 0.0) {
            return -INF;
        }
        if ($x == 0.0) {
            return -$lambda;
        }
        if (is_infinite($x)) {
            return -INF;
        }

        return -self::stirlingError($x) - self::devianceTerm($x, $lambda) - 0.5 * log(2.0 * M_PI * $x);
    }

    /** ln of the binomial density at real 0 ≤ x ≤ n, with q = 1 - p supplied. */
    public static function lnBinomialRaw(float $x, float $n, float $p, float $q): float
    {
        if ($x < 0.0 || $x > $n) {
            return -INF;
        }

        return self::lnBinomialRawSplit($x, $n - $x, $p, $q);
    }

    /** As lnBinomialRaw with n = x + m. */
    public static function lnBinomialRawSplit(float $x, float $m, float $p, float $q): float
    {
        $n = $x + $m;
        if ($p == 0.0) {
            return $x == 0.0 ? 0.0 : -INF;
        }
        if ($q == 0.0) {
            return $m == 0.0 ? 0.0 : -INF;
        }
        if ($x == 0.0) {
            if ($n == 0.0) {
                return 0.0;
            }

            return $p < 0.1 ? -self::devianceTerm($n, $n * $q) - $n * $p : $n * log($q);
        }
        if ($m == 0.0) {
            return $q < 0.1 ? -self::devianceTerm($n, $n * $p) - $n * $q : $n * log($p);
        }
        $lc = self::stirlingError($n)
            - self::stirlingError($x)
            - self::stirlingError($m)
            - self::devianceTerm($x, $n * $p)
            - self::devianceTerm($m, $n * $q);
        $lf = log(2.0 * M_PI) + log($x) + log($m / $n);

        return $lc - 0.5 * $lf;
    }

    /** ln(1 + e^t) without overflow. */
    public static function softplus(float $t): float
    {
        return ($t > 0.0 ? $t : 0.0) + log1p(exp(-abs($t)));
    }

    /** ln(1 - e^-t) for t ≥ 0 (Mächler's log1mexp). */
    public static function lnOneMinusExpNeg(float $t): float
    {
        if ($t <= M_LN2) {
            return log(-expm1(-$t));
        }

        return log1p(-exp(-$t));
    }
}
