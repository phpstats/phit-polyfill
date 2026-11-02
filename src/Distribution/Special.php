<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/**
 * Numerical pieces cdflib does not provide: log-scale tails beyond the double
 * range and logistic-scale helpers, mirroring phit's stats-core
 * distributions/special.rs. Loader's saddle-point densities and the incomplete
 * gamma and beta continued fractions are in Saddle, re-exported here as
 * special.rs re-exports saddle.rs.
 *
 * @internal
 */
final class Special
{
    /**
     * Probabilities below this are recomputed on the log scale instead of taking
     * ln(p), so log tails stay accurate after p itself underflows.
     */
    public const float TINY = 1e-280;

    public const float LN_SQRT_2PI = 0.9189385332046728;

    /**
     * Replace a probability below TINY by exp(ln p) from a log-scale tail, which
     * stays accurate where a transformed argument has underflowed (the `refine`
     * of continuous.rs and discrete.rs).
     *
     * @param  callable(): float  $lnValue
     */
    public static function refine(float $value, callable $lnValue): float
    {
        if ($value >= self::TINY) {
            return $value;
        }
        $ln = $lnValue();

        return is_nan($ln) ? $value : exp($ln);
    }

    /**
     * ln(lower) from independently computed complementary probabilities.
     * `$extreme` supplies ln(lower) once `lower` is too small for ln(lower).
     *
     * @param  callable(): float  $extreme
     */
    public static function logTail(float $lower, float $upper, callable $extreme): float
    {
        if ($lower > 0.5) {
            return log1p(-$upper);
        }
        if ($lower >= self::TINY) {
            return log($lower);
        }

        return $extreme();
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
            $denominator = $x + fdiv((float) $k, $denominator);
        }

        return -0.5 * $x * $x - self::LN_SQRT_2PI - log($denominator);
    }

    /** ln(1 + e^t) without overflow. */
    public static function softplus(float $t): float
    {
        // Rust's f64::max returns the non-NaN operand, so NaN gives 0.0 here.
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

    /** ln(λ^x e^-λ / Γ(x + 1)); see Saddle::lnPoissonRaw. */
    public static function lnPoissonRaw(float $x, float $lambda): float
    {
        return Saddle::lnPoissonRaw($x, $lambda);
    }

    /** Binomial log probability; see Saddle::lnBinomialRaw. */
    public static function lnBinomialRaw(float $x, float $n, float $p, float $q): float
    {
        return Saddle::lnBinomialRaw($x, $n, $p, $q);
    }

    /** Binomial log probability with m = n - x given; see Saddle::lnBinomialRawSplit. */
    public static function lnBinomialRawSplit(float $x, float $m, float $p, float $q): float
    {
        return Saddle::lnBinomialRawSplit($x, $m, $p, $q);
    }

    /** ln P(a, x); see Saddle::lnGammaLowerSeries. */
    public static function lnGammaLowerSeries(float $a, float $x): float
    {
        return Saddle::lnGammaLowerSeries($a, $x);
    }

    /** ln Q(a, x); see Saddle::lnGammaUpperFraction. */
    public static function lnGammaUpperFraction(float $a, float $x): float
    {
        return Saddle::lnGammaUpperFraction($a, $x);
    }

    /** ln I_x(a, b); see Saddle::lnBetaLowerFraction. */
    public static function lnBetaLowerFraction(float $a, float $b, float $x, float $lnX, float $lnY): float
    {
        return Saddle::lnBetaLowerFraction($a, $b, $x, $lnX, $lnY);
    }
}
