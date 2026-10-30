<?php

declare(strict_types=1);

namespace Phit\Polyfill;

/**
 * Descriptive statistics, ported from the phit extension (stats-core
 * descriptive.rs plus the PHP glue). null and NAN are missing values;
 * missing: 'raise' rejects them, 'drop' removes them (pairwise for two arrays).
 * Infinite values are always errors. Sums are Neumaier-compensated and second
 * moments use the corrected two-pass algorithm.
 */
final class Descriptive
{
    /** @param array<int|float|null> $values */
    public static function mean(array $values, string $missing = 'raise'): float
    {
        $v = self::values($values, $missing);
        self::check($v, 1, 'mean');

        return self::computed(self::plainMean($v));
    }

    /** @param array<int|float|null> $values */
    public static function var(array $values, string $type = 'sample', string $missing = 'raise'): float
    {
        $v = self::values($values, $missing);

        return self::variance($v, self::population($type));
    }

    /** @param array<int|float|null> $values */
    public static function sd(array $values, string $type = 'sample', string $missing = 'raise'): float
    {
        $v = self::values($values, $missing);

        return sqrt(self::variance($v, self::population($type)));
    }

    /** @param array<int|float|null> $values */
    public static function median(array $values, string $missing = 'raise'): float
    {
        return self::quantiles(self::values($values, $missing), [0.5])[0];
    }

    /** @param array<int|float|null> $values */
    public static function quantile(array $values, float $p, string $missing = 'raise'): float
    {
        return self::quantiles(self::values($values, $missing), [$p])[0];
    }

    /**
     * @param array<int|float|null> $values
     * @param array<int|float|null> $probabilities
     * @return list<float>
     */
    public static function quantileArray(array $values, array $probabilities, string $missing = 'raise'): array
    {
        $p = [];
        foreach ($probabilities as $probability) {
            if (! \is_int($probability) && ! \is_float($probability)) {
                throw self::error(ErrorCode::INVALID_INPUT, 'probabilities must be numbers');
            }
            $p[] = (float) $probability;
        }

        return self::quantiles(self::values($values, $missing), $p);
    }

    /**
     * @param array<int|float|null> $x
     * @param array<int|float|null> $y
     */
    public static function cov(array $x, array $y, string $type = 'sample', string $missing = 'raise'): float
    {
        [$a, $b] = self::pairs($x, $y, $missing);
        $population = self::population($type);
        self::checkPair($a, $b, 'covariance');
        $n = (float) \count($a);

        return self::computed(self::crossProducts($a, $b)[0] / ($population ? $n : $n - 1.0));
    }

    /**
     * @param array<int|float|null> $x
     * @param array<int|float|null> $y
     */
    public static function cor(array $x, array $y, string $missing = 'raise'): float
    {
        [$a, $b] = self::pairs($x, $y, $missing);
        self::checkPair($a, $b, 'correlation');
        [$sxy, $sxx, $syy] = self::crossProducts($a, $b);
        if ($sxx == 0.0 || $syy == 0.0) {
            return NAN;
        }
        $r = self::computed(($sxy / sqrt($sxx)) / sqrt($syy));

        return min(1.0, max(-1.0, $r));
    }

    /**
     * @param array<int|float|null> $values
     * @param array<int|float|null> $weights
     */
    public static function weightedMean(array $values, array $weights, string $missing = 'raise'): float
    {
        [$v, $w] = self::weightedPairs($values, $weights, $missing);
        $total = self::checkWeights($v, $w);

        return self::computed(self::plainWeightedMean($v, $w, $total));
    }

    /**
     * @param array<int|float|null> $values
     * @param array<int|float|null> $weights
     */
    public static function weightedVar(array $values, array $weights, string $type = 'frequency', string $missing = 'raise'): float
    {
        if (! \in_array($type, ['frequency', 'reliability', 'population'], true)) {
            throw self::error(ErrorCode::INVALID_INPUT, "type must be 'frequency', 'reliability' or 'population'");
        }
        [$v, $w] = self::weightedPairs($values, $weights, $missing);
        $total = self::checkWeights($v, $w);
        $m = self::plainWeightedMean($v, $w, $total);
        $sq = [];
        $lin = [];
        foreach ($v as $i => $value) {
            $sq[] = $w[$i] * ($value - $m) * ($value - $m);
            $lin[] = $w[$i] * ($value - $m);
        }
        $s = self::sum($lin);
        $ss = self::computed(self::sum($sq) - $s * $s / $total);
        $ss = $ss < 0.0 ? 0.0 : $ss;
        $denominator = self::computed(match ($type) {
            'frequency' => $total - 1.0,
            'reliability' => $total - self::sum(array_map(static fn (float $x): float => $x * $x, $w)) / $total,
            'population' => $total,
        });
        if ($denominator <= 0.0) {
            throw self::error(ErrorCode::INVALID_INPUT, $type === 'frequency'
                ? 'frequency weights must sum to more than 1'
                : 'reliability weights need at least two positive weights');
        }

        return self::computed($ss / $denominator);
    }

    // ---- PHP glue -------------------------------------------------------

    private static function error(int $code, string $message): \Exception
    {
        return ErrorCode::exception($code, $message);
    }

    /** A float, or null for a missing value. */
    private static function element(mixed $value): ?float
    {
        if ($value === null) {
            return null;
        }
        if (! \is_int($value) && ! \is_float($value)) {
            throw self::error(ErrorCode::INVALID_INPUT, 'values must be integers, floats, or null');
        }
        $value = (float) $value;
        if (is_nan($value)) {
            return null;
        }
        if (is_infinite($value)) {
            throw self::error(ErrorCode::NONFINITE, 'infinite values are not allowed');
        }

        return $value;
    }

    /** @internal Also Smoothing's missing-value policy. */
    public static function dropMissing(string $missing): bool
    {
        return match ($missing) {
            'drop' => true,
            'raise' => false,
            default => throw self::error(ErrorCode::INVALID_INPUT, "missing must be 'raise' or 'drop'"),
        };
    }

    private static function missingError(): \Exception
    {
        return self::error(ErrorCode::INVALID_INPUT, "missing values found with missing='raise'");
    }

    /**
     * @internal Also Smoothing's missing-value policy.
     *
     * @return list<float>
     */
    public static function values(array $input, string $missing): array
    {
        $drop = self::dropMissing($missing);
        $out = [];
        foreach ($input as $value) {
            $v = self::element($value);
            if ($v !== null) {
                $out[] = $v;
            } elseif (! $drop) {
                throw self::missingError();
            }
        }

        return $out;
    }

    /**
     * @internal Also Smoothing's missing-value policy.
     *
     * @return array{list<float>, list<float>}
     */
    public static function pairs(array $x, array $y, string $missing): array
    {
        $drop = self::dropMissing($missing);
        if (\count($x) !== \count($y)) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'arrays must have the same length');
        }
        $a = [];
        $b = [];
        $ys = array_values($y);
        foreach (array_values($x) as $i => $u) {
            $u = self::element($u);
            $v = self::element($ys[$i]);
            if ($u !== null && $v !== null) {
                $a[] = $u;
                $b[] = $v;
            } elseif (! $drop) {
                throw self::missingError();
            }
        }

        return [$a, $b];
    }

    /** @return array{list<float>, list<float>} */
    private static function weightedPairs(array $values, array $weights, string $missing): array
    {
        foreach ($weights as $weight) {
            $w = self::element($weight);
            if ($w !== null && $w < 0.0) {
                throw self::error(ErrorCode::INVALID_INPUT, 'weights must be nonnegative');
            }
        }

        return self::pairs($values, $weights, $missing);
    }

    private static function population(string $type): bool
    {
        return match ($type) {
            'sample' => false,
            'population' => true,
            default => throw self::error(ErrorCode::INVALID_INPUT, "type must be 'sample' or 'population'"),
        };
    }

    // ---- numerical core -------------------------------------------------

    private static function computed(float $value): float
    {
        if (! is_finite($value)) {
            throw self::error(ErrorCode::NONFINITE, 'descriptive calculation exceeded floating-point range');
        }

        return $value;
    }

    /** @param list<float> $values */
    private static function check(array $values, int $minimum, string $what): void
    {
        foreach ($values as $v) {
            if (! is_finite($v)) {
                throw self::error(ErrorCode::NONFINITE, 'values must be finite');
            }
        }
        if (\count($values) < $minimum) {
            $count = $minimum === 1 ? 'one value' : 'two values';
            throw self::error(ErrorCode::INVALID_INPUT, "$what requires at least $count");
        }
    }

    /**
     * @param list<float> $x
     * @param list<float> $y
     */
    private static function checkPair(array $x, array $y, string $what): void
    {
        if (\count($x) !== \count($y)) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'x and y must have the same length');
        }
        self::check($x, 2, $what);
        self::check($y, 2, $what);
    }

    /** Neumaier compensated sum. @param list<float> $values */
    private static function sum(array $values): float
    {
        $total = 0.0;
        $compensation = 0.0;
        foreach ($values as $v) {
            $t = $total + $v;
            if (abs($total) >= abs($v)) {
                $compensation += ($total - $t) + $v;
            } else {
                $compensation += ($v - $t) + $total;
            }
            $total = $t;
        }

        return $total + $compensation;
    }

    /** @param list<float> $values */
    private static function plainMean(array $values): float
    {
        $n = (float) \count($values);
        $m = self::sum($values) / $n;

        return $m + self::sum(array_map(static fn (float $v): float => $v - $m, $values)) / $n;
    }

    /** @param list<float> $values */
    private static function variance(array $values, bool $population): float
    {
        self::check($values, $population ? 1 : 2, 'variance');
        $n = (float) \count($values);
        $m = self::plainMean($values);
        $ss = self::sum(array_map(static fn (float $v): float => ($v - $m) * ($v - $m), $values));
        $s = self::sum(array_map(static fn (float $v): float => $v - $m, $values));
        $ss = self::computed($ss - $s * $s / $n);
        $ss = $ss < 0.0 ? 0.0 : $ss;

        return self::computed($ss / ($population ? $n : $n - 1.0));
    }

    /**
     * @param list<float> $x
     * @param list<float> $y
     * @return array{float, float, float}
     */
    private static function crossProducts(array $x, array $y): array
    {
        $mx = self::plainMean($x);
        $my = self::plainMean($y);
        $n = (float) \count($x);
        $dx = array_map(static fn (float $v): float => $v - $mx, $x);
        $dy = array_map(static fn (float $v): float => $v - $my, $y);
        $sx = self::sum($dx);
        $sy = self::sum($dy);
        $xy = [];
        $xx = [];
        $yy = [];
        foreach ($dx as $i => $a) {
            $b = $dy[$i];
            $xy[] = $a * $b;
            $xx[] = $a * $a;
            $yy[] = $b * $b;
        }
        $sxy = self::computed(self::sum($xy) - $sx * $sy / $n);
        $sxx = self::computed(self::sum($xx) - $sx * $sx / $n);
        $syy = self::computed(self::sum($yy) - $sy * $sy / $n);

        return [$sxy, $sxx < 0.0 ? 0.0 : $sxx, $syy < 0.0 ? 0.0 : $syy];
    }

    /**
     * Type 7 quantiles from one sort.
     *
     * @param list<float> $values
     * @param list<float> $probabilities
     * @return list<float>
     */
    private static function quantiles(array $values, array $probabilities): array
    {
        self::check($values, 1, 'quantile');
        foreach ($probabilities as $p) {
            if (is_nan($p)) {
                throw self::error(ErrorCode::NONFINITE, 'probabilities must not be NaN');
            }
            if ($p < 0.0 || $p > 1.0) {
                throw self::error(ErrorCode::INVALID_INPUT, 'probabilities must be in [0, 1]');
            }
        }
        sort($values, SORT_NUMERIC);
        $n = \count($values);
        $out = [];
        foreach ($probabilities as $p) {
            $index = 1.0 + ($n - 1) * $p;
            $lo = floor($index);
            $h = $index - $lo;
            $lo = (int) $lo - 1;
            $low = $values[$lo];
            $out[] = ($h > 0.0 && $values[$lo + 1] != $low)
                ? (1.0 - $h) * $low + $h * $values[$lo + 1]
                : $low;
        }

        return $out;
    }

    /**
     * @param list<float> $values
     * @param list<float> $weights
     */
    private static function checkWeights(array $values, array $weights): float
    {
        if (\count($values) !== \count($weights)) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'values and weights must have the same length');
        }
        self::check($values, 1, 'weighted statistics');
        foreach ($weights as $w) {
            if (! is_finite($w)) {
                throw self::error(ErrorCode::NONFINITE, 'weights must be finite');
            }
        }
        foreach ($weights as $w) {
            if ($w < 0.0) {
                throw self::error(ErrorCode::INVALID_INPUT, 'weights must be nonnegative');
            }
        }
        $total = self::computed(self::sum($weights));
        if ($total <= 0.0) {
            throw self::error(ErrorCode::INVALID_INPUT, 'weights must have a positive sum');
        }

        return $total;
    }

    /**
     * @param list<float> $values
     * @param list<float> $weights
     */
    private static function plainWeightedMean(array $values, array $weights, float $total): float
    {
        $products = [];
        foreach ($values as $i => $v) {
            $products[] = $v * $weights[$i];
        }
        $m = self::sum($products) / $total;
        $residuals = [];
        foreach ($values as $i => $v) {
            $residuals[] = $weights[$i] * ($v - $m);
        }

        return $m + self::sum($residuals) / $total;
    }
}
