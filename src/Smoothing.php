<?php

declare(strict_types=1);

namespace Phit\Polyfill;

/**
 * Kernel density estimation, rule-of-thumb bandwidths, interpolation and the
 * empirical CDF, ported from the phit extension (stats-core smoothing.rs plus
 * the PHP glue), which follows R's stats::density, bw.nrd0, bw.nrd, approx
 * and ecdf. Data arrays follow the descriptive functions: null and NAN are
 * missing, rejected by missing: 'raise' or removed by 'drop'; infinite values
 * are errors. Evaluation points (xout, q) may be null, NAN or infinite; a null
 * or NAN point gives NAN, as R gives NA.
 */
final class Smoothing
{
    private const array KERNELS = ['gaussian', 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine', 'optcosine'];

    private const float FRAC_1_SQRT_2PI = 0.398942280401432677939946059934;

    /**
     * Phit\stat_density
     *
     * @param array<int|float|null> $x
     * @param array<int|float|null>|null $weights
     * @return array{x: list<float>, y: list<float>, bw: float, n: int}
     */
    public static function density(array $x, mixed $bw = 'nrd0', float $adjust = 1.0, string $kernel = 'gaussian', int $n = 512, ?float $from = null, ?float $to = null, float $cut = 3.0, ?array $weights = null, string $missing = 'raise'): array
    {
        if (! \in_array($kernel, self::KERNELS, true)) {
            throw self::invalid("kernel must be 'gaussian', 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine' or 'optcosine'");
        }
        if ($n < 1) {
            throw self::invalid('n must be at least 1');
        }
        if ($weights === null) {
            $x = Descriptive::values($x, $missing);
        } else {
            // As R's na.rm: weights that summed to one are renormalized after
            // dropping, so the estimate stays a density.
            $total = 0.0;
            foreach ($weights as $w) {
                if (\is_int($w) || \is_float($w)) {
                    $total += (float) $w;
                }
            }
            [$x, $weights] = Descriptive::pairs($x, $weights, $missing);
            if (Descriptive::dropMissing($missing) && abs($total - 1.0) < 1.5e-8) {
                $kept = 0.0;
                foreach ($weights as $w) {
                    $kept += $w;
                }
                foreach ($weights as $i => $w) {
                    $weights[$i] = $w / $kept;
                }
            }
        }
        $rule = self::bandwidthRule($bw);

        return self::estimate($x, $rule, $adjust, $kernel, $n, $from, $to, $cut, $weights);
    }

    /**
     * Phit\stat_bw_nrd0
     *
     * @param array<int|float|null> $x
     */
    public static function bwNrd0(array $x, string $missing = 'raise'): float
    {
        return self::nrd0(Descriptive::values($x, $missing));
    }

    /**
     * Phit\stat_bw_nrd
     *
     * @param array<int|float|null> $x
     */
    public static function bwNrd(array $x, string $missing = 'raise'): float
    {
        return self::nrd(Descriptive::values($x, $missing));
    }

    /**
     * Phit\stat_approx
     *
     * @param array<int|float|null> $x
     * @param array<int|float|null> $y
     * @param array<int|float|null> $xout
     * @return list<float>
     */
    public static function approx(array $x, array $y, array $xout, string $method = 'linear', int $rule = 1, string $ties = 'mean', float $f = 0.0, string $missing = 'raise'): array
    {
        $constant = match ($method) {
            'linear' => false,
            'constant' => true,
            default => throw self::invalid("method must be 'linear' or 'constant'"),
        };
        if ($rule !== 1 && $rule !== 2) {
            throw self::invalid('rule must be 1 or 2');
        }
        if (! \in_array($ties, ['mean', 'min', 'max'], true)) {
            throw self::invalid("ties must be 'mean', 'min' or 'max'");
        }
        [$x, $y] = Descriptive::pairs($x, $y, $missing);
        $points = self::points($xout, 'xout');
        if ($constant && (is_nan($f) || $f < 0.0 || $f > 1.0)) {
            throw self::invalid('f must be between 0 and 1');
        }

        // A stable sort by Rust's total order, which puts -0.0 before 0.0.
        $order = array_keys($x);
        usort($order, static fn (int $a, int $b): int => self::totalCompare($x[$a], $x[$b]));
        $ux = [];
        $uy = [];
        $count = \count($order);
        for ($start = 0; $start < $count;) {
            $end = $start + 1;
            while ($end < $count && $x[$order[$end]] == $x[$order[$start]]) {
                $end++;
            }
            $group = [];
            for ($k = $start; $k < $end; $k++) {
                $group[] = $y[$order[$k]];
            }
            $ux[] = $x[$order[$start]];
            $uy[] = match ($ties) {
                'mean' => Descriptive::mean($group),
                'min' => array_reduce($group, self::min(...), INF),
                'max' => array_reduce($group, self::max(...), -INF),
            };
            $start = $end;
        }
        if ($ux === []) {
            throw self::invalid('approx requires at least one value');
        }
        if (\count($ux) === 1 && ! $constant) {
            throw self::invalid('linear interpolation requires at least two distinct x values');
        }
        [$low, $high] = $rule === 1 ? [NAN, NAN] : [$uy[0], $uy[\count($uy) - 1]];
        $out = [];
        foreach ($points as $v) {
            $out[] = self::interpolate($ux, $uy, $v, $constant ? $f : null, $low, $high);
        }

        return $out;
    }

    /**
     * Phit\stat_ecdf
     *
     * @param array<int|float|null> $x
     * @param array<int|float|null> $q
     * @return list<float>
     */
    public static function ecdf(array $x, array $q, string $missing = 'raise'): array
    {
        $sorted = Descriptive::values($x, $missing);
        $points = self::points($q, 'q');
        if ($sorted === []) {
            throw self::invalid('ecdf requires at least one value');
        }
        sort($sorted, SORT_NUMERIC);
        $count = \count($sorted);
        $n = (float) $count;
        $out = [];
        foreach ($points as $v) {
            if (is_nan($v)) {
                $out[] = NAN;
                continue;
            }
            // The number of sorted values at or below $v.
            [$lo, $hi] = [0, $count];
            while ($lo < $hi) {
                $mid = ($lo + $hi) >> 1;
                if ($sorted[$mid] <= $v) {
                    $lo = $mid + 1;
                } else {
                    $hi = $mid;
                }
            }
            $out[] = $lo / $n;
        }

        return $out;
    }

    // ---- PHP glue -------------------------------------------------------

    private static function invalid(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    /** 'nrd0', 'nrd' or a fixed bandwidth. */
    private static function bandwidthRule(mixed $bw): string|float
    {
        if (\is_string($bw)) {
            return match (strtolower($bw)) {
                'nrd0' => 'nrd0',
                'nrd' => 'nrd',
                default => throw self::invalid("bw must be 'nrd0', 'nrd' or a positive number"),
            };
        }
        if (\is_int($bw) || \is_float($bw)) {
            return (float) $bw;
        }

        throw new \TypeError('Phit\stat_density(): Argument #2 ($bw) must be of type string|int|float, '.get_debug_type($bw).' given');
    }

    /**
     * Evaluation points: numbers, with null read as NAN.
     *
     * @return list<float>
     */
    private static function points(array $input, string $what): array
    {
        $out = [];
        foreach ($input as $v) {
            if ($v === null) {
                $out[] = NAN;
            } elseif (\is_int($v) || \is_float($v)) {
                $out[] = (float) $v;
            } else {
                throw self::invalid("$what must be integers, floats, or null");
            }
        }

        return $out;
    }

    private static function totalCompare(float $a, float $b): int
    {
        if ($a == $b && $a == 0.0) {
            return self::negative($b) <=> self::negative($a);
        }

        return $a <=> $b;
    }

    private static function negative(float $zero): bool
    {
        return fdiv(1.0, $zero) < 0.0;
    }

    /** Rust's f64::min on aarch64 (fminnm): -0.0 is below 0.0. */
    private static function min(float $a, float $b): float
    {
        if ($a == $b && $a == 0.0) {
            return self::negative($a) ? $a : $b;
        }

        return $b < $a ? $b : $a;
    }

    /** Rust's f64::max on aarch64 (fmaxnm): 0.0 is above -0.0. */
    private static function max(float $a, float $b): float
    {
        if ($a == $b && $a == 0.0) {
            return self::negative($a) ? $b : $a;
        }

        return $b > $a ? $b : $a;
    }

    // ---- numerical core -------------------------------------------------

    /**
     * R's seq.int(from, to, length.out = n): exact ends, the first half
     * stepped from $from and the second half from $to.
     *
     * @return list<float>
     */
    private static function seq(float $from, float $to, int $n): array
    {
        $out = array_fill(0, $n, 0.0);
        if ($n > 0) {
            $out[0] = $from;
        }
        if ($n > 1) {
            $out[$n - 1] = $to;
        }
        if ($n > 2) {
            $by = ($to - $from) / ($n - 1);
            $half = intdiv($n, 2);
            for ($i = 1; $i < $n - 1; $i++) {
                $out[$i] = $i < $half ? $from + $i * $by : $to - ($n - 1 - $i) * $by;
            }
        }

        return $out;
    }

    /** @param list<float> $x */
    private static function checkBandwidthInput(array $x): void
    {
        self::checkFinite($x, 'values');
        if (\count($x) < 2) {
            throw self::invalid('bandwidth selection requires at least two values');
        }
    }

    /** @param list<float> $values */
    private static function checkFinite(array $values, string $what): void
    {
        foreach ($values as $v) {
            if (! is_finite($v)) {
                throw ErrorCode::exception(ErrorCode::NONFINITE, "$what must be finite");
            }
        }
    }

    /**
     * R's bw.nrd0: 0.9 * min(sd, IQR / 1.34) * n^(-1/5), with R's fallbacks
     * when that is zero.
     *
     * @param list<float> $x
     */
    private static function nrd0(array $x): float
    {
        self::checkBandwidthInput($x);
        $hi = Descriptive::sd($x);
        $q = Descriptive::quantileArray($x, [0.25, 0.75]);
        $lo = self::min($hi, ($q[1] - $q[0]) / 1.34);
        if ($lo == 0.0) {
            $lo = match (true) {
                $hi != 0.0 => $hi,
                $x[0] != 0.0 => abs($x[0]),
                default => 1.0,
            };
        }

        return 0.9 * $lo * ((float) \count($x)) ** -0.2;
    }

    /**
     * R's bw.nrd: 1.06 * min(sd, IQR / 1.34) * n^(-1/5).
     *
     * @param list<float> $x
     */
    private static function nrd(array $x): float
    {
        self::checkBandwidthInput($x);
        $q = Descriptive::quantileArray($x, [0.25, 0.75]);
        $h = ($q[1] - $q[0]) / 1.34;

        return 1.06 * self::min(Descriptive::sd($x), $h) * ((float) \count($x)) ** (-1.0 / 5.0);
    }

    /** A kernel of R's density(), scaled to standard deviation $bw. */
    private static function kernel(string $kernel, float $x, float $bw): float
    {
        $ax = abs($x);
        switch ($kernel) {
            case 'gaussian':
                return self::dnorm($x, $bw);
            case 'rectangular':
                $a = $bw * sqrt(3.0);

                return $ax < $a ? 0.5 / $a : 0.0;
            case 'triangular':
                $a = $bw * sqrt(6.0);

                return $ax < $a ? (1.0 - $ax / $a) / $a : 0.0;
            case 'epanechnikov':
                $a = $bw * sqrt(5.0);
                $t = $ax / $a;

                return $ax < $a ? 3.0 / 4.0 * (1.0 - $t * $t) / $a : 0.0;
            case 'biweight':
                $a = $bw * sqrt(7.0);
                $t = $ax / $a;
                $u = 1.0 - $t * $t;

                return $ax < $a ? 15.0 / 16.0 * ($u * $u) / $a : 0.0;
            case 'cosine':
                $a = $bw / sqrt(1.0 / 3.0 - 2.0 / (M_PI * M_PI));

                return $ax < $a ? (1.0 + cos(M_PI * $x / $a)) / (2.0 * $a) : 0.0;
            default:
                $a = $bw / sqrt(1.0 - 8.0 / (M_PI * M_PI));

                return $ax < $a ? M_PI / 4.0 * cos(M_PI * $x / (2.0 * $a)) / $a : 0.0;
        }
    }

    /**
     * R's dnorm(x, 0, sigma), including its split evaluation of exp(-x²/2) in
     * the tail, so kernel heights agree with R to the last bit.
     */
    private static function dnorm(float $x, float $sigma): float
    {
        $x = abs($x / $sigma);
        if (! is_finite($x) || $x >= 2.0 * sqrt(PHP_FLOAT_MAX)) {
            return 0.0;
        }
        if ($x < 5.0) {
            return self::FRAC_1_SQRT_2PI * exp(-0.5 * $x * $x) / $sigma;
        }
        // sqrt(-2 ln 2 (DBL_MIN_EXP + 1 - DBL_MANT_DIG))
        if ($x > sqrt(-2.0 * M_LN2 * (-1021 + 1 - 53))) {
            return 0.0;
        }
        $x1 = self::roundTiesEven($x * 65536.0) / 65536.0;
        $x2 = $x - $x1;

        return self::FRAC_1_SQRT_2PI / $sigma * (exp(-0.5 * $x1 * $x1) * exp((-0.5 * $x2 - $x1) * $x2));
    }

    /** Rust's round_ties_even, exact for every double, unlike PHP 8.3's rounding. */
    private static function roundTiesEven(float $v): float
    {
        $f = floor($v);
        $d = $v - $f;
        if ($d > 0.5) {
            return $f + 1.0;
        }
        if ($d < 0.5) {
            return $f;
        }

        return fmod($f, 2.0) == 0.0 ? $f : $f + 1.0;
    }

    /**
     * R's C BinDist: linear binning of weighted points onto n grid points
     * from lo to hi, returned zero-padded to length 2n.
     *
     * @param list<float> $x
     * @param list<float> $w
     * @return list<float>
     */
    private static function binDist(array $x, array $w, float $lo, float $hi, int $n): array
    {
        $y = array_fill(0, 2 * $n, 0.0);
        $ixmax = $n - 2;
        $delta = ($hi - $lo) / ($n - 1);
        foreach ($x as $k => $xi) {
            $wi = $w[$k];
            $pos = ($xi - $lo) / $delta;
            $ix = floor($pos);
            $fx = $pos - $ix;
            if ($ix >= 0.0 && $ix <= $ixmax) {
                $i = (int) $ix;
                $y[$i] += $wi * (1.0 - $fx);
                $y[$i + 1] += $wi * $fx;
            } elseif ($ix == -1.0) {
                $y[0] += $wi * $fx;
            } elseif ($ix == $ixmax + 1) {
                $y[$ixmax + 1] += $wi * (1.0 - $fx);
            }
        }

        return $y;
    }

    /**
     * In-place iterative radix-2 FFT; $inverse uses the positive exponent and
     * does not scale, like R's fft(inverse = TRUE).
     *
     * @param list<float> $re
     * @param list<float> $im
     */
    private static function fft(array &$re, array &$im, bool $inverse): void
    {
        $n = \count($re);
        $j = 0;
        for ($i = 1; $i < $n; $i++) {
            $bit = $n >> 1;
            while (($j & $bit) !== 0) {
                $j ^= $bit;
                $bit >>= 1;
            }
            $j |= $bit;
            if ($i < $j) {
                [$re[$i], $re[$j]] = [$re[$j], $re[$i]];
                [$im[$i], $im[$j]] = [$im[$j], $im[$i]];
            }
        }
        $sign = $inverse ? 1.0 : -1.0;
        for ($len = 2; $len <= $n; $len <<= 1) {
            $half = $len >> 1;
            $step = $sign * 2.0 * M_PI / $len;
            for ($k = 0; $k < $half; $k++) {
                // Twiddles computed directly rather than by recurrence, as in
                // the extension.
                $angle = $step * $k;
                $s = sin($angle);
                $c = cos($angle);
                for ($a = $k; $a < $n; $a += $len) {
                    $b = $a + $half;
                    $tr = $re[$b] * $c - $im[$b] * $s;
                    $ti = $re[$b] * $s + $im[$b] * $c;
                    $re[$b] = $re[$a] - $tr;
                    $im[$b] = $im[$a] - $ti;
                    $re[$a] += $tr;
                    $im[$a] += $ti;
                }
            }
        }
    }

    /**
     * R's density() for values the glue has already checked for missing ones.
     *
     * @param list<float> $x
     * @param list<float>|null $weights
     * @return array{x: list<float>, y: list<float>, bw: float, n: int}
     */
    private static function estimate(array $x, string|float $rule, float $adjust, string $kernel, int $nUser, ?float $from, ?float $to, float $cut, ?array $weights): array
    {
        self::checkFinite($x, 'values');
        if ($x === []) {
            throw self::invalid('density requires at least one value');
        }
        $nx = \count($x);
        if ($weights === null) {
            $weights = array_fill(0, $nx, 1.0 / $nx);
        } else {
            if (\count($weights) !== $nx) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'values and weights must have the same length');
            }
            self::checkFinite($weights, 'weights');
            foreach ($weights as $w) {
                if ($w < 0.0) {
                    throw self::invalid('weights must be nonnegative');
                }
            }
        }
        // At least 512, else the next power of two.
        $n = 512;
        while ($n < $nUser) {
            $n <<= 1;
        }
        $bw = match ($rule) {
            'nrd0' => self::nrd0($x),
            'nrd' => self::nrd($x),
            default => $rule,
        };
        if (! is_finite($bw) || ! is_finite($adjust) || ! is_finite($cut)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'bw, adjust and cut must be finite');
        }
        $bw = $adjust * $bw;
        if ($bw <= 0.0) {
            throw self::invalid('the bandwidth must be positive');
        }
        $min = INF;
        $max = -INF;
        foreach ($x as $v) {
            $min = self::min($min, $v);
            $max = self::max($max, $v);
        }
        $from ??= $min - $cut * $bw;
        $to ??= $max + $cut * $bw;
        if (! is_finite($from) || ! is_finite($to)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'from and to must be finite');
        }
        if ($from >= $to) {
            throw self::invalid('from must be less than to');
        }
        $lo = $from - 4.0 * $bw;
        $up = $to + 4.0 * $bw;

        $yRe = self::binDist($x, $weights, $lo, $up, $n);
        $kords = self::seq(0.0, (float) (2 * $n - 1) / ($n - 1) * ($up - $lo), 2 * $n);
        for ($j = $n + 1; $j < 2 * $n; $j++) {
            $kords[$j] = -$kords[2 * $n - $j];
        }
        $kRe = [];
        foreach ($kords as $k) {
            $kRe[] = self::kernel($kernel, $k, $bw);
        }
        $yIm = array_fill(0, 2 * $n, 0.0);
        $kIm = $yIm;
        self::fft($yRe, $yIm, false);
        self::fft($kRe, $kIm, false);
        // fft(y) * Conj(fft(kords))
        for ($i = 0; $i < 2 * $n; $i++) {
            [$a, $b, $c, $d] = [$yRe[$i], $yIm[$i], $kRe[$i], -$kIm[$i]];
            $yRe[$i] = $a * $c - $b * $d;
            $yIm[$i] = $a * $d + $b * $c;
        }
        self::fft($yRe, $yIm, true);
        $scale = (float) (2 * $n);
        $heights = [];
        for ($i = 0; $i < $n; $i++) {
            $h = $yRe[$i] / $scale;
            $heights[] = $h > 0.0 ? $h : 0.0;
        }
        $grid = self::seq($lo, $up, $n);
        $outX = self::seq($from, $to, $nUser);
        $outY = [];
        foreach ($outX as $v) {
            $outY[] = self::interpolate($grid, $heights, $v, null, NAN, NAN);
        }

        return ['x' => $outX, 'y' => $outY, 'bw' => $bw, 'n' => $nx];
    }

    /**
     * R's C approx1 for sorted unique x. $f is null for linear interpolation,
     * otherwise the right-hand weight of the step function.
     *
     * @param list<float> $x
     * @param list<float> $y
     */
    private static function interpolate(array $x, array $y, float $v, ?float $f, float $low, float $high): float
    {
        if (is_nan($v)) {
            return NAN;
        }
        $i = 0;
        $j = \count($x) - 1;
        if ($v < $x[$i]) {
            return $low;
        }
        if ($v > $x[$j]) {
            return $high;
        }
        while ($i + 1 < $j) {
            $ij = ($i + $j) >> 1;
            if ($v < $x[$ij]) {
                $j = $ij;
            } else {
                $i = $ij;
            }
        }
        if ($v == $x[$j]) {
            return $y[$j];
        }
        if ($v == $x[$i]) {
            return $y[$i];
        }
        if ($f === null) {
            return $y[$i] + ($y[$j] - $y[$i]) * (($v - $x[$i]) / ($x[$j] - $x[$i]));
        }
        $f1 = 1.0 - $f;

        return ($f1 != 0.0 ? $y[$i] * $f1 : 0.0) + ($f != 0.0 ? $y[$j] * $f : 0.0);
    }
}
