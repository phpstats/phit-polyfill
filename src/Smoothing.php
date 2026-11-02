<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Frontend;
use Phit\Polyfill\Internal\Fma;

/**
 * Kernel density estimation, interpolation and the empirical CDF, as R's
 * stats::density, stats::approx and stats::ecdf, following the phit
 * extension: stats-core smoothing.rs, the routines of numeric.rs (linear
 * binning, interpolation and the normal density, written from the clean-room
 * specification specs/cleanroom/smoothing.md) and the PHP glue in
 * php-extension smoothing.rs.
 *
 * Data arrays follow the descriptive functions: null and NAN are missing,
 * rejected by missing='raise' (the default) or removed by 'drop'; infinite
 * values are errors. Evaluation points (xout, q) may be null, NAN or infinite;
 * a null or NAN point gives NAN, as R gives NA.
 */
final class Smoothing
{
    private const array KERNELS = ['gaussian', 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine', 'optcosine'];

    /** 1 / sqrt(2π) rounded to double. */
    private const float RECIP_SQRT_TAU = 0.3989422804014327;

    /**
     * Kernel density estimate, as R's density(): x and y on an n-point grid,
     * the bandwidth used and the number of observations.
     *
     * @param array<int|float|null> $x
     * @param ?array<int|float|null> $weights
     * @return array{x: list<float>, y: list<float>, bw: float, n: int}
     */
    public static function density(
        array $x,
        mixed $bw = 'nrd0',
        float $adjust = 1.0,
        string $kernel = 'gaussian',
        int $n = 512,
        ?float $from = null,
        ?float $to = null,
        float $cut = 3.0,
        ?array $weights = null,
        string $missing = 'raise',
    ): array {
        $bw = self::bwArgument($bw);
        Frontend::stringArgument($kernel, 'kernel');
        Frontend::stringArgument($missing, 'missing');
        if (! \in_array($kernel, self::KERNELS, true)) {
            throw self::invalid("kernel must be 'gaussian', 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine' or 'optcosine'");
        }
        if ($n < 1) {
            throw self::invalid('n must be at least 1');
        }
        if ($weights === null) {
            $values = Descriptive::values($x, $missing);
            $w = null;
        } else {
            // As R's na.rm: weights that summed to one are renormalized after
            // dropping, so the estimate stays a density.
            $total = -0.0;
            foreach ($weights as $v) {
                if (\is_float($v) || \is_int($v)) {
                    $total += (float) $v;
                }
            }
            [$values, $w] = Descriptive::pairs($x, $weights, $missing);
            if (Descriptive::dropMissing($missing) && abs($total - 1.0) < 1.5e-8) {
                $kept = -0.0;
                foreach ($w as $v) {
                    $kept += $v;
                }
                foreach ($w as $i => $v) {
                    $w[$i] = fdiv($v, $kept);
                }
            }
        }
        $bandwidth = self::bandwidth($bw);

        return self::estimate($values, $bandwidth, $adjust, $kernel, $n, $from, $to, $cut, $w);
    }

    /**
     * R's bw.nrd0: Silverman's rule-of-thumb bandwidth.
     *
     * @param array<int|float|null> $x
     */
    public static function bwNrd0(array $x, string $missing = 'raise'): float
    {
        Frontend::stringArgument($missing, 'missing');

        return self::nrd0(Descriptive::values($x, $missing));
    }

    /**
     * R's bw.nrd: Scott's variation of the rule-of-thumb bandwidth.
     *
     * @param array<int|float|null> $x
     */
    public static function bwNrd(array $x, string $missing = 'raise'): float
    {
        Frontend::stringArgument($missing, 'missing');

        return self::nrd(Descriptive::values($x, $missing));
    }

    /**
     * Interpolate (x, y) at xout, as R's approx().
     *
     * @param array<int|float|null> $x
     * @param array<int|float|null> $y
     * @param array<int|float|null> $xout
     * @return list<float>
     */
    public static function approx(
        array $x,
        array $y,
        array $xout,
        string $method = 'linear',
        int $rule = 1,
        string $ties = 'mean',
        float $f = 0.0,
        string $missing = 'raise',
    ): array {
        Frontend::stringArgument($method, 'method');
        Frontend::stringArgument($ties, 'ties');
        Frontend::stringArgument($missing, 'missing');
        $constant = match ($method) {
            'linear' => false,
            'constant' => true,
            default => throw self::invalid("method must be 'linear' or 'constant'"),
        };
        $nearest = match ($rule) {
            1 => false,
            2 => true,
            default => throw self::invalid('rule must be 1 or 2'),
        };
        if (! \in_array($ties, ['mean', 'min', 'max'], true)) {
            throw self::invalid("ties must be 'mean', 'min' or 'max'");
        }
        [$xs, $ys] = Descriptive::pairs($x, $y, $missing);
        $points = self::points($xout, 'xout');

        return self::interpolateAll($xs, $ys, $points, $constant, $f, $nearest, $ties);
    }

    /**
     * The empirical CDF of x at each q: the proportion of x at or below q.
     *
     * @param array<int|float|null> $x
     * @param array<int|float|null> $q
     * @return list<float>
     */
    public static function ecdf(array $x, array $q, string $missing = 'raise'): array
    {
        Frontend::stringArgument($missing, 'missing');
        $values = Descriptive::values($x, $missing);
        $points = self::points($q, 'q');
        self::checkFinite($values, 'values');
        if ($values === []) {
            throw self::invalid('ecdf requires at least one value');
        }
        $sorted = $values;
        usort($sorted, Descriptive::totalOrder(...));
        $n = (float) \count($sorted);
        $out = [];
        foreach ($points as $v) {
            if (is_nan($v)) {
                $out[] = NAN;

                continue;
            }
            // The number of sorted values at or below v.
            $lo = 0;
            $hi = \count($sorted);
            while ($lo < $hi) {
                $mid = ($lo + $hi) >> 1;
                if ($sorted[$mid] <= $v) {
                    $lo = $mid + 1;
                } else {
                    $hi = $mid;
                }
            }
            $out[] = (float) $lo / $n;
        }

        return $out;
    }

    // ---- PHP glue -------------------------------------------------------

    private static function error(int $code, string $message): \Exception
    {
        return ErrorCode::exception($code, $message);
    }

    private static function invalid(string $message): \Exception
    {
        return self::error(ErrorCode::INVALID_INPUT, $message);
    }

    /**
     * The bw argument as the extension reads it: a rule name or a number.
     */
    private static function bwArgument(mixed $bw): string|float
    {
        if ($bw === null) {
            throw new \TypeError('Argument `$bw` must not be null');
        }
        if (\is_string($bw) && mb_check_encoding($bw, 'UTF-8')) {
            return $bw;
        }
        if (\is_float($bw) || \is_int($bw)) {
            return (float) $bw;
        }
        throw new \Exception('Invalid value given for argument `bw`.', 0);
    }

    /** 'nrd0', 'nrd' or a fixed bandwidth. */
    private static function bandwidth(string|float $bw): string|float
    {
        if (\is_float($bw)) {
            return $bw;
        }

        return match (strtolower($bw)) {
            'nrd0' => 'nrd0',
            'nrd' => 'nrd',
            default => throw self::invalid("bw must be 'nrd0', 'nrd' or a positive number"),
        };
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
            } elseif (\is_float($v) || \is_int($v)) {
                $out[] = (float) $v;
            } else {
                throw self::invalid("$what must be integers, floats, or null");
            }
        }

        return $out;
    }

    // ---- numerical core (stats-core smoothing.rs) ---------------------------

    /** @param list<float> $values */
    private static function checkFinite(array $values, string $what): void
    {
        foreach ($values as $v) {
            if (! is_finite($v)) {
                throw self::error(ErrorCode::NONFINITE, "$what must be finite");
            }
        }
    }

    /** Rust's f64::min on aarch64 (FMINNM): a NAN operand yields the other, and -0.0 is below 0.0. */
    private static function min(float $a, float $b): float
    {
        if (is_nan($a)) {
            return $b;
        }
        if (is_nan($b)) {
            return $a;
        }
        if ($a == $b) {
            return fdiv(1.0, $a) < 0.0 ? $a : $b;
        }

        return $b < $a ? $b : $a;
    }

    /** Rust's f64::max on aarch64 (FMAXNM): a NAN operand yields the other, and 0.0 is above -0.0. */
    private static function max(float $a, float $b): float
    {
        if (is_nan($a)) {
            return $b;
        }
        if (is_nan($b)) {
            return $a;
        }
        if ($a == $b) {
            return fdiv(1.0, $a) > 0.0 ? $a : $b;
        }

        return $b > $a ? $b : $a;
    }

    /**
     * R's seq.int(from, to, length.out = n): exact ends, the first half
     * stepped from `from` and the second half from `to`.
     *
     * @return list<float>
     */
    private static function seq(float $from, float $to, int $n): array
    {
        $out = $n > 0 ? array_fill(0, $n, 0.0) : [];
        if ($n > 0) {
            $out[0] = $from;
        }
        if ($n > 1) {
            $out[$n - 1] = $to;
        }
        if ($n > 2) {
            $by = ($to - $from) / (float) ($n - 1);
            $half = intdiv($n, 2);
            for ($i = 1; $i < $n - 1; $i++) {
                $out[$i] = $i < $half
                    ? $from + (float) $i * $by
                    : $to - (float) ($n - 1 - $i) * $by;
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

    /**
     * R's bw.nrd0: Silverman's rule of thumb, 0.9 * min(sd, IQR / 1.34) * n^(-1/5),
     * with R's fallbacks when that is zero.
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
            if ($hi != 0.0) {
                $lo = $hi;
            } elseif ($x[0] != 0.0) {
                $lo = abs($x[0]);
            } else {
                $lo = 1.0;
            }
        }

        return 0.9 * $lo * pow((float) \count($x), -0.2);
    }

    /**
     * R's bw.nrd: Scott's variation, 1.06 * min(sd, IQR / 1.34) * n^(-1/5).
     *
     * @param list<float> $x
     */
    private static function nrd(array $x): float
    {
        self::checkBandwidthInput($x);
        $q = Descriptive::quantileArray($x, [0.25, 0.75]);
        $h = ($q[1] - $q[0]) / 1.34;

        return 1.06 * self::min(Descriptive::sd($x), $h) * pow((float) \count($x), -1.0 / 5.0);
    }

    /** One of R's smoothing kernels at x, scaled to standard deviation bw. */
    private static function kernel(string $kernel, float $x, float $bw): float
    {
        $ax = abs($x);
        switch ($kernel) {
            case 'gaussian':
                return fdiv(self::standardNormalDensity(fdiv($x, $bw)), $bw);
            case 'rectangular':
                $a = $bw * sqrt(3.0);
                $value = fdiv(0.5, $a);
                break;
            case 'triangular':
                $a = $bw * sqrt(6.0);
                $value = fdiv(1.0 - fdiv($ax, $a), $a);
                break;
            case 'epanechnikov':
                $a = $bw * sqrt(5.0);
                $t = fdiv($ax, $a);
                $value = fdiv(0.75 * (1.0 - $t * $t), $a);
                break;
            case 'biweight':
                $a = $bw * sqrt(7.0);
                $t = fdiv($ax, $a);
                $u = 1.0 - $t * $t;
                $value = fdiv(0.9375 * ($u * $u), $a);
                break;
            case 'cosine':
                $a = fdiv($bw, sqrt(1.0 / 3.0 - 2.0 / (M_PI * M_PI)));
                $value = fdiv(1.0 + cos(fdiv(M_PI * $x, $a)), 2.0 * $a);
                break;
            default: // optcosine
                $a = fdiv($bw, sqrt(1.0 - 8.0 / (M_PI * M_PI)));
                $value = fdiv(M_PI / 4.0 * cos(fdiv(M_PI * $x, 2.0 * $a)), $a);
                break;
        }

        return $ax < $a ? $value : 0.0;
    }

    /**
     * In-place iterative radix-2 FFT; inverse uses the positive exponent and
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
            $step = $sign * 2.0 * M_PI / (float) $len;
            for ($k = 0; $k < $half; $k++) {
                $angle = $step * (float) $k;
                $s = sin($angle);
                $c = cos($angle);
                for ($start = 0; $start < $n; $start += $len) {
                    $a = $start + $k;
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
     * R's density() for finite x.
     *
     * @param list<float> $x
     * @param ?list<float> $weights
     * @return array{x: list<float>, y: list<float>, bw: float, n: int}
     */
    private static function estimate(
        array $x,
        string|float $bandwidth,
        float $adjust,
        string $kernel,
        int $nUser,
        ?float $from,
        ?float $to,
        float $cut,
        ?array $weights,
    ): array {
        self::checkFinite($x, 'values');
        if ($x === []) {
            throw self::invalid('density requires at least one value');
        }
        $nx = \count($x);
        if ($weights === null) {
            $weights = array_fill(0, $nx, 1.0 / (float) $nx);
        } else {
            if (\count($weights) !== $nx) {
                throw self::error(ErrorCode::DIMENSION_MISMATCH, 'values and weights must have the same length');
            }
            self::checkFinite($weights, 'weights');
            foreach ($weights as $w) {
                if ($w < 0.0) {
                    throw self::invalid('weights must be nonnegative');
                }
            }
        }
        if ($nUser === 0) {
            throw self::invalid('n must be at least 1');
        }
        $n = 512;
        if ($nUser > 512) {
            $n = 1;
            while ($n < $nUser) {
                $n <<= 1;
            }
        }
        $bw = match ($bandwidth) {
            'nrd0' => self::nrd0($x),
            'nrd' => self::nrd($x),
            default => $bandwidth,
        };
        if (! is_finite($bw) || ! is_finite($adjust) || ! is_finite($cut)) {
            throw self::error(ErrorCode::NONFINITE, 'bw, adjust and cut must be finite');
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
            throw self::error(ErrorCode::NONFINITE, 'from and to must be finite');
        }
        if ($from >= $to) {
            throw self::invalid('from must be less than to');
        }
        $lo = $from - 4.0 * $bw;
        $up = $to + 4.0 * $bw;

        $yRe = self::linearBinning($x, $weights, $lo, $up, $n);
        $kords = self::seq(0.0, (float) (2 * $n - 1) / (float) ($n - 1) * ($up - $lo), 2 * $n);
        for ($j = $n + 1; $j < 2 * $n; $j++) {
            $kords[$j] = -$kords[2 * $n - $j];
        }
        $kRe = [];
        foreach ($kords as $k) {
            $kRe[] = self::kernel($kernel, $k, $bw);
        }
        $yIm = array_fill(0, 2 * $n, 0.0);
        $kIm = array_fill(0, 2 * $n, 0.0);
        self::fft($yRe, $yIm, false);
        self::fft($kRe, $kIm, false);
        // fft(y) * Conj(fft(kords))
        for ($i = 0; $i < 2 * $n; $i++) {
            $a = $yRe[$i];
            $b = $yIm[$i];
            $c = $kRe[$i];
            $d = -$kIm[$i];
            $yRe[$i] = $a * $c - $b * $d;
            $yIm[$i] = $a * $d + $b * $c;
        }
        self::fft($yRe, $yIm, true);
        $scale = (float) (2 * $n);
        $heights = [];
        for ($i = 0; $i < $n; $i++) {
            $heights[] = self::max(fdiv($yRe[$i], $scale), 0.0);
        }
        $grid = self::seq($lo, $up, $n);
        $outX = self::seq($from, $to, $nUser);
        $outY = [];
        foreach ($outX as $v) {
            $outY[] = self::interpolate($grid, $heights, $v, false, 0.0, NAN, NAN);
        }

        return ['x' => $outX, 'y' => $outY, 'bw' => $bw, 'n' => $nx];
    }

    /**
     * R's approx(x, y, xout): sorts by x, collapses ties, then interpolates.
     *
     * @param list<float> $x
     * @param list<float> $y
     * @param list<float> $xout
     * @return list<float>
     */
    private static function interpolateAll(array $x, array $y, array $xout, bool $constant, float $f, bool $nearest, string $ties): array
    {
        if (\count($x) !== \count($y)) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'x and y must have the same length');
        }
        self::checkFinite($x, 'x');
        self::checkFinite($y, 'y');
        if ($constant && (is_nan($f) || $f < 0.0 || $f > 1.0)) {
            throw self::invalid('f must be between 0 and 1');
        }
        // A stable sort of the indices by x in Rust's total order.
        $order = array_keys($x);
        usort($order, static fn (int $a, int $b): int => Descriptive::totalOrder($x[$a], $x[$b]));
        $ux = [];
        $uy = [];
        $count = \count($order);
        $start = 0;
        while ($start < $count) {
            $end = $start + 1;
            while ($end < $count && $x[$order[$end]] == $x[$order[$start]]) {
                $end++;
            }
            $group = [];
            for ($i = $start; $i < $end; $i++) {
                $group[] = $y[$order[$i]];
            }
            $ux[] = $x[$order[$start]];
            if ($ties === 'mean') {
                $uy[] = Descriptive::mean($group);
            } elseif ($ties === 'min') {
                $m = INF;
                foreach ($group as $v) {
                    $m = self::min($m, $v);
                }
                $uy[] = $m;
            } else {
                $m = -INF;
                foreach ($group as $v) {
                    $m = self::max($m, $v);
                }
                $uy[] = $m;
            }
            $start = $end;
        }
        if ($ux === []) {
            throw self::invalid('approx requires at least one value');
        }
        if (\count($ux) === 1 && ! $constant) {
            throw self::invalid('linear interpolation requires at least two distinct x values');
        }
        [$low, $high] = $nearest ? [$uy[0], $uy[\count($uy) - 1]] : [NAN, NAN];
        $out = [];
        foreach ($xout as $v) {
            $out[] = self::interpolate($ux, $uy, $v, $constant, $f, $low, $high);
        }

        return $out;
    }

    // ---- numeric.rs -----------------------------------------------------

    /**
     * Linear binning of weighted points onto n equally spaced grid points from
     * lo to hi (n >= 2, lo < hi), returned zero-padded to length 2n.
     *
     * @param list<float> $x
     * @param list<float> $w
     * @return list<float>
     */
    private static function linearBinning(array $x, array $w, float $lo, float $hi, int $n): array
    {
        $bins = array_fill(0, 2 * $n, 0.0);
        $step = ($hi - $lo) / ((float) $n - 1.0);
        $last = (float) $n - 1.0;
        foreach ($x as $i => $xi) {
            $wi = $w[$i];
            if (! is_finite($xi)) {
                continue;
            }
            $t = fdiv($xi - $lo, $step);
            $left = floor($t);
            // Both target bins lie outside 0..n-1: nothing to add.
            if (is_nan($left) || $left < -1.0 || $left > $last) {
                continue;
            }
            $frac = $t - $left;
            $k = (int) $left;
            if ($k >= 0) {
                $bins[$k] += $wi * (1.0 - $frac);
            }
            $k1 = $k + 1;
            if ($k1 >= 0 && $k1 < $n) {
                $bins[$k1] += $wi * $frac;
            }
        }

        return $bins;
    }

    /**
     * Interpolation at v from points (x[i], y[i]) with x strictly increasing.
     * Values left of x[0] are low; values right of the last x are high.
     * Linear when constant is false, otherwise R's constant method with weight f.
     *
     * @param list<float> $x
     * @param list<float> $y
     */
    private static function interpolate(array $x, array $y, float $v, bool $constant, float $f, float $low, float $high): float
    {
        if (is_nan($v)) {
            return NAN;
        }
        $m = \count($x);
        if ($v < $x[0]) {
            return $low;
        }
        if ($v > $x[$m - 1]) {
            return $high;
        }
        // First index whose x is >= v.
        $lo = 0;
        $hi = $m;
        while ($lo < $hi) {
            $mid = ($lo + $hi) >> 1;
            if ($x[$mid] < $v) {
                $lo = $mid + 1;
            } else {
                $hi = $mid;
            }
        }
        $pos = $lo;
        if ($x[$pos] == $v) {
            return $y[$pos];
        }
        $a = $pos - 1;
        $b = $pos;
        if (! $constant) {
            // y0 + (y1 - y0) * t with the multiply-add fused, as the extension does.
            return Fma::fma($y[$b] - $y[$a], fdiv($v - $x[$a], $x[$b] - $x[$a]), $y[$a]);
        }
        $leftWeight = 1.0 - $f;
        if ($leftWeight != 0.0 && $f != 0.0) {
            return $y[$a] * $leftWeight + $y[$b] * $f;
        }
        if ($leftWeight != 0.0) {
            return $y[$a] * $leftWeight;
        }
        if ($f != 0.0) {
            return $y[$b] * $f;
        }

        return 0.0;
    }

    /**
     * The standard normal density φ(z) = exp(-z²/2) / sqrt(2π), with z²/2
     * computed so that its rounding error does not reach the exponential.
     *
     * @internal
     */
    public static function standardNormalDensity(float $z): float
    {
        if (is_nan($z)) {
            return NAN;
        }
        $a = abs($z);
        // Beyond 40 the true value (below 1e-347) rounds to zero; this also covers infinities.
        if ($a > 40.0) {
            return 0.0;
        }
        // Split a = head + tail with head on a 2^-16 grid, so head has at most 22 significant bits and head² is exact.
        $scale = 65536.0;
        $head = floor($a * $scale) / $scale;
        $tail = $a - $head; // exact
        $halfHeadSq = 0.5 * ($head * $head); // exact
        $rest = $head * $tail + 0.5 * ($tail * $tail);
        $smallFactor = exp(-$rest);
        if ($halfHeadSq < 700.0) {
            return exp(-$halfHeadSq) * $smallFactor * self::RECIP_SQRT_TAU;
        }
        // Deep tail: keep every intermediate in the normal range, then scale
        // down by an exact power of two so the only rounding into the
        // subnormal range happens once, at the end.
        $shift = 600.0;
        $base = exp(-$shift) * 2.0 ** 800;

        return ($base * exp(-($halfHeadSq - $shift)) * $smallFactor * self::RECIP_SQRT_TAU) * 2.0 ** -800;
    }
}
