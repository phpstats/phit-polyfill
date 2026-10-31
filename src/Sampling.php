<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Distribution\Support;
use Phit\Polyfill\Internal\ChaCha20;

/**
 * Uniform, categorical and permutation sampling, ported from the phit
 * extension (stats-core distributions/rng.rs plus the PHP glue). Every draw
 * takes one Open01 uniform from the same freshly seeded ChaCha20 stream as
 * the Phit\dist_*_rng families, and the checks run in the extension's order.
 */
final class Sampling
{
    /** Phit\dist_uniform_rng: R's runif(), min + (max - min) * u. */
    public static function uniform(int $n, int $seed, float $min = 0.0, float $max = 1.0): array
    {
        self::draws($n, $seed);
        Support::finite('min', $min);
        Support::finite('max', $max);
        if ($min > $max) {
            throw Support::invalid('min must not exceed max');
        }

        return Support::sample($n, $seed, static fn (float $u): float => $min + ($max - $min) * $u);
    }

    /**
     * Phit\dist_categorical_rng: n 0-based indices into prob, drawn with
     * those weights, as R's sample(length(prob), n, replace = TRUE, prob).
     *
     * @return list<int>
     */
    public static function categorical(int $n, int $seed, array $prob): array
    {
        self::draws($n, $seed);
        $weights = self::weights($prob, 'prob');
        self::count($n, self::TOO_MANY);
        $cumulative = self::cumulative($weights, 'prob');
        $rng = new ChaCha20($seed);
        $out = [];
        for ($i = 0; $i < $n; $i++) {
            $out[] = self::category($cumulative, $rng->open01());
        }

        return $out;
    }

    /**
     * Phit\dist_categorical_rows_rng: one 0-based index per row of prob,
     * drawn with that row's weights; row i takes the i-th uniform.
     *
     * @return list<int>
     */
    public static function categoricalRows(int $seed, array $prob): array
    {
        self::draws(0, $seed);
        if (! array_is_list($prob)) {
            throw Support::invalid('prob must be a list of rows');
        }
        $rows = [];
        foreach ($prob as $i => $row) {
            if (! \is_array($row)) {
                throw Support::invalid('prob must be a list of rows');
            }
            $rows[] = self::weights($row, "prob row $i");
        }
        self::count(\count($rows), 'prob must not have more than 1,000,000 rows');
        foreach ($rows as $i => $row) {
            $rows[$i] = self::cumulative($row, "prob row $i");
        }
        $rng = new ChaCha20($seed);
        $out = [];
        foreach ($rows as $row) {
            $out[] = self::category($row, $rng->open01());
        }

        return $out;
    }

    /**
     * Phit\permutation_rng: size distinct 0-based indices from 0..n-1 in
     * random order, as R's sample(n, size), by a partial Fisher-Yates shuffle.
     *
     * @return list<int>
     */
    public static function permutation(int $n, int $seed, ?int $size = null): array
    {
        self::draws($n, $seed);
        $size ??= $n;
        if ($size < 0) {
            throw Support::invalid('size must be nonnegative');
        }
        self::count($n, self::TOO_MANY);
        if ($size > $n) {
            throw Support::invalid('size must not exceed n');
        }
        $out = $n === 0 ? [] : range(0, $n - 1);
        $rng = new ChaCha20($seed);
        for ($i = 0; $i < $size; $i++) {
            $remaining = $n - $i;
            $j = $i + min((int) ($rng->open01() * $remaining), $remaining - 1);
            [$out[$i], $out[$j]] = [$out[$j], $out[$i]];
        }

        return \array_slice($out, 0, $size);
    }

    private const string TOO_MANY = 'n must not exceed 1,000,000 samples per call';

    private static function draws(int $n, int $seed): void
    {
        if ($n < 0) {
            throw Support::invalid('n must be nonnegative');
        }
        if ($seed < 0) {
            throw Support::invalid('seed must be nonnegative');
        }
    }

    private static function count(int $n, string $message): void
    {
        if ($n > Support::MAX_RNG_DRAWS) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, $message);
        }
    }

    /**
     * A list of ints or floats, since draws are positions in it.
     *
     * @return list<float>
     */
    private static function weights(array $values, string $name): array
    {
        if (! array_is_list($values)) {
            throw Support::invalid("$name must be a list");
        }
        $out = [];
        foreach ($values as $value) {
            if (\is_float($value)) {
                $out[] = $value;
            } elseif (\is_int($value)) {
                $out[] = (float) $value;
            } else {
                throw Support::invalid('array elements must be integers or floats');
            }
        }

        return $out;
    }

    /**
     * Running sums of the weights, validated as R's sample() validates prob:
     * finite, nonnegative and with a positive sum.
     *
     * @param list<float> $prob
     * @return non-empty-list<float>
     */
    private static function cumulative(array $prob, string $name): array
    {
        if ($prob === []) {
            throw Support::invalid("$name must not be empty");
        }
        $total = 0.0;
        $out = [];
        foreach ($prob as $w) {
            if (! is_finite($w)) {
                throw ErrorCode::exception(ErrorCode::NONFINITE, "$name probabilities must be finite");
            }
            if ($w < 0.0) {
                throw Support::invalid("$name probabilities must be nonnegative");
            }
            $total += $w;
            $out[] = $total;
        }
        if (! is_finite($total)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, "$name probabilities must have a finite sum");
        }
        if ($total <= 0.0) {
            throw Support::invalid("$name probabilities must have a positive sum");
        }

        return $out;
    }

    /**
     * The first k with u * total < cumulative[k]. If rounding leaves no such
     * k, the last category with positive weight, never a zero weight.
     *
     * @param non-empty-list<float> $cumulative
     */
    private static function category(array $cumulative, float $u): int
    {
        $last = \count($cumulative) - 1;
        $total = $cumulative[$last];
        $target = $u * $total;
        foreach ($cumulative as $k => $c) {
            if ($target < $c) {
                return $k;
            }
        }
        $k = $last;
        while ($k > 0 && $cumulative[$k - 1] == $total) {
            $k--;
        }

        return $k;
    }
}
