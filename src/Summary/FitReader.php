<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

use Phit\Polyfill\ErrorCode;

/**
 * Reads a fit array as the extension's printer.rs reads the PHP hash table:
 * a float key may hold an int, an int key may hold an integral float, a
 * string must be valid UTF-8, and a null list element reads as NaN.
 */
final class FitReader
{
    public static function notAFit(string $detail): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, "fit must be a result of a Phit model function: $detail");
    }

    public static function float(array $fit, string $key): ?float
    {
        $v = $fit[$key] ?? null;

        return is_float($v) || is_int($v) ? (float) $v : null;
    }

    public static function requiredFloat(array $fit, string $key): float
    {
        return self::float($fit, $key) ?? throw self::notAFit("missing $key");
    }

    public static function long(array $fit, string $key): ?int
    {
        $v = $fit[$key] ?? null;
        if (is_int($v)) {
            return $v;
        }

        return is_float($v) && is_finite($v) && fmod($v, 1.0) == 0.0 ? (int) $v : null;
    }

    public static function text(array $fit, string $key): ?string
    {
        return self::utf8($fit[$key] ?? null);
    }

    public static function dfResidual(array $fit): int
    {
        return self::long($fit, 'df_residual') ?? throw self::notAFit('missing df_residual');
    }

    public static function nObservations(array $fit): int
    {
        return self::long($fit, 'n_observations') ?? throw self::notAFit('missing n_observations');
    }

    public static function intercept(array $fit): ?bool
    {
        $v = $fit['intercept'] ?? null;

        return is_bool($v) ? $v : null;
    }

    /** The coefficient table of any fit: names from `columns`, or x1, x2, ... for a matrix fit. */
    public static function coefficients(array $fit): Coefficients
    {
        $estimates = self::floats($fit, 'coefficients');
        if (array_key_exists('t_statistics', $fit)) {
            [$statistics, $statistic] = [self::floats($fit, 't_statistics'), 't'];
        } else {
            [$statistics, $statistic] = [self::floats($fit, 'z_statistics'), 'z'];
        }
        $names = self::strings($fit, 'columns') ?? self::numbered('x', count($estimates));
        $stdErrors = self::floats($fit, 'std_errors');
        $pValues = self::floats($fit, 'p_values');
        $n = count($estimates);
        if (count($names) !== $n || count($stdErrors) !== $n || count($statistics) !== $n || count($pValues) !== $n) {
            throw self::notAFit('coefficient lists differ in length');
        }

        return new Coefficients($names, $estimates, $stdErrors, $statistics, $pValues, $statistic);
    }

    /** @return list<RandomTerm> */
    public static function randomTerms(array $fit): array
    {
        $oneIntercept = fn (string $group): array => [new RandomTerm(
            $group,
            ['(Intercept)'],
            [[self::requiredFloat($fit, 'random_intercept_variance')]],
            false,
            self::long($fit, 'n_groups') ?? throw self::notAFit('missing n_groups'),
            self::randomValues($fit, 'group_labels', 'random_effects', 'random_effects_cond_var', 1, true),
        )];
        $random = $fit['random'] ?? null;
        if (! is_array($random)) {
            return $oneIntercept('group');
        }
        $terms = [];
        $i = 0;
        foreach ($random as $term) {
            if (! is_array($term)) {
                throw self::notAFit('random terms must be arrays');
            }
            $group = self::text($term, 'group') ?? 'g'.($i + 1);
            $i++;
            if (! array_key_exists('covariance', $term)) {
                // A binomial model's single (1 | g) term keeps its variance at the top level.
                return $oneIntercept($group);
            }
            $covariance = self::matrix($term['covariance'], 'covariance');
            $k = count($covariance);
            if ($k === 0 || array_filter($covariance, fn (array $row): bool => count($row) !== $k) !== []) {
                throw self::notAFit("a random term's covariance must be square");
            }
            $columns = self::strings($term, 'columns') ?? self::numbered('z', $k);
            if (count($columns) !== $k) {
                throw self::notAFit("a random term's columns and covariance differ in size");
            }
            $labels = $term['group_labels'] ?? null;
            $correlated = $term['correlated'] ?? null;
            $terms[] = new RandomTerm(
                $group,
                $columns,
                $covariance,
                is_bool($correlated) ? $correlated : true,
                is_array($labels) ? count($labels) : 0,
                self::randomValues($term, 'group_labels', 'effects', 'effects_cond_var', $k, false),
            );
        }

        return $terms;
    }

    /**
     * A term's conditional modes and standard errors, when the fit has the labels, the effects, and their
     * variances: `effects` rows with `effects_cond_var` k×k matrices, or ($flat) one value per group with one
     * variance per group.
     */
    private static function randomValues(array $source, string $labelsKey, string $effectsKey, string $condVarKey, int $k, bool $flat): ?RandomValues
    {
        $levels = self::labels($source[$labelsKey] ?? null);
        if ($levels === null || ! array_key_exists($effectsKey, $source) || ! array_key_exists($condVarKey, $source)) {
            return null;
        }
        if ($flat) {
            $e = self::floats($source, $effectsKey);
            $v = self::floats($source, $condVarKey);
            if (count($v) !== count($e)) {
                throw self::notAFit('random effects and their variances differ in length');
            }
            $estimates = array_map(fn (float $x): array => [$x], $e);
            $stdErrors = array_map(fn (float $x): array => [sqrt($x)], $v);
        } else {
            $estimates = self::matrix($source[$effectsKey], $effectsKey);
            $matrices = $source[$condVarKey];
            if (! is_array($matrices)) {
                throw self::notAFit("missing $condVarKey");
            }
            $v = [];
            foreach ($matrices as $m) {
                $v[] = self::matrix($m, $condVarKey);
            }
            $square = fn (array $m): bool => count($m) === $k && array_filter($m, fn (array $r): bool => count($r) !== $k) === [];
            if (count($v) !== count($estimates) || array_filter($v, fn (array $m): bool => ! $square($m)) !== []) {
                throw self::notAFit('random effects and their variances differ in size');
            }
            $stdErrors = [];
            foreach ($v as $m) {
                $se = [];
                for ($j = 0; $j < $k; $j++) {
                    $se[] = sqrt($m[$j][$j]);
                }
                $stdErrors[] = $se;
            }
        }
        if (count($estimates) !== count($levels) || array_filter($estimates, fn (array $r): bool => count($r) !== $k) !== []) {
            throw self::notAFit('random effects and group labels differ in size');
        }

        return new RandomValues($levels, $estimates, $stdErrors);
    }

    /**
     * Group labels as text, as lme4 prints integer and string labels alike: a float as Rust's Display writes it.
     *
     * @return ?list<string>
     */
    private static function labels(mixed $values): ?array
    {
        if (! is_array($values)) {
            return null;
        }
        $out = [];
        foreach ($values as $v) {
            $text = match (true) {
                is_int($v) => (string) $v,
                is_float($v) => self::rustDisplay($v),
                default => self::utf8($v),
            };
            if ($text === null) {
                return null;
            }
            $out[] = $text;
        }

        return $out;
    }

    /** Rust's Display for f64: the shortest round-trip digits in plain decimal notation, never an exponent. */
    private static function rustDisplay(float $x): string
    {
        if (is_nan($x)) {
            return 'NaN';
        }
        if (is_infinite($x)) {
            return $x > 0 ? 'inf' : '-inf';
        }
        $sign = $x < 0.0 || ($x == 0.0 && fdiv(1.0, $x) < 0.0) ? '-' : '';
        if ($x == 0.0) {
            return $sign.'0';
        }
        // The shortest significant digits that parse back to |x|, and the decimal exponent of the first one.
        for ($precision = 0; $precision < 17; $precision++) {
            $text = sprintf("%.{$precision}e", abs($x));
            if ((float) $text === abs($x)) {
                break;
            }
        }
        [$mantissa, $exponent] = explode('e', $text);
        $digits = str_replace('.', '', $mantissa);
        $point = (int) $exponent + 1;
        // sprintf breaks an exact tie to even; Rust's shortest form takes the larger neighbour, which also parses
        // back to |x| (936733262222091.25 prints as 936733262222091.3).
        [$exactMantissa, $exactExponent] = explode('e', sprintf('%.52e', abs($x)));
        $exact = str_replace('.', '', $exactMantissa);
        $kept = strlen($digits);
        if ((int) $exactExponent + 1 === $point && rtrim(substr($exact, $kept), '0') === '5') {
            $digits = substr($exact, 0, $kept);
            for ($i = $kept - 1; $i >= 0 && $digits[$i] === '9'; $i--) {
                $digits[$i] = '0';
            }
            if ($i < 0) {
                $digits = '1'.$digits;
                $point++;
            } else {
                $digits[$i] = (string) ((int) $digits[$i] + 1);
            }
        }
        $digits = rtrim($digits, '0');
        if ($point <= 0) {
            return $sign.'0.'.str_repeat('0', -$point).$digits;
        }
        if ($point >= strlen($digits)) {
            return $sign.$digits.str_repeat('0', $point - strlen($digits));
        }

        return $sign.substr($digits, 0, $point).'.'.substr($digits, $point);
    }

    /** @return list<list<float>> */
    private static function matrix(mixed $rows, string $key): array
    {
        if (! is_array($rows)) {
            throw self::notAFit("missing $key");
        }
        $out = [];
        foreach ($rows as $row) {
            $out[] = self::floats([$key => $row], $key);
        }

        return $out;
    }

    private static function utf8(mixed $v): ?string
    {
        return is_string($v) && preg_match('//u', $v) === 1 ? $v : null;
    }

    /** @return list<float> */
    private static function floats(array $fit, string $key): array
    {
        $values = $fit[$key] ?? null;
        if (! is_array($values)) {
            throw self::notAFit("missing $key");
        }
        $out = [];
        foreach ($values as $v) {
            $out[] = match (true) {
                $v === null => NAN,
                is_float($v), is_int($v) => (float) $v,
                default => throw self::notAFit("$key must be numbers"),
            };
        }

        return $out;
    }

    /** @return ?list<string> */
    private static function strings(array $fit, string $key): ?array
    {
        $values = $fit[$key] ?? null;
        if (! is_array($values)) {
            return null;
        }
        $out = [];
        foreach ($values as $v) {
            $text = self::utf8($v);
            if ($text === null) {
                return null;
            }
            $out[] = $text;
        }

        return $out;
    }

    /** @return list<string> {$prefix}1 to {$prefix}{$n} */
    private static function numbered(string $prefix, int $n): array
    {
        $names = [];
        for ($i = 1; $i <= $n; $i++) {
            $names[] = $prefix.$i;
        }

        return $names;
    }
}
