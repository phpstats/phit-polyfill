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
            $rows = $term['covariance'];
            if (! is_array($rows)) {
                throw self::notAFit('missing covariance');
            }
            $covariance = [];
            foreach ($rows as $row) {
                $covariance[] = self::floats(['covariance' => $row], 'covariance');
            }
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
            );
        }

        return $terms;
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
