<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\Distribution\Support;
use Phit\Polyfill\ErrorCode;

/**
 * Validation and listwise deletion shared by the regression functions, in the
 * same order as the extension so the same input raises the same code.
 *
 * @internal
 */
final class ModelData
{
    /**
     * @param  list<list<float>>  $x  complete rows, intercept column prepended
     * @param  list<float>  $y
     * @param  list<float>  $weights
     * @param  list<float>  $offset
     * @param  list<int>  $retained
     * @param  list<int>  $dropped
     * @param  list<int|string>  $groups  group labels of the complete rows, when groups were given
     */
    private function __construct(
        public array $x,
        public array $y,
        public array $weights,
        public array $offset,
        public int $nInput,
        public array $retained,
        public array $dropped,
        public array $groups = [],
    ) {}

    /** Sample accounting only, for results whose rows were prepared elsewhere (formulas). */
    public static function accounting(int $nInput, array $retained, array $dropped): self
    {
        return new self([], [], [], [], $nInput, $retained, $dropped);
    }

    /** A PHP value as a number; null and NaN are missing, INF and non-numbers are errors. */
    public static function number(mixed $value): ?float
    {
        if ($value === null) {
            return null;
        }
        if (is_int($value)) {
            return (float) $value;
        }
        if (! is_float($value)) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'inputs must contain integers, floats, or null');
        }
        if (is_nan($value)) {
            return null;
        }
        if (is_infinite($value)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'infinite input values are not allowed');
        }

        return $value;
    }

    /** @return list<?float> */
    private static function vector(?array $values, int $n, float $default): array
    {
        if ($values === null) {
            return array_fill(0, $n, $default);
        }
        $out = [];
        foreach ($values as $value) {
            $out[] = self::number($value);
        }
        if (count($out) !== $n) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'optional vectors must have one value per original observation');
        }

        return $out;
    }

    /**
     * Group identifiers: integers or strings, with null marking a missing label.
     *
     * @return list<int|string|null>
     */
    public static function labels(array $values): array
    {
        $out = [];
        foreach ($values as $value) {
            if ($value === null || is_int($value) || is_string($value)) {
                $out[] = $value;
            } else {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'group IDs must be integers, strings, or null');
            }
        }

        return $out;
    }

    /**
     * Dense indices of group labels in order of first appearance. Integer 5 and
     * string "5" are different groups, as in the extension.
     *
     * @param  list<int|string>  $labels
     * @return array{0: list<int>, 1: list<int|string>} indices and distinct labels
     */
    public static function groupIndices(array $labels): array
    {
        $seen = [];
        $distinct = [];
        $indices = [];
        foreach ($labels as $label) {
            $key = (is_int($label) ? 'i:' : 's:').$label;
            if (! isset($seen[$key])) {
                $seen[$key] = count($distinct);
                $distinct[] = $label;
            }
            $indices[] = $seen[$key];
        }

        return [$indices, $distinct];
    }

    /** @param 'continuous'|'binary'|'count' $outcome */
    public static function prepare(
        array $x,
        array $y,
        string $missing,
        string $outcome,
        ?array $weights,
        ?array $offset,
        bool $intercept,
        ?array $groups = null,
    ): self {
        if ($missing !== 'drop' && $missing !== 'raise') {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "missing must be 'drop' or 'raise'");
        }
        // Rows are converted with the intercept already in front; null marks a missing value.
        $rows = [];
        $complete = [];
        foreach ($x as $row) {
            if (! is_array($row)) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each x row must be an array');
            }
            $values = $intercept ? [1.0] : [];
            $full = true;
            foreach ($row as $value) {
                // Fast path for finite floats; everything else goes through number().
                if (is_float($value) && is_finite($value)) {
                    $values[] = $value;
                } else {
                    $value = self::number($value);
                    $values[] = $value;
                    $full = $full && $value !== null;
                }
            }
            $rows[] = $values;
            $complete[] = $full;
        }
        $outcomes = [];
        foreach ($y as $value) {
            $value = self::number($value);
            if ($value !== null && $outcome === 'binary' && (is_nan($value) || $value < 0.0 || $value > 1.0)) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binary/grouped-binomial y must be in [0,1]');
            }
            if ($value !== null && $outcome === 'count' && ($value < 0.0 || fmod($value, 1.0) != 0.0)) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'Poisson y must contain nonnegative integer counts');
            }
            $outcomes[] = $value;
        }
        $labels = $groups === null ? null : self::labels($groups);
        $n = count($rows);
        if ($n === 0 || (! $intercept && count($rows[0]) === 0)) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'x must have at least one row; without an intercept, at least one column is required');
        }
        if ($labels !== null && count($labels) !== $n) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'groups must have one identifier per observation');
        }
        if (count($outcomes) !== $n) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'x and y must have the same number of observations');
        }
        $width = count($rows[0]);  // includes the intercept column, the same for every row
        foreach ($rows as $row) {
            if (count($row) !== $width) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'x must be rectangular');
            }
        }
        $w = self::vector($weights, $n, 1.0);
        $o = self::vector($offset, $n, 0.0);
        $binary = $outcome === 'binary';
        foreach ($w as $value) {
            if ($value !== null && ($value <= 0.0 || ($binary && fmod($value, 1.0) != 0.0))) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'weights must be positive; binary weights must be integer trial counts');
            }
        }
        if ($binary) {
            foreach ($w as $value) {
                if ($value !== null && $value > Limits::MAX_TRIALS) {
                    throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'binomial trial weights must not exceed 1000000');
                }
            }
            foreach ($outcomes as $i => $value) {
                if ($value !== null && $w[$i] !== null) {
                    $events = $value * $w[$i];
                    if (abs($events - Support::round($events)) > 1e-8) {
                        throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binary/grouped-binomial outcomes require integer trial counts and successes');
                    }
                }
            }
        }
        $cx = $cy = $cw = $co = $cg = $retained = $dropped = [];
        foreach ($rows as $i => $row) {
            if (! $complete[$i] || $outcomes[$i] === null || ($labels !== null && $labels[$i] === null) || $w[$i] === null || $o[$i] === null) {
                $dropped[] = $i;

                continue;
            }
            if ($labels !== null) {
                $cg[] = $labels[$i];
            }
            $cx[] = $row;
            $cy[] = $outcomes[$i];
            $cw[] = $w[$i];
            $co[] = $o[$i];
            $retained[] = $i;
        }
        if ($missing === 'raise' && $dropped !== []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "missing values found with missing='raise'");
        }

        return new self($cx, $cy, $cw, $co, $n, $retained, $dropped, $cg);
    }

    /**
     * Add diagnostics and sample accounting, and raise the listwise-deletion warning.
     *
     * @param  array<string, mixed>  $output
     * @return array<string, mixed>
     */
    public function finish(array $output, string $function, string $solver): array
    {
        $diagnostics = ['converged' => true, 'solver' => $solver];
        foreach (['iterations', 'evaluations', 'gradient_max', 'boundary', 'singular'] as $key) {
            if (array_key_exists($key, $output)) {
                $diagnostics[$key] = $output[$key];
            }
        }
        $warnings = [];
        $dropped = count($this->dropped);
        if ($dropped > 0) {
            $warnings[] = "listwise deletion removed $dropped of {$this->nInput} original rows";
        }
        if (($output['singular'] ?? null) === true) {
            $warnings[] = 'random-effect covariance is singular or near singular';
        }
        if (in_array($function, ['Phit\fit_mixed_logistic_matrix', 'Phit\fit_mixed_logistic'], true) && ($output['quadrature_points'] ?? null) === 1) {
            $warnings[] = 'Laplace integration is approximate; assess quadrature sensitivity for rare events or small groups';
        }
        $diagnostics['warnings'] = $warnings;
        $output['diagnostics'] = $diagnostics;
        $output['n_input'] = $this->nInput;
        $output['n_observations'] = count($this->retained);
        $output['n_dropped'] = $dropped;
        $output['retained_rows'] = $this->retained;
        $output['dropped_rows'] = $this->dropped;
        if ($dropped > 0) {
            // Userland code cannot raise E_WARNING itself; E_USER_WARNING is the closest
            // equivalent. The "fn(): " prefix matches what PHP adds to native warnings:
            // the name the caller used, which differs through an alias.
            $caller = CalledAs::name() ?? $function;
            trigger_error(
                "{$caller}(): dropped $dropped of {$this->nInput} observations with missing values (listwise deletion)",
                E_USER_WARNING,
            );
        }

        return $output;
    }
}
