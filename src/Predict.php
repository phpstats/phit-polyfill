<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Distribution\Normal;
use Phit\Polyfill\Distribution\StudentT;
use Phit\Polyfill\Formula\Data;
use Phit\Polyfill\Formula\Design;
use Phit\Polyfill\Formula\Frontend;

/**
 * Phit\predict: predictions for new data from a formula fit, as R's
 * predict() does for lm, glm, and merMod objects (the extension's
 * formula.rs). The fit's formula is rebuilt over the new data with the
 * fitted factor levels and poly() bases from its predvars.
 */
final class Predict
{
    /** Phit\predict */
    public static function predict(
        mixed $fit,
        mixed $newdata,
        string $type = 'link',
        bool $seFit = false,
        string $interval = 'none',
        float $level = 0.95,
        bool $random = true,
        bool $allowNewLevels = false,
    ): array {
        Frontend::stringArgument($type, 'type');
        Frontend::stringArgument($interval, 'interval');
        $fit = Objects\Wrap::data($fit) ?? throw self::malformed('not a fit object or array');
        $formula = $fit['structural_formula'] ?? $fit['formula'] ?? null;
        if (! is_string($formula)) {
            throw self::malformed('missing formula');
        }
        $mixed = array_key_exists('random', $fit);
        $ridge = ($fit['model'] ?? null) === 'ridge';
        if ($ridge && $interval === 'prediction') throw ErrorCode::exception(1007, 'ridge prediction intervals are unavailable');
        $link = $fit['link'] ?? null;
        if (is_string($link)) {
            if (! in_array($link, ['identity', 'logit', 'probit', 'log'], true)) {
                throw self::malformed('unknown link');
            }
        } else {
            $link = 'identity';
        }
        $isLm = array_key_exists('family', $fit) === false && $mixed === false;
        $response = match ($type) {
            'link' => false,
            'response' => true,
            default => throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "type must be 'link' or 'response'"),
        };
        $intervalKind = match (true) {
            $interval === 'none' => null,
            ($interval === 'confidence' || $interval === 'prediction') && $isLm => $interval === 'prediction',
            $interval === 'confidence' || $interval === 'prediction' => throw ErrorCode::exception(
                ErrorCode::UNSUPPORTED,
                'intervals are available for Phit\fit_ols fits; use seFit for GLMs',
            ),
            default => throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "interval must be 'none', 'confidence' or 'prediction'"),
        };
        if (is_nan($level) || $level <= 0.0 || $level >= 1.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'level must be strictly between 0 and 1');
        }
        $uncertainty = $seFit || $intervalKind !== null;
        if ($mixed && $uncertainty) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'standard errors and intervals are not available for mixed-model predictions');
        }
        $coefficients = self::floatList($fit['coefficients'] ?? null, 'coefficients');
        $columns = self::stringList($fit['columns'] ?? null, 'columns');
        $rows = $fit['covariance'] ?? null;
        if (! is_array($rows)) {
            throw self::malformed('missing covariance');
        }
        $covariance = [];
        foreach ($rows as $row) {
            $covariance[] = self::floatList($ridge ? array_map(fn ($v) => $v ?? 0.0, $row) : $row, 'covariance');
        }
        $p = count($columns);
        $shapes = count($coefficients) === $p && count($covariance) === $p;
        foreach ($covariance as $row) {
            $shapes = $shapes && count($row) === $p;
        }
        if ($shapes === false) {
            throw self::malformed('coefficients, columns and covariance disagree');
        }
        $predvars = self::predvars($fit);
        $frame = Data::from($newdata);
        $design = Design::prediction($formula, $frame, $predvars, $mixed && $random);
        if ($design->columns !== $columns) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, sprintf(
                'newdata produce columns [%s] but the fit has [%s]',
                implode(', ', $design->columns),
                implode(', ', $columns),
            ));
        }
        $effects = $mixed && $random ? self::fittedEffects($fit) : null;
        $n = $frame->rows;
        $eta = array_fill(0, $n, null);
        $se = array_fill(0, $n, null);
        $newLevels = [];
        foreach ($design->retainedRows as $k => $row) {
            $x = $design->x[$k];
            $value = -0.0; // Rust's f64 sum starts from -0.0
            foreach ($x as $j => $a) {
                if (array_key_exists($j, $coefficients)) {
                    $value += $a * $coefficients[$j];
                }
            }
            if ($design->offset !== null) {
                $value += $design->offset[$k];
            }
            if ($effects !== null) {
                foreach ($design->random as $t => $term) {
                    if (! array_key_exists($t, $effects)) {
                        break;
                    }
                    $key = (string) $term['groups'][$k];
                    $b = $effects[$t][$key] ?? null;
                    if ($b !== null) {
                        $sum = -0.0;
                        foreach ($term['z'][$k] as $j => $z) {
                            if (array_key_exists($j, $b)) {
                                $sum += $z * $b[$j];
                            }
                        }
                        $value += $sum;
                    } elseif ($allowNewLevels === false) {
                        $label = "{$term['group']}=$key";
                        if (! in_array($label, $newLevels, true)) {
                            $newLevels[] = $label;
                        }
                    }
                }
            }
            $eta[$row] = $value;
            if ($uncertainty) {
                $quadratic = self::quadForm($x, $covariance);
                // f64::max: a NaN form becomes 0.
                $unavailable = false;
                if ($ridge) foreach ($fit['at_bound'] as $j=>$bound) if ($bound && $x[$j] != 0.0) $unavailable = true;
                $se[$row] = $unavailable ? null : sqrt(is_nan($quadratic) || $quadratic < 0.0 ? 0.0 : $quadratic);
            }
        }
        if ($newLevels !== []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'new levels detected in newdata: '.implode(', ', $newLevels).'; pass allowNewLevels: true to predict them at the population level');
        }
        $fitted = $eta;
        if ($response) {
            foreach ($fitted as $i => $e) {
                if ($e !== null) {
                    $fitted[$i] = self::inverseLink($link, $e);
                    if ($se[$i] !== null) {
                        $se[$i] = $se[$i] * abs(self::derivative($link, $e));
                    }
                }
            }
        }
        $out = ['fit' => $fitted];
        if ($seFit) {
            $out['se_fit'] = $se;
        }
        $dfValue=$fit['inference_df']??$fit['df_residual']??null;
        $df = is_numeric($dfValue) ? $dfValue : null;
        if ($isLm && $uncertainty) {
            $out['df'] = $df;
        }
        if ($mixed === false && $uncertainty) {
            // R's residual.scale: sigma for Gaussian fits, 1 for binomial and Poisson.
            $out['residual_scale'] = $link === 'identity'
                ? (is_float($fit['sigma'] ?? null) ? $fit['sigma'] : null)
                : 1.0;
        }
        if ($intervalKind !== null) {
            if ($df === null) {
                throw self::malformed('missing df_residual');
            }
            $sigma = $fit['sigma'] ?? null;
            if (! is_float($sigma)) {
                throw self::malformed('missing sigma');
            }
            $critical = StudentT::create((float) $df)->quantileAt((1.0 + $level) / 2.0);
            if ($intervalKind) {
                try {
                    $weights = self::floatList($fit['weights'] ?? null, 'weights');
                } catch (\Exception) {
                    $weights = [];
                }
                foreach ($weights as $w) {
                    if ($w != 1.0) {
                        // Userland code cannot raise E_WARNING; the text matches the extension's.
                        trigger_error('Phit\predict(): Phit\predict assumes constant prediction variance even though the model fit is weighted', E_USER_WARNING);
                        break;
                    }
                }
            }
            $lower = array_fill(0, $n, null);
            $upper = array_fill(0, $n, null);
            for ($i = 0; $i < $n; $i++) {
                if ($fitted[$i] !== null && $se[$i] !== null) {
                    $s = $se[$i];
                    $width = $critical * ($intervalKind ? sqrt($s * $s + $sigma * $sigma) : $s);
                    $lower[$i] = $fitted[$i] - $width;
                    $upper[$i] = $fitted[$i] + $width;
                }
            }
            $out['lower'] = $lower;
            $out['upper'] = $upper;
        }
        $out['level'] = $intervalKind !== null ? $level : null;
        $out['n_input'] = $n;
        $out['n_predicted'] = count($design->retainedRows);
        $out['missing_rows'] = $design->droppedRows;

        return $out;
    }

    private static function malformed(string $detail): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, "fit must be a result of Phit\\fit_ols, Phit\\fit_generalized_linear, Phit\\fit_mixed_gaussian or Phit\\fit_mixed_logistic: $detail");
    }

    /** @return list<float> */
    private static function floatList(mixed $value, string $what): array
    {
        if (! is_array($value)) {
            throw self::malformed("missing $what");
        }
        $out = [];
        foreach ($value as $v) {
            if (is_float($v)) {
                $out[] = $v;
            } elseif (is_int($v)) {
                $out[] = (float) $v;
            } else {
                throw self::malformed("$what must be numbers");
            }
        }

        return $out;
    }

    /** @return list<string> */
    private static function stringList(mixed $value, string $what): array
    {
        if (! is_array($value)) {
            throw self::malformed("missing $what");
        }
        $out = [];
        foreach ($value as $v) {
            if (! is_string($v)) {
                throw self::malformed("$what must be strings");
            }
            $out[] = $v;
        }

        return $out;
    }

    /** The fit's factor levels and poly() coefficients. */
    private static function predvars(array $fit): array
    {
        $table = $fit['predvars'] ?? null;
        if (! is_array($table)) {
            throw self::malformed('missing predvars');
        }
        $predvars = ['levels' => [], 'poly' => []];
        if (is_array($table['levels'] ?? null)) {
            foreach ($table['levels'] as $name => $values) {
                $predvars['levels'][(string) $name] = self::stringList($values, 'levels');
            }
        }
        if (is_array($table['poly'] ?? null)) {
            foreach ($table['poly'] as $text => $coefs) {
                if (! is_array($coefs)) {
                    throw self::malformed('poly coefficients must be arrays');
                }
                $predvars['poly'][(string) $text] = [
                    'alpha' => self::floatList($coefs['alpha'] ?? null, 'poly alpha'),
                    'norm2' => self::floatList($coefs['norm2'] ?? null, 'poly norm2'),
                ];
            }
        }

        return $predvars;
    }

    private static function labelKey(mixed $value): ?string
    {
        return is_int($value) || is_string($value) ? (string) $value : null;
    }

    /**
     * Fitted random effects by term: group label => coefficient vector.
     *
     * @return list<array<string, list<float>>>
     */
    private static function fittedEffects(array $fit): array
    {
        $malformed = fn () => self::malformed('random effects are incomplete');
        if (! is_array($fit['random'] ?? null)) {
            throw $malformed();
        }
        if (is_array($fit['random_effects'] ?? null)) {
            // Phit\fit_mixed_logistic: one random intercept with top-level effects and labels.
            $labels = $fit['group_labels'] ?? null;
            if (! is_array($labels)) {
                throw $malformed();
            }
            $map = [];
            $effects = array_values($fit['random_effects']);
            foreach (array_values($labels) as $i => $label) {
                if (! array_key_exists($i, $effects)) {
                    break;
                }
                if (! is_float($effects[$i])) {
                    throw $malformed();
                }
                $key = self::labelKey($label) ?? throw $malformed();
                $map[$key] = [$effects[$i]];
            }

            return [$map];
        }
        $out = [];
        foreach ($fit['random'] as $term) {
            if (! is_array($term) || ! is_array($term['group_labels'] ?? null) || ! is_array($term['effects'] ?? null)) {
                throw $malformed();
            }
            $map = [];
            $effects = array_values($term['effects']);
            foreach (array_values($term['group_labels']) as $i => $label) {
                if (! array_key_exists($i, $effects)) {
                    break;
                }
                $key = self::labelKey($label) ?? throw $malformed();
                $map[$key] = self::floatList($effects[$i], 'random effects');
            }
            $out[] = $map;
        }

        return $out;
    }

    private static function quadForm(array $x, array $covariance): float
    {
        $total = 0.0;
        foreach ($x as $i => $xi) {
            foreach ($x as $j => $xj) {
                $total += $xi * $covariance[$i][$j] * $xj;
            }
        }

        return $total;
    }

    private static function inverseLink(string $link, float $eta): float
    {
        switch ($link) {
            case 'logit':
                if ($eta >= 0.0) {
                    return 1.0 / (1.0 + exp(-$eta));
                }
                $e = exp($eta);

                return $e / (1.0 + $e);
            case 'probit':
                return Normal::create(0.0, 1.0)->cdfAt($eta);
            case 'log':
                return exp($eta);
            default:
                return $eta;
        }
    }

    private static function derivative(string $link, float $eta): float
    {
        switch ($link) {
            case 'logit':
                $mu = self::inverseLink($link, $eta);

                return $mu * (1.0 - $mu);
            case 'probit':
                return Normal::create(0.0, 1.0)->pdfAt($eta);
            case 'log':
                return exp($eta);
            default:
                return 1.0;
        }
    }
}
