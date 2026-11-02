<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Internal\LinearMixed;
use Phit\Polyfill\Internal\LogisticRandomIntercept;
use Phit\Polyfill\Internal\Limits;
use Phit\Polyfill\Internal\ModelData;
use Phit\Polyfill\Internal\RandomIntercept;

/**
 * Mixed models: Phit\fit_mixed_gaussian_intercept_matrix, Phit\fit_mixed_gaussian_matrix, and Phit\fit_mixed_logistic_matrix (the
 * extension's lib.rs glue around stats-core mixed.rs, lmm.rs, and
 * logistic_mixed.rs).
 */
final class MixedModels
{
    /** Phit\fit_mixed_gaussian_intercept_matrix */
    public static function lmm(array $x, array $y, array $groups, string $missing = 'drop', bool $intercept = true): array
    {
        $data = ModelData::prepare($x, $y, $missing, 'continuous', null, null, $intercept, $groups);
        [$indices, $labels] = ModelData::groupIndices($data->groups);
        $result = RandomIntercept::fit($data->x, $data->y, $indices);
        $output = [
            'coefficients' => $result['coefficients'],
            'covariance' => $result['covariance'],
            'std_errors' => $result['std_errors'],
            'z_statistics' => $result['z_statistics'],
            'p_values' => $result['p_values'],
            'conf_int' => $result['conf_int'],
            'fitted' => $result['fitted'],
            'residuals' => $result['residuals'],
            'log_likelihood' => $result['log_likelihood'],
            'aic' => $result['aic'],
            'df_residual' => $result['df_residual'],
            'evaluations' => $result['evaluations'],
            'fixed_fitted' => $result['fixed_fitted'],
            'marginal_residuals' => $result['marginal_residuals'],
            'random_effects' => $result['random_effects'],
            'random_effects_cond_var' => $result['random_effects_cond_var'],
            'group_sizes' => $result['group_sizes'],
            'n_groups' => count($labels),
            'group_labels' => $labels,
            'residual_variance' => $result['residual_variance'],
            'random_intercept_variance' => $result['random_intercept_variance'],
            'sigma' => $result['sigma'],
            'icc' => $result['icc'],
            'boundary' => $result['boundary'],
            'method' => 'ML',
            'converged' => true,
        ];

        return $data->finish($output, 'Phit\fit_mixed_gaussian_intercept_matrix', 'profile_variance_search');
    }

    /** Phit\fit_mixed_gaussian_matrix */
    public static function lmmFit(
        array $x,
        array $y,
        array $random,
        string $missing = 'drop',
        bool $reml = false,
        ?array $weights = null,
        ?array $offset = null,
        ?array $start = null,
        int $maxEvaluations = 20000,
        float $tolerance = 1e-7,
        string $optimizer = 'bobyqa',
        float $objectiveTolerance = 1e-12,
        bool $intercept = true,
    ): array {
        if ($missing !== 'drop' && $missing !== 'raise') {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "missing must be 'drop' or 'raise'");
        }
        $rows = self::numericRows($x);
        $n = count($rows);
        $rectangular = true;
        foreach ($rows as $row) {
            $rectangular = $rectangular && count($row) === count($rows[0]);
        }
        if ($n === 0 || ($intercept === false && $rows[0] === []) || $rectangular === false) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'x must be nonempty and rectangular');
        }
        $outcomes = self::numericVector($y, $n, 0.0);
        $w = self::numericVector($weights, $n, 1.0);
        $o = self::numericVector($offset, $n, 0.0);
        foreach ($w as $value) {
            if ($value !== null && $value <= 0.0) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'weights must be positive');
            }
        }
        $terms = [];
        foreach ($random as $term) {
            if (! is_array($term)) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each random term must be an array');
            }
            if (! array_key_exists('z', $term) || ! is_array($term['z'])) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each random term needs z');
            }
            if (! array_key_exists('groups', $term) || ! is_array($term['groups'])) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each random term needs groups');
            }
            $z = self::numericRows($term['z']);
            $groups = ModelData::labels($term['groups']);
            $zRectangular = true;
            foreach ($z as $row) {
                $zRectangular = $zRectangular && count($row) === count($z[0]);
            }
            if (count($z) !== $n || count($groups) !== $n || $z[0] === [] || $zRectangular === false) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'each random z must be rectangular with one row and group per original observation');
            }
            $correlated = true;
            if (array_key_exists('correlated', $term)) {
                if (! is_bool($term['correlated'])) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'correlated must be a boolean');
                }
                $correlated = $term['correlated'];
            }
            $terms[] = [$z, $groups, $correlated];
        }
        if ($terms === []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'supply at least one random term');
        }
        $retained = [];
        $dropped = [];
        for ($i = 0; $i < $n; $i++) {
            $complete = $outcomes[$i] !== null && $w[$i] !== null && $o[$i] !== null && ! in_array(null, $rows[$i], true);
            foreach ($terms as [$z, $groups]) {
                $complete = $complete && $groups[$i] !== null && ! in_array(null, $z[$i], true);
            }
            if ($complete) {
                $retained[] = $i;
            } else {
                $dropped[] = $i;
            }
        }
        if ($missing === 'raise' && $dropped !== []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "missing values found with missing='raise'");
        }
        $coreTerms = [];
        $termLabels = [];
        foreach ($terms as [$z, $groups, $correlated]) {
            [$indices, $labels] = ModelData::groupIndices(array_map(fn (int $i) => $groups[$i], $retained));
            $termLabels[] = $labels;
            $coreTerms[] = [
                'z' => array_map(fn (int $i) => $z[$i], $retained),
                'groups' => $indices,
                'correlated' => $correlated,
            ];
        }
        $startValues = null;
        if ($start !== null) {
            $startValues = [];
            foreach ($start as $value) {
                $value = ModelData::number($value);
                if ($value === null) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'start cannot contain missing values');
                }
                $startValues[] = $value;
            }
        }
        $method = match ($optimizer) {
            'bobyqa', 'nelder_mead' => $optimizer,
            default => throw ErrorCode::exception(ErrorCode::UNSUPPORTED, "optimizer must be 'bobyqa' or 'nelder_mead'; R wrapper names are not runtime aliases"),
        };
        if ($maxEvaluations < 0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'maxEvaluations must be positive');
        }
        $design = [];
        foreach ($retained as $i) {
            $design[] = $intercept ? [1.0, ...$rows[$i]] : $rows[$i];
        }
        $result = LinearMixed::fit(
            $design,
            array_map(fn (int $i) => $outcomes[$i], $retained),
            $coreTerms,
            [
                'reml' => $reml,
                'optimizer' => $method,
                'objective_tolerance' => $objectiveTolerance,
                'weights' => array_map(fn (int $i) => $w[$i], $retained),
                'offset' => array_map(fn (int $i) => $o[$i], $retained),
                'start' => $startValues,
                'max_evaluations' => $maxEvaluations,
                'tolerance' => $tolerance,
            ],
        );

        return ModelData::accounting($n, $retained, $dropped)
            ->finish(self::lmmOutput($result, $reml, $termLabels), 'Phit\fit_mixed_gaussian_matrix', $method);
    }

    /**
     * The Phit\fit_mixed_gaussian_matrix / Phit\fit_mixed_gaussian result table.
     *
     * @param  list<list<int|string>>  $termLabels
     */
    public static function lmmOutput(array $result, bool $reml, array $termLabels): array
    {
        $terms = [];
        foreach ($result['terms'] as $t => $term) {
            $terms[] = [
                'covariance' => $term['covariance'],
                'effects' => $term['effects'],
                'effects_cond_var' => $term['effects_cond_var'],
                'group_labels' => $termLabels[$t],
                'group_sizes' => $term['group_sizes'],
            ];
        }

        return [
            'coefficients' => $result['coefficients'],
            'covariance' => $result['covariance'],
            'std_errors' => $result['std_errors'],
            'z_statistics' => $result['z_statistics'],
            'p_values' => $result['p_values'],
            'conf_int' => $result['conf_int'],
            'fixed_fitted' => $result['fixed_fitted'],
            'fitted' => $result['fitted'],
            'residuals' => $result['residuals'],
            'marginal_residuals' => $result['marginal_residuals'],
            'residual_variance' => $result['residual_variance'],
            'sigma' => $result['sigma'],
            'log_likelihood' => $result['log_likelihood'],
            'aic' => $result['aic'],
            'bic' => $result['bic'],
            'df_residual' => $result['df_residual'],
            'n_parameters' => $result['n_parameters'],
            'theta' => $result['theta'],
            'singular' => $result['singular'],
            'method' => $reml ? 'REML' : 'ML',
            'optimizer' => $result['optimizer'],
            'optimizer_implementation' => $result['optimizer'] === 'bobyqa' ? 'bobyqa 0.3.2 (Rust implementation of PRIMA BOBYQA)' : 'argmin 0.11.0 / Nelder-Mead',
            'converged' => true,
            'evaluations' => $result['evaluations'],
            'random' => $terms,
        ];
    }

    /** Phit\fit_mixed_logistic_matrix */
    public static function logisticMixed(
        array $x,
        array $y,
        array $groups,
        string $missing = 'drop',
        ?array $weights = null,
        ?array $offset = null,
        int $quadraturePoints = 1,
        int $maxEvaluations = 20000,
        float $tolerance = 1e-7,
        bool $intercept = true,
    ): array {
        if ($weights !== null) {
            foreach ($weights as $value) {
                $value = ModelData::number($value);
                if ($value !== null && $value > Limits::MAX_TRIALS) {
                    throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'trial counts must not exceed 1000000');
                }
            }
        }
        $data = ModelData::prepare($x, $y, $missing, 'binary', $weights, $offset, $intercept, $groups);
        [$indices, $labels] = ModelData::groupIndices($data->groups);
        if ($quadraturePoints < 0) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'quadraturePoints must be 1..61');
        }
        if ($maxEvaluations < 0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'maxEvaluations must be positive');
        }
        $r = LogisticRandomIntercept::fit($data->x, $data->y, $indices, [
            'weights' => $data->weights,
            'offset' => $data->offset,
            'quadrature_points' => $quadraturePoints,
            'max_evaluations' => $maxEvaluations,
            'tolerance' => $tolerance,
        ]);

        return $data->finish(self::logisticMixedOutput($r, $labels), 'Phit\fit_mixed_logistic_matrix', 'bobyqa');
    }

    /**
     * The Phit\fit_mixed_logistic_matrix / Phit\fit_mixed_logistic result table.
     *
     * @param  list<int|string>  $labels
     */
    public static function logisticMixedOutput(array $r, array $labels): array
    {
        return [
            'coefficients' => $r['coefficients'],
            'covariance' => $r['covariance'],
            'std_errors' => $r['std_errors'],
            'z_statistics' => $r['z_statistics'],
            'p_values' => $r['p_values'],
            'conf_int' => $r['conf_int'],
            'fixed_fitted' => $r['fixed_fitted'],
            'fitted' => $r['fitted'],
            'linear_predictor' => $r['linear_predictor'],
            'residuals' => $r['residuals'],
            'random_effects' => $r['random_effects'],
            'random_effects_cond_var' => $r['random_effects_cond_var'],
            'group_sizes' => $r['group_sizes'],
            'group_labels' => $labels,
            'n_groups' => count($labels),
            'random_intercept_variance' => $r['random_intercept_variance'],
            'random_intercept_sd' => $r['random_intercept_sd'],
            'log_likelihood' => $r['log_likelihood'],
            'neg2_log_likelihood' => $r['neg2_log_likelihood'],
            'aic' => $r['aic'],
            'bic' => $r['bic'],
            'n_parameters' => $r['n_parameters'],
            'df_residual' => $r['df_residual'],
            'quadrature_points' => $r['quadrature_points'],
            'integration' => $r['quadrature_points'] === 1 ? 'laplace' : 'adaptive_gauss_hermite',
            'evaluations' => $r['evaluations'],
            'gradient_max' => $r['gradient_max'],
            'boundary' => $r['boundary'],
            'singular' => $r['singular'],
            'family' => 'binomial',
            'link' => 'logit',
            'method' => 'ML',
            'optimizer' => 'bobyqa',
            'converged' => true,
        ];
    }

    /**
     * Rows of numbers, null, or NaN (missing); each row must be an array.
     *
     * @return list<list<?float>>
     */
    private static function numericRows(array $values): array
    {
        $out = [];
        foreach ($values as $row) {
            if (! is_array($row)) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each matrix row must be an array');
            }
            $converted = [];
            foreach ($row as $value) {
                $converted[] = ModelData::number($value);
            }
            $out[] = $converted;
        }

        return $out;
    }

    /** @return list<?float> */
    private static function numericVector(?array $values, int $n, float $default): array
    {
        if ($values === null) {
            return array_fill(0, $n, $default);
        }
        $out = [];
        foreach ($values as $value) {
            $out[] = ModelData::number($value);
        }
        if (count($out) !== $n) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'optional vectors must have one value per original observation');
        }

        return $out;
    }
}
