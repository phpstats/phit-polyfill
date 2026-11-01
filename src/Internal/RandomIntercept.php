<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\ErrorCode;

/**
 * Gaussian random-intercept maximum likelihood with one grouping factor
 * (stats-core mixed.rs, behind Phit\fit_mixed_gaussian_intercept_matrix). The variance ratio rho = tau/sigma^2
 * is profiled out on t = ln(1 + rho): a fixed grid of 64 steps, bisection of
 * every negative-to-positive score crossing, and the lowest objective wins.
 *
 * @internal
 */
final class RandomIntercept
{
    /** @var list<list<float>> R of the within-group QR, as columns */
    private array $within;

    private int $withinRows;

    /** @var list<list<float>> sqrt(m_g) times each group's [X, y] mean, $means[$g][$j] */
    private array $means;

    /** @var list<int> */
    private array $sizes;

    private function __construct(private int $n, private int $p) {}

    /**
     * @param  list<list<float>>  $x  rows of the fixed design
     * @param  list<float>  $y
     * @param  list<int>  $groups  dense group index per row
     * @return array<string, mixed>
     */
    public static function fit(array $x, array $y, array $groups): array
    {
        $n = count($x);
        if ($n === 0 || $x[0] === []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'x must have at least one row and column');
        }
        $p = count($x[0]);
        $rectangular = true;
        foreach ($x as $row) {
            $rectangular = $rectangular && count($row) === $p;
        }
        if (count($y) !== $n || count($groups) !== $n || $rectangular === false) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'x must be rectangular with one y and group per row');
        }
        if ($n <= $p) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'mixed model requires n > p');
        }
        $members = [];
        $indices = [];
        $map = [];
        foreach ($groups as $i => $id) {
            if (! isset($map[$id])) {
                $map[$id] = count($members);
                $members[] = [];
            }
            $members[$map[$id]][] = $i;
            $indices[] = $map[$id];
        }
        if (count($members) < 2 || count($members) === $n) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'need at least two groups and at least one repeated group to identify variances');
        }
        $design = [];
        for ($j = 0; $j < $p; $j++) {
            $design[$j] = array_column($x, $j);
        }
        $scales = [];
        foreach ($design as $column) {
            $scale = sqrt(LeastSquares::dot($column, $column, 0, $n));
            if (is_nan($scale) || is_infinite($scale) || $scale == 0.0) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'fixed design is rank deficient or exceeds floating-point range');
            }
            $scales[] = $scale;
        }

        $data = new self($n, $p);
        $data->prepare($design, $scales, $y, $members);
        $atZero = $data->profile(0.0, false);
        // Catch numerically perfect fixed fits: the Gaussian ML scale has no positive optimum.
        $yNorm = LeastSquares::dot($y, $y, 0, $n);
        if ($atZero['variance'] <= PHP_FLOAT_EPSILON * PHP_FLOAT_EPSILON * ($n * $n) * $yNorm / $n) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'numerically perfect fit: no positive residual-variance MLE');
        }
        $bestT = 0.0;
        $best = $atZero['objective'];
        $previousT = 0.0;
        $previousDerivative = $atZero['derivative'];
        $evaluations = 1;
        for ($k = 1; $k <= 64; $k++) {
            $t = $k * 0.5;
            $value = $data->profile($t, false);
            $evaluations++;
            if ($previousDerivative < 0.0 && $value['derivative'] >= 0.0) {
                $lower = $previousT;
                $upper = $t;
                for ($i = 0; $i < 40; $i++) {
                    $middle = ($lower + $upper) / 2.0;
                    $candidate = $data->profile($middle, false);
                    $evaluations++;
                    if ($candidate['derivative'] < 0.0) {
                        $lower = $middle;
                    } else {
                        $upper = $middle;
                    }
                }
                $candidateT = ($lower + $upper) / 2.0;
                $candidate = $data->profile($candidateT, false);
                $evaluations++;
                if ($candidate['objective'] < $best) {
                    $best = $candidate['objective'];
                    $bestT = $candidateT;
                }
            }
            if ($k === 64 && ($value['derivative'] <= 0.0 || $value['objective'] < $best)) {
                throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'mixed-model optimum exceeds variance-ratio search range; a zero residual-variance boundary may be present');
            }
            $previousT = $t;
            $previousDerivative = $value['derivative'];
        }
        $final = $data->profile($bestT, true);
        $evaluations++;
        $rho = expm1($bestT);
        $variance = $final['variance'];
        $tau = $variance * $rho;
        $beta = [];
        for ($j = 0; $j < $p; $j++) {
            $beta[] = $final['beta'][$j] / $scales[$j];
        }
        $information = $final['information'];
        $covariance = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $covariance[$i][$j] = $variance * $information[$i][$j] / $scales[$i] / $scales[$j];
            }
        }
        $stdErrors = [];
        $z = [];
        $pValues = [];
        $confInt = [];
        for ($j = 0; $j < $p; $j++) {
            $se = sqrt($covariance[$j][$j]);
            $stdErrors[] = $se;
            $z[] = $beta[$j] / $se;
            $pValues[] = 2.0 * Statrs::normalSf(abs($beta[$j] / $se));
            $confInt[] = [$beta[$j] - Statrs::Z_975 * $se, $beta[$j] + Statrs::Z_975 * $se];
        }
        // design * beta, accumulated column by column as nalgebra's gemv does.
        $fixedFitted = [];
        for ($i = 0; $i < $n; $i++) {
            $fixedFitted[$i] = $design[0][$i] * $beta[0];
        }
        for ($j = 1; $j < $p; $j++) {
            for ($i = 0; $i < $n; $i++) {
                $fixedFitted[$i] = $design[$j][$i] * $beta[$j] + $fixedFitted[$i];
            }
        }
        $marginal = [];
        for ($i = 0; $i < $n; $i++) {
            $marginal[] = $y[$i] - $fixedFitted[$i];
        }
        $randomEffects = [];
        $condVar = [];
        foreach ($members as $rows) {
            $sum = 0.0;
            foreach ($rows as $i) {
                $sum += $marginal[$i];
            }
            $size = count($rows);
            $shrinkage = ($size * $rho) / (1.0 + $size * $rho);
            $randomEffects[] = $shrinkage * ($sum / $size);
            $condVar[] = $tau / (1.0 + $size * $rho);
        }
        $fitted = [];
        $residuals = [];
        for ($i = 0; $i < $n; $i++) {
            $fitted[] = $fixedFitted[$i] + $randomEffects[$indices[$i]];
            $residuals[] = $y[$i] - $fitted[$i];
        }
        $all = [$beta, ...$covariance, $stdErrors, $z, $pValues, ...$confInt, $randomEffects, $condVar, $fitted, $residuals, $marginal];
        $finite = is_finite($tau);
        foreach ($all as $values) {
            foreach ($values as $v) {
                $finite = $finite && is_finite($v);
            }
        }
        if ($finite === false) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'mixed-model inference exceeds floating-point range');
        }

        return [
            'coefficients' => $beta,
            'covariance' => $covariance,
            'std_errors' => $stdErrors,
            'z_statistics' => $z,
            'p_values' => $pValues,
            'conf_int' => $confInt,
            'fixed_fitted' => $fixedFitted,
            'fitted' => $fitted,
            'residuals' => $residuals,
            'marginal_residuals' => $marginal,
            'random_effects' => $randomEffects,
            'random_effects_cond_var' => $condVar,
            'group_sizes' => array_map('count', $members),
            'residual_variance' => $variance,
            'random_intercept_variance' => $tau,
            'sigma' => sqrt($variance),
            'icc' => $rho / (1.0 + $rho),
            'log_likelihood' => -$final['objective'],
            'aic' => 2.0 * ($p + 2) + 2.0 * $final['objective'],
            'df_residual' => $n - $p,
            'evaluations' => $evaluations,
            'boundary' => $bestT == 0.0,
        ];
    }

    /**
     * Within-group QR and scaled group means of [X / scales, y].
     *
     * @param  list<list<float>>  $design  columns
     * @param  list<list<int>>  $members
     */
    private function prepare(array $design, array $scales, array $y, array $members): void
    {
        $p = $this->p;
        $centered = [];
        for ($j = 0; $j < $p; $j++) {
            $column = [];
            foreach ($design[$j] as $v) {
                $column[] = $v / $scales[$j];
            }
            $centered[$j] = $column;
        }
        $centered[$p] = $y;
        $means = [];
        foreach ($members as $g => $rows) {
            $size = count($rows);
            for ($j = 0; $j <= $p; $j++) {
                // Deviations from one observation keep precision under a large common offset.
                $origin = $centered[$j][$rows[0]];
                $sum = 0.0;
                foreach ($rows as $i) {
                    $sum += $centered[$j][$i] - $origin;
                }
                $mean = $origin + $sum / $size;
                $means[$g][$j] = $mean * sqrt($size);
                foreach ($rows as $i) {
                    $centered[$j][$i] -= $mean;
                }
            }
        }
        [$this->within] = LeastSquares::householder($centered, $this->n, $p + 1);
        $this->withinRows = $p + 1;
        $this->means = $means;
        $this->sizes = array_map('count', $members);
    }

    /**
     * The ML profile at t = ln(1 + rho): beta, variance, objective, and its
     * derivative in t; with $inference, also the scaled (A'A)^-1.
     *
     * @return array{beta: list<float>, variance: float, objective: float, derivative: float, information: ?list<list<float>>}
     */
    private function profile(float $t, bool $inference): array
    {
        $rho = expm1($t);
        $n = $this->n;
        $p = $this->p;
        $rows = $this->withinRows + count($this->sizes);
        // The compact (p+1+G)-row system, as columns.
        $a = [];
        for ($j = 0; $j <= $p; $j++) {
            $column = [];
            for ($i = 0; $i < $this->withinRows; $i++) {
                $column[] = $this->within[$j][$i];
            }
            foreach ($this->sizes as $g => $size) {
                $column[] = $this->means[$g][$j] / sqrt(1.0 + $size * $rho);
            }
            $a[$j] = $column;
        }
        $b = $a[$p];
        unset($a[$p]);
        $svd = new Svd($a, true, true);
        $singular = $svd->singularValues;
        $tolerance = PHP_FLOAT_EPSILON * max($n, $p) * $singular[0];
        foreach ($singular as $s) {
            if (is_nan($s) || is_infinite($s) || $s <= $tolerance) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'mixed-model fixed design/information is numerically rank deficient');
            }
        }
        $beta = $svd->solveVector($b, $tolerance);
        if (count($beta) !== $p) {
            throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'least-squares factorization failed');
        }
        // b - A beta, with A beta accumulated column by column (gemv).
        $fitted = [];
        for ($i = 0; $i < $rows; $i++) {
            $fitted[$i] = $a[0][$i] * $beta[0];
        }
        for ($j = 1; $j < $p; $j++) {
            for ($i = 0; $i < $rows; $i++) {
                $fitted[$i] = $a[$j][$i] * $beta[$j] + $fitted[$i];
            }
        }
        $residual = [];
        for ($i = 0; $i < $rows; $i++) {
            $residual[] = $b[$i] - $fitted[$i];
        }
        $q = LeastSquares::dot($residual, $residual, 0, $rows);
        if (is_nan($q) || is_infinite($q) || $q <= 0.0) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'mixed-model residual variance is zero or exceeds floating-point range');
        }
        $variance = $q / $n; // ML scale, not RSS/(n-p), and not REML.
        $logdet = 0.0;
        foreach ($this->sizes as $m) {
            $logdet += log1p($m * $rho);
        }
        $objective = 0.5 * ($n * (log(2.0 * M_PI) + log($variance) + 1.0) + $logdet);
        // Envelope theorem: q' from the group residual sums of the cached means.
        $qDerivative = 0.0;
        foreach ($this->sizes as $g => $m) {
            $sum = 0.0;
            for ($j = 0; $j < $p; $j++) {
                $sum += $this->means[$g][$j] * $beta[$j];
            }
            $meanResidual = $this->means[$g][$p] - $sum;
            $ratio = $meanResidual / (1.0 + $m * $rho);
            $qDerivative += $m * ($ratio * $ratio);
        }
        $qDerivative = -$qDerivative;
        $logdetDerivative = 0.0;
        foreach ($this->sizes as $m) {
            $logdetDerivative += $m / (1.0 + $m * $rho);
        }
        $derivative = 0.5 * (1.0 + $rho) * ($n * $qDerivative / $q + $logdetDerivative);
        if (is_nan($objective) || is_infinite($objective) || is_nan($derivative) || is_infinite($derivative)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'mixed-model likelihood exceeds floating-point range');
        }
        $information = null;
        if ($inference) {
            $vT = $svd->vT;
            $information = [];
            for ($i = 0; $i < $p; $i++) {
                for ($j = 0; $j < $p; $j++) {
                    $sum = 0.0;
                    for ($k = 0; $k < $p; $k++) {
                        $sum += $vT[$i][$k] * $vT[$j][$k] / ($singular[$k] * $singular[$k]);
                    }
                    $information[$i][$j] = $sum;
                }
            }
        }

        return ['beta' => $beta, 'variance' => $variance, 'objective' => $objective, 'derivative' => $derivative, 'information' => $information];
    }
}
