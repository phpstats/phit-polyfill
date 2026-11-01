<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Internal\LeastSquares;
use Phit\Polyfill\Internal\ModelData;
use Phit\Polyfill\Internal\Statrs;

/**
 * Ordinary and weighted least squares with known offsets (Phit\fit_ols_matrix).
 */
final class Ols
{
    public static function fit(
        array $x,
        array $y,
        string $missing = 'drop',
        ?array $weights = null,
        ?array $offset = null,
        bool $intercept = true,
    ): array {
        $data = ModelData::prepare($x, $y, $missing, 'continuous', $weights, $offset, $intercept);
        $r = self::fitWithData($data->x, $data->y, $data->weights, $data->offset);

        return $data->finish([
            'coefficients' => $r['coefficients'],
            'covariance' => $r['covariance'],
            'std_errors' => $r['std_errors'],
            't_statistics' => $r['t_statistics'],
            'p_values' => $r['p_values'],
            'conf_int' => $r['conf_int'],
            'fitted' => $r['fitted'],
            'residuals' => $r['residuals'],
            'rss' => $r['rss'],
            'residual_variance' => $r['residual_variance'],
            'sigma' => $r['sigma'],
            'r_squared' => $r['r_squared'],
            'adj_r_squared' => $r['adj_r_squared'],
            'df_residual' => $r['df_residual'],
            'log_likelihood' => $r['log_likelihood'],
        ], 'Phit\fit_ols_matrix', 'scaled_svd');
    }

    private static function error(int $code, string $message): \Exception
    {
        return ErrorCode::exception($code, $message);
    }

    /**
     * Weighted Gaussian fit with known additive offsets. Prior weights specify
     * Var(e_i) = sigma^2 / w_i; residuals stay on the response scale.
     * Shared with the formula interface (Phit\fit_ols, Phit\fit_generalized_linear gaussian).
     *
     * @internal
     */
    public static function fitWithData(array $x, array $y, array $weights, array $offset): array
    {
        $n = count($x);
        foreach ($weights as $w) {
            if ($w <= 0.0) {
                throw self::error(ErrorCode::INVALID_INPUT, 'weights must be positive (integer trials for binary outcomes)');
            }
        }
        if ($n === 0 || count($y) !== $n || $x[0] === []) {
            throw self::error(ErrorCode::NONFINITE, 'x and y must be finite, rectangular and aligned');
        }
        $plain = true;
        foreach ($weights as $i => $w) {
            if ($w != 1.0 || $offset[$i] != 0.0) {
                $plain = false;
                break;
            }
        }
        if ($plain) {
            return self::fitOls($x, $y);
        }
        $p = count($x[0]);
        $intercept = self::hasIntercept($x, $p);
        // Householder QR is unstable when row scales differ by many orders of
        // magnitude unless the heavy rows come first (Powell and Reid), so the
        // whitened rows are factored in decreasing-weight order, as in the
        // extension. fitted and residuals are rebuilt in the original order below.
        $order = array_keys($weights);
        usort($order, fn (int $a, int $b): int => $weights[$b] <=> $weights[$a]);
        $wx = [];
        $wy = [];
        foreach ($order as $i) {
            $row = $x[$i];
            $root = sqrt($weights[$i]);
            $scaled = [];
            foreach ($row as $value) {
                $scaled[] = $value * $root;
            }
            $wx[] = $scaled;
            $wy[] = ($y[$i] - $offset[$i]) * $root;
        }
        $result = self::fitOls($wx, $wy);
        $fitted = [];
        $residuals = [];
        foreach ($x as $i => $row) {
            $sum = -0.0; // Rust's f64 sum starts from -0.0
            foreach ($row as $j => $value) {
                $sum += $value * $result['coefficients'][$j];
            }
            $fitted[$i] = $offset[$i] + $sum;
            $residuals[$i] = $y[$i] - $fitted[$i];
        }
        $result['fitted'] = $fitted;
        $result['residuals'] = $residuals;
        // R summary.lm removes the offset from fitted values before calculating MSS.
        $mean = 0.0;
        if ($intercept) {
            $num = 0.0;
            foreach ($fitted as $i => $f) {
                $num += ($f - $offset[$i]) * $weights[$i];
            }
            $mean = $num / array_sum($weights);
        }
        $mss = 0.0;
        foreach ($fitted as $i => $f) {
            $d = $f - $offset[$i] - $mean;
            $mss += $weights[$i] * ($d * $d);
        }
        if (! is_finite($mss)) {
            throw self::error(ErrorCode::NONFINITE, 'weighted fitted variation exceeds floating-point range; rescale the outcome');
        }
        if ($p === 1 && $intercept) {
            $result['r_squared'] = 0.0;
            $result['adj_r_squared'] = 0.0;
        } else {
            $result['r_squared'] = fdiv($mss, $mss + $result['rss']);
            $result['adj_r_squared'] = 1.0 - (1.0 - $result['r_squared']) * ($n - (int) $intercept) / ($n - $p);
        }
        $logWeights = 0.0;
        foreach ($weights as $w) {
            $logWeights += log($w);
        }
        $result['log_likelihood'] += 0.5 * $logWeights;
        foreach ([...$fitted, ...$residuals] as $value) {
            if (! is_finite($value)) {
                throw self::error(ErrorCode::NONFINITE, 'unweighted response-scale results exceed floating-point range');
            }
        }

        return $result;
    }

    /** @param list<list<float>> $x */
    private static function hasIntercept(array $x, int $p): bool
    {
        for ($j = 0; $j < $p; $j++) {
            $ones = true;
            foreach ($x as $row) {
                if ($row[$j] != 1.0) {
                    $ones = false;
                    break;
                }
            }
            if ($ones) {
                return true;
            }
        }

        return false;
    }

    /** Full-rank fit; an intercept must be supplied as an exact column of ones. */
    private static function fitOls(array $x, array $y): array
    {
        $n = count($x);
        if ($n === 0 || $x[0] === []) {
            throw self::error(ErrorCode::INVALID_INPUT, 'x must have at least one row and column');
        }
        $p = count($x[0]);
        foreach ($x as $row) {
            foreach ($row as $value) {
                if (! is_finite($value)) {
                    throw self::error(ErrorCode::NONFINITE, 'x and y must contain only finite numbers');
                }
            }
        }
        foreach ($y as $value) {
            if (! is_finite($value)) {
                throw self::error(ErrorCode::NONFINITE, 'x and y must contain only finite numbers');
            }
        }
        if ($n <= $p) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'OLS requires n > p for residual degrees of freedom');
        }
        $intercept = self::hasIntercept($x, $p);
        // Scale columns before factorizing so rank detection is invariant to predictor units.
        $columns = [];
        $scales = [];
        for ($j = 0; $j < $p; $j++) {
            $column = array_column($x, $j);
            // nalgebra's column norm, summed the same way.
            $scale = sqrt(LeastSquares::dot($column, $column, 0, $n));
            if (! is_finite($scale) || $scale == 0.0) {
                throw self::error(ErrorCode::RANK_DEFICIENCY, 'design is rank deficient or its scale exceeds floating-point range');
            }
            foreach ($column as $i => $value) {
                $column[$i] = $value / $scale;
            }
            $columns[] = $column;
            $scales[] = $scale;
        }
        $ls = new LeastSquares($columns, $n, $p, $y);
        if ($ls->rankDeficient()) {
            throw self::error(ErrorCode::RANK_DEFICIENCY, 'design is rank deficient at numerical precision');
        }
        $solution = $ls->solution()
            ?? throw self::error(ErrorCode::FACTORIZATION_FAILURE, 'least-squares factorization failed');
        $beta = [];
        for ($j = 0; $j < $p; $j++) {
            $beta[$j] = $solution[$j] / $scales[$j];
        }
        $fitted = [];
        $residuals = [];
        $nonzeroResidual = false;
        foreach ($x as $i => $row) {
            // design * beta as nalgebra's gemv: the first product starts the sum.
            $f = $row[0] * $beta[0];
            for ($j = 1; $j < $p; $j++) {
                $f = $row[$j] * $beta[$j] + $f;
            }
            $fitted[$i] = $f;
            $residual = $y[$i] - $f;
            $residuals[$i] = $residual;
            $nonzeroResidual = $nonzeroResidual || $residual != 0.0;
        }
        // nalgebra's norm_squared: the 8-accumulator dot product.
        $rss = LeastSquares::dot($residuals, $residuals, 0, $n);
        $df = $n - $p;
        $variance = $rss / $df;
        if (($rss == 0.0 && $nonzeroResidual) || ($rss > 0.0 && $variance == 0.0)) {
            throw self::error(ErrorCode::NONFINITE, 'nonzero residual variance underflows floating-point range; rescale the outcome');
        }
        $covariance = $ls->inverseGram($scales, $variance)
            ?? throw self::error(ErrorCode::FACTORIZATION_FAILURE, 'triangular factor could not be inverted');
        $critical = Statrs::tCritical95((float) $df);
        $stdErrors = [];
        $tStatistics = [];
        $pValues = [];
        $confInt = [];
        for ($j = 0; $j < $p; $j++) {
            $se = sqrt($covariance[$j][$j]);
            $t = fdiv($beta[$j], $se);
            $stdErrors[$j] = $se;
            $tStatistics[$j] = $t;
            $pValues[$j] = match (true) {
                is_nan($t) => NAN,
                is_infinite($t) => 0.0,
                default => 2.0 * Statrs::tSf(abs($t), (float) $df),
            };
            $confInt[$j] = [$beta[$j] - $critical * $se, $beta[$j] + $critical * $se];
        }
        $mean = $intercept ? array_sum($y) / $n : 0.0;
        $tss = 0.0;
        foreach ($y as $value) {
            $tss += ($value - $mean) * ($value - $mean);
        }
        if (! is_finite($tss)) {
            throw self::error(ErrorCode::NONFINITE, 'total outcome variation exceeds floating-point range; rescale the outcome');
        }
        if ($p === 1 && $intercept) {
            $rSquared = 0.0;
            $adjRSquared = 0.0;
        } else {
            $rSquared = 1.0 - fdiv($rss, $tss);
            $adjRSquared = 1.0 - (1.0 - $rSquared) * ($n - (int) $intercept) / $df;
        }
        // Gaussian likelihood uses the ML variance RSS/n, unlike the inference variance RSS/(n-p).
        $logLikelihood = -0.5 * $n * (log(2.0 * M_PI) + log($rss / $n) + 1.0);
        $finite = is_finite($rss);
        foreach ([...$beta, ...$fitted, ...$residuals] as $value) {
            $finite = $finite && is_finite($value);
        }
        foreach ($covariance as $row) {
            foreach ($row as $value) {
                $finite = $finite && is_finite($value);
            }
        }
        if (! $finite) {
            throw self::error(ErrorCode::NONFINITE, 'calculation exceeds floating-point range');
        }

        return [
            'coefficients' => $beta,
            'covariance' => $covariance,
            'std_errors' => $stdErrors,
            't_statistics' => $tStatistics,
            'p_values' => $pValues,
            'conf_int' => $confInt,
            'fitted' => $fitted,
            'residuals' => $residuals,
            'rss' => $rss,
            'residual_variance' => $variance,
            'sigma' => sqrt($variance),
            'r_squared' => $rSquared,
            'adj_r_squared' => $adjRSquared,
            'df_residual' => $df,
            'log_likelihood' => $logLikelihood,
        ];
    }
}
