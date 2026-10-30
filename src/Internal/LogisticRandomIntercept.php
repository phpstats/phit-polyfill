<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\Distribution\Support;
use Phit\Polyfill\ErrorCode;
use Phit\Polyfill\Internal\Family\LogisticFamily;

/**
 * Binomial logit with one Gaussian random intercept, by marginal maximum
 * likelihood with the Laplace approximation or adaptive Gauss-Hermite
 * quadrature (stats-core logistic_mixed.rs, behind Phit\fit_logistic_mixed and
 * Phit\fit_glmer). BOBYQA maximizes over the scaled fixed effects and tau.
 *
 * @internal
 */
final class LogisticRandomIntercept
{
    /** @var list<list<float>> scaled design rows */
    private array $x;

    private int $p;

    /** @var list<float> */
    private array $y;

    /** @var list<float> */
    private array $weights;

    /** @var list<float> */
    private array $offset;

    /** @var list<list<int>> */
    private array $members;

    /** @var list<array{0: float, 1: float}> nodes and weights */
    private array $quadrature;

    private float $constant;

    /**
     * @param  list<list<float>>  $x
     * @param  list<float>  $y
     * @param  list<int>  $groups
     * @param  array{weights: list<float>, offset: list<float>, quadrature_points: int, max_evaluations: int, tolerance: float}  $options
     * @return array<string, mixed>
     */
    public static function fit(array $x, array $y, array $groups, array $options): array
    {
        $n = count($x);
        if ($n === 0 || $x[0] === []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'x must have rows and columns');
        }
        $p = count($x[0]);
        $rectangular = true;
        foreach ($x as $row) {
            $rectangular = $rectangular && count($row) === $p;
        }
        if ($n <= $p || count($y) !== $n || count($groups) !== $n || $rectangular === false) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'require rectangular finite x, proportions y in [0,1], aligned groups and n > p');
        }
        foreach ($y as $v) {
            if (is_nan($v) || $v < 0.0 || $v > 1.0) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binomial y must be proportions in [0,1]');
            }
        }
        $points = $options['quadrature_points'];
        if ($points < 1 || $points > 61) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'quadraturePoints must be 1 (Laplace) or 2..61 (adaptive quadrature)');
        }
        $tolerance = $options['tolerance'];
        if (is_nan($tolerance) || is_infinite($tolerance) || $tolerance <= 0.0 || $tolerance > 0.1) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'tolerance must be finite in (0,0.1]');
        }
        [$weights, $offset] = self::observations($n, $options['weights'], $options['offset']);
        foreach ($y as $i => $v) {
            $events = $v * $weights[$i];
            if ($weights[$i] > Limits::MAX_TRIALS || abs($events - Support::round($events)) > 1e-8) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'require integer trials <= 1000000 and integer success counts');
            }
        }
        // Canonical integer successes: tiny proportions must not become fictitious events.
        $allZero = true;
        $allOne = true;
        foreach ($y as $i => $v) {
            $y[$i] = Support::round($v * $weights[$i]) / $weights[$i];
            $allZero = $allZero && $y[$i] == 0.0;
            $allOne = $allOne && $y[$i] == 1.0;
        }
        if ($allZero || $allOne) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'require successes and failures for a finite unpenalized fit');
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
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'require at least two groups and repeated observations');
        }
        // RMS scaling keeps the optimizer's coefficient units manageable.
        $scales = [];
        for ($j = 0; $j < $p; $j++) {
            $column = array_column($x, $j);
            $scales[] = sqrt(LeastSquares::dot($column, $column, 0, $n)) / sqrt($n);
        }
        foreach ($scales as $s) {
            if (is_nan($s) || is_infinite($s) || $s == 0.0) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'fixed design is rank deficient or has invalid scale');
            }
        }
        $scaled = [];
        foreach ($x as $i => $row) {
            foreach ($row as $j => $v) {
                $scaled[$i][$j] = $v / $scales[$j];
            }
        }
        $columns = [];
        for ($j = 0; $j < $p; $j++) {
            $columns[$j] = array_column($scaled, $j);
        }
        $singular = (new Svd($columns, false, false))->singularValues;
        $rankTolerance = PHP_FLOAT_EPSILON * max($n, $p) * $singular[0];
        foreach ($singular as $s) {
            if ($s <= $rankTolerance) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'fixed design is numerically rank deficient');
            }
        }
        $constant = 0.0;
        foreach ($y as $i => $v) {
            if ($v == 0.0 || $v == 1.0) {
                $term = 0.0;
            } else {
                $w = $weights[$i];
                $k = $v * $w;
                $term = Statrs::lnGamma($w + 1.0) - Statrs::lnGamma($k + 1.0) - Statrs::lnGamma($w - $k + 1.0);
            }
            $constant += $term;
        }
        if (is_nan($constant) || is_infinite($constant)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'binomial likelihood constant exceeds floating-point range');
        }
        // Separation in the fixed design also prevents a finite unpenalized GLMM MLE.
        $glm = (new LogisticFamily)->fitData($x, $y, $weights, $offset);

        $data = new self;
        $data->x = $scaled;
        $data->p = $p;
        $data->y = $y;
        $data->weights = $weights;
        $data->offset = $offset;
        $data->members = $members;
        $data->quadrature = self::normalQuadrature($points);
        $data->constant = $constant;

        $parameters = [];
        for ($j = 0; $j < $p; $j++) {
            $parameters[] = $glm['coefficients'][$j] * $scales[$j];
        }
        $parameters[] = 1.0;
        $lower = array_fill(0, $p + 1, -100.0);
        $lower[$p] = 0.0;
        $upper = array_fill(0, $p + 1, 100.0);
        $upper[$p] = 20.0;
        foreach ($parameters as $j => $v) {
            if ($v <= $lower[$j] || $v >= $upper[$j]) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'initial fixed effects exceed mixed-model optimizer bounds; rescale/center predictors');
            }
        }
        $config = new Bobyqa\Config($p + 1);
        $config->rhoBegin = 0.1;
        $config->rhoEnd = $tolerance;
        $config->maxFun = $options['max_evaluations'];
        try {
            $optimizer = new Bobyqa($p + 1, $config);
        } catch (Bobyqa\BobyqaException $e) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'invalid binomial mixed control: '.$e->getMessage());
        }
        $outcome = $optimizer->minimize(function (array $t) use ($data): float {
            try {
                return -$data->logLikelihood($t);
            } catch (\Exception $e) {
                return INF;
            }
        }, $parameters, $lower, $upper);
        if ($outcome->status !== Bobyqa\Status::Converged || $outcome->foundFinite() === false) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binomial mixed optimizer failed: '.$outcome->status->label());
        }
        foreach ($parameters as $j => $v) {
            if ($j < $p && abs($v) >= 99.9) {
                throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'binomial mixed fit reaches search bound; separation/weak identification may be present');
            }
        }
        if ($parameters[$p] >= 19.9) {
            throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'binomial mixed fit reaches search bound; separation/weak identification may be present');
        }
        $likelihood = $data->logLikelihood($parameters);
        if ($glm['log_likelihood'] >= $likelihood
            || ($parameters[$p] < 1e-4 && $likelihood - $glm['log_likelihood'] <= 64.0 * PHP_FLOAT_EPSILON * (1.0 + abs($likelihood)))) {
            for ($j = 0; $j < $p; $j++) {
                $parameters[$j] = $glm['coefficients'][$j] * $scales[$j];
            }
            $parameters[$p] = 0.0;
            $likelihood = $glm['log_likelihood'];
        }
        $boundary = $parameters[$p] == 0.0;
        $tau = $parameters[$p];
        if ($boundary) {
            $covariance = $glm['covariance'];
        } else {
            $coarse = self::information($data, $parameters, 1e-3);
            $fine = self::information($data, $parameters, 5e-4);
            $m = $p + 1;
            $hessian = [];
            for ($i = 0; $i < $m; $i++) {
                for ($j = 0; $j < $m; $j++) {
                    $hessian[$i][$j] = (4.0 * $fine[$i][$j] - $coarse[$i][$j]) / 3.0;
                }
            }
            $inverse = self::choleskyInverse($hessian);
            if ($inverse === null) {
                throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'marginal observed information is not positive definite; inference is unreliable');
            }
            $covariance = [];
            for ($i = 0; $i < $p; $i++) {
                for ($j = 0; $j < $p; $j++) {
                    $covariance[$i][$j] = $inverse[$i][$j] / $scales[$i] / $scales[$j];
                }
            }
        }
        $beta = [];
        for ($j = 0; $j < $p; $j++) {
            $beta[] = $parameters[$j] / $scales[$j];
        }
        $eta = $data->eta($parameters);
        $effects = [];
        $condVar = [];
        foreach ($data->members as $rows) {
            [$u, $precision] = $data->mode($rows, $eta, $tau);
            $effects[] = $tau * $u;
            $condVar[] = $tau * $tau / $precision;
        }
        $linearPredictor = [];
        $fitted = [];
        for ($i = 0; $i < $n; $i++) {
            $linearPredictor[] = $eta[$i] + $effects[$indices[$i]];
            $fitted[] = self::sigmoid($linearPredictor[$i]);
        }
        $stdErrors = [];
        $z = [];
        for ($j = 0; $j < $p; $j++) {
            $stdErrors[] = sqrt($covariance[$j][$j]);
            $z[] = $beta[$j] / $stdErrors[$j];
        }
        $gradientMax = 0.0;
        for ($j = 0; $j <= $p; $j++) {
            if ($boundary && $j === $p) {
                continue;
            }
            $h = 1e-5 * (1.0 + abs($parameters[$j]));
            $a = $parameters;
            $b = $parameters;
            $a[$j] += $h;
            $b[$j] -= $h;
            $gradient = abs(($data->logLikelihood($a) - $data->logLikelihood($b)) / (2.0 * $h));
            // f64::max: a NaN argument loses.
            if (is_nan($gradientMax) || $gradient > $gradientMax) {
                $gradientMax = $gradient;
            }
        }
        if (is_nan($gradientMax) || is_infinite($gradientMax) || $gradientMax > 1e-3) {
            throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'binomial mixed optimizer stopped with a nonstationary likelihood score');
        }
        $pValues = [];
        $confInt = [];
        for ($j = 0; $j < $p; $j++) {
            $pValues[] = 2.0 * Statrs::normalSf(abs($z[$j]));
            $confInt[] = [$beta[$j] - Statrs::Z_975 * $stdErrors[$j], $beta[$j] + Statrs::Z_975 * $stdErrors[$j]];
        }
        $fixedFitted = [];
        $residuals = [];
        for ($i = 0; $i < $n; $i++) {
            $fixedFitted[] = self::sigmoid($eta[$i]);
            $residuals[] = $y[$i] - $fitted[$i];
        }
        $nParameters = $p + 1;
        $finite = true;
        foreach ([$beta, ...$covariance, $stdErrors, $z, $fitted, $effects, $condVar] as $values) {
            foreach ($values as $v) {
                $finite = $finite && is_finite($v);
            }
        }
        if ($finite === false) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'binomial mixed results exceed floating-point range');
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
            'linear_predictor' => $linearPredictor,
            'residuals' => $residuals,
            'random_effects' => $effects,
            'random_effects_cond_var' => $condVar,
            'group_sizes' => array_map('count', $data->members),
            'random_intercept_variance' => $tau * $tau,
            'random_intercept_sd' => $tau,
            'log_likelihood' => $likelihood,
            'neg2_log_likelihood' => -2.0 * $likelihood,
            'aic' => 2.0 * $nParameters - 2.0 * $likelihood,
            'bic' => $nParameters * log($n) - 2.0 * $likelihood,
            'n_parameters' => $nParameters,
            'df_residual' => $n - $nParameters,
            'quadrature_points' => $points,
            'evaluations' => $outcome->nEval,
            'gradient_max' => $gradientMax,
            'boundary' => $boundary,
            'singular' => $tau < 1e-4,
        ];
    }

    /** stats-core observations::prepare for binary outcomes. */
    private static function observations(int $n, array $weights, array $offset): array
    {
        $weights = $weights === [] ? array_fill(0, $n, 1.0) : $weights;
        $offset = $offset === [] ? array_fill(0, $n, 0.0) : $offset;
        if (count($weights) !== $n || count($offset) !== $n) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'weights and offset must have one value per row');
        }
        foreach ([...$weights, ...$offset] as $v) {
            if (is_nan($v) || is_infinite($v)) {
                throw ErrorCode::exception(ErrorCode::NONFINITE, 'weights and offset must be finite');
            }
        }
        foreach ($weights as $v) {
            if ($v <= 0.0 || fmod($v, 1.0) != 0.0) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'weights must be positive (integer trials for binary outcomes)');
            }
        }
        foreach ($weights as $v) {
            if ($v > Limits::MAX_TRIALS) {
                throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'binomial trial weights must not exceed 1000000; larger-count likelihood arithmetic is not verified');
            }
        }

        return [$weights, $offset];
    }

    private static function sigmoid(float $v): float
    {
        if ($v >= 0.0) {
            return 1.0 / (1.0 + exp(-$v));
        }
        $e = exp($v);

        return $e / (1.0 + $e);
    }

    private static function softplus(float $v): float
    {
        return ($v > 0.0 ? $v : 0.0) + log1p(exp(-abs($v)));
    }

    /**
     * Standard-normal quadrature from the Hermite Jacobi matrix (Golub-Welsch),
     * with Christoffel weights from the orthonormal polynomial sum.
     *
     * @return list<array{0: float, 1: float}>
     */
    private static function normalQuadrature(int $points): array
    {
        $jacobi = array_fill(0, $points, array_fill(0, $points, 0.0));
        for ($j = 1; $j < $points; $j++) {
            $jacobi[$j - 1][$j] = sqrt($j);
            $jacobi[$j][$j - 1] = sqrt($j);
        }
        $eigenvalues = (new SymmetricEigen($jacobi))->eigenvalues;
        $result = [];
        foreach ($eigenvalues as $z) {
            $previous = 0.0;
            $current = 1.0;
            $sum = 1.0;
            for ($k = 1; $k < $points; $k++) {
                $next = ($z * $current - sqrt($k - 1) * $previous) / sqrt($k);
                $sum += $next * $next;
                $previous = $current;
                $current = $next;
            }
            $result[] = [$z, 1.0 / $sum];
        }
        // total_cmp ordering; the nodes are finite and distinct.
        usort($result, fn (array $a, array $b) => $a[0] <=> $b[0]);

        return $result;
    }

    /** @return list<float> */
    private function eta(array $parameters): array
    {
        $eta = [];
        foreach ($this->x as $i => $row) {
            $sum = 0.0;
            for ($j = 0; $j < $this->p; $j++) {
                $sum += $row[$j] * $parameters[$j];
            }
            $eta[] = $this->offset[$i] + $sum;
        }

        return $eta;
    }

    private function groupKernel(array $rows, array $eta, float $tau, float $u): float
    {
        $sum = 0.0;
        foreach ($rows as $i) {
            $e = $eta[$i] + $tau * $u;
            $sum += $this->weights[$i] * ($this->y[$i] * self::softplus(-$e) + (1.0 - $this->y[$i]) * self::softplus($e));
        }

        return -0.5 * $u * $u - $sum;
    }

    /** @return array{0: float, 1: float} */
    private function scorePrecision(array $rows, array $eta, float $tau, float $u): array
    {
        $score = -$u;
        $precision = 1.0;
        foreach ($rows as $i) {
            $e = $eta[$i] + $tau * $u;
            $score += $tau * $this->weights[$i] * ($this->y[$i] * self::sigmoid(-$e) - (1.0 - $this->y[$i]) * self::sigmoid($e));
            $tail = exp(-abs($e));
            $denominator = 1.0 + $tail;
            $precision += $tau * $tau * $this->weights[$i] * $tail / ($denominator * $denominator);
        }

        return [$score, $precision];
    }

    /** @return array{0: float, 1: float} conditional mode u and the posterior precision there */
    private function mode(array $rows, array $eta, float $tau): array
    {
        if ($tau == 0.0) {
            return [0.0, 1.0];
        }
        // The log posterior is strictly concave; Newton is safeguarded by a score bracket.
        $total = 0.0;
        foreach ($rows as $i) {
            $total += $this->weights[$i];
        }
        $bound = abs($tau) * $total + 1.0;
        $lower = -$bound;
        $upper = $bound;
        $u = 0.0;
        for ($iteration = 0; $iteration < 100; $iteration++) {
            [$score, $precision] = $this->scorePrecision($rows, $eta, $tau, $u);
            if (is_nan($score) || is_infinite($score) || is_nan($precision) || is_infinite($precision)) {
                throw ErrorCode::exception(ErrorCode::NONFINITE, 'binomial mixed inner calculation overflowed');
            }
            if (abs($score) <= 1e-11 * (1.0 + abs($u)) || abs($score) / $precision <= 1e-12 * (1.0 + abs($u))) {
                return [$u, $precision];
            }
            if ($score > 0.0) {
                $lower = $u;
            } else {
                $upper = $u;
            }
            $trial = $u + $score / $precision;
            $u = $trial > $lower && $trial < $upper ? $trial : 0.5 * ($lower + $upper);
        }
        throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'binomial mixed conditional mode did not converge');
    }

    private function logLikelihood(array $parameters): float
    {
        $tau = $parameters[$this->p];
        $eta = $this->eta($parameters);
        foreach ($eta as $e) {
            if (is_nan($e) || is_infinite($e)) {
                throw ErrorCode::exception(ErrorCode::NONFINITE, 'linear predictor exceeds floating-point range');
            }
        }
        $ll = $this->constant;
        foreach ($this->members as $rows) {
            [$mode, $precision] = $this->mode($rows, $eta, $tau);
            $root = sqrt($precision);
            $logs = [];
            $maximum = -INF;
            foreach ($this->quadrature as [$z, $w]) {
                $value = log($w) + $this->groupKernel($rows, $eta, $tau, $mode + $z / $root) + 0.5 * $z * $z;
                $logs[] = $value;
                // f64::max: NaN loses.
                if ($value > $maximum || is_nan($maximum)) {
                    $maximum = $value;
                }
            }
            $sum = 0.0;
            foreach ($logs as $v) {
                $sum += exp($v - $maximum);
            }
            // Laplace is exactly the one-node case: kernel(mode) - log(h)/2.
            $ll += $maximum + log($sum) - 0.5 * log($precision);
        }
        if (is_nan($ll) || is_infinite($ll)) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binomial mixed likelihood is not finite');
        }

        return $ll;
    }

    /**
     * Central-difference observed information of the marginal likelihood.
     *
     * @return list<list<float>>
     */
    private static function information(self $data, array $parameters, float $step): array
    {
        $m = count($parameters);
        $f0 = -$data->logLikelihood($parameters);
        $steps = [];
        foreach ($parameters as $v) {
            $steps[] = $step * (1.0 + abs($v));
        }
        $h = array_fill(0, $m, array_fill(0, $m, 0.0));
        for ($i = 0; $i < $m; $i++) {
            $a = $parameters;
            $b = $parameters;
            $a[$i] += $steps[$i];
            $b[$i] -= $steps[$i];
            $fa = -$data->logLikelihood($a);
            $h[$i][$i] = ($fa - $data->logLikelihood($b) - 2.0 * $f0) / ($steps[$i] * $steps[$i]);
            for ($j = 0; $j < $i; $j++) {
                $pp = $parameters;
                $pm = $parameters;
                $mp = $parameters;
                $mm = $parameters;
                $pp[$i] += $steps[$i];
                $pp[$j] += $steps[$j];
                $pm[$i] += $steps[$i];
                $pm[$j] -= $steps[$j];
                $mp[$i] -= $steps[$i];
                $mp[$j] += $steps[$j];
                $mm[$i] -= $steps[$i];
                $mm[$j] -= $steps[$j];
                $value = (-$data->logLikelihood($pp) + $data->logLikelihood($pm) + $data->logLikelihood($mp) - $data->logLikelihood($mm))
                    / (4.0 * $steps[$i] * $steps[$j]);
                $h[$i][$j] = $value;
                $h[$j][$i] = $value;
            }
        }

        return $h;
    }

    /**
     * nalgebra Cholesky::new(h).solve(I): L from the lower triangle, then
     * L y = e_j forward and L' x = y backward with dot products.
     *
     * @return list<list<float>>|null
     */
    private static function choleskyInverse(array $h): ?array
    {
        $m = count($h);
        // Column-major working copy: $c[$col][$row].
        $c = [];
        for ($j = 0; $j < $m; $j++) {
            $c[$j] = array_column($h, $j);
        }
        for ($j = 0; $j < $m; $j++) {
            for ($k = 0; $k < $j; $k++) {
                $factor = -$c[$k][$j];
                for ($i = $j; $i < $m; $i++) {
                    $c[$j][$i] = $factor * $c[$k][$i] + $c[$j][$i];
                }
            }
            $diagonal = $c[$j][$j];
            if ($diagonal == 0.0 || is_nan($diagonal) || $diagonal < 0.0) {
                return null;
            }
            $root = sqrt($diagonal);
            $c[$j][$j] = $root;
            for ($i = $j + 1; $i < $m; $i++) {
                $c[$j][$i] /= $root;
            }
        }
        $inverse = [];
        for ($col = 0; $col < $m; $col++) {
            $b = array_fill(0, $m, 0.0);
            $b[$col] = 1.0;
            // L y = b: forward, subtracting each solved coefficient times its column.
            for ($i = 0; $i < $m; $i++) {
                $coeff = $b[$i] / $c[$i][$i];
                $b[$i] = $coeff;
                for ($r = $i + 1; $r < $m; $r++) {
                    $b[$r] = -$coeff * $c[$i][$r] + $b[$r];
                }
            }
            // L' x = y: backward, with dot products of the column below the diagonal.
            for ($i = $m - 1; $i >= 0; $i--) {
                $below = array_slice($c[$i], $i + 1);
                $tail = array_slice($b, $i + 1);
                $dot = LeastSquares::dot($below, $tail, 0, count($below));
                $b[$i] = ($b[$i] - $dot) / $c[$i][$i];
            }
            $inverse[$col] = $b;
        }
        // $inverse holds columns; the inverse is symmetric up to rounding, so return rows.
        $rows = [];
        for ($i = 0; $i < $m; $i++) {
            for ($j = 0; $j < $m; $j++) {
                $rows[$i][$j] = $inverse[$j][$i];
            }
        }

        return $rows;
    }
}
