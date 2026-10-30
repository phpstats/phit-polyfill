<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\ErrorCode;

/**
 * Gaussian mixed models by square-root penalized least squares (stats-core
 * lmm.rs, behind Phit\fit_mixed and Phit\fit_lmer): explicit random-design terms,
 * a dense QR summary of the weighted [Z, X, y - offset], and the relative
 * covariance factors theta optimized by BOBYQA or Nelder-Mead.
 *
 * Matrices are column-major ($m[$col][$row]) unless named as rows.
 *
 * @internal
 */
final class LinearMixed
{
    /** @var list<list<float>> R of the augmented design, columns of min(n, q+p+1) rows */
    private array $r;

    private int $rRows;

    /** @var list<array{k: int, base: int, ids: list<int>, sizes: list<int>, indices: list<int>, parameters: list<array{0: int, 1: int, 2: int}>}> */
    private array $terms;

    private int $n;

    private int $p;

    private int $q;

    /** @var list<float> */
    private array $scales;

    private float $logWeights;

    private bool $reml;

    /**
     * @param  list<list<float>>  $x  rows of the fixed design
     * @param  list<float>  $y
     * @param  list<array{z: list<list<float>>, groups: list<int>, correlated: bool}>  $terms
     * @param  array{reml: bool, optimizer: string, objective_tolerance: float, weights: list<float>, offset: list<float>, start: ?list<float>, max_evaluations: int, tolerance: float}  $options
     * @return array<string, mixed>
     */
    public static function fit(array $x, array $y, array $terms, array $options): array
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
        if ($n <= $p || count($y) !== $n || $rectangular === false) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'require rectangular finite x, finite y, and n > p');
        }
        if ($terms === []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'supply at least one random term');
        }
        $weights = $options['weights'] === [] ? array_fill(0, $n, 1.0) : $options['weights'];
        $offset = $options['offset'] === [] ? array_fill(0, $n, 0.0) : $options['offset'];
        $valid = count($weights) === $n && count($offset) === $n;
        foreach ($weights as $v) {
            $valid = $valid && is_finite($v) && $v > 0.0;
        }
        foreach ($offset as $v) {
            $valid = $valid && is_finite($v);
        }
        if ($valid === false) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'weights must be positive finite and offset finite, one value per row');
        }
        $tolerance = $options['tolerance'];
        if (is_nan($tolerance) || is_infinite($tolerance) || $tolerance <= 0.0 || $tolerance > 0.1) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'tolerance must be finite in (0, 0.1]');
        }
        $metadata = [];
        $q = 0;
        $nt = 0;
        $lower = [];
        $initial = [];
        foreach ($terms as $term) {
            if (count($term['z']) !== $n || count($term['groups']) !== $n || $term['z'][0] === []) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each random term needs n nonempty z rows and n groups');
            }
            $k = count($term['z'][0]);
            $zValid = true;
            foreach ($term['z'] as $row) {
                $zValid = $zValid && count($row) === $k;
                foreach ($row as $v) {
                    $zValid = $zValid && is_finite($v);
                }
            }
            if ($zValid === false) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'random z must be rectangular and finite');
            }
            $ids = [];
            $sizes = [];
            $map = [];
            $indices = [];
            foreach ($term['groups'] as $id) {
                if (! isset($map[$id])) {
                    $map[$id] = count($ids);
                    $ids[] = $id;
                    $sizes[] = 0;
                }
                $sizes[$map[$id]]++;
                $indices[] = $map[$id];
            }
            if (count($ids) < 2 || count($ids) === $n) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each term needs >=2 groups and repeated observations');
            }
            if ($term['correlated'] && $k * count($ids) >= $n) {
                throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'correlated random term requires more observations than random coefficients (lme4 default identifiability check)');
            }
            $parameters = [];
            for ($j = 0; $j < $k; $j++) {
                for ($i = $j; $i < $k; $i++) {
                    if ($i === $j || $term['correlated']) {
                        $parameters[] = [$i, $j, $nt];
                        $nt++;
                        $lower[] = $i === $j ? 0.0 : -1e4;
                        $initial[] = $i === $j ? 1.0 : 0.0;
                    }
                }
            }
            $metadata[] = ['k' => $k, 'base' => $q, 'ids' => $ids, 'sizes' => $sizes, 'indices' => $indices, 'parameters' => $parameters];
            $q += $k * count($map);
        }
        // A dense square-root summary, not lme4's sparse Cholesky: bounded explicitly.
        if ($q > 512 || $n * ($q + $p + 1) > 50_000_000) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'general mixed-model dense design limit exceeded (q <= 512, n*(q+p+1) <= 50000000); use fit_mixed_intercept for large random-intercept fits');
        }
        $scales = [];
        for ($j = 0; $j < $p; $j++) {
            $sum = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $sum += $weights[$i] * ($x[$i][$j] * $x[$i][$j]);
            }
            $scale = sqrt($sum);
            if (is_nan($scale) || is_infinite($scale) || $scale == 0.0) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'fixed design is rank deficient or has invalid scale');
            }
            $scales[] = $scale;
        }
        $width = $q + $p + 1;
        $augmented = array_fill(0, $width, array_fill(0, $n, 0.0));
        for ($i = 0; $i < $n; $i++) {
            $w = sqrt($weights[$i]);
            foreach ($terms as $t => $term) {
                $meta = $metadata[$t];
                for ($j = 0; $j < $meta['k']; $j++) {
                    $augmented[$meta['base'] + $meta['indices'][$i] * $meta['k'] + $j][$i] = $w * $term['z'][$i][$j];
                }
            }
            for ($j = 0; $j < $p; $j++) {
                $augmented[$q + $j][$i] = $w * $x[$i][$j] / $scales[$j];
            }
            $augmented[$q + $p][$i] = $w * ($y[$i] - $offset[$i]);
        }
        foreach ($augmented as $column) {
            foreach ($column as $v) {
                if (is_nan($v) || is_infinite($v)) {
                    throw ErrorCode::exception(ErrorCode::NONFINITE, 'weighted design exceeds floating-point range');
                }
            }
        }
        $data = new self;
        [$data->r] = LeastSquares::householder($augmented, $n, $width);
        $data->rRows = min($n, $width);
        $data->terms = $metadata;
        $data->n = $n;
        $data->p = $p;
        $data->q = $q;
        $data->scales = $scales;
        $logWeights = 0.0;
        foreach ($weights as $w) {
            $logWeights += log($w);
        }
        $data->logWeights = $logWeights;
        $data->reml = $options['reml'];

        // Validate rank and positive variance before any objective reaches the optimizer.
        $zero = array_fill(0, $nt, 0.0);
        $zeroFit = $data->profile($zero);
        $outcomeScale = 0.0;
        for ($i = 0; $i < $n; $i++) {
            $d = $y[$i] - $offset[$i];
            $outcomeScale += $weights[$i] * ($d * $d);
        }
        if ($zeroFit['variance'] <= PHP_FLOAT_EPSILON * PHP_FLOAT_EPSILON * ($n * $n) * $outcomeScale / $n) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'numerically perfect fit: no positive residual-variance optimum');
        }
        $theta = $options['start'] ?? $initial;
        $startValid = count($theta) === $nt;
        foreach ($theta as $j => $v) {
            $startValid = $startValid && is_finite($v) && $v >= $lower[$j] && $v <= 1e4;
        }
        if ($startValid === false) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'start must match theta length and bounds');
        }
        $objectiveTolerance = $options['objective_tolerance'];
        if (is_nan($objectiveTolerance) || is_infinite($objectiveTolerance) || $objectiveTolerance <= 0.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'objectiveTolerance must be positive finite');
        }
        $objective = function (array $t) use ($data): float {
            try {
                return $data->profile($t)['objective'];
            } catch (\Exception $e) {
                return INF;
            }
        };
        if ($options['optimizer'] === 'bobyqa') {
            $config = new Bobyqa\Config($nt);
            $config->rhoBegin = 0.1;
            $config->rhoEnd = $tolerance;
            $config->maxFun = $options['max_evaluations'];
            try {
                $optimizer = new Bobyqa($nt, $config);
            } catch (Bobyqa\BobyqaException $e) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'invalid BOBYQA control: '.$e->getMessage());
            }
            $outcome = $optimizer->minimize($objective, $theta, $lower, array_fill(0, $nt, 1e4));
            if ($outcome->status !== Bobyqa\Status::Converged || $outcome->foundFinite() === false) {
                throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'mixed-model BOBYQA failed: '.$outcome->status->label());
            }
            $evaluations = $outcome->nEval;
        } else {
            $evaluations = 0;
            $budget = $options['max_evaluations'];
            $simplex = [$theta];
            for ($j = 0; $j < $nt; $j++) {
                $vertex = $theta;
                $vertex[$j] = $vertex[$j] + 0.1 > 1e4 ? $vertex[$j] - 0.1 : $vertex[$j] + 0.1;
                $simplex[] = $vertex;
            }
            try {
                $solver = NelderMead::new($simplex)->withSdTolerance($objectiveTolerance);
            } catch (NelderMeadError $e) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'invalid Nelder-Mead control: '.$e->getMessage());
            }
            $cost = function (array $t) use ($data, $lower, &$evaluations, $budget): float {
                if ($evaluations >= $budget) {
                    throw new NelderMeadError('Nelder-Mead evaluation budget exhausted');
                }
                $evaluations++;
                foreach ($t as $j => $v) {
                    if (is_nan($v) || is_infinite($v) || $v < $lower[$j] || $v > 1e4) {
                        return INF;
                    }
                }
                try {
                    return $data->profile($t)['objective'];
                } catch (\Exception $e) {
                    return INF;
                }
            };
            try {
                $result = $solver->run($cost, $budget);
            } catch (NelderMeadError $e) {
                throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'mixed-model Nelder-Mead failed: '.$e->getMessage());
            }
            if ($result->termination !== NelderMeadTermination::SolverConverged || is_nan($result->bestCost) || is_infinite($result->bestCost)) {
                throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'mixed-model Nelder-Mead did not converge');
            }
            if ($result->bestParam === null) {
                throw ErrorCode::exception(ErrorCode::INTERNAL, 'Nelder-Mead returned no parameters');
            }
            $theta = $result->bestParam;
        }
        foreach ($theta as $v) {
            if (abs($v) >= 9999.0) {
                throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'optimum reaches theta search bound');
            }
        }

        return $data->finish($x, $y, $terms, $offset, $theta, $zero, $zeroFit, $nt, $evaluations, $options);
    }

    /** @return array<string, mixed> */
    private function finish(array $x, array $y, array $terms, array $offset, array $theta, array $zero, array $zeroFit, int $nt, int $evaluations, array $options): array
    {
        $n = $this->n;
        $p = $this->p;
        $q = $this->q;
        // The all-zero covariance is an admissible boundary, not a log-variance limit.
        $fit = $this->profile($theta);
        if ($zeroFit['objective'] < $fit['objective']) {
            $theta = $zero;
            $fit = $zeroFit;
        }
        $beta = [];
        for ($j = 0; $j < $p; $j++) {
            $beta[] = $fit['solution'][$q + $j] / $this->scales[$j];
        }
        $inverse = self::solveUpperIdentity($fit['fixedR']);
        if ($inverse === null) {
            throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'cannot solve fixed-effect covariance factor');
        }
        $scaledCov = Dense::multiply($inverse, Dense::transpose($inverse));
        $covariance = [];
        $stdErrors = [];
        $z = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $covariance[$i][$j] = $fit['variance'] * $scaledCov[$j][$i] / $this->scales[$i] / $this->scales[$j];
            }
        }
        for ($j = 0; $j < $p; $j++) {
            $stdErrors[] = sqrt($covariance[$j][$j]);
            $z[] = $beta[$j] / $stdErrors[$j];
        }
        // effects = Lambda * u, as a matrix-vector product (gemv).
        $u = array_slice($fit['solution'], 0, $q);
        $effects = Dense::multiplyVector($fit['lambda'], $u);
        $fixed = [];
        for ($i = 0; $i < $n; $i++) {
            $sum = 0.0;
            for ($j = 0; $j < $p; $j++) {
                $sum += $x[$i][$j] * $beta[$j];
            }
            $fixed[] = $offset[$i] + $sum;
        }
        $fitted = [];
        for ($i = 0; $i < $n; $i++) {
            $random = 0.0;
            foreach ($terms as $t => $term) {
                $m = $this->terms[$t];
                $sum = 0.0;
                for ($j = 0; $j < $m['k']; $j++) {
                    $sum += $term['z'][$i][$j] * $effects[$m['base'] + $m['indices'][$i] * $m['k'] + $j];
                }
                $random += $sum;
            }
            $fitted[] = $fixed[$i] + $random;
        }
        // lme4 condVar: sigma^2 Lambda (L L')^-1 Lambda'; with R_11'R_11 = L L',
        // solve R_11' S = Lambda' so the product is S'S.
        $spread = self::solveTransposedUpper($fit['randomR'], Dense::transpose($fit['lambda']));
        if ($spread === null) {
            throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'cannot solve random-effect conditional covariance factor');
        }
        $singular = false;
        $termResults = [];
        foreach ($this->terms as $m) {
            $k = $m['k'];
            $l = array_fill(0, $k, array_fill(0, $k, 0.0));
            foreach ($m['parameters'] as [$i, $j, $t]) {
                $l[$j][$i] = $theta[$t];
            }
            for ($j = 0; $j < $k; $j++) {
                $singular = $singular || $l[$j][$j] < 1e-4;
            }
            $scaledL = [];
            foreach ($l as $j => $column) {
                foreach ($column as $i => $v) {
                    $scaledL[$j][$i] = $fit['variance'] * $v;
                }
            }
            $termCovariance = Dense::multiply($scaledL, Dense::transpose($l));
            $groupEffects = [];
            $condVar = [];
            foreach ($m['ids'] as $g => $id) {
                $vector = [];
                for ($j = 0; $j < $k; $j++) {
                    $vector[] = $effects[$m['base'] + $g * $k + $j];
                }
                $groupEffects[] = $vector;
                $base = $m['base'] + $g * $k;
                $block = [];
                for ($i = 0; $i < $k; $i++) {
                    for ($j = 0; $j < $k; $j++) {
                        $a = $spread[$base + $i];
                        $b = $spread[$base + $j];
                        $block[$i][$j] = $fit['variance'] * LeastSquares::dot($a, $b, 0, count($a));
                    }
                }
                $condVar[] = $block;
            }
            $termResults[] = [
                'covariance' => Dense::transpose($termCovariance),
                'effects' => $groupEffects,
                'effects_cond_var' => $condVar,
                'group_sizes' => $m['sizes'],
            ];
        }
        $nParameters = $p + $nt + 1;
        $pValues = [];
        $confInt = [];
        for ($j = 0; $j < $p; $j++) {
            $pValues[] = 2.0 * Statrs::normalSf(abs($z[$j]));
            $confInt[] = [$beta[$j] - Statrs::Z_975 * sqrt($covariance[$j][$j]), $beta[$j] + Statrs::Z_975 * sqrt($covariance[$j][$j])];
        }
        $residuals = [];
        $marginal = [];
        for ($i = 0; $i < $n; $i++) {
            $residuals[] = $y[$i] - $fitted[$i];
            $marginal[] = $y[$i] - $fixed[$i];
        }
        $finite = true;
        $check = [$beta, ...$covariance, $fitted, $stdErrors, $z];
        foreach ($termResults as $term) {
            foreach ($term['effects_cond_var'] as $block) {
                foreach ($block as $row) {
                    $check[] = $row;
                }
            }
        }
        foreach ($check as $values) {
            foreach ($values as $v) {
                $finite = $finite && is_finite($v);
            }
        }
        if ($finite === false) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'mixed-model output exceeds floating-point range');
        }

        return [
            'coefficients' => $beta,
            'covariance' => $covariance,
            'std_errors' => $stdErrors,
            'z_statistics' => $z,
            'p_values' => $pValues,
            'conf_int' => $confInt,
            'fixed_fitted' => $fixed,
            'fitted' => $fitted,
            'residuals' => $residuals,
            'marginal_residuals' => $marginal,
            'terms' => $termResults,
            'theta' => $theta,
            'residual_variance' => $fit['variance'],
            'sigma' => sqrt($fit['variance']),
            'log_likelihood' => -$fit['objective'],
            'aic' => 2.0 * $nParameters + 2.0 * $fit['objective'],
            'bic' => $nParameters * log($n) + 2.0 * $fit['objective'],
            'df_residual' => $n - $nParameters,
            'n_parameters' => $nParameters,
            'evaluations' => $evaluations + 2,
            'optimizer' => $options['optimizer'],
            'singular' => $singular,
        ];
    }

    /**
     * nalgebra's solve_upper_triangular(I): columns of R^-1 by back-substitution,
     * or null for a zero diagonal.
     *
     * @param  list<list<float>>  $r  upper triangular, columns
     * @return list<list<float>>|null
     */
    private static function solveUpperIdentity(array $r): ?array
    {
        $m = count($r);
        $out = [];
        for ($c = 0; $c < $m; $c++) {
            $b = array_fill(0, $m, 0.0);
            $b[$c] = 1.0;
            for ($i = $m - 1; $i >= 0; $i--) {
                if ($r[$i][$i] == 0.0) {
                    return null;
                }
                $coeff = $b[$i] / $r[$i][$i];
                $b[$i] = $coeff;
                for ($k = 0; $k < $i; $k++) {
                    $b[$k] = -$coeff * $r[$i][$k] + $b[$k];
                }
            }
            $out[] = $b;
        }

        return $out;
    }

    /**
     * Solve R' S = B with R upper triangular, as nalgebra's
     * R.transpose().solve_lower_triangular(B): forward substitution with
     * column updates, L(r, i) = R(i, r).
     *
     * @param  list<list<float>>  $r  columns
     * @param  list<list<float>>  $b  columns
     * @return list<list<float>>|null
     */
    private static function solveTransposedUpper(array $r, array $b): ?array
    {
        $m = count($r);
        foreach ($b as $c => $column) {
            for ($i = 0; $i < $m; $i++) {
                $diagonal = $r[$i][$i];
                if ($diagonal == 0.0) {
                    return null;
                }
                $coeff = $column[$i] / $diagonal;
                $column[$i] = $coeff;
                for ($k = $i + 1; $k < $m; $k++) {
                    $column[$k] = -$coeff * $r[$k][$i] + $column[$k];
                }
            }
            $b[$c] = $column;
        }

        return $b;
    }

    /** Lambda: the block-diagonal relative covariance factor for theta (q×q, columns). */
    private function lambda(array $theta): array
    {
        $l = array_fill(0, $this->q, array_fill(0, $this->q, 0.0));
        foreach ($this->terms as $term) {
            foreach ($term['ids'] as $g => $id) {
                $offset = $term['base'] + $g * $term['k'];
                foreach ($term['parameters'] as [$i, $j, $t]) {
                    $l[$offset + $j][$offset + $i] = $theta[$t];
                }
            }
        }

        return $l;
    }

    /**
     * The profiled objective at theta, from the QR of [R_Z Lambda, R_X, R_y; I, 0, 0].
     *
     * @return array{solution: list<float>, randomR: list<list<float>>, fixedR: list<list<float>>, variance: float, objective: float, lambda: list<list<float>>}
     */
    private function profile(array $theta): array
    {
        $q = $this->q;
        $p = $this->p;
        $nr = $this->rRows;
        $lambda = $this->lambda($theta);
        $u = Dense::multiply(array_slice($this->r, 0, $q), $lambda, $nr);
        $rows = $nr + $q;
        $a = [];
        for ($j = 0; $j < $q; $j++) {
            $column = array_fill(0, $rows, 0.0);
            for ($i = 0; $i < $nr; $i++) {
                $column[$i] = $u[$j][$i];
            }
            $column[$nr + $j] = 1.0;
            $a[] = $column;
        }
        for ($j = 0; $j < $p; $j++) {
            $column = array_fill(0, $rows, 0.0);
            for ($i = 0; $i < $nr; $i++) {
                $column[$i] = $this->r[$q + $j][$i];
            }
            $a[] = $column;
        }
        $b = array_fill(0, $rows, 0.0);
        for ($i = 0; $i < $nr; $i++) {
            $b[$i] = $this->r[$q + $p][$i];
        }
        $m = $q + $p;
        [$r, $qtb] = LeastSquares::householder($a, $rows, $m, $b);
        // Solve the leading (q+p)×(q+p) triangle by nalgebra's back-substitution.
        $solution = array_slice($qtb, 0, $m);
        for ($i = $m - 1; $i >= 0; $i--) {
            $diagonal = $r[$i][$i];
            if ($diagonal == 0.0) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'mixed-model fixed design is rank deficient');
            }
            $coeff = $solution[$i] / $diagonal;
            $solution[$i] = $coeff;
            for ($k = 0; $k < $i; $k++) {
                $solution[$k] = -$coeff * $r[$i][$k] + $solution[$k];
            }
        }
        // QR keeps residual directions; never form y'y - beta'X'y.
        $fitted = Dense::multiplyVector($a, $solution);
        $residual = [];
        for ($i = 0; $i < $rows; $i++) {
            $residual[] = $b[$i] - $fitted[$i];
        }
        $rss = LeastSquares::dot($residual, $residual, 0, $rows);
        $df = $this->n - ($this->reml ? $p : 0);
        $variance = $rss / $df;
        if (is_nan($variance) || is_infinite($variance) || $variance <= 0.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'no positive finite residual variance');
        }
        $fixedR = [];
        for ($j = 0; $j < $p; $j++) {
            $fixedR[] = array_slice($r[$q + $j], $q, $p);
        }
        $singular = (new Svd($fixedR, false, false))->singularValues;
        $cutoff = PHP_FLOAT_EPSILON * max($this->n, $p) * $singular[0];
        foreach ($singular as $s) {
            if (is_nan($s) || is_infinite($s) || $s <= $cutoff) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'mixed-model fixed design/information is numerically rank deficient');
            }
        }
        // The random block of R satisfies R_U'R_U = I + U'U.
        $sum = 0.0;
        for ($j = 0; $j < $q; $j++) {
            $sum += log(abs($r[$j][$j]));
        }
        $logdetV = 2.0 * $sum - $this->logWeights;
        $logdetX = 0.0;
        if ($this->reml) {
            $fixedSum = 0.0;
            for ($j = 0; $j < $p; $j++) {
                $fixedSum += log(abs($fixedR[$j][$j]));
            }
            $scaleSum = 0.0;
            foreach ($this->scales as $s) {
                $scaleSum += log($s);
            }
            $logdetX = 2.0 * ($fixedSum + $scaleSum);
        }
        $objective = 0.5 * ($df * (log(2.0 * M_PI) + log($variance) + 1.0) + $logdetV + $logdetX);
        if (is_nan($objective) || is_infinite($objective)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'mixed likelihood exceeds floating-point range');
        }
        $randomR = [];
        for ($j = 0; $j < $q; $j++) {
            $randomR[] = array_slice($r[$j], 0, $q);
        }

        return ['solution' => $solution, 'randomR' => $randomR, 'fixedR' => $fixedR, 'variance' => $variance, 'objective' => $objective, 'lambda' => $lambda];
    }
}
