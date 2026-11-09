<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Frontend;
use Phit\Polyfill\Internal\LeastSquares;
use Phit\Polyfill\Internal\Statrs;

/** Gaussian ridge and grouped CV, from the clean-room ridge specification. */
final class Ridge
{
    public static function fit(string $formula, mixed $data, float $lambda, array $penaltyFactor = [], bool $standardize = true, array $bounds = [], string $missing = 'drop', mixed $weights = null): array
    {
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        $options = self::options($design, $lambda, $penaltyFactor, $standardize, $bounds);
        $result = self::estimate($design->x, Frontend::numericResponse($design, 'Phit\\fit_ridge'), $w, $options);
        return self::finish($result, $design);
    }

    public static function cv(string $formula, mixed $data, ?array $lambda = null, ?array $foldid = null, int $nfolds = 10, int $seed = 0, array $penaltyFactor = [], bool $standardize = true, array $bounds = [], string $missing = 'drop', mixed $weights = null): array
    {
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        $options = self::options($design, 0.0, $penaltyFactor, $standardize, $bounds);
        $y = Frontend::numericResponse($design, 'Phit\\cv_ridge');
        $n = count($y);
        if ($lambda === null) {
            $options['lambda'] = 1.0;
            $probe = self::estimate($design->x, $y, $w, $options);
            $upper = max($probe['response_scale'], PHP_FLOAT_MIN) * 100.0;
            $lambda = array_map(fn ($j) => $upper * 10 ** (-4.0 * $j / 99.0), range(0, 99));
        }
        $lambda = array_values($lambda);
        if ($lambda === []) throw self::invalid('lambda must be a nonempty strictly decreasing sequence of finite nonnegative values');
        foreach ($lambda as $j => $value) {
            if ((!is_int($value) && !is_float($value)) || !is_finite((float)$value) || $value < 0 || ($j > 0 && $lambda[$j-1] <= $value)) throw self::invalid('lambda must be a nonempty strictly decreasing sequence of finite nonnegative values');
            $lambda[$j] = (float)$value;
        }
        if ($foldid !== null) {
            if (count($foldid) !== $design->nInput) throw ErrorCode::exception(1002, 'foldid must have one label per original row');
            $original = array_values($foldid);
            foreach ($original as $v) if (!is_int($v)) throw self::invalid('foldid must contain integer labels');
            $foldid = array_map(fn ($i) => $original[$i], $design->retainedRows);
        } else {
            if ($seed < 0 || $nfolds < 3 || $nfolds > $n) throw self::invalid('seed must be nonnegative and nfolds between 3 and the retained row count');
            $foldid = array_fill(0, $n, 0);
            foreach (Sampling::permutation($n, $seed) as $j => $i) $foldid[$i] = $j % $nfolds;
        }
        $folds = array_values(array_unique($foldid)); sort($folds);
        if (count($folds) < 3) throw self::invalid('cv_ridge requires at least three nonempty folds');
        $grouped = intdiv($n, count($folds)) >= 3;
        $losses = $foldWeights = [];
        foreach ($folds as $fold) {
            $tx = $ty = $tw = $test = [];
            foreach ($foldid as $i => $label) {
                if ($label === $fold) $test[] = $i;
                else { $tx[] = $design->x[$i]; $ty[] = $y[$i]; $tw[] = $w[$i]; }
            }
            $sum = array_sum(array_map(fn ($i) => $w[$i], $test));
            $foldLoss = [];
            $observationLoss = array_fill(0, count($test), []);
            foreach ($lambda as $l) {
                $options['lambda'] = $l;
                $fit = self::estimate($tx, $ty, $tw, $options);
                $error = 0.0;
                foreach ($test as $j => $i) {
                    $predicted = self::dot($design->x[$i], $fit['coefficients']);
                    $loss = ($y[$i] - $predicted) ** 2;
                    $observationLoss[$j][] = $loss;
                    $error += $w[$i] * $loss;
                }
                $foldLoss[] = $error / $sum;
            }
            if ($grouped) { $losses[] = $foldLoss; $foldWeights[] = $sum; }
            else foreach ($test as $j => $i) { $losses[] = $observationLoss[$j]; $foldWeights[] = $w[$i]; }
        }
        $sum = array_sum($foldWeights); $cvm = $cvsd = [];
        foreach ($lambda as $j => $l) {
            $mean = 0.0;
            foreach ($losses as $k => $loss) $mean += $foldWeights[$k] * $loss[$j];
            $mean /= $sum; $variance = 0.0;
            foreach ($losses as $k => $loss) $variance += $foldWeights[$k] * ($loss[$j] - $mean) ** 2;
            $cvm[] = $mean; $cvsd[] = sqrt($variance / $sum / (count($losses) - 1));
        }
        $indexMin = array_search(min($cvm), $cvm, true); $index1se = 0;
        foreach ($cvm as $j => $mean) if ($mean <= $cvm[$indexMin] + $cvsd[$indexMin]) { $index1se = $j; break; }
        $options['lambda'] = $lambda[$indexMin];
        $result = self::estimate($design->x, $y, $w, $options);
        return ['lambda'=>$lambda, 'cvm'=>$cvm, 'cvsd'=>$cvsd,
            'cvup'=>array_map(fn ($m,$s) => $m+$s, $cvm, $cvsd), 'cvlo'=>array_map(fn ($m,$s) => $m-$s, $cvm, $cvsd),
            'lambda_min'=>$lambda[$indexMin], 'lambda_1se'=>$lambda[$index1se], 'index_min'=>$indexMin, 'index_1se'=>$index1se,
            'fit'=>self::finish($result, $design)];
    }

    private static function invalid(string $message): \Exception { return ErrorCode::exception(1001, $message); }

    private static function matches(string $pattern, string $name): bool
    {
        $expression = str_replace(['\\*','\\?'], ['.*','.'], preg_quote($pattern, '~'));
        return preg_match('~^'.$expression.'$~uD', $name) === 1;
    }

    private static function options(object $design, float $lambda, array $penaltyFactor, bool $standardize, array $bounds): array
    {
        if ($design->random !== [] || $design->offset !== null) throw ErrorCode::exception(1007, 'ridge does not support random effects or offsets');
        $intercept = array_search('(Intercept)', $design->columns, true);
        $intercept = $intercept === false ? null : $intercept;
        $pf = array_fill(0, count($design->columns), 1.0);
        if ($intercept !== null) $pf[$intercept] = 0.0;
        $limits = array_fill(0, count($pf), [null, null]);
        foreach (['penalty'=>$penaltyFactor, 'bounds'=>$bounds] as $type=>$map) foreach ($map as $pattern=>$value) {
            $matched = [];
            foreach ($design->columns as $j=>$name) if ($j !== $intercept && self::matches((string)$pattern, $name)) $matched[] = $j;
            if ($matched === []) throw self::invalid("pattern '$pattern' matches no predictor columns");
            if ($type === 'penalty') {
                if ((!is_int($value) && !is_float($value)) || !is_finite((float)$value)) throw self::invalid('ridge options must be finite numbers');
                foreach ($matched as $j) $pf[$j] = (float)$value;
            } else {
                if (!is_array($value) || count($value) !== 2) throw self::invalid('each bound must be [lower|null, upper|null]');
                $value = array_values($value);
                foreach ($value as $v) if ($v !== null && ((!is_int($v) && !is_float($v)) || !is_finite((float)$v))) throw self::invalid('ridge options must be finite numbers');
                foreach ($matched as $j) $limits[$j] = array_map(fn ($v) => $v === null ? null : (float)$v, $value);
            }
        }
        return ['lambda'=>$lambda, 'penalty_factor'=>$pf, 'standardize'=>$standardize, 'bounds'=>$limits, 'intercept'=>$intercept];
    }

    private static function dot(array $a, array $b): float
    {
        $sum = 0.0;
        foreach ($a as $i=>$v) $sum += $v*$b[$i];
        return $sum;
    }
    private static function multiply(array $a, array $b): array
    {
        if ($a === []) return [];
        $columns = count($b[0] ?? []); $out = [];
        foreach ($a as $row) {
            $result = array_fill(0, $columns, 0.0);
            foreach ($row as $k=>$value) foreach ($result as $j=>$_) $result[$j] += $value*$b[$k][$j];
            $out[] = $result;
        }
        return $out;
    }
    private static function transpose(array $a, int $columns): array
    {
        $out = [];
        for ($j=0; $j<$columns; $j++) $out[] = array_column($a, $j);
        return $out;
    }
    private static function select(array $a, array $indices): array
    {
        return array_map(fn ($row) => array_map(fn ($j) => $row[$j], $indices), $a);
    }
    /** Column-scaled augmented QR with a unit penalized covariance. */
    private static function solve(array $a, array $y): array
    {
        $n = count($a); $p = count($a[0]);
        if ($p === 0) return [[],[]];
        $columns = $scales = [];
        for ($j=0; $j<$p; $j++) {
            $column = array_column($a, $j); $scale = sqrt(self::dot($column, $column));
            if (!is_finite($scale) || $scale == 0.0) throw ErrorCode::exception(1004, 'unpenalized ridge design is rank deficient');
            $scales[] = $scale; $columns[] = array_map(fn ($v) => $v/$scale, $column);
        }
        $qr = new LeastSquares($columns, $n, $p, $y);
        if ($qr->rankDeficient()) throw ErrorCode::exception(1004, 'unpenalized ridge design is rank deficient at numerical precision');
        return [array_map(fn ($b,$s) => $b/$s, $qr->solution(), $scales), $qr->inverseGram($scales)];
    }

    private static function estimate(array $x, array $y, array $weights, array $options): array
    {
        $n = count($x); $p = count($x[0]); $intercept = $options['intercept'];
        $lambda = $options['lambda']; $pf = $options['penalty_factor']; $bounds = $options['bounds'];
        if (!is_finite($lambda) || $lambda < 0) throw self::invalid('lambda and penalty factors must be finite and nonnegative');
        foreach ($pf as $factor) if (!is_finite($factor) || $factor < 0) throw self::invalid('lambda and penalty factors must be finite and nonnegative');
        foreach ($bounds as [$lo,$hi]) if ($lo !== null && $hi !== null && $lo > $hi) throw self::invalid('bounds must be finite or null, with lower <= upper');
        $maxWeight = max($weights); $sum = array_sum(array_map(fn ($w) => $w/$maxWeight, $weights));
        $q = array_map(fn ($w) => ($w/$maxWeight)/$sum, $weights);
        $centers = $means = $scales = $constant = [];
        for ($j=0; $j<$p; $j++) {
            $mean = 0.0; foreach ($x as $i=>$row) $mean += $q[$i]*$row[$j];
            $centers[$j] = $j === $intercept ? 0.0 : $mean;
            $means[$j] = $intercept !== null ? $centers[$j] : 0.0;
            $variance = 0.0; foreach ($x as $i=>$row) $variance += $q[$i]*($row[$j]-$centers[$j])**2;
            $constant[$j] = $j !== $intercept && (min(array_column($x,$j)) === max(array_column($x,$j)) || $variance == 0.0);
            $scales[$j] = $options['standardize'] && $j !== $intercept && !$constant[$j] ? sqrt($variance) : 1.0;
        }
        $ymean = $intercept !== null ? self::dot($q, $y) : 0.0; $variation = 0.0;
        foreach ($y as $i=>$value) $variation += $q[$i]*($value-$ymean)**2;
        $sy = sqrt($variation);
        if (!is_finite($sy) || $sy == 0.0) throw self::invalid('ridge requires nonzero finite response variation (RMS without an intercept)');
        $effective = $lambda/$sy; $sumFactor = array_sum($pf); $predictors = $p-(int)($intercept !== null);
        $pf = array_map(fn ($v) => $sumFactor > 0 ? $v*$predictors/$sumFactor : 0.0, $pf);
        $design = $augmented = $outcome = [];
        foreach ($x as $i=>$row) {
            $z = []; foreach ($row as $j=>$value) $z[] = $constant[$j] ? 0.0 : ($value-$means[$j])/$scales[$j];
            $design[] = $z; $augmented[] = array_map(fn ($v) => $v*sqrt($q[$i]), $z);
            $outcome[] = sqrt($q[$i])*($y[$i]-$ymean);
        }
        for ($j=0; $j<$p; $j++) { $row = array_fill(0,$p,0.0); $row[$j] = sqrt($effective*$pf[$j]); $augmented[] = $row; $outcome[] = 0.0; }
        $lower = $upper = $available = []; $theta = array_fill(0,$p,0.0); $bounded = false;
        foreach ($bounds as $j=>[$lo,$hi]) {
            $lower[$j] = $lo === null ? -INF : $lo*$scales[$j]; $upper[$j] = $hi === null ? INF : $hi*$scales[$j];
            $bounded = $bounded || $lo !== null || $hi !== null;
            if ($constant[$j] && ($lower[$j] > 0 || $upper[$j] < 0)) throw self::invalid('bounds for constant predictors must include zero');
            if ($lower[$j] === $upper[$j]) $theta[$j] = $lower[$j];
            elseif (!$constant[$j]) $available[] = $j;
        }
        $adjusted = array_map(fn ($row,$y) => $y-self::dot($row,$theta), $augmented, $outcome);
        [$initial,$initialInverse] = self::solve(self::select($augmented,$available), $adjusted);
        foreach ($available as $h=>$j) $theta[$j] = max($lower[$j],min($upper[$j],$initial[$h]));
        $iterations = 0;
        if ($bounded) {
            $residual = array_map(fn ($row,$y) => $y-self::dot($row,$theta), $augmented, $outcome);
            $columns = self::transpose($augmented, $p);
            $diagonals = array_map(fn ($col) => self::dot($col,$col), $columns);
            $tolerance = 1e-12*(1+sqrt(self::dot($outcome,$outcome))); $converged = false;
            for ($iteration=0; $iteration<100000; $iteration++) {
                foreach ($available as $j) {
                    $next = max($lower[$j],min($upper[$j],$theta[$j]+self::dot($columns[$j],$residual)/$diagonals[$j]));
                    $change = $next-$theta[$j]; $theta[$j] = $next;
                    if ($change != 0) foreach ($residual as $i=>$value) $residual[$i] -= $columns[$j][$i]*$change;
                }
                $iterations = $iteration+1; $violation = 0.0;
                foreach ($available as $j) {
                    $gradient = -self::dot($columns[$j],$residual);
                    $g = $theta[$j] <= $lower[$j] ? min(0.0,$gradient) : ($theta[$j] >= $upper[$j] ? max(0.0,$gradient) : $gradient);
                    $violation = max($violation,abs($g)/sqrt($diagonals[$j]));
                }
                if ($violation <= $tolerance) {$converged = true; break;}
            }
            if (!$converged) throw ErrorCode::exception(1006, 'bounded ridge did not satisfy its KKT tolerance');
        }
        $atBound = $free = [];
        for ($j=0; $j<$p; $j++) {
            $atBound[$j] = !$constant[$j] && ($theta[$j] === $lower[$j] || $theta[$j] === $upper[$j]);
            if (!$constant[$j] && !$atBound[$j]) $free[] = $j;
        }
        $inverse = $free === $available ? $initialInverse : self::solve(self::select($augmented,$free),$outcome)[1];
        $weightedFree = [];
        foreach ($design as $i=>$row) $weightedFree[] = array_map(fn ($j) => sqrt($q[$i])*$row[$j], $free);
        $gram = self::multiply(self::transpose($weightedFree,count($free)),$weightedFree);
        $df = 0.0;
        foreach ($inverse as $i=>$row) foreach ($row as $j=>$value) $df += $value*$gram[$j][$i];
        $dfResidual = $n-$df;
        if ($dfResidual <= PHP_FLOAT_EPSILON*$n) throw ErrorCode::exception(1002,'ridge requires positive effective residual degrees of freedom');
        $beta = array_map(fn ($t,$s) => $t/$s, $theta,$scales);
        if ($intercept !== null) { $correction = 0.0; foreach ($beta as $j=>$b) if ($j !== $intercept) $correction += $means[$j]*$b; $beta[$intercept] += $ymean-$correction; }
        $fitted = array_map(fn ($row) => self::dot($row,$beta), $x);
        $residuals = array_map(fn ($a,$b) => $a-$b, $y,$fitted); $rss = 0.0;
        foreach ($residuals as $i=>$r) $rss += $q[$i]*$n*$r**2;
        $variance = $rss/$dfResidual;
        $transform = array_fill(0,$p,array_fill(0,count($free),0.0));
        foreach ($free as $h=>$j) {
            $transform[$j][$h] = 1/$scales[$j];
            if ($intercept !== null && $j !== $intercept) $transform[$intercept][$h] = -$means[$j]/$scales[$j];
        }
        $unit = self::multiply(self::multiply($inverse,$gram),$inverse);
        $covariance = count($free) === 0 ? array_fill(0,$p,array_fill(0,$p,0.0)) : self::multiply(self::multiply($transform,$unit),self::transpose($transform,count($free)));
        $se = $ci = $statistics = $pValues = []; $critical = Statrs::tCritical95($dfResidual);
        foreach ($covariance as $i=>&$row) foreach ($row as $j=>&$value) $value = $atBound[$i] || $atBound[$j] ? null : $value*$variance/$n;
        unset($row,$value);
        foreach ($beta as $j=>$b) {
            $s = $atBound[$j] ? null : sqrt(max(0.0,$covariance[$j][$j])); $se[] = $s;
            $stat = $s === null ? null : fdiv($b,$s); $statistics[] = $stat;
            $pValues[] = $stat === null ? null : (is_nan($stat) ? NAN : 2*Statrs::tSf(abs($stat),$dfResidual));
            $ci[] = $s === null ? null : [$b-$critical*$s,$b+$critical*$s];
        }
        $convertible = $options['standardize'] && $intercept !== null && !$bounded && min($weights) === max($weights);
        foreach ($pf as $j=>$v) if ($j !== $intercept && $v != 1.0) $convertible = false;
        foreach ([...$beta,...$fitted,...$residuals,$rss,$variance] as $value) if (!is_finite($value)) throw ErrorCode::exception(1003,'ridge results exceed floating-point range');
        return ['coefficients'=>$beta,'covariance'=>$covariance,'std_errors'=>$se,'conf_int'=>$ci,'t_statistics'=>$statistics,'p_values'=>$pValues,
            'fitted'=>$fitted,'residuals'=>$residuals,'rss'=>$rss,'sigma'=>sqrt($variance),'df'=>$df,'df_residual'=>$dfResidual,
            'r_squared'=>1-$rss/($n*$sy**2),'penalty_factor'=>$pf,'at_bound'=>$atBound,'lambda'=>$lambda,'lambda_effective'=>$effective,
            'response_scale'=>$sy,'lambda_lm_ridge'=>$convertible ? $n*$effective : null,'iterations'=>$iterations,
            'standardize'=>$options['standardize'],'bounds'=>$bounds];
    }
    private static function finish(array $result, object $design): array
    {
        $result = Frontend::addNames($result,$design);
        return Frontend::finish($result,$design,'Phit\\fit_ridge','augmented_qr_bounded_coordinate_descent');
    }
}
