<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Frontend;
use Phit\Polyfill\Internal\LeastSquares;
use Phit\Polyfill\Internal\Statrs;
use Phit\Polyfill\Distribution\FisherF;
use Phit\Polyfill\Distribution\ChiSquare;

/** Two-stage least squares, implemented from specs/cleanroom/iv.md. */
final class Iv
{
    /** Split only top-level bars, preserving random-effect bars and quoted names. */
    private static function split(string $formula): array
    {
        $depth = 0;
        $quote = null;
        $escaped = false;
        $bar = null;
        $tilde = null;
        for ($i = 0; $i < strlen($formula); $i++) {
            $c = $formula[$i];
            if ($quote !== null) {
                if ($escaped) $escaped = false;
                elseif ($c === '\\') $escaped = true;
                elseif ($c === $quote) $quote = null;
                continue;
            }
            if (in_array($c, ['`', "'", '"'], true)) $quote = $c;
            elseif ($c === '(') $depth++;
            elseif ($c === ')') $depth--;
            elseif ($c === '~' && $depth === 0 && $bar !== null) {
                if ($tilde !== null) throw ErrorCode::exception(1001, 'IV formula needs one ~ after |');
                $tilde = $i;
            }
            elseif ($c === '|' && $depth === 0) {
                if ($bar !== null) throw ErrorCode::exception(1001, 'IV formula requires exactly one top-level |');
                $bar = $i;
            }
            if ($depth < 0) throw ErrorCode::exception(1001, 'unbalanced IV formula parentheses');
        }
        if ($depth !== 0 || $quote !== null) throw ErrorCode::exception(1001, 'unbalanced IV formula');
        if ($bar === null) throw ErrorCode::exception(1001, "IV formula must be 'y ~ exog | endog ~ instruments'");
        if ($tilde === null) throw ErrorCode::exception(1001, 'IV formula needs endog ~ instruments after |');
        $endog = trim(substr($formula, $bar + 1, $tilde - $bar - 1));
        $instruments = trim(substr($formula, $tilde + 1));
        if ($endog === '' || $instruments === '') throw ErrorCode::exception(1001, 'IV formula needs endog ~ instruments after |');
        return [trim(substr($formula, 0, $bar)), $endog, $instruments];
    }

    public static function fit(string $formula, mixed $data, string $missing = 'drop', mixed $weights = null, mixed $vcov = 'iid'): array
    {
        Frontend::stringArgument($formula, 'formula');
        if ($vcov !== 'iid') throw ErrorCode::exception(1007, "fit_iv currently supports only vcov: 'iid'");
        [$left, $endog, $instruments] = self::split($formula);
        $structuralFormula = "$left + $endog";
        [$a, $w] = Frontend::prepare($structuralFormula, $data, 'drop', $weights, true);
        [$b] = Frontend::prepare("~ $instruments", $data, 'drop', $weights, true);
        Frontend::rejectRandom($a, 'Phit\\fit_iv');
        Frontend::rejectRandom($b, 'Phit\\fit_iv');
        $shared = array_fill(0, $a->nInput, null);
        foreach ($a->retainedRows as $i => $row) {
            if (in_array($row, $b->retainedRows, true)) $shared[$row] = $w[$i];
        }
        [$design, $w] = Frontend::prepare($structuralFormula, $data, $missing, $shared, true);
        [$exog] = Frontend::prepare($left, $data, $missing, $shared, true, true);
        [$excluded] = Frontend::prepare("$left + ($instruments) + ".($exog->intercept ? '1' : '0'), $data, $missing, $shared, true);
        foreach ([$design, $exog, $excluded] as $d) {
            Frontend::rejectRandom($d, 'Phit\\fit_iv');
            if ($d->offset !== null) throw ErrorCode::exception(1007, 'fit_iv does not support offsets');
        }
        if ($design->intercept !== $exog->intercept) throw ErrorCode::exception(1001, 'specify intercept removal in the exogenous formula part');
        foreach ($exog->columns as $name) {
            if (!in_array($name, $design->columns, true) || !in_array($name, $excluded->columns, true)) throw ErrorCode::exception(1001, 'endogenous and instrument parts must retain every exogenous column');
        }
        $endogenous = array_keys(array_diff($design->columns, $exog->columns));
        if ($endogenous === []) throw ErrorCode::exception(1001, 'IV formula must add at least one endogenous design column');
        $excludedIndices = [];
        foreach ($excluded->columns as $j => $name) {
            if ($name !== '(Intercept)' && !in_array($name, $exog->columns, true)) $excludedIndices[] = $j;
        }
        $instrumentNames = array_map(fn ($j) => $excluded->columns[$j], $excludedIndices);
        $zNames = $excluded->columns;
        $z = $excluded->x;
        $output = self::estimate($design->x, Frontend::numericResponse($design, 'Phit\\fit_iv'), $z, $endogenous, $w, $design->intercept);
        $stages = [];
        foreach ($endogenous as $j => $index) {
            $stage = $output['first_stage'][$j];
            $stage['coefficients'] = array_combine($zNames, $stage['coefficients']);
            $stages[$design->columns[$index]] = $stage;
        }
        $output['first_stage'] = $stages;
        $output = Frontend::addNames($output, $design);
        $output['formula'] = $formula;
        $output['structural_formula'] = $structuralFormula;
        $output['endogenous'] = array_map(fn ($j) => $design->columns[$j], $endogenous);
        $output['instruments'] = $instrumentNames;
        $output['vcov'] = 'iid';
        return Frontend::finish($output, $design, 'Phit\\fit_iv', 'scaled_qr_2sls');
    }

    /** Scaled QR solve, unit covariance, fitted values, and residuals. */
    private static function leastSquares(array $x, array $y): array
    {
        $n = count($x);
        $p = count($x[0]);
        if ($p === 0) return [[], [], array_fill(0, $n, 0.0), $y];
        $columns = $scales = [];
        for ($j = 0; $j < $p; $j++) {
            $column = array_column($x, $j);
            $scale = sqrt(LeastSquares::dot($column, $column, 0, $n));
            if (!is_finite($scale) || $scale == 0.0) throw ErrorCode::exception(1004, 'IV design is rank deficient');
            $scales[] = $scale;
            $columns[] = array_map(fn ($v) => $v / $scale, $column);
        }
        $qr = new LeastSquares($columns, $n, $p, $y);
        if ($qr->rankDeficient()) throw ErrorCode::exception(1004, 'IV design is rank deficient at numerical precision');
        $beta = array_map(fn ($b, $s) => $b / $s, $qr->solution(), $scales);
        $fit = array_map(fn ($row) => self::dot($row, $beta), $x);
        return [$beta, $qr->inverseGram($scales), $fit, self::subtract($y, $fit)];
    }
    private static function dot(array $a, array $b): float
    {
        $sum = 0.0;
        foreach ($a as $i => $v) $sum += $v * $b[$i];
        return $sum;
    }
    private static function subtract(array $a, array $b): array
    {
        return array_map(fn ($a, $b) => $a - $b, $a, $b);
    }
    private static function select(array $x, array $indices): array
    {
        return array_map(fn ($row) => array_map(fn ($j) => $row[$j], $indices), $x);
    }
    private static function fTest(float $f, int $df1, int $df2): array
    {
        return ['statistic' => $f, 'df1' => $df1, 'df2' => $df2, 'p_value' => is_nan($f) ? NAN : FisherF::sf($f, (float)$df1, (float)$df2)];
    }
    private static function estimate(array $x, array $y, array $z, array $endogenous, array $weights, bool $intercept): array
    {
        $n = count($x); $k = count($x[0]); $l = count($z[0]); $m = count($endogenous);
        if ($l < $k) throw ErrorCode::exception(1004, 'under-identified IV model: fewer instruments than regressors');
        if ($n <= $l || $n <= $k) throw ErrorCode::exception(1002, 'IV requires positive instrument and structural residual degrees of freedom');
        $wx = $wz = $wy = [];
        foreach ($x as $i => $row) {
            $root = sqrt($weights[$i]);
            $wx[] = array_map(fn ($v) => $v * $root, $row);
            $wz[] = array_map(fn ($v) => $v * $root, $z[$i]);
            $wy[] = $y[$i] * $root;
        }
        $projected = array_fill(0, $n, []);
        $firstBeta = [];
        for ($j = 0; $j < $k; $j++) {
            [$b, , $fit] = self::leastSquares($wz, array_column($wx, $j));
            $firstBeta[$j] = $b;
            foreach ($fit as $i => $v) $projected[$i][] = $v;
        }
        [$beta, $unit] = self::leastSquares($projected, $wy);
        $fitted = array_map(fn ($row) => self::dot($row, $beta), $x);
        $residuals = self::subtract($y, $fitted);
        $rss = 0.0; $mean = 0.0;
        foreach ($y as $i => $v) { $rss += $weights[$i] * $residuals[$i] ** 2; $mean += $weights[$i] * $v; }
        $mean = $intercept ? $mean / array_sum($weights) : 0.0;
        $tss = 0.0;
        foreach ($y as $i => $v) $tss += $weights[$i] * ($v - $mean) ** 2;
        $df = $n - $k; $variance = $rss / $df;
        $covariance = array_map(fn ($r) => array_map(fn ($v) => $v * $variance, $r), $unit);
        $se = $t = $p = $ci = [];
        $critical = Statrs::tCritical95((float)$df);
        foreach ($beta as $j => $b) {
            $se[] = $s = sqrt($covariance[$j][$j]);
            $t[] = $stat = fdiv($b, $s);
            $p[] = is_nan($stat) ? NAN : 2 * Statrs::tSf(abs($stat), (float)$df);
            $ci[] = [$b - $critical * $s, $b + $critical * $s];
        }
        foreach ([...$beta, ...$fitted, ...$residuals, $rss, $tss, ...array_merge(...$covariance)] as $v) {
            if (!is_finite($v)) throw ErrorCode::exception(1003, 'IV calculation exceeds floating-point range');
        }
        $exog = array_values(array_diff(range(0, $k - 1), $endogenous));
        $c = self::select($wx, $exog);
        $d = self::select($wx, $endogenous);
        $h = self::select($projected, $endogenous);
        $v = array_map(self::subtract(...), $d, $h);
        $q = $l - count($exog); $first = $warnings = [];
        foreach ($endogenous as $j => $index) {
            [, , , $restricted] = self::leastSquares($c, array_column($d, $j));
            $marginal = self::dot($restricted, $restricted);
            $vj = array_column($v, $j); $error = self::dot($vj, $vj);
            $test = self::fTest(fdiv(max(0.0, $marginal - $error) / $q, $error / ($n - $l)), $q, $n - $l);
            $sw = $m === 1 ? $test : null; $shea = null;
            $others = array_values(array_diff(range(0, $m - 1), [$j]));
            if ($m > 1) {
                $otherH = self::select($h, $others);
                $control = array_map(fn ($a, $b) => [...$a, ...$b], $c, $otherH);
                try {
                    [$delta, , , $conditional] = self::leastSquares($control, array_column($d, $j));
                    [, , $signal] = self::leastSquares($wz, $conditional);
                    $correction = array_slice($delta, count($exog));
                    $noise = self::subtract($vj, array_map(fn ($row) => self::dot($row, $correction), self::select($v, $others)));
                    $f = fdiv(self::dot($signal, $signal) / ($q - $m + 1), self::dot($noise, $noise) / ($n - $l));
                    if (is_finite($f)) $sw = self::fTest($f, $q - $m + 1, $n - $l);
                    $actual = array_map(fn ($a, $b) => [...$a, ...$b], $c, self::select($d, $others));
                    [, , , $actualR] = self::leastSquares($actual, array_column($d, $j));
                    [, , , $fittedR] = self::leastSquares($control, array_column($h, $j));
                    $ratio = fdiv(self::dot($fittedR, $fittedR), self::dot($actualR, $actualR));
                    if (is_finite($ratio)) $shea = max(0.0, min(1.0, $ratio));
                } catch (\Exception $e) { if ($e->getCode() !== 1004) throw $e; }
            }
            if ($sw === null) $warnings[] = "Sanderson-Windmeijer F unavailable for endogenous column $index: singular conditional regression or undefined variance";
            $first[] = ['f_statistic' => $test['statistic'], 'df1' => $q, 'df2' => $n - $l, 'p_value' => $test['p_value'],
                'partial_r_squared' => max(0.0, min(1.0, 1 - fdiv($error, $marginal))), 'shea_partial_r_squared' => $shea,
                'sw_f_statistic' => $sw['statistic'] ?? null, 'sw_df1' => $q - $m + 1, 'sw_df2' => $n - $l,
                'sw_p_value' => $sw['p_value'] ?? null, 'coefficients' => $firstBeta[$index]];
        }
        $wu = null;
        if ($n > $k + $m) {
            try {
                [, , , $r] = self::leastSquares($wx, $wy);
                $augmented = array_map(fn ($a, $b) => [...$a, ...$b], $wx, $v);
                [, , , $u] = self::leastSquares($augmented, $wy);
                $wu = self::fTest(fdiv(max(0.0, self::dot($r, $r) - self::dot($u, $u)) / $m, self::dot($u, $u) / ($n - $k - $m)), $m, $n - $k - $m);
            } catch (\Exception $e) { if ($e->getCode() !== 1004) throw $e; }
        }
        if ($wu === null) $warnings[] = 'Wu-Hausman test unavailable: insufficient degrees of freedom or singular augmented regression';
        $sargan = null;
        if ($l > $k) {
            $wr = array_map(fn ($r, $w) => $r * sqrt($w), $residuals, $weights);
            [, , , $r] = self::leastSquares($wz, $wr);
            $stat = $n * max(0.0, min(1.0, 1 - fdiv(self::dot($r, $r), $rss)));
            $sargan = ['statistic' => $stat, 'df' => $l - $k, 'p_value' => ChiSquare::sf($stat, (float)($l - $k))];
        }
        $r2 = 1 - fdiv($rss, $tss);
        return ['coefficients' => $beta, 'covariance' => $covariance, 'std_errors' => $se, 't_statistics' => $t, 'p_values' => $p,
            'conf_int' => $ci, 'fitted' => $fitted, 'residuals' => $residuals, 'rss' => $rss, 'residual_variance' => $variance,
            'sigma' => sqrt($variance), 'r_squared' => $r2, 'adj_r_squared' => 1 - (1 - $r2) * ($n - (int)$intercept) / $df,
            'df_residual' => $df, 'log_likelihood' => -.5 * $n * (log(2 * M_PI) + log($rss / $n) + 1) + .5 * array_sum(array_map(log(...), $weights)),
            'first_stage' => $first, 'wu_hausman' => $wu ?? ['statistic' => null, 'df1' => null, 'df2' => null, 'p_value' => null],
            'sargan' => $sargan, 'iv_warnings' => $warnings];
    }
}
