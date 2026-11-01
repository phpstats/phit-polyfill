<?php

// Polyfill vs native phit timings for Phit\fit_ols_matrix and Phit\fit_logistic_matrix.
//
// Run with the extension and JIT enabled, for example:
//   php85 -d extension=.ext/libphit.dylib -d opcache.enable_cli=1 \
//         -d opcache.jit=tracing -d opcache.jit_buffer_size=64M benchmarks/run.php
// Options: --sizes=10000,100000 --p=10 --repeat=5 --seed=5 --json=benchmarks/results.json
//
// With --json the results are also written as JSON: run context, then per case
// the timings of each implementation (median and minimum seconds) and whether
// the polyfill's estimates agree with the extension's.

declare(strict_types=1);

require __DIR__.'/../vendor/autoload.php';

use Phit\Polyfill\Logistic;
use Phit\Polyfill\Ols;

$options = getopt('', ['sizes::', 'p::', 'repeat::', 'seed::', 'json::']);
$sizes = array_map(intval(...), explode(',', $options['sizes'] ?? '10000,100000'));
$p = (int) ($options['p'] ?? 10);
$repeat = (int) ($options['repeat'] ?? 5);
$seed = (int) ($options['seed'] ?? 5);
$json = $options['json'] ?? null;
$native = extension_loaded('phit');
$opcache = function_exists('opcache_get_status') ? (opcache_get_status(false) ?: []) : [];
$jit = (bool) ($opcache['jit']['on'] ?? false);

function gaussian(): float
{
    $u = (mt_rand() + 1) / (mt_getrandmax() + 2.0);
    $v = mt_rand() / (mt_getrandmax() + 1.0);

    return sqrt(-2.0 * log($u)) * cos(2.0 * M_PI * $v);
}

/**
 * Run $f $repeat times after one warm-up call.
 *
 * @return array{result: mixed, median_seconds: float, min_seconds: float}
 */
function measure(callable $f, int $repeat): array
{
    $result = $f();
    $times = [];
    for ($i = 0; $i < $repeat; $i++) {
        $start = hrtime(true);
        $f();
        $times[] = (hrtime(true) - $start) / 1e9;
    }
    sort($times);
    $middle = intdiv(count($times), 2);
    $median = count($times) % 2 ? $times[$middle] : ($times[$middle - 1] + $times[$middle]) / 2;

    return ['result' => $result, 'median_seconds' => $median, 'min_seconds' => $times[0]];
}

/**
 * Plain-PHP baseline: coefficients from the normal equations (X'X) b = X'y by
 * Cholesky, without standard errors, rank checks, or any other output.
 */
function normal_equations(array $x, array $y): array
{
    $p = count($x[0]) + 1;
    $xtx = array_fill(0, $p, array_fill(0, $p, 0.0));
    $xty = array_fill(0, $p, 0.0);
    foreach ($x as $i => $row) {
        $r = [1.0, ...$row];
        for ($a = 0; $a < $p; $a++) {
            $ra = $r[$a];
            $xty[$a] += $ra * $y[$i];
            for ($b = 0; $b <= $a; $b++) {
                $xtx[$a][$b] += $ra * $r[$b];
            }
        }
    }
    $l = array_fill(0, $p, array_fill(0, $p, 0.0));
    for ($a = 0; $a < $p; $a++) {
        for ($b = 0; $b <= $a; $b++) {
            $sum = $xtx[$a][$b];
            for ($k = 0; $k < $b; $k++) {
                $sum -= $l[$a][$k] * $l[$b][$k];
            }
            $l[$a][$b] = $a === $b ? sqrt($sum) : $sum / $l[$b][$b];
        }
    }
    $z = [];
    for ($a = 0; $a < $p; $a++) {
        $sum = $xty[$a];
        for ($k = 0; $k < $a; $k++) {
            $sum -= $l[$a][$k] * $z[$k];
        }
        $z[$a] = $sum / $l[$a][$a];
    }
    $beta = array_fill(0, $p, 0.0);
    for ($a = $p - 1; $a >= 0; $a--) {
        $sum = $z[$a];
        for ($k = $a + 1; $k < $p; $k++) {
            $sum -= $l[$k][$a] * $beta[$k];
        }
        $beta[$a] = $sum / $l[$a][$a];
    }

    return ['coefficients' => $beta];
}

/** Largest coefficient difference in units of the extension's standard errors. */
function coefficient_difference(array $polyfill, array $native): float
{
    $max = 0.0;
    foreach ($native['coefficients'] as $j => $b) {
        $max = max($max, abs($polyfill['coefficients'][$j] - $b) / (abs($b) + $native['std_errors'][$j]));
    }

    return $max;
}

$cases = [];
printf("PHP %s, JIT %s, extension %s, p = %d predictors + intercept, median of %d\n\n",
    PHP_VERSION, $jit ? 'on' : 'off', $native ? 'loaded' : 'not loaded', $p, $repeat);
printf("%-10s %9s %14s %14s %8s\n", 'model', 'n', 'polyfill ms', 'extension ms', 'ratio');

foreach ($sizes as $n) {
    mt_srand($seed + $n);
    $x = [];
    $continuous = [];
    $binary = [];
    for ($i = 0; $i < $n; $i++) {
        $row = [];
        $eta = 0.3;
        for ($j = 0; $j < $p; $j++) {
            $value = gaussian();
            $row[] = $value;
            $eta += $value * (($j % 3) - 1) * 0.2;
        }
        $x[] = $row;
        $continuous[] = $eta + gaussian();
        $binary[] = (mt_rand() / mt_getrandmax()) < 1 / (1 + exp(-$eta)) ? 1 : 0;
    }
    $models = [
        ['ols', Ols::fit(...), 'Phit\fit_ols_matrix', $continuous],
        ['logistic', Logistic::fit(...), 'Phit\fit_logistic_matrix', $binary],
    ];
    foreach ($models as [$name, $polyfill, $fn, $y]) {
        $case = ['model' => $name, 'function' => $fn, 'n' => $n, 'p' => $p, 'implementations' => []];
        try {
            $ours = measure(fn () => $polyfill($x, $y), $repeat);
            $theirs = $native ? measure(fn () => $fn($x, $y), $repeat) : null;
        } catch (Exception $e) {
            printf("%-10s %9d   failed with code %d: %s\n", $name, $n, $e->getCode(), $e->getMessage());

            continue;
        }
        $case['implementations']['polyfill'] = ['median_seconds' => $ours['median_seconds'], 'min_seconds' => $ours['min_seconds']];
        if ($theirs !== null) {
            $case['implementations']['extension'] = ['median_seconds' => $theirs['median_seconds'], 'min_seconds' => $theirs['min_seconds']];
            $difference = coefficient_difference($ours['result'], $theirs['result']);
            $case['agreement'] = [
                'max_coefficient_difference_in_standard_errors' => $difference,
                'same_result_keys' => array_keys($ours['result']) === array_keys($theirs['result']),
                'agrees' => $difference <= 1e-7 && array_keys($ours['result']) === array_keys($theirs['result']),
            ];
        }
        if ($name === 'ols') {
            $plain = measure(fn () => normal_equations($x, $y), $repeat);
            $case['implementations']['plain_php_normal_equations'] = [
                'median_seconds' => $plain['median_seconds'],
                'min_seconds' => $plain['min_seconds'],
                'note' => 'coefficients only, by Cholesky on X\'X; no inference or rank checks',
            ];
        }
        $cases[] = $case;
        printf("%-10s %9d %14.1f %14.1f %8s\n", $name, $n, 1e3 * $ours['median_seconds'],
            $theirs ? 1e3 * $theirs['median_seconds'] : NAN,
            $theirs ? sprintf('%.1fx', $ours['median_seconds'] / $theirs['median_seconds']) : '-');
    }
}

if ($json !== null) {
    $cpu = trim((string) @shell_exec('sysctl -n machdep.cpu.brand_string 2>/dev/null'))
        ?: trim((string) @shell_exec("grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2"));
    $commit = trim((string) @shell_exec('git -C '.escapeshellarg(dirname(__DIR__)).' rev-parse --short HEAD 2>/dev/null'));
    $reference = parse_ini_file(dirname(__DIR__).'/tests/fixtures/REFERENCE') ?: [];
    $document = [
        'recorded_at' => gmdate('Y-m-d\TH:i:s\Z'),
        'cpu' => $cpu ?: null,
        'platform' => PHP_OS_FAMILY.' '.php_uname('r').' '.php_uname('m'),
        'php_version' => PHP_VERSION,
        'opcache' => [
            'enabled' => (bool) ($opcache['opcache_enabled'] ?? false),
            'jit' => $jit,
            'jit_mode' => ini_get('opcache.jit') ?: null,
            'jit_buffer_size' => ini_get('opcache.jit_buffer_size') ?: null,
        ],
        'polyfill_commit' => $commit ?: null,
        'extension_reference_commit' => $reference['commit'] ?? null,
        'extension_loaded' => $native,
        'design' => "$p standard normal predictors plus an intercept; seed $seed",
        'repeat' => $repeat,
        'warmup_calls' => 1,
        'cases' => $cases,
    ];
    file_put_contents($json, json_encode($document, JSON_PRETTY_PRINT | JSON_UNESCAPED_SLASHES | JSON_PRESERVE_ZERO_FRACTION)."\n");
    echo "\nwrote $json\n";
}
