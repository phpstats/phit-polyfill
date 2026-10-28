<?php

// Polyfill vs native phit timings for the package site: one file per area,
// in the site's benchmark format, under site/benchmarks/.
//
// Run with the extension, OPcache, and no JIT, as the README recommends:
//   php85 -d extension=.ext/libphit.dylib -d opcache.enable_cli=1 \
//         -d opcache.jit=disable benchmarks/site.php
// Options: --repeat=7 --seed=5 --out=site/benchmarks
//
// Each case is called once to warm up and then timed --repeat times. "value"
// is the median wall time in milliseconds and "reps" holds every timing.

declare(strict_types=1);

require __DIR__.'/../vendor/autoload.php';

use Phit\Polyfill\Polyfill;

$options = getopt('', ['repeat::', 'seed::', 'out::']);
$repeat = (int) ($options['repeat'] ?? 7);
$seed = (int) ($options['seed'] ?? 5);
$out = $options['out'] ?? dirname(__DIR__).'/site/benchmarks';
if (! extension_loaded('phit')) {
    fwrite(STDERR, "load the phit extension to compare with it\n");
    exit(1);
}

$gaussian = function (): float {
    $u = (mt_rand() + 1) / (mt_getrandmax() + 2.0);
    $v = mt_rand() / (mt_getrandmax() + 1.0);

    return sqrt(-2.0 * log($u)) * cos(2.0 * M_PI * $v);
};

/** @return list<float> milliseconds per call, after one warm-up call */
$time = function (callable $f) use ($repeat): array {
    $f();
    $times = [];
    for ($i = 0; $i < $repeat; $i++) {
        $start = hrtime(true);
        $f();
        $times[] = round((hrtime(true) - $start) / 1e6, 4);
    }

    return $times;
};

$median = function (array $values): float {
    sort($values);
    $middle = intdiv(count($values), 2);

    return count($values) % 2 ? $values[$middle] : ($values[$middle - 1] + $values[$middle]) / 2;
};

$rows = function (string $case, string $fn, array $arguments) use ($time, $median): array {
    $polyfill = Polyfill::implementation($fn);
    $quiet = fn (callable $f) => function () use ($f, $arguments) {
        set_error_handler(fn () => true);
        try {
            return $f(...$arguments);
        } finally {
            restore_error_handler();
        }
    };
    $result = [];
    foreach (['polyfill' => $polyfill, 'phit extension' => $fn] as $implementation => $f) {
        $reps = $time($quiet($f));
        $result[] = ['case' => $case, 'implementation' => $implementation, 'value' => $median($reps), 'reps' => $reps];
        fprintf(STDERR, "%-48s %-15s %10.2f ms\n", $case, $implementation, $median($reps));
    }

    return $result;
};

$cpu = trim((string) @shell_exec('sysctl -n machdep.cpu.brand_string 2>/dev/null'))
    ?: trim((string) @shell_exec("grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2"));
$opcache = function_exists('opcache_get_status') ? (opcache_get_status(false) ?: []) : [];
$jit = ($opcache['jit']['on'] ?? false) ? 'JIT '.ini_get('opcache.jit') : 'no JIT';
$commit = trim((string) @shell_exec('git -C '.escapeshellarg(dirname(__DIR__)).' rev-parse --short HEAD 2>/dev/null'));
$reference = parse_ini_file(dirname(__DIR__).'/tests/fixtures/REFERENCE') ?: [];
$environment = [
    'cpu' => $cpu ?: php_uname('m'),
    'php' => PHP_VERSION.', OPcache with '.$jit,
    'date' => gmdate('Y-m-d'),
    'commit' => $commit.' (phit '.substr($reference['commit'] ?? '', 0, 7).')',
];
$write = function (string $name, array $document) use ($out, $environment): void {
    $document = [
        'title' => $document['title'],
        'description' => $document['description'],
        'measure' => 'median wall time per call',
        'unit' => 'ms',
        'environment' => $environment,
        'rows' => $document['rows'],
    ];
    @mkdir($out, 0777, true);
    file_put_contents("$out/$name.json", json_encode($document, JSON_PRETTY_PRINT | JSON_UNESCAPED_SLASHES | JSON_UNESCAPED_UNICODE | JSON_PRESERVE_ZERO_FRACTION)."\n");
};

// Models: 10 standard normal predictors plus an intercept.
mt_srand($seed);
$n = 10_000;
$x = $continuous = $binary = $counts = [];
$data = ['y' => [], 'b' => []];
for ($i = 0; $i < $n; $i++) {
    $row = [];
    $eta = 0.3;
    for ($j = 0; $j < 10; $j++) {
        $row[] = $value = $gaussian();
        $data['x'.($j + 1)][] = $value;
        $eta += $value * (($j % 3) - 1) * 0.2;
    }
    $x[] = $row;
    $continuous[] = $data['y'][] = $eta + $gaussian();
    $binary[] = $data['b'][] = (mt_rand() / mt_getrandmax()) < 1 / (1 + exp(-$eta)) ? 1 : 0;
    $counts[] = (int) floor(-log((mt_rand() + 1) / (mt_getrandmax() + 2.0)) * exp(0.5 + 0.3 * $eta));
}
$rhs = implode(' + ', array_map(fn ($j) => "x$j", range(1, 10)));

// Mixed models: 2,000 rows in 50 groups with two predictors.
$m = 2_000;
$mixed = ['y' => [], 'b' => [], 'x1' => [], 'x2' => [], 'g' => []];
$effects = array_map(fn () => 0.8 * $gaussian(), range(1, 50));
$mx = $labels = $my = $mb = [];
for ($i = 0; $i < $m; $i++) {
    $g = $i % 50;
    $x1 = $gaussian();
    $x2 = $gaussian();
    $eta = 0.2 + 0.5 * $x1 - 0.3 * $x2 + $effects[$g];
    $mx[] = [$x1, $x2];
    $labels[] = 'g'.$g;
    $my[] = $eta + $gaussian();
    $mb[] = (mt_rand() / mt_getrandmax()) < 1 / (1 + exp(-$eta)) ? 1 : 0;
    $mixed['x1'][] = $x1;
    $mixed['x2'][] = $x2;
    $mixed['g'][] = 'g'.$g;
    $mixed['y'][] = end($my);
    $mixed['b'][] = end($mb);
}

$write('models', [
    'title' => 'Models: polyfill and native phit',
    'description' => 'Time per fit for each model function, with the polyfill and the phit extension in the same process. '
        .'Regression designs have 10,000 rows and 10 standard normal predictors plus an intercept. '
        ."Mixed models have 2,000 rows in 50 groups with two predictors. Each case runs once to warm up and then $repeat times.",
    'rows' => array_merge(
        $rows('fit_ols, n = 10,000', 'Phit\fit_ols', [$x, $continuous]),
        $rows('fit_logistic, n = 10,000', 'Phit\fit_logistic', [$x, $binary]),
        $rows('fit_probit, n = 10,000', 'Phit\fit_probit', [$x, $binary]),
        $rows('fit_poisson, n = 10,000', 'Phit\fit_poisson', [$x, $counts]),
        $rows('fit_lm, n = 10,000', 'Phit\fit_lm', ["y ~ $rhs", $data]),
        $rows('fit_glm (binomial), n = 10,000', 'Phit\fit_glm', ["b ~ $rhs", $data, 'binomial']),
        $rows('fit_mixed_intercept, n = 2,000', 'Phit\fit_mixed_intercept', [$mx, $my, $labels]),
        $rows('fit_lmer, n = 2,000', 'Phit\fit_lmer', ['y ~ x1 + x2 + (1 | g)', $mixed]),
        $rows('fit_logistic_mixed, n = 2,000', 'Phit\fit_logistic_mixed', [$mx, $mb, $labels]),
        $rows('fit_glmer, n = 2,000', 'Phit\fit_glmer', ['b ~ x1 + x2 + (1 | g)', $mixed]),
    ),
]);

// Distributions: the cdf and quantile of 10,000 points per family.
$families = [
    'normal' => [[0.0, 1.0], fn () => 3 * $gaussian()],
    't' => [[5.0], fn () => 3 * $gaussian()],
    'chisquare' => [[4.0], fn () => 10 * mt_rand() / mt_getrandmax()],
    'f' => [[5.0, 12.0], fn () => 4 * mt_rand() / mt_getrandmax()],
    'exponential' => [[1.5], fn () => 3 * mt_rand() / mt_getrandmax()],
    'logistic' => [[0.0, 1.0], fn () => 3 * $gaussian()],
    'lognormal' => [[0.0, 1.0], fn () => 5 * mt_rand() / mt_getrandmax()],
    'gamma' => [[2.5, 1.5], fn () => 10 * mt_rand() / mt_getrandmax()],
    'beta' => [[2.0, 5.0], fn () => mt_rand() / mt_getrandmax()],
    'binomial' => [[40, 0.3], fn () => (float) mt_rand(0, 40)],
    'poisson' => [[6.0], fn () => (float) mt_rand(0, 20)],
    'nbinom' => [[3.0, 0.4], fn () => (float) mt_rand(0, 20)],
    'nbinom_mu' => [[3.0, 4.5], fn () => (float) mt_rand(0, 20)],
];
$distributionRows = [];
foreach ($families as $family => [$params, $draw]) {
    mt_srand($seed + crc32($family));
    $points = array_map(fn () => $draw(), range(1, 10_000));
    $probabilities = array_map(fn () => (mt_rand() + 1) / (mt_getrandmax() + 2.0), range(1, 10_000));
    $distributionRows = array_merge(
        $distributionRows,
        $rows("$family cdf, 10,000 points", "Phit\\dist_{$family}_cdf_array", [$points, ...$params]),
        $rows("$family quantile, 10,000 points", "Phit\\dist_{$family}_quantile_array", [$probabilities, ...$params]),
    );
}
$write('distributions', [
    'title' => 'Distributions: polyfill and native phit',
    'description' => 'Time per call of each family\'s cdf_array and quantile_array on 10,000 points, with the polyfill and the phit extension in the same process. '
        ."Each case runs once to warm up and then $repeat times.",
    'rows' => $distributionRows,
]);
