<?php

// When the polyfill is fast enough: time a workload with the phit extension
// and with the polyfill, and compare it with a time budget per request.

use Phit\Polyfill\Polyfill;

return [
    'inputs' => [
        'workload' => ['type' => 'radio', 'label' => 'Workload', 'options' => [
            'pvalues' => 'Two-sided p-values for 10,000 t statistics',
            'quantiles' => '10,000 gamma quantiles',
            'logistic' => 'Logistic regression, 5,000 rows, 5 predictors',
            'glmer' => 'Random-intercept logistic model, 1,000 rows in 40 groups',
        ], 'default' => 'logistic'],
        'budget' => ['type' => 'number', 'label' => 'Time budget per request (ms)', 'min' => 10, 'max' => 2000, 'step' => 10, 'default' => 100],
    ],
    'run' => function (array $input, string $data): array {
        $workload = (string) $input['workload'];
        $budget = (float) $input['budget'];
        mt_srand(11);
        $gaussian = function (): float {
            $u = (mt_rand() + 1) / (mt_getrandmax() + 2.0);
            $v = mt_rand() / (mt_getrandmax() + 1.0);

            return sqrt(-2.0 * log($u)) * cos(2.0 * M_PI * $v);
        };
        $uniform = fn (): float => (mt_rand() + 1) / (mt_getrandmax() + 2.0);

        // [function, arguments, code shown on the page]
        [$fn, $arguments, $code] = match ($workload) {
            'pvalues' => ['Phit\dist_t_sf_array', [array_map(fn () => abs(2.5 * $gaussian()), range(1, 10_000)), 24.0],
                "use function Phit\\dist_t_sf_array;\n\n// Two-sided p-values: 2 · P(T > |t|) with 24 degrees of freedom.\n\$p = array_map(fn (\$v) => 2 * \$v, dist_t_sf_array(array_map(abs(...), \$t), df: 24));"],
            'quantiles' => ['Phit\dist_gamma_quantile_array', [array_map(fn () => $uniform(), range(1, 10_000)), 2.5, 1.5],
                "use function Phit\\dist_gamma_quantile_array;\n\n\$q = dist_gamma_quantile_array(\$probabilities, shape: 2.5, scale: 1.5);"],
            'logistic' => (function () use ($gaussian) {
                $x = $y = [];
                for ($i = 0; $i < 5_000; $i++) {
                    $row = array_map(fn () => $gaussian(), range(1, 5));
                    $eta = 0.2 + 0.5 * $row[0] - 0.4 * $row[1] + 0.3 * $row[2];
                    $x[] = $row;
                    $y[] = mt_rand() / mt_getrandmax() < 1 / (1 + exp(-$eta)) ? 1 : 0;
                }

                return ['Phit\fit_logistic', [$x, $y], "use function Phit\\fit_logistic;\n\n\$fit = fit_logistic(\$x, \$y);"];
            })(),
            'glmer' => (function () use ($gaussian) {
                $effects = array_map(fn () => 0.8 * $gaussian(), range(1, 40));
                $frame = ['passed' => [], 'hours' => [], 'school' => []];
                for ($i = 0; $i < 1_000; $i++) {
                    $hours = 2.0 + $gaussian();
                    $eta = -1.0 + 0.6 * $hours + $effects[$i % 40];
                    $frame['passed'][] = mt_rand() / mt_getrandmax() < 1 / (1 + exp(-$eta)) ? 1 : 0;
                    $frame['hours'][] = $hours;
                    $frame['school'][] = 's'.($i % 40);
                }

                return ['Phit\fit_glmer', ['passed ~ hours + (1 | school)', $frame],
                    "use function Phit\\fit_glmer;\n\n\$fit = fit_glmer('passed ~ hours + (1 | school)', \$frame);"];
            })(),
        };

        $time = function (callable $f) use ($arguments): float {
            $f(...$arguments);
            $times = [];
            for ($i = 0; $i < 3; $i++) {
                $start = hrtime(true);
                $f(...$arguments);
                $times[] = (hrtime(true) - $start) / 1e6;
            }
            sort($times);

            return $times[1];
        };
        $polyfillMs = $time(Polyfill::implementation($fn));
        $extensionMs = extension_loaded('phit') ? $time($fn) : null;
        $fits = $polyfillMs <= $budget;

        $blocks = [
            ['type' => 'text', 'content' => 'The phit extension is the faster choice wherever it can be installed. The '
                .'polyfill is for hosts where it cannot: shared hosting, managed platforms that allow no custom '
                .'extensions, or a development machine without a Rust toolchain. Because the polyfill returns the same '
                .'results bit for bit, the choice is only about speed. This page times one workload with each, '
                .'taking the median of three runs after a warm-up call, on the server that renders it.'],
            ['type' => 'code', 'chain' => $code],
            ['type' => 'metrics', 'items' => array_values(array_filter([
                ['label' => 'Polyfill (ms)', 'value' => $polyfillMs, 'format' => '%.1f'],
                $extensionMs !== null ? ['label' => 'Extension (ms)', 'value' => $extensionMs, 'format' => '%.1f'] : null,
                $extensionMs !== null ? ['label' => 'Polyfill / extension', 'value' => $polyfillMs / max($extensionMs, 1e-9), 'format' => '%.1f×'] : null,
                ['label' => 'Share of the budget used by the polyfill', 'value' => 100 * $polyfillMs / $budget, 'format' => '%.0f%%'],
            ]))],
            ['type' => 'text', 'content' => $fits
                ? sprintf('The polyfill completes this workload in %.1f ms, within the %.0f ms budget. On a host without the extension, the polyfill is a workable choice for it.', $polyfillMs, $budget)
                : sprintf('The polyfill needs %.1f ms, more than the %.0f ms budget. Install the extension for this workload, or run it outside the request, for example in a queued job.', $polyfillMs, $budget)],
            ['type' => 'text', 'content' => 'A library can require the polyfill and suggest the extension, so that it '
                .'works everywhere and runs faster where phit is installed. When the extension is loaded, the polyfill '
                .'defines nothing and adds no overhead to the calls.'],
            ['type' => 'code', 'chain' => "{\n    \"require\": {\n        \"phpstats/phit-polyfill\": \"^0.1\"\n    },\n    \"suggest\": {\n        \"ext-phit\": \"Native implementation, about 10 to 20 times faster\"\n    }\n}"],
            ['type' => 'note', 'content' => 'On PHP 8.3, run the polyfill with OPcache and without the JIT, which is '
                .'PHP 8.3\'s default. On PHP 8.4 and 8.5, avoid the function JIT with hot-counter triggers '
                .'(`opcache.jit=1235`). See "PHP versions and the JIT" in the docs.'],
        ];
        if ($extensionMs === null) {
            array_splice($blocks, 3, 0, [['type' => 'note', 'content' => 'The phit extension is not loaded on this server, so only the polyfill was timed.']]);
        }

        return $blocks;
    },
];
