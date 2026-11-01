<?php

// One model, fitted by the phit extension and by the polyfill in the same
// process, compared bit for bit.

use Phit\Polyfill\Polyfill;

return [
    'inputs' => [
        'model' => ['type' => 'radio', 'label' => 'Model', 'options' => [
            'fit_ols_matrix' => 'Linear (fit_ols_matrix)',
            'fit_logistic_matrix' => 'Logistic (fit_logistic_matrix)',
            'fit_poisson_matrix' => 'Poisson (fit_poisson_matrix)',
        ], 'default' => 'fit_logistic_matrix'],
        'rows' => ['type' => 'select', 'label' => 'Rows', 'options' => [500 => '500', 5000 => '5,000', 20000 => '20,000'], 'default' => 5000],
        'seed' => ['type' => 'number', 'label' => 'Seed', 'min' => 1, 'max' => 1000, 'step' => 1, 'default' => 7],
    ],
    'run' => function (array $input, string $data): array {
        $model = (string) $input['model'];
        $n = (int) $input['rows'];
        mt_srand((int) $input['seed']);

        // Simulated study: hours of study and prior score predict the outcome.
        $gaussian = function (): float {
            $u = (mt_rand() + 1) / (mt_getrandmax() + 2.0);
            $v = mt_rand() / (mt_getrandmax() + 1.0);

            return sqrt(-2.0 * log($u)) * cos(2.0 * M_PI * $v);
        };
        $x = $y = [];
        for ($i = 0; $i < $n; $i++) {
            $hours = 2.0 + 1.5 * $gaussian();
            $prior = $gaussian();
            $eta = -0.4 + 0.35 * $hours + 0.6 * $prior;
            $x[] = [$hours, $prior];
            $y[] = match ($model) {
                'fit_ols_matrix' => 60.0 + 4.0 * $eta + 5.0 * $gaussian(),
                'fit_logistic_matrix' => mt_rand() / mt_getrandmax() < 1 / (1 + exp(-$eta)) ? 1 : 0,
                'fit_poisson_matrix' => (function () use ($eta): int {
                    // Knuth's method for a Poisson draw with mean exp(eta / 2).
                    $limit = exp(-exp($eta / 2));
                    $k = 0;
                    $p = 1.0;
                    do {
                        $k++;
                        $p *= mt_rand() / mt_getrandmax();
                    } while ($p > $limit);

                    return $k - 1;
                })(),
            };
        }

        $time = function (callable $fit) use ($x, $y): array {
            $start = hrtime(true);
            $result = $fit($x, $y);

            return [$result, (hrtime(true) - $start) / 1e6];
        };
        $native = 'Phit\\'.$model;
        $loaded = extension_loaded('phit');
        [$polyfill, $polyfillMs] = $time(Polyfill::implementation($native));
        [$extension, $extensionMs] = $loaded ? $time($native) : [null, null];

        $blocks = [
            ['type' => 'text', 'content' => 'Code written against phit calls `Phit\\'.$model.'()`. With the extension '
                .'loaded, that is the native function. Without it, the polyfill defines the same function. This page '
                .'runs both in one process: the extension through its function name and the polyfill through '
                .'`Polyfill::implementation()`, which returns the polyfill even when the extension is loaded.'],
            ['type' => 'code', 'chain' => <<<PHP
                use Phit\\Polyfill\\Polyfill;
                use function Phit\\{$model};

                \$native = {$model}(\$x, \$y);                                   // the extension
                \$polyfill = Polyfill::implementation('Phit\\{$model}')(\$x, \$y); // the polyfill

                serialize(\$native) === serialize(\$polyfill);                 // bit for bit
                PHP],
        ];
        if (! $loaded) {
            return [...$blocks, ['type' => 'note', 'content' => 'The phit extension is not loaded on this server, so only the polyfill ran.'],
                ['type' => 'model', 'title' => 'Polyfill fit', 'result' => $polyfill]];
        }

        // Compare every number in the two results, bit for bit.
        $leaves = function (array $result): array {
            $out = [];
            array_walk_recursive($result, function ($v) use (&$out) {
                $out[] = $v;
            });

            return $out;
        };
        $a = $leaves($polyfill);
        $b = $leaves($extension);
        $numbers = $differing = 0;
        foreach ($b as $i => $v) {
            if (is_float($v) || is_int($v)) {
                $numbers++;
                $same = $a[$i] === $v || (is_float($v) && is_nan($v) && is_nan($a[$i]));
                $differing += $same ? 0 : 1;
            }
        }
        $identical = serialize($polyfill) === serialize($extension);

        $rows = [];
        foreach ($extension['coefficients'] as $j => $coefficient) {
            $rows[] = [
                ['(Intercept)', 'hours', 'prior'][$j],
                $coefficient,
                bin2hex(pack('E', $coefficient)),
                bin2hex(pack('E', $polyfill['coefficients'][$j])),
            ];
        }

        return [
            ...$blocks,
            ['type' => 'metrics', 'items' => [
                ['label' => 'Numbers compared', 'value' => $numbers, 'format' => '%d'],
                ['label' => 'Numbers that differ', 'value' => $differing, 'format' => '%d'],
                ['label' => 'Extension (ms)', 'value' => $extensionMs, 'format' => '%.1f'],
                ['label' => 'Polyfill (ms)', 'value' => $polyfillMs, 'format' => '%.1f'],
            ]],
            ['type' => 'table', 'title' => 'Coefficients and their IEEE 754 bit patterns', 'columns' => ['term', 'estimate', 'extension bits', 'polyfill bits'],
                'rows' => $rows, 'digits' => 6],
            ['type' => 'text', 'content' => $identical
                ? 'The two results are identical: the same keys in the same order, and every number equal bit for bit. '
                    .sprintf('The polyfill took %.1f times as long as the extension on this fit.', $polyfillMs / max($extensionMs, 1e-9))
                : 'The two results differ. Please report this with the model, row count, and seed.'],
            ['type' => 'model', 'title' => 'Fit (extension)', 'result' => $extension],
        ];
    },
];
