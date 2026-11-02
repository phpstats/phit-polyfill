<?php

declare(strict_types=1);

// Shared by OptimizeTest and OptimizeParityTest, so that either runs on its own.

/** The PHP twins of the R functions in phit's tests/generate-optimize-fixtures.R. */
function optimizeFunctions(): array
{
    return [
        'parabola' => fn (float $x) => ($x - 1 / 3) * ($x - 1 / 3),
        'sine' => fn (float $x) => sin($x),
        'kink' => fn (float $x) => abs($x - 0.7),
        'flat' => fn (float $x) => 1.0,
        'linear' => fn (float $x) => $x,
        'decay' => fn (float $x) => exp(-$x),
        'double_well' => fn (float $x) => ($x * $x - 1) * ($x * $x - 1) + 0.3 * $x,
        'wall' => fn (float $x) => $x < 0.45 ? INF : ($x - 0.5) * ($x - 0.5),
        'far' => fn (float $x) => ($x - 1e6) * ($x - 1e6),
        'hump' => fn (float $x) => $x * exp(-$x),
        'sqrt2' => fn (float $x) => $x * $x - 2,
        'fixed_point' => fn (float $x) => cos($x) - $x,
        'exp5' => fn (float $x) => exp($x) - 5,
        'at_one' => fn (float $x) => $x - 1,
        'at_three' => fn (float $x) => $x - 3,
        'step' => fn (float $x) => $x < 0.3 ? -1.0 : 1.0,
        'cubic' => fn (float $x) => $x * $x * $x - 2 * $x - 5,
        'log' => fn (float $x) => log($x),
        'cliff' => fn (float $x) => $x < 0.5 ? -INF : $x - 0.6,
    ];
}
