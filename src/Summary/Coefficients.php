<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/** Coefficient estimates with their standard errors, test statistics, and p-values. */
final readonly class Coefficients
{
    /**
     * @param  list<string>  $names
     * @param  list<float>  $estimates
     * @param  list<float>  $stdErrors
     * @param  list<float>  $statistics
     * @param  list<float>  $pValues
     * @param  't'|'z'  $statistic
     */
    public function __construct(
        public array $names,
        public array $estimates,
        public array $stdErrors,
        public array $statistics,
        public array $pValues,
        public string $statistic,
    ) {}
}
