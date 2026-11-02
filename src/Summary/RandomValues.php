<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/**
 * One random term's conditional modes and their conditional standard errors,
 * as lme4's ranef(fit, condVar = TRUE): one row per group, one column per term
 * column.
 */
final readonly class RandomValues
{
    /**
     * @param  list<string>  $levels
     * @param  list<list<float>>  $estimates
     * @param  list<list<float>>  $stdErrors
     */
    public function __construct(
        public array $levels,
        public array $estimates,
        public array $stdErrors,
    ) {}
}
