<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/** One random-effects term: its grouping factor, column names, covariance, and number of groups. */
final readonly class RandomTerm
{
    /**
     * @param  list<string>  $columns
     * @param  list<list<float>>  $covariance
     */
    public function __construct(
        public string $group,
        public array $columns,
        public array $covariance,
        public bool $correlated,
        public int $nGroups,
    ) {}
}
