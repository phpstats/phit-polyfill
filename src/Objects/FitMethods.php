<?php

declare(strict_types=1);

namespace Phit\Polyfill\Objects;

use Phit\Polyfill\Internal\Jit;
use Phit\Polyfill\Predict;
use Phit\Polyfill\Summary;

/**
 * The methods of the four fit classes: the result array of a model function,
 * read-only, printing its summary.
 *
 * @internal
 */
trait FitMethods
{
    use HeldArray;

    public function __toString(): string
    {
        return $this->summary();
    }

    public function summary(): string
    {
        Jit::refuse();

        return Summary::summary($this->data);
    }

    public function tidy(bool $confInt = false, float $level = 0.95, string $effects = 'fixed'): mixed
    {
        Jit::refuse();

        return Wrap::table(Summary::tidy($this->data, $confInt, $level, $effects));
    }

    public function predict(
        mixed $newdata,
        string $type = 'link',
        bool $seFit = false,
        string $interval = 'none',
        float $level = 0.95,
        bool $random = true,
        bool $allowNewLevels = false,
    ): array {
        Jit::refuse();

        return Predict::predict($this->data, $newdata, $type, $seFit, $interval, $level, $random, $allowNewLevels);
    }

    /** A Jupyter display bundle: the summary as text/plain and the array as application/json. */
    public function toMimeBundle(): array
    {
        return ['text/plain' => $this->summary(), 'application/json' => $this->data];
    }

    public function count(): int
    {
        return count($this->data);
    }

    private static function readOnly(): string
    {
        return 'A Phit fit is read-only; call toArray() for an array copy';
    }
}
