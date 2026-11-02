<?php

declare(strict_types=1);

namespace Phit\Polyfill\Objects;

use Phit\Polyfill\Summary\Format;

/**
 * The methods of Phit\Table: named columns, read-only, printed as R's
 * print.data.frame().
 *
 * @internal
 */
trait TableMethods
{
    use HeldArray;

    public function __toString(): string
    {
        return Format::dataFrame($this->data);
    }

    /** @return array<string|int, mixed> */
    public function toColumns(): array
    {
        return $this->data;
    }

    /** A Jupyter display bundle: the printed table as text/plain and the columns as application/json. */
    public function toMimeBundle(): array
    {
        return ['text/plain' => (string) $this, 'application/json' => $this->data];
    }

    /** The number of rows: the length of the first column, or 0. */
    public function count(): int
    {
        $first = reset($this->data);

        return is_array($first) ? count($first) : 0;
    }

    private static function readOnly(): string
    {
        return 'A Phit table is read-only; call toArray() for an array copy';
    }
}
