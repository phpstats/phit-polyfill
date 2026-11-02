<?php

declare(strict_types=1);

namespace Phit\Polyfill\Objects;

/**
 * Reading, iteration, and refused writes over a held array, for the fit and
 * table objects: the extension's fit.rs.
 *
 * @internal
 */
trait HeldArray
{
    /** @var array<string|int, mixed> */
    private array $data = [];

    private int $position = 0;

    public function __construct()
    {
        throw new \Exception('You cannot instantiate this class from PHP.');
    }

    /** @return array<string|int, mixed> */
    public function toArray(): array
    {
        return $this->data;
    }

    /** @return array<string|int, mixed> */
    public function jsonSerialize(): array
    {
        return $this->data;
    }

    /** @return array<string|int, mixed> */
    public function __debugInfo(): array
    {
        return $this->data;
    }

    public function offsetExists(mixed $offset): bool
    {
        $key = self::keyText($offset);

        return $key !== null && isset($this->data[$key]);
    }

    /** `$x[$key]` on the held array, with PHP's warning for a missing key. */
    public function offsetGet(mixed $offset): mixed
    {
        $key = self::keyText($offset);
        if ($key !== null && array_key_exists($key, $this->data)) {
            return $this->data[$key];
        }
        // Userland code cannot raise E_WARNING; the prefix is the one PHP gives the extension's warning.
        trigger_error(static::class.'::offsetGet(): Undefined array key "'.($key ?? '').'"', E_USER_WARNING);

        return null;
    }

    public function offsetSet(mixed $offset, mixed $value): void
    {
        throw new \LogicException(self::readOnly());
    }

    public function offsetUnset(mixed $offset): void
    {
        throw new \LogicException(self::readOnly());
    }

    public function current(): mixed
    {
        return array_values($this->data)[$this->position] ?? null;
    }

    public function key(): ?string
    {
        $key = array_keys($this->data)[$this->position] ?? null;

        return $key === null ? null : (string) $key;
    }

    public function next(): void
    {
        $this->position++;
    }

    public function rewind(): void
    {
        $this->position = 0;
    }

    public function valid(): bool
    {
        return $this->position < count($this->data);
    }

    /** A string key, or an int key as its decimal string; other offsets find nothing. */
    private static function keyText(mixed $offset): ?string
    {
        return match (true) {
            is_string($offset) => $offset,
            is_int($offset) => (string) $offset,
            default => null,
        };
    }

    abstract private static function readOnly(): string;
}
