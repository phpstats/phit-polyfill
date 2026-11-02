<?php

declare(strict_types=1);

namespace Phit;

/** The table of named columns that tidy() returns, printed as R's print.data.frame(). */
final class Table implements \Stringable, \ArrayAccess, \Iterator, \Countable, \JsonSerializable
{
    use \Phit\Polyfill\Objects\TableMethods;
}
