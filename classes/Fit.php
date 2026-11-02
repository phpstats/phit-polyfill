<?php

declare(strict_types=1);

namespace Phit;

/**
 * A model fit: every model function's result implements it. Loaded by
 * bootstrap.php only when the phit extension is not.
 */
interface Fit extends \ArrayAccess, \Iterator, \Countable, \JsonSerializable, \Stringable {}
