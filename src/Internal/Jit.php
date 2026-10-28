<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\ErrorCode;

/**
 * Whether OPcache's JIT is compiling this process. The mixed-model optimizers
 * refuse to run under it: PHP 8.4.26's tracing JIT moved BOBYQA to a minimum
 * 18% away from the true one without an error (docs/notes/function-jit.md).
 *
 * @internal
 */
final class Jit
{
    public static function compiling(): bool
    {
        if (! function_exists('opcache_get_status')) {
            return false;
        }
        $status = @opcache_get_status(false);
        if (is_array($status)) {
            return ($status['jit']['enabled'] ?? false) && ($status['jit']['on'] ?? false);
        }
        // opcache.restrict_api hides the status from this script, so judge
        // from the settings. opcache.jit alone is not enough: it keeps its
        // value when OPcache itself is off and nothing compiles.
        $enabled = (bool) ini_get('opcache.enable') && (PHP_SAPI !== 'cli' || (bool) ini_get('opcache.enable_cli'));
        $mode = strtolower(trim((string) ini_get('opcache.jit')));

        return $enabled
            && ! in_array($mode, ['', '0', 'off', 'disable', 'false', 'no'], true)
            && (int) ini_get('opcache.jit_buffer_size') > 0;
    }

    /** Throw before a mixed-model optimizer runs under the JIT. */
    public static function refuseMixedFit(): void
    {
        if (self::compiling()) {
            throw ErrorCode::exception(
                ErrorCode::UNSUPPORTED,
                "mixed models do not run under OPcache's JIT, which can move their optimizer to a wrong minimum without an error; set opcache.jit=disable or install the phit extension",
            );
        }
    }
}
