<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\ErrorCode;

/**
 * Whether OPcache's JIT is compiling this process. Every Phit\ function the
 * polyfill defines refuses to run under it, on every PHP version: the JIT has
 * computed wrong results in the polyfill without an error on PHP 8.3, 8.4 and
 * 8.5 (docs/notes/function-jit.md, docs/notes/php83-jit.md).
 *
 * @internal
 */
final class Jit
{
    /** The OPcache settings refuse() last judged, and its verdict. */
    private static ?string $settings = null;

    private static bool $refusing = false;

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

    /**
     * Throw before a Phit\ function runs under the JIT. The verdict is cached
     * until one of the settings that can change at run time changes.
     */
    public static function refuse(): void
    {
        $settings = ini_get('opcache.enable').'|'.ini_get('opcache.enable_cli').'|'.ini_get('opcache.jit');
        if ($settings !== self::$settings) {
            self::$refusing = self::compiling();
            self::$settings = $settings;
        }
        if (self::$refusing) {
            throw ErrorCode::exception(
                ErrorCode::UNSUPPORTED,
                "the phit polyfill does not run under OPcache's JIT, which has computed wrong results in it without an error; set opcache.jit=disable or install the phit extension",
            );
        }
    }
}
