<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Jit;

// Every Phit\ function the polyfill defines refuses to run while a JIT
// compiles the process (docs/notes/function-jit.md). The tests that run a
// child PHP under the tracing JIT check that in every job; the in-process
// test checks it in the CI jobs that enable a JIT.

/**
 * Run $code in a child PHP under the tracing JIT, after the autoloader, and
 * return what it printed, or null when no JIT compiles the child or the
 * extension defines the functions there.
 */
function underJit(string $code): ?string
{
    $script = tempnam(sys_get_temp_dir(), 'phit-jit');
    file_put_contents($script, '<?php require '.var_export(dirname(__DIR__, 2).'/vendor/autoload.php', true).';'
        .' if (! \Phit\Polyfill\Internal\Jit::compiling() || extension_loaded("phit")) { echo "skip"; exit; }'
        .' function code(callable $f): string { try { $f(); return "ok"; } catch (\Throwable $e) { return (string) $e->getCode(); } }'
        .$code);
    $command = escapeshellarg(PHP_BINARY).' -d opcache.enable_cli=1 -d opcache.jit=tracing -d opcache.jit_buffer_size=64M '.escapeshellarg($script);
    $output = (string) shell_exec($command.' 2>&1');
    unlink($script);

    return $output === 'skip' ? null : $output;
}

it('refuses every kind of Phit\ function under the JIT', function () {
    $output = underJit('
        echo implode(" ", [
            code(fn () => Phit\stat_mean([1, 2, 3])),
            code(fn () => Phit\dist_beta_cdf(0.5, 1e5, 1e5)),
            code(fn () => Phit\dist_normal_rng(3, 42)),
            code(fn () => Phit\fit_ols([[1], [2], [3], [4]], [1, 3, 2, 4])),
            code(fn () => Phit\fit_lmer("y ~ x + (1 | g)", ["y" => [1, 2, 3, 4], "x" => [1, 2, 3, 4], "g" => [1, 1, 2, 2]])),
            code(fn () => Phit\fit_mixed_intercept([[1], [2], [3], [4]], [1, 3, 2, 4], [1, 1, 2, 2])),
            code(fn () => Phit\stat_mean([])),
        ]);');
    if ($output === null) {
        test()->markTestSkipped('OPcache cannot start a JIT in a child PHP here');
    }

    expect($output)->toBe('1007 1007 1007 1007 1007 1007 1007');
});

it('runs once ini_set switches the JIT off, and refuses again when it comes back', function () {
    $output = underJit('
        $off = ini_set("opcache.jit", "off");
        $a = code(fn () => Phit\stat_mean([1, 2, 3]));
        ini_set("opcache.jit", $off);
        echo $a, " ", code(fn () => Phit\stat_mean([1, 2, 3]));');
    if ($output === null) {
        test()->markTestSkipped('OPcache cannot start a JIT in a child PHP here');
    }

    expect($output)->toBe('ok 1007');
});

it('refuses in this process when a JIT compiles it', function () {
    if (! Jit::compiling() || extension_loaded('phit')) {
        test()->markTestSkipped('no JIT is compiling this process, or the extension defines the functions');
    }

    expect(thrownCode(fn () => Phit\stat_mean([1.0, 2.0])))->toBe(1007);
});

it('guards every function in the generated bootstrap', function () {
    $bootstrap = (string) file_get_contents(dirname(__DIR__, 2).'/bootstrap.php');

    expect(substr_count($bootstrap, "\\Phit\\Polyfill\\Internal\\Jit::refuse();\n"))
        ->toBe(preg_match_all('/^    function /m', $bootstrap))
        ->toBeGreaterThan(200);
});
