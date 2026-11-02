<?php
// Each variant of the series sum in its own required file; hashes of the exact bits over the table.
declare(strict_types=1);
require __DIR__.'/../vendor/autoload.php';
$names = ['orig', 'fresh', 'hoisted', 'both', 'whileLoop', 'saddleCopy'];
foreach ($names as $n) { require __DIR__."/F_$n.php"; }
$rows = json_decode(file_get_contents($argv[1]), true)['ln_gamma_lower'];
$sumOf = (new ReflectionMethod(\Phit\Polyfill\Distribution\Saddle::class, 'lowerSeriesSum'))->getClosure();
for ($rep = 0; $rep < 3; $rep++) {
    $lines = [];
    foreach ([...$names, 'real'] as $n) {
        $bits = '';
        foreach ($rows as $r) {
            [$a, $x] = [(float) $r[0], (float) $r[1]];
            $bits .= pack('E', $n === 'real' ? $sumOf($a, $x) : ("Probe\\F_$n::$n")($a, $x));
        }
        $lines[] = sprintf("  %-10s %s", $n, substr(md5($bits), 0, 12));
    }
}
echo implode("\n", $lines), "\n";
