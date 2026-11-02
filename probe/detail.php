<?php
// Per-row components of Saddle::lnGammaLowerSeries. With "write", records them; with "compare", reports differences.
declare(strict_types=1);
require __DIR__.'/../vendor/autoload.php';
use Phit\Polyfill\Distribution\Saddle;
use Phit\Polyfill\Internal\Statrs;
$rows = json_decode(file_get_contents($argv[2]), true)['ln_gamma_lower'];
$private = fn (string $m) => (new ReflectionMethod(Saddle::class, $m))->getClosure();
$sumOf = $private('lowerSeriesSum');
$prefactor = $private('lnGammaPrefactor');
$out = [];
for ($rep = 0; $rep < 3; $rep++) {
    foreach ($rows as $i => $r) {
        [$a, $x] = [(float) $r[0], (float) $r[1]];
        $out[$i] = [
            'a' => $a, 'x' => $x,
            'real' => Saddle::lnGammaLowerSeries($a, $x),
            'sum' => $sumOf($a, $x),
            'prefactor' => $prefactor($a, $x),
            'poisson' => Saddle::lnPoissonRaw($a, $x),
            'lngamma' => Statrs::lnGamma($a + 1.0),
            'logsum' => log($sumOf($a, $x)),
        ];
    }
}
$hex = fn (array $row) => array_map(fn ($v) => bin2hex(pack('E', $v)), $row);
if ($argv[1] === 'write') {
    file_put_contents('/tmp/detail.json', json_encode(array_map($hex, $out)));
    exit;
}
$ref = json_decode(file_get_contents('/tmp/detail.json'), true);
$counts = [];
$shown = 0;
foreach ($out as $i => $row) {
    $h = $hex($row);
    $bad = array_keys(array_filter($h, fn ($v, $k) => $v !== $ref[$i][$k], ARRAY_FILTER_USE_BOTH));
    foreach ($bad as $k) { $counts[$k] = ($counts[$k] ?? 0) + 1; }
    if ($bad !== [] && $shown++ < 6) {
        printf("  row %d a=%.17g x=%.17g differs in %s\n", $i, $row['a'], $row['x'], implode(', ', $bad));
        foreach ($bad as $k) { printf("    %-9s jit %.17g  off %.17g\n", $k, $row[$k], unpack('E', hex2bin($ref[$i][$k]))[1]); }
    }
}
echo '  rows differing per component: ', json_encode($counts), "\n";
