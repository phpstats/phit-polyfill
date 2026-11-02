<?php
// Variants of Saddle::lnGammaLowerSeries's loop; prints, per variant, a hash of the exact bits of the
// series sums over the ln_gamma_lower table, and how many rows differ from the original computed
// in the same process. Compare the hashes across JIT settings.
declare(strict_types=1);
require __DIR__.'/../vendor/autoload.php';
final class V
{
    private const float EPS = PHP_FLOAT_EPSILON;
    private const int MAX_ITER = 1_000_000;
    public static function orig(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a; $converged = false;
        for ($i = 0; $i < self::MAX_ITER; $i++) {
            $denom += 1.0;
            $term *= $x / $denom;
            $y = $term - $comp;
            $t = $sum + $y;
            $comp = ($t - $sum) - $y;
            $sum = $t;
            if ($term <= $sum * self::EPS * 0.5) { $converged = true; break; }
        }
        return $converged ? $sum : NAN;
    }
    // B: fresh temporaries for each step of the compensation.
    public static function fresh(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a; $converged = false;
        for ($i = 0; $i < self::MAX_ITER; $i++) {
            $denom += 1.0;
            $ratio = $x / $denom;
            $term *= $ratio;
            $y = $term - $comp;
            $t = $sum + $y;
            $grown = $t - $sum;
            $comp = $grown - $y;
            $sum = $t;
            if ($term <= $sum * self::EPS * 0.5) { $converged = true; break; }
        }
        return $converged ? $sum : NAN;
    }
    // C: constants hoisted into locals, return from inside the loop.
    public static function hoisted(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a;
        $halfEps = self::EPS * 0.5; $cap = self::MAX_ITER;
        for ($i = 0; $i < $cap; $i++) {
            $denom += 1.0;
            $term *= $x / $denom;
            $y = $term - $comp;
            $t = $sum + $y;
            $comp = ($t - $sum) - $y;
            $sum = $t;
            if ($term <= $sum * $halfEps) { return $sum; }
        }
        return NAN;
    }
    // D: B and C together.
    public static function both(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a;
        $halfEps = self::EPS * 0.5; $cap = self::MAX_ITER;
        for ($i = 0; $i < $cap; $i++) {
            $denom += 1.0;
            $ratio = $x / $denom;
            $term *= $ratio;
            $y = $term - $comp;
            $t = $sum + $y;
            $grown = $t - $sum;
            $comp = $grown - $y;
            $sum = $t;
            if ($term <= $sum * $halfEps) { return $sum; }
        }
        return NAN;
    }
    // E: a while loop over the state, with no counter variable reused as a float.
    public static function whileLoop(float $a, float $x): float
    {
        $term = 1.0; $sum = 1.0; $comp = 0.0; $denom = $a; $left = self::MAX_ITER;
        while ($left-- > 0) {
            $denom += 1.0;
            $term *= $x / $denom;
            $y = $term - $comp;
            $t = $sum + $y;
            $comp = ($t - $sum) - $y;
            $sum = $t;
            if ($term <= $sum * self::EPS * 0.5) { return $sum; }
        }
        return NAN;
    }
}
$rows = json_decode(file_get_contents($argv[1]), true)['ln_gamma_lower'];
$cell = fn ($c) => is_string($c) ? (float) $c : (float) $c;
$jit = function_exists('opcache_get_status') ? ((opcache_get_status(false)['jit']['on'] ?? false) ? ini_get('opcache.jit') : 'off') : 'off';
echo "jit: $jit, PHP ", PHP_VERSION, ' ', php_uname('m'), "\n";
for ($rep = 0; $rep < 3; $rep++) { // repeat so hot counters trigger
    $ref = [];
    foreach (['orig', 'fresh', 'hoisted', 'both', 'whileLoop'] as $v) {
        $bits = '';
        $diff = 0;
        foreach ($rows as $i => $r) {
            $s = V::$v($cell($r[0]), $cell($r[1]));
            $b = pack('E', $s);
            $bits .= $b;
            if ($v === 'orig') { $ref[$i] = $b; } elseif ($b !== $ref[$i]) { $diff++; }
        }
        if ($rep === 2) printf("  %-10s %s  differs from orig in %d rows\n", $v, substr(md5($bits), 0, 12), $diff);
    }
    $bits = '';
    foreach ($rows as $r) { $bits .= pack('E', \Phit\Polyfill\Distribution\Saddle::lnGammaLowerSeries($cell($r[0]), $cell($r[1]))); }
    if ($rep === 2) printf("  %-10s %s  (Saddle::lnGammaLowerSeries, ln of the full value)\n", 'real', substr(md5($bits), 0, 12));
}
