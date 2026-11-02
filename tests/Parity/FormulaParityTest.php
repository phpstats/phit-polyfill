<?php

declare(strict_types=1);

// Polyfill vs native extension for Phit\model_matrix, Phit\fit_ols and Phit\fit_generalized_linear:
// the same keys in the same order, identical design matrices, names and row
// accounting, fit values within a scale-aware tolerance, the same warning text
// and the same exception codes. Runs only when phit is loaded.

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

/** Run $f, capturing the result or exception code plus warning messages. */
function formulaOutcome(callable $f, array $arguments): array
{
    $warnings = [];
    set_error_handler(function (int $severity, string $message) use (&$warnings): bool {
        $warnings[] = $message;

        return true;
    });
    try {
        return ['result' => nativeValue($f(...$arguments)), 'warnings' => $warnings];
    } catch (\Throwable $e) {
        return ['code' => $e->getCode(), 'class' => $e::class, 'message' => $e->getMessage(), 'warnings' => $warnings];
    } finally {
        restore_error_handler();
    }
}

/** @return list<mixed> */
function formulaFlatten(mixed $value): array
{
    if (! is_array($value)) {
        return [$value];
    }

    return array_merge([], ...array_map(formulaFlatten(...), array_values($value)));
}

/**
 * Largest scaled difference for a numeric fit key, judged against the natural
 * scale of each quantity as in ModelsParityTest: coefficients and intervals
 * against their standard errors, statistics against 1, p-values on the log
 * scale, fitted values and residuals against max |y|, everything else against
 * 1e-3 of the key's largest magnitude.
 */
function formulaFitDifference(array $polyfill, array $native, string $key): float
{
    $a = formulaFlatten($polyfill[$key]);
    $b = formulaFlatten($native[$key]);
    if (count($a) !== count($b)) {
        return INF;
    }
    $largest = 0.0;
    foreach ($b as $value) {
        if (is_float($value) && is_finite($value)) {
            $largest = max($largest, abs($value));
        }
    }
    $se = $native['std_errors'];
    $ymax = max(array_map(fn ($f, $r) => abs($f + $r), $native['fitted'], $native['residuals']));
    $max = 0.0;
    foreach ($b as $i => $value) {
        $x = $a[$i];
        if (! is_int($value) && ! is_float($value)) {
            $max = max($max, $x === $value ? 0.0 : INF);

            continue;
        }
        $x = (float) $x;
        $y = (float) $value;
        if (is_nan($x) || is_nan($y) || is_infinite($x) || is_infinite($y)) {
            $max = max($max, ($x === $y || (is_nan($x) && is_nan($y))) ? 0.0 : INF);

            continue;
        }
        if ($key === 'p_values' && $x > 0.0 && $y > 0.0) {
            $max = max($max, abs(log($x) - log($y)) / max(1.0, abs(log($y))));

            continue;
        }
        $floor = match ($key) {
            'coefficients', 'named' => $se[$i],
            'conf_int' => $se[intdiv($i, 2)],
            't_statistics', 'z_statistics', 'r_squared', 'adj_r_squared' => 1.0,
            'fitted', 'residuals' => $ymax,
            default => 1e-3 * $largest,
        };
        $max = max($max, abs($x - $y) / (abs($y) + $floor + PHP_FLOAT_MIN));
    }

    return $max;
}

function recordFormulaParity(string $fn, string $group, float $difference): void
{
    $GLOBALS['formula_parity'][$fn][$group] = max($GLOBALS['formula_parity'][$fn][$group] ?? 0.0, $difference);
}

/** Exact keys/order and design; fits within $tolerance; identical exception codes. */
function assertFormulaParity(string $fn, array $arguments, string $group, float $tolerance = 1e-9): void
{
    $label = "$fn(".json_encode($arguments[0] ?? null).") [$group]";
    $polyfill = formulaOutcome(\Phit\Polyfill\Polyfill::implementation($fn), $arguments);
    $native = formulaOutcome($fn, $arguments);
    parityLog('formula', $fn, $polyfill, $native, isset($native['code']) || isset($polyfill['code']));
    if (isset($native['code']) || isset($polyfill['code'])) {
        expect($polyfill['code'] ?? 'result')->toBe($native['code'] ?? 'result', "$label: native '".($native['message'] ?? '')."', polyfill '".($polyfill['message'] ?? '')."'");
        // Messages are not part of the contract, but the port keeps them identical.
        expect($polyfill['class'] ?? null)->toBe($native['class'] ?? null, $label)
            ->and($polyfill['message'] ?? null)->toBe($native['message'] ?? null, $label);
        $GLOBALS['formula_parity_errors'][$fn] = ($GLOBALS['formula_parity_errors'][$fn] ?? 0) + 1;

        return;
    }
    $a = $polyfill['result'];
    $b = $native['result'];
    expect(array_keys($a))->toBe(array_keys($b), $label)
        ->and($polyfill['warnings'])->toBe($native['warnings'], $label);
    if ($fn === 'Phit\model_matrix') {
        // Designs are built by the same arithmetic, so they agree exactly.
        expect($a)->toBe($b, $label);
        recordFormulaParity($fn, $group, 0.0);

        return;
    }
    if ($fn !== 'Phit\fit_ols' && ($b['link'] ?? 'identity') !== 'identity') {
        // Iterative fits stop once the full step is below 1e-8 relative.
        $tolerance = max($tolerance, 1e-7);
    }
    foreach (array_keys($b) as $key) {
        $va = $a[$key];
        $vb = $b[$key];
        if ($key === 'iterations') {
            expect(abs($va - $vb))->toBeLessThanOrEqual(1, "$label iterations");

            continue;
        }
        if ($key === 'diagnostics') {
            unset($va['iterations'], $vb['iterations']);
            expect($va)->toBe($vb, $label);

            continue;
        }
        if ($key === 'named') {
            expect(array_keys($va))->toBe(array_keys($vb), $label);
        }
        $numeric = array_filter(formulaFlatten($vb), fn ($v) => is_float($v));
        if ($numeric === [] || in_array($key, ['df_residual', 'n_input', 'n_observations', 'n_dropped', 'retained_rows', 'dropped_rows', 'assign'], true)) {
            expect($va)->toBe($vb, "$label $key");

            continue;
        }
        $difference = formulaFitDifference($a, $b, $key);
        recordFormulaParity($fn, $group, $difference);
        expect($difference)->toBeLessThanOrEqual($tolerance, "$label $key differs by $difference");
    }
}

/** A frame with every column kind, missing values, unused and reordered levels. */
function formulaFrame(): array
{
    mt_srand(2024);
    $n = 30;
    $u = fn () => mt_rand() / (mt_getrandmax() + 1.0);
    $data = [];
    $text = ['b', 'a', 'B', 'é', 'a', 'Zed', 'b', 'a b'];
    for ($i = 0; $i < $n; $i++) {
        $x = round(4 * $u() - 2, 3);
        $z = mt_rand(1, 9);
        $data['x'][] = $x;
        $data['z'][] = $z;
        $data['f'][] = $text[$i % count($text)];
        $data['gcode'][] = $i % 3 === 0 ? 2 : $i % 2;
        $data['ocode'][] = $i % 4;
        $data['l'][] = $i === 5 ? null : ($i % 3 === 1);
        $data['bin'][] = (int) ($u() < 1 / (1 + exp(-$x)));
        $data['cnt'][] = mt_rand(0, 6) + ($x > 0 ? 2 : 0);
        $data['expo'][] = round(0.5 + 3 * $u(), 2);
        $data['w'][] = round(0.5 + 2 * $u(), 2);
        $data['trials'][] = mt_rand(3, 8);
        $data['neg'][] = $i % 7 === 3 ? -1.0 : 1.0 + $u();
        $data['zero'][] = $i === 4 ? 0 : 2;
        $data['y'][] = $i === 9 ? null : ($i === 17 ? NAN : round(1.0 + 0.8 * $x - 0.1 * $z + $u(), 4));
        $data['k'][] = [10, 2, 2.5, 10, -1][$i % 5];
        $data['grp'][] = $i % 6;
    }
    $data['succ'] = array_map(fn ($t, $b) => min($t, $b * 2 + $t % 3), $data['trials'], $data['bin']);
    $data['g'] = ['levels' => ['hi', 'lo', 'mid', 'unused'], 'codes' => $data['gcode'], 'ordered' => false];
    $data['o'] = ['levels' => ['low', 'mid', 'high', 'top'], 'codes' => $data['ocode'], 'ordered' => true];
    $data['o']['codes'][3] = null;
    unset($data['gcode'], $data['ocode']);
    $data['name with space'] = array_map(fn ($v) => $v * 2.0, $data['x']);
    $data['é'] = array_map(fn ($v) => $v + 1.0, $data['x']);
    $data['meta'] = array_fill(0, $n, [1]);
    $data['mixed'] = array_merge([1, 'a', true], array_fill(0, $n - 3, 2));

    return $data;
}

afterAll(function () {
    $path = getenv('PARITY_REPORT');
    if ($path && isset($GLOBALS['formula_parity'])) {
        $existing = is_file($path) ? json_decode(file_get_contents($path), true) : [];
        file_put_contents($path, json_encode(array_replace_recursive($existing, $GLOBALS['formula_parity']), JSON_PRETTY_PRINT));
    }
});

it('agrees with the extension on the R reference cases', function () {
    $fixture = loadFixture('formula-r.json');
    foreach ($fixture['cases'] as $case) {
        $data = $fixture['datasets'][$case['dataset']];
        [$fn, $arguments] = match ($case['function']) {
            'model_matrix' => ['Phit\model_matrix', [$case['formula'], $data]],
            'lm' => ['Phit\fit_ols', [$case['formula'], $data, 'drop', $case['weights']]],
            'glm' => ['Phit\fit_generalized_linear', [$case['formula'], $data, $case['family'], $case['link'], 'drop', $case['weights']]],
        };
        assertFormulaParity($fn, $arguments, 'fixtures');
        if ($case['function'] !== 'model_matrix') {
            // The design and the fit also agree when the formula keeps unused levels.
            assertFormulaParity('Phit\model_matrix', [$case['formula'], $data], 'fixtures');
        }
    }
});

it('agrees on every supported syntax form and transform', function () {
    $data = formulaFrame();
    $formulas = [
        'y ~ x', 'y ~ x + z', 'y~x+z', "y ~ x\n\t+ z", 'y ~ x + x', 'y ~ x:x', 'y ~ +x', 'y ~ x - 1', 'y ~ x + 0',
        'y ~ 0 + x', 'y ~ -1 + x', 'y ~ x + 0 + 1', 'y ~ x - 0', 'y ~ 1', 'y ~ 0 + 1', 'y ~ -x', 'y ~ x - x',
        'y ~ x * z', 'y ~ x:z', 'y ~ z:x + x', 'y ~ (x + z + f)^2', 'y ~ (x + z + l)^3', 'y ~ (x + z)^1', 'y ~ (x + z)**2',
        'y ~ x * z - x:z', 'y ~ f / x', 'y ~ f / (x + z)', 'y ~ x %in% f', 'y ~ x + z %in% g', 'y ~ (x + z) %in% f',
        'y ~ f', 'y ~ g', 'y ~ o', 'y ~ l', 'y ~ 0 + f', 'y ~ 0 + g', 'y ~ 0 + o', 'y ~ f - 1', 'y ~ f:g', 'y ~ f * g',
        'y ~ g:f + f', 'y ~ f + g + 0', 'y ~ f * g - 1', 'y ~ x:f', 'y ~ x:f - 1', 'y ~ f:x - 1', 'y ~ o:x', 'y ~ l:g',
        'y ~ x * l * o', 'y ~ f:g:l', 'y ~ 0 + x:g + f', 'y ~ .', 'y ~ . - x', 'log(y) ~ x + .', '~ x + g', '~ 0 + g',
        'y ~ `name with space`', 'y ~ é + x', 'y ~ TRUE', 'y ~ T + x', 'y ~ log(x + 3)', 'y ~ log(expo, 2)',
        'y ~ log(expo, base = 10)', 'y ~ log(expo, base = -2)', 'y ~ log(expo, base = (2))', 'y ~ log2(expo) + log10(expo)',
        'y ~ log1p(expo) + expm1(x)', 'y ~ exp(x) + sqrt(expo) + abs(x)', 'y ~ sin(x) + cos(x) + tan(x)',
        'y ~ I(x^2) + I(x^-1) + I(2^x) + I(x**3)', 'y ~ I(x * 100000) + I(x / 3) + I(x * 0.1) + I(x - 19723)',
        'y ~ I(-x) + I(+x)', 'y ~ I((x + 1)^2)', 'y ~ I(x > 0)', 'y ~ I(x >= 0) + I(z <= 3)', 'y ~ I(x == z) + I(x != z)',
        'y ~ I(x < 0 & z > 3)', 'y ~ I(x < 0 && z > 3)', 'y ~ I(x < 0 | z > 3)', 'y ~ I(!l)', 'y ~ I(!x)', 'y ~ !l',
        'y ~ I(f == "a")', 'y ~ I("a" != f)', "y ~ I(g == 'lo')", 'y ~ I(f == "\\"q\\"")', 'y ~ I(l & x)',
        'y ~ I(1^x) + I(x^0)', 'y ~ I(TRUE)', 'y ~ poly(x, 3)', 'y ~ poly(x, degree = 2)', 'y ~ poly(x, 2, raw = TRUE)',
        'y ~ poly(x, raw = T, 3)', 'y ~ poly(x)', 'y ~ poly(z, 4) + x', 'y ~ poly(x, 2):f', 'y ~ poly(x = x, 2)',
        'y ~ factor(z)', 'y ~ as.factor(z) + x', 'y ~ ordered(z)', 'y ~ factor(k)', 'y ~ factor(l)', 'y ~ factor(g)',
        'y ~ factor(f)', 'y ~ ordered(f)', 'y ~ factor(o)', 'y ~ as.numeric(g) + as.numeric(x) + as.double(l)',
        'y ~ as.numeric(o)', 'y ~ x + offset(z)', 'y ~ offset(log(expo)) + x + offset(z)', 'y ~ offset(x) + offset(x)',
        'y ~ x + offset(l)', 'cbind(succ, trials - succ) ~ x', 'cbind(succ, fail = trials - succ) ~ x', 'f ~ x', 'g ~ x',
        'l ~ x', 'o ~ x', 'factor(z) ~ x', 'I(y > 1) ~ x', 'y ~ log(neg)', 'y ~ sqrt(neg) + x', 'y ~ x + 2L', 'y ~ 1L + x',
        'y ~ x + (1 | f)', 'y ~ x + (1 + x | g)', 'y ~ x + (0 + x || f)', 'y ~ x + (1 | f/g)', 'y ~ x + (1 | f:g)',
        'y ~ x + (x | g:f) + (1 | grp)', 'y ~ x + (1 | z)', 'y ~ x + (1 | k)', 'y ~ x + (1 | l)', 'y ~ x + (f | g)',
        'y ~ x + 1 | g', 'y ~ x - (1 | g)', 'y ~ (1 | g) + 0 + x', 'y ~ x + (1 | (f))', 'y ~ x + (1 | factor(z))',
        'y ~ x + (1 | x + z)', 'y ~ x + (offset(z) | g)', 'y ~ x + ((1 | f) | g)', 'y ~ x + (0 | g)', 'y ~ x + (1 | "a")',
        'y ~ x + (1 | cbind(x, z))',
    ];
    foreach ($formulas as $formula) {
        assertFormulaParity('Phit\model_matrix', [$formula, $data], 'syntax');
        assertFormulaParity('Phit\model_matrix', [$formula, $data, 'raise'], 'syntax');
        assertFormulaParity('Phit\fit_ols', [$formula, $data], 'syntax');
        assertFormulaParity('Phit\fit_generalized_linear', [$formula, $data], 'syntax');
    }
});

it('agrees on malformed formulas', function () {
    $data = formulaFrame();
    $formulas = [
        '', ' ', '~', 'y ~', 'y ~ x +', 'y ~ (x', 'y ~ x)', 'y ~ ((x + z)', 'y + x', 'y', 'y ~ x ~ z', '~ ~ x',
        'y ~ nope', 'y ~ x + nope', 'nope ~ x', 'y ~ foo(x)', 'y ~ ns(x, 3)', 'y ~ log()', 'y ~ log(x, 2, 3)',
        'y ~ log(x = expo, 2)', 'y ~ log(expo, base = z)', 'y ~ log2(x, 2)', 'y ~ log2(y = x)', 'y ~ I()', 'y ~ I(x, z)',
        'y ~ offset()', 'y ~ factor(f, levels = 2)', 'y ~ factor()', 'y ~ factor(cbind(x, z))', 'y ~ as.numeric(f)',
        'y ~ as.numeric(x, z)', 'y ~ poly(x, 0)', 'y ~ poly(x, 26)', 'y ~ poly(x, 2.5)', 'y ~ poly(x, z)',
        'y ~ poly(x, degree = 2, raw = 1)', 'y ~ poly(x, z, 2)', 'y ~ poly(x, coefs = 2)', 'y ~ poly(f, 2)',
        'y ~ poly(l, 1)', 'y ~ poly(y, 2)', 'y ~ poly(x, 40)', 'y ~ poly(grp, 6)', 'y ~ cbind()', 'cbind(cnt, x, z) ~ x',
        'cbind(x) ~ z', 'y ~ cbind(x, z)', 'y ~ I(f + 1)', 'y ~ I(x == "a")', 'y ~ I("a" == "b")', 'y ~ I(x || z)',
        'y ~ I(x %in% z)', 'y ~ I(x:z)', 'y ~ I(x ~ z)', 'y ~ I(!f)', 'y ~ I(-f)', 'y ~ x > 1', 'y ~ x & z',
        'y ~ "a"', 'y ~ x + "a"', 'y ~ 2', 'y ~ x + 3', 'y ~ x^2', 'y ~ (x + z)^0', 'y ~ (x + z)^1.5', 'y ~ (x + z)^z',
        'y ~ x | g', 'y ~ log(x | g)', 'y ~ x %o% z', 'y ~ x % z', 'y ~ x # comment', 'y ~ x $ z', 'y ~ x; z',
        'y ~ x [1]', 'y ~ "unterminated', 'y ~ `unterminated', 'y ~ ``', 'y ~ 1.2.3', 'y ~ x + .5x', 'y ~ log(zero)',
        'y ~ exp(x * 1000)', 'y ~ I(x / 0)', 'y ~ log(x - 10) + x', 'y ~ y + x', 'y ~ x + log(y)', 'log(y) ~ y',
        'y ~ x - y', 'y ~ 0', 'y ~ -1', 'y ~ x - x - 1', 'y ~ meta', 'y ~ mixed', 'y ~ x + meta', 'meta ~ x',
        'y ~ .', 'y ~ f + Zed', 'y ~ f,x', 'y ~ f(x)=1', 'y ~ log(x,)', 'y ~ log(,x)', 'y ~ x + ', '(y ~ x)',
        "y ~ x\xff", 'y ~ (1 | g) + (1 | f)', '~ (1 | g)', 'y ~ x:(1 | g)',
    ];
    foreach ($formulas as $formula) {
        assertFormulaParity('Phit\model_matrix', [$formula, $data], 'malformed');
        assertFormulaParity('Phit\fit_ols', [$formula, $data], 'malformed');
        assertFormulaParity('Phit\fit_generalized_linear', [$formula, $data, 'binomial'], 'malformed');
    }
});

it('agrees on level ordering, collation and unused levels', function () {
    $strings = ['b', 'a', 'B', 'é', 'a', 'Zed', 'b', '', '10', '9', 'e', 'É', "a\0b", 'ä', 'a'];
    $n = count($strings);
    $y = array_map(fn ($i) => sin($i) + $i / 10, range(0, $n - 1));
    $data = ['y' => $y, 'x' => range(1, $n), 's' => $strings,
        'd' => ['levels' => ['z', 'y', 'x', 'w'], 'codes' => array_map(fn ($i) => $i % 3, range(0, $n - 1)), 'ordered' => false],
        'od' => ['levels' => ['c', 'a', 'b'], 'codes' => array_map(fn ($i) => $i % 3, range(0, $n - 1)), 'ordered' => true],
        'num' => array_map(fn ($i) => [0.1 + 0.2, 0.3, -0.0, 0.0, 1e5, 1e-5, 123456.7, -2][$i % 8], range(0, $n - 1)),
        'bad' => array_merge(["\xff"], array_fill(0, $n - 1, 'ok')),
    ];
    foreach (['y ~ s', 'y ~ 0 + s', 'y ~ d', 'y ~ 0 + d', 'y ~ od', 'y ~ factor(num)', 'y ~ ordered(num)', 'y ~ factor(d)',
        'y ~ s:d', 'y ~ d * od', 'y ~ factor(s)', 'y ~ bad', 'y ~ x + (1 | s)', 'y ~ x + (1 | num)', 'y ~ x + (1 | d:s)',
        'y ~ I(s == "é")', 's ~ x', 'd ~ x'] as $formula) {
        assertFormulaParity('Phit\model_matrix', [$formula, $data], 'levels');
        assertFormulaParity('Phit\fit_ols', [$formula, $data], 'levels');
    }
    // Rows that remove a level: kept by model_matrix, dropped by the fits.
    $sparse = $data;
    $sparse['y'][0] = null;
    $sparse['y'][3] = null;
    foreach (['y ~ s + x', 'y ~ 0 + s', 'y ~ d + x', 'y ~ od + x', 'y ~ factor(s)'] as $formula) {
        assertFormulaParity('Phit\model_matrix', [$formula, $sparse], 'levels');
        assertFormulaParity('Phit\fit_ols', [$formula, $sparse], 'levels');
        assertFormulaParity('Phit\fit_generalized_linear', [$formula, $sparse], 'levels');
    }
});

it('agrees on missing values under drop and raise, and on weights', function () {
    $data = formulaFrame();
    foreach (['drop', 'raise', 'keep', 'DROP', ''] as $missing) {
        foreach (['y ~ x', 'y ~ x + l', 'y ~ o + x', 'z ~ x', 'y ~ log(neg)'] as $formula) {
            assertFormulaParity('Phit\model_matrix', [$formula, $data, $missing], 'missing');
            assertFormulaParity('Phit\fit_ols', [$formula, $data, $missing], 'missing');
            assertFormulaParity('Phit\fit_generalized_linear', [$formula, $data, 'gaussian', '', $missing], 'missing');
        }
    }
    $n = count($data['x']);
    $weights = [
        null, 'w', 'expo', 'z', 'x', 'f', 'l', 'g', 'meta', 'nope', '', array_fill(0, $n, 1.0), $data['w'],
        array_merge([null], array_slice($data['w'], 1)), array_merge([NAN], array_slice($data['w'], 1)),
        array_merge([0], array_slice($data['w'], 1)), array_merge([-1.0], array_slice($data['w'], 1)),
        array_merge([INF], array_slice($data['w'], 1)), array_merge(['1'], array_slice($data['w'], 1)),
        array_merge([true], array_slice($data['w'], 1)), [1, 2, 3], 5, 1.5, true, "\xff", new ArrayObject([1]),
        array_combine(array_map(fn ($i) => "k$i", range(1, $n)), $data['w']), $data['trials'],
    ];
    foreach ($weights as $w) {
        assertFormulaParity('Phit\fit_ols', ['y ~ x + f', $data, 'drop', $w], 'weights');
        assertFormulaParity('Phit\fit_ols', ['y ~ x + f', $data, 'raise', $w], 'weights');
        assertFormulaParity('Phit\fit_generalized_linear', ['y ~ x + offset(z)', $data, 'gaussian', '', 'drop', $w], 'weights');
        assertFormulaParity('Phit\fit_generalized_linear', ['cnt ~ x', $data, 'poisson', '', 'drop', $w], 'weights');
        assertFormulaParity('Phit\fit_generalized_linear', ['bin ~ x', $data, 'binomial', '', 'drop', $w], 'weights');
        assertFormulaParity('Phit\fit_generalized_linear', ['cbind(succ, trials - succ) ~ x', $data, 'binomial', 'probit', 'drop', $w], 'weights');
    }
});

it('agrees on every family and link, and rejects the others', function () {
    $data = formulaFrame();
    $combinations = [
        ['gaussian', ''], ['gaussian', 'identity'], ['gaussian', 'log'], ['gaussian', 'logit'],
        ['binomial', ''], ['binomial', 'logit'], ['binomial', 'probit'], ['binomial', 'cloglog'], ['binomial', 'identity'],
        ['poisson', ''], ['poisson', 'log'], ['poisson', 'identity'], ['poisson', 'sqrt'],
        ['gamma', ''], ['Gamma', 'inverse'], ['', ''], ['Binomial', 'logit'], ['quasipoisson', 'log'], ["binomial\xff", ''],
        ['binomial', "logit\xff"],
    ];
    $formulas = [
        'y ~ x + f', 'y ~ 0 + f + x', 'y ~ x + offset(z)', 'bin ~ x + z', 'bin ~ x * z', 'l ~ x', 'f ~ x', 'g ~ x', 'o ~ x',
        'I(f == "a") ~ x', 'cbind(succ, trials - succ) ~ x', 'cbind(succ, x) ~ z', 'cbind(trials - succ - 9, succ) ~ x',
        'I(succ / trials) ~ x', 'cnt ~ x + g', 'cnt ~ x + offset(log(expo))', 'cnt ~ 0 + g', 'I(cnt + 0.5) ~ x', 'I(-cnt) ~ x',
        'I(cnt * 0) ~ x', 'I(x > 10) ~ z', 'z ~ x', '~ x', 'y ~ x + (1 | g)', 'cbind(succ, trials - succ, z) ~ x',
    ];
    foreach ($combinations as [$family, $link]) {
        foreach ($formulas as $formula) {
            assertFormulaParity('Phit\fit_generalized_linear', [$formula, $data, $family, $link], 'families');
        }
        assertFormulaParity('Phit\fit_generalized_linear', ['I(succ / trials) ~ x', $data, $family, $link, 'drop', 'trials'], 'families');
        assertFormulaParity('Phit\fit_generalized_linear', ['I(succ / trials) ~ x', $data, $family, $link, 'drop', 'w'], 'families');
    }
});

final class FormulaParitySource
{
    public function __construct(private mixed $c) {}

    public function toModelData(): mixed
    {
        return $this->c;
    }
}

final class FormulaParityColumns
{
    public function __construct(private array $c) {}

    public function toColumns(): array
    {
        return $this->c;
    }

    public function toModelData(int $required): array
    {
        return $this->c;
    }
}

final class FormulaParityPrivate
{
    public function __construct(private array $c) {}

    private function toColumns(): array
    {
        return $this->c;
    }
}

final class FormulaParityMagic
{
    public function __call(string $name, array $arguments): array
    {
        return ['y' => [1.0, 2.0, 4.0], 'x' => [1, 2, 3]];
    }
}

it('agrees on data shapes', function () {
    $base = ['y' => [1.5, 2.1, 2.9, 4.2, 4.8, 6.1], 'x' => [1, 2, 3, 4, 5, 6], 'g' => ['a', 'b', 'a', 'b', 'c', 'c']];
    // A column slot that is a PHP reference counts as a nested value in the extension.
    $referenced = $base;
    $alias = &$referenced['x'][2];
    $referencedNull = $base;
    $referencedNull['x'][1] = null;
    $aliasNull = &$referencedNull['x'][1];
    $shapes = [
        $base, [], [[1, 2], [3, 4]], [['y' => 1, 'x' => 2], ['y' => 2, 'x' => 3]], 5, 'y', null, true,
        ['y' => $base['y'], 0 => $base['x']], ['y' => $base['y'], 'x' => 7], ['y' => $base['y'], 'x' => null],
        ['y' => $base['y'], 'x' => []], ['y' => [], 'x' => []], ['y' => $base['y'], 'x' => [1, 2, 3]],
        ['y' => [1, 2, 3], 'x' => $base['x']], ['y' => $base['y'], 'x' => $base['x'], 'q' => [1, 2]],
        ['y' => $base['y'], 'x' => [1, 2, INF, 4, 5, 6]], ['y' => $base['y'], 'x' => $base['x'], 'q' => [1, 2, -INF, 4, 5, 6]],
        ['y' => $base['y'], 'x' => [1, 2, NAN, 4, 5, 6]], ['y' => $base['y'], 'x' => ['1', '2', '3', '4', '5', '6']],
        ['y' => $base['y'], 'x' => [true, false, true, true, false, true]], ['y' => $base['y'], 'x' => [1, 2, 3, 4, 5, '6']],
        ['y' => $base['y'], 'x' => [1, 2, 3, 4, 5, true]], ['y' => $base['y'], 'x' => [1, 2, 3, 4, 5, [6]]],
        ['y' => $base['y'], 'x' => [1, 2, 3, 4, 5, new stdClass]], ['y' => $base['y'], 'x' => [5 => 1, 4 => 2, 3 => 3, 2 => 4, 1 => 5, 0 => 6]],
        ['y' => $base['y'], 'x' => [10 => 1, 2, 3, 4, 5, 6]], ['y' => $base['y'], 'x' => ['a' => 1, 'b' => 2, 'c' => 3]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 0, 1, 0, null], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 0, 1, 0, 1]]],
        ['y' => $base['y'], 'x' => ['levels' => 'a', 'codes' => [0, 1, 0, 1, 0, 1], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 0, 1, 0, 1.0], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 0, 1, 0, -1], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 0, 1, 0, '1'], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => 'x', 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 0, 1, 0, 1], 'ordered' => 1]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', "\xff"], 'codes' => [0, 1, 0, 1, 0, 1], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'A'], 'codes' => [0, 1, 0, 1, 0, 1], 'ordered' => true]],
        ['y' => $base['y'], 'x' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 0], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['levels' => [], 'codes' => [null, null, null, null, null, null], 'ordered' => false]],
        ['y' => $base['y'], 'x' => ['codes' => [0, 1, 0, 1, 0, 1], 'ordered' => false, 'levels' => ['a', 'b']]],
        ['y' => $base['y'], 'x' => [0 => 1, 'levels' => ['a'], 'ordered' => false]],
        ['y' => $base['y'], 'x' => $base['x'], 'x2' => [1, 2, 3, 4, 5, ["\xff"]]],
        ['y' => $base['y'], 'x' => $base['x'], 'g' => ['a', "\xff", 'a', 'b', 'c', 'c']],
        new FormulaParitySource($base), new FormulaParitySource('nope'), new FormulaParitySource([[1, 2]]),
        new FormulaParityColumns($base), new FormulaParityPrivate($base), new FormulaParityMagic, new ArrayObject($base),
        new stdClass, $referenced, $referencedNull,
    ];
    foreach ($shapes as $data) {
        foreach (['y ~ x', 'y ~ x + g', 'y ~ .'] as $formula) {
            assertFormulaParity('Phit\model_matrix', [$formula, $data], 'shapes');
            assertFormulaParity('Phit\fit_ols', [$formula, $data], 'shapes');
        }
    }
});

it('agrees on random designs with factors, transforms, offsets and weights', function () {
    mt_srand(99);
    $u = fn () => mt_rand() / (mt_getrandmax() + 1.0);
    for ($trial = 0; $trial < 30; $trial++) {
        $n = [15, 40, 120][$trial % 3];
        $data = [];
        for ($i = 0; $i < $n; $i++) {
            $x1 = $u() * 4 - 2;
            $x2 = $u() * 10;
            $f = ['p', 'q', 'r', 's'][mt_rand(0, 3)];
            $eta = 0.3 + 0.5 * $x1 - 0.05 * $x2 + ($f === 'q' ? 0.4 : 0.0);
            $data['x1'][] = $trial % 5 === 1 && $i === 2 ? null : $x1;
            $data['x2'][] = $x2;
            $data['f'][] = $f;
            $data['h'][] = $u() < 0.5;
            $data['y'][] = $eta + $u() - 0.5;
            $data['b'][] = (int) ($u() < 1 / (1 + exp(-$eta)));
            $data['c'][] = (int) floor(exp($eta) * 2 * $u() * 2);
            $data['wt'][] = (float) mt_rand(1, 4);
        }
        $weights = $trial % 2 ? 'wt' : null;
        foreach (['y ~ x1 * f + log(x2 + 1)', 'y ~ poly(x2, 2) + f:h', 'y ~ (x1 + f + h)^2 + offset(x2 / 10)'] as $formula) {
            assertFormulaParity('Phit\model_matrix', [$formula, $data], 'random');
            assertFormulaParity('Phit\fit_ols', [$formula, $data, 'drop', $weights], 'random');
            assertFormulaParity('Phit\fit_generalized_linear', [$formula, $data, 'gaussian', '', 'drop', $weights], 'random');
        }
        assertFormulaParity('Phit\fit_generalized_linear', ['b ~ x1 + f + I(x2^2)', $data, 'binomial', $trial % 3 ? 'logit' : 'probit', 'drop', $weights], 'random');
        assertFormulaParity('Phit\fit_generalized_linear', ['c ~ x1 + h + offset(log(x2 + 1))', $data, 'poisson', '', 'drop', $weights], 'random');
    }
});
