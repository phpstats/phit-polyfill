<?php

declare(strict_types=1);

// Port of the extension's tests/formula.php: Phit\model_matrix, Phit\fit_ols and
// Phit\fit_generalized_linear against R's model.matrix()/lm()/glm() (tests/fixtures/formula-r.json),
// plus the error codes and data conventions. The polyfill raises E_USER_WARNING
// where the extension raises E_WARNING; the text is the same.

/** Compare as tests/formula.php's fcompare: strings/bools exactly, numbers within abs + rel*|expected|. */
function formulaCompare(mixed $actual, mixed $expected, string $label, float $abs, float $rel): void
{
    if (is_array($expected)) {
        expect(is_array($actual) && count($actual) === count($expected))->toBeTrue("$label shape");
        foreach ($expected as $i => $value) {
            formulaCompare($actual[$i] ?? null, $value, "$label/$i", $abs, $rel);
        }

        return;
    }
    if (is_string($expected) || is_bool($expected)) {
        expect($actual)->toBe($expected, $label);

        return;
    }
    expect(is_int($actual) || is_float($actual))->toBeTrue("$label: not numeric");
    expect(abs($actual - $expected) <= $abs + $rel * abs($expected))->toBeTrue("$label: $actual != $expected");
}

/** Assert the exception code and that the message contains $fragment. */
function formulaThrows(callable $call, int $code, string $fragment, string $label): void
{
    try {
        warnings($call);
    } catch (\Throwable $e) {
        expect($e->getCode())->toBe($code, "$label: {$e->getMessage()}")
            ->and(str_contains($e->getMessage(), $fragment))->toBeTrue("$label: message '{$e->getMessage()}' lacks '$fragment'");

        return;
    }
    throw new \RuntimeException("$label: expected exception $code");
}

dataset('formula cases', function () {
    foreach (loadFixture('formula-r.json')['cases'] as $case) {
        yield $case['function'].': '.$case['formula'].($case['weights'] ? " [weights {$case['weights']}]" : '') => [$case];
    }
});

it('matches R on the formula reference cases', function (array $case) {
    $fixture = loadFixture('formula-r.json');
    $data = $fixture['datasets'][$case['dataset']];
    $expected = $case['expected'];
    $warnings = warnings(fn () => match ($case['function']) {
        'model_matrix' => pf('Phit\model_matrix', $case['formula'], $data),
        'lm' => pf('Phit\fit_ols', $case['formula'], $data, weights: $case['weights']),
        'glm' => pf('Phit\fit_generalized_linear', $case['formula'], $data, $case['family'], $case['link'], weights: $case['weights']),
    }, $result);
    $first = $data[array_key_first($data)];
    $dropped = count($first['codes'] ?? $first) - count($expected['retained_rows']);
    expect($warnings)->toHaveCount($dropped > 0 ? 1 : 0);
    if ($dropped > 0) {
        expect($warnings[0][0])->toBe(E_USER_WARNING)
            ->and($warnings[0][1])->toContain("dropped $dropped of");
    }
    expect($result['retained_rows'])->toBe($expected['retained_rows'])
        ->and($result['n_dropped'])->toBe($dropped)
        ->and($result['columns'])->toBe($expected['columns'])
        ->and($result['terms'])->toBe($expected['terms']);
    [$abs, $rel] = $case['function'] === 'model_matrix' ? [1e-12, 1e-12] : [1e-9, 1e-7];
    if (str_contains($case['formula'], 'poly(area_m2, 2)')) {
        [$abs, $rel] = [1e-9, 1e-7];
    }
    foreach ($expected as $key => $value) {
        if (in_array($key, ['columns', 'terms', 'retained_rows'], true)) {
            continue;
        }
        expect($result)->toHaveKey($key);
        formulaCompare($result[$key], $value, "{$case['name']}/$key", $abs, $rel);
    }
    if ($case['function'] !== 'model_matrix') {
        expect(array_keys($result['named']))->toBe($expected['columns'])
            ->and(array_values($result['named']))->toBe($result['coefficients'])
            ->and($result['diagnostics']['converged'])->toBeTrue();
    }
})->with('formula cases');

it('covers every fixture function', function () {
    $counts = array_count_values(array_column(loadFixture('formula-r.json')['cases'], 'function'));
    expect($counts['model_matrix'])->toBeGreaterThanOrEqual(30)
        ->and($counts['lm'])->toBeGreaterThanOrEqual(12)
        ->and($counts['glm'])->toBeGreaterThanOrEqual(12)
        ->and(array_sum($counts))->toBe(62);
});

function formulaBase(): array
{
    return ['y' => [1.5, 2.1, 2.9, 4.2, 4.8, 6.1], 'x' => [1, 2, 3, 4, 5, 6], 'g' => ['a', 'b', 'a', 'b', 'c', 'c']];
}

it('raises the extension error codes', function () {
    $syn = loadFixture('formula-r.json')['datasets']['synthetic'];
    $base = formulaBase();
    $cases = [
        ['unknown column', fn () => pf('Phit\fit_ols', 'y ~ nope', $base), 1001, "unknown column 'nope'"],
        ['length mismatch', fn () => pf('Phit\fit_ols', 'y ~ x', $base + ['z' => [1, 2]]), 1002, "column 'z'"],
        ['list column', fn () => pf('Phit\fit_ols', 'y ~ x + z', $base + ['z' => [[1], [2], [3], [4], [5], [6]]]), 1001, 'list-columns'],
        ['mixed column', fn () => pf('Phit\fit_ols', 'y ~ x + z', $base + ['z' => [1, 'a', 2, 3, 4, 5]]), 1001, 'mix'],
        ['descriptor keys', fn () => pf('Phit\fit_ols', 'y ~ f', $base + ['f' => ['levels' => ['a'], 'codes' => [0, 0, 0, 0, 0, 0], 'ordered' => false, 'x' => 1]]), 1001, 'factor descriptor'],
        ['descriptor code', fn () => pf('Phit\fit_ols', 'y ~ f', $base + ['f' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 2, 0, 1, 0], 'ordered' => false]]), 1001, '0-based'],
        ['descriptor levels', fn () => pf('Phit\fit_ols', 'y ~ f', $base + ['f' => ['levels' => [1, 2], 'codes' => [0, 1, 1, 0, 1, 0], 'ordered' => false]]), 1001, 'strings'],
        ['descriptor ordered', fn () => pf('Phit\fit_ols', 'y ~ f', $base + ['f' => ['levels' => ['a', 'b'], 'codes' => [0, 1, 1, 0, 1, 0], 'ordered' => 'no']]), 1001, 'boolean'],
        ['duplicate levels', fn () => pf('Phit\fit_ols', 'y ~ f', $base + ['f' => ['levels' => ['a', 'a'], 'codes' => [0, 1, 1, 0, 1, 0], 'ordered' => false]]), 1001, 'duplicate levels'],
        ['row list data', fn () => pf('Phit\fit_ols', 'y ~ x', [['y' => 1, 'x' => 2]]), 1001, 'associative array of columns'],
        ['scalar data', fn () => pf('Phit\fit_ols', 'y ~ x', 5), 1001, 'toModelData'],
        ['raise', fn () => pf('Phit\fit_ols', 'y ~ x1', $syn, 'raise'), 1001, "missing='raise'"],
        ['model matrix raise', fn () => pf('Phit\model_matrix', 'y ~ x1', $syn, 'raise'), 1001, "missing='raise'"],
        ['policy', fn () => pf('Phit\fit_ols', 'y ~ x', $base, 'keep'), 1001, "'drop' or 'raise'"],
        ['syntax', fn () => pf('Phit\fit_ols', 'y ~ x +', $base), 1001, 'formula'],
        ['no tilde', fn () => pf('Phit\fit_ols', 'y + x', $base), 1001, "'~'"],
        ['bare bar', fn () => pf('Phit\model_matrix', 'y ~ x + log(x | g)', $base), 1001, ''],
        ['unsupported function', fn () => pf('Phit\fit_ols', 'y ~ ns(x, 3)', $base), 1007, "'ns()'"],
        ['random in lm', fn () => pf('Phit\fit_ols', 'y ~ x + (1 | g)', $base), 1007, 'Phit\fit_mixed_gaussian'],
        ['random in glm', fn () => pf('Phit\fit_generalized_linear', 'y ~ x + (1 | g)', $base), 1007, 'Phit\fit_mixed_gaussian_matrix'],
        ['family', fn () => pf('Phit\fit_generalized_linear', 'y ~ x', $base, 'gamma'), 1007, "family 'gamma'"],
        ['link', fn () => pf('Phit\fit_generalized_linear', 'y ~ x', $base, 'poisson', 'identity'), 1007, "link 'identity'"],
        ['rank deficiency', fn () => pf('Phit\fit_ols', 'y ~ x + I(2 * x)', $base), 1004, ''],
        ['unused level with intercept-free factor', fn () => pf('Phit\fit_ols', 'y ~ 0 + g + I(g == "a")', $base), 1004, ''],
        ['binomial range', fn () => pf('Phit\fit_generalized_linear', 'y ~ x', $base, 'binomial'), 1001, '[0,1]'],
        ['multi-level response in lm', fn () => pf('Phit\fit_ols', 'g ~ x', $base), 1001, 'numeric or logical response'],
        ['poisson counts', fn () => pf('Phit\fit_generalized_linear', 'y ~ x', $base, 'poisson'), 1001, 'nonnegative integer'],
        ['one-sided lm', fn () => pf('Phit\fit_ols', '~ x', $base), 1001, 'response'],
        ['empty design', fn () => pf('Phit\fit_ols', 'y ~ 0', $base), 1001, 'no predictor columns'],
        ['log zero', fn () => pf('Phit\fit_ols', 'y ~ log(x - 1)', $base), 1003, 'infinite'],
        ['infinite data', fn () => pf('Phit\fit_ols', 'y ~ x', ['y' => [1, 2, 3], 'x' => [1, INF, 3]]), 1003, 'infinite'],
        ['single level contrast', fn () => pf('Phit\fit_ols', 'y ~ x + f', $base + ['f' => ['u', 'u', 'u', 'u', 'u', 'u']]), 1001, '2 or more levels'],
        ['weights length', fn () => pf('Phit\fit_ols', 'y ~ x', $base, weights: [1, 2]), 1002, 'one value per'],
        ['weights sign', fn () => pf('Phit\fit_ols', 'y ~ x', $base, weights: [1, 2, 0, 1, 1, 1]), 1001, 'positive'],
        ['weights column', fn () => pf('Phit\fit_ols', 'y ~ x', $base, weights: 'g'), 1001, 'numeric'],
        ['weights unknown', fn () => pf('Phit\fit_ols', 'y ~ x', $base, weights: 'nope'), 1001, 'unknown weights column'],
        ['orthogonal poly NA', fn () => pf('Phit\fit_ols', 'y ~ poly(x2, 2)', $syn), 1001, 'raw = TRUE'],
        ['binomial weights', fn () => pf('Phit\fit_generalized_linear', 'I(x / 10) ~ y', $base, 'binomial', weights: [1.5, 1, 1, 1, 1, 1]), 1001, 'integer'],
    ];
    foreach ($cases as [$label, $call, $code, $fragment]) {
        formulaThrows($call, $code, $fragment, $label);
    }
});

it('validates only the columns a formula references', function () {
    $base = formulaBase();
    $withList = $base + ['meta' => [[1], [2], [3], [4], [5], [6]], 'mixed' => [1, 'a', true, null, 2, 3]];
    expect(pf('Phit\fit_ols', 'y ~ x', $withList)['coefficients'])->toBe(pf('Phit\fit_ols', 'y ~ x', $base)['coefficients']);
});

it('accepts weights by column name or aligned array, and deletes null weights', function () {
    $base = formulaBase();
    $wdata = $base + ['w' => [1, 2, 1, 3, 1, 2]];
    expect(pf('Phit\fit_ols', 'y ~ x', $wdata, weights: 'w')['coefficients'])
        ->toBeClose(pf('Phit\fit_ols', 'y ~ x', $base, weights: [1, 2, 1, 3, 1, 2])['coefficients'], 1e-14, 1e-14);
    $warn = warnings(fn () => pf('Phit\fit_ols', 'y ~ x', $base, weights: [1, 2, null, 3, 1, 2]), $wnull);
    expect($wnull['dropped_rows'])->toBe([2])
        ->and($warn)->toHaveCount(1)
        ->and($warn[0][1])->toBe('Phit\fit_ols(): dropped 1 of 6 observations with missing values (listwise deletion)');
    expect(pf('Phit\fit_ols_matrix', array_map(fn ($v) => [$v], $base['x']), $base['y'])['coefficients'])
        ->toBe(pf('Phit\fit_ols', 'y ~ x', $base)['coefficients']);
});

it('deletes rows only over the variables the formula uses', function () {
    $base = formulaBase();
    $syn = loadFixture('formula-r.json')['datasets']['synthetic'];
    $warn = warnings(fn () => pf('Phit\fit_ols', 'y ~ x', $base + ['z' => array_fill(0, 6, null)]), $lm);
    expect($lm['n_dropped'])->toBe(0)->and($warn)->toBe([]);
    $warn = warnings(fn () => pf('Phit\model_matrix', 'y ~ x1 + f', $syn), $mm);
    expect($mm['n_input'])->toBe(90)
        ->and($mm['n_observations'])->toBe(count($mm['x']))
        ->and($mm['n_observations'] + $mm['n_dropped'])->toBe(90)
        ->and($warn)->toHaveCount(1)
        ->and($warn[0][0])->toBe(E_USER_WARNING)
        ->and($warn[0][1])->toStartWith('Phit\model_matrix(): dropped ');
    set_error_handler(function (int $severity, string $message): never {
        throw new ErrorException($message, 0, $severity);
    });
    $promoted = false;
    try {
        pf('Phit\fit_ols', 'y ~ x1', $syn);
    } catch (ErrorException) {
        $promoted = true;
    } finally {
        restore_error_handler();
    }
    expect($promoted)->toBeTrue();
});

it('fits logical, factor and character binomial responses with names', function () {
    $bdata = ['ok' => [true, false, true, true, false, true, false, false], 'x' => [1, 2, 3, 4, 5, 6, 7, 8],
        'res' => ['levels' => ['fail', 'pass'], 'codes' => [1, 0, 1, 1, 0, 1, 0, 0], 'ordered' => false],
        'txt' => ['pass', 'fail', 'pass', 'pass', 'fail', 'pass', 'fail', 'fail']];
    $logical = pf('Phit\fit_generalized_linear', 'ok ~ x', $bdata, 'binomial');
    expect(pf('Phit\fit_generalized_linear', 'res ~ x', $bdata, 'binomial')['coefficients'])->toBeClose($logical['coefficients'], 1e-12, 1e-12)
        ->and(pf('Phit\fit_generalized_linear', 'txt ~ x', $bdata, 'binomial')['coefficients'])->toBeClose($logical['coefficients'], 1e-12, 1e-12)
        ->and(pf('Phit\fit_logistic_matrix', array_map(fn ($v) => [$v], $bdata['x']), [1, 0, 1, 1, 0, 1, 0, 0])['coefficients'])->toBe($logical['coefficients'])
        ->and($logical['family'])->toBe('binomial')
        ->and($logical['link'])->toBe('logit')
        ->and($logical['diagnostics']['solver'])->toBe('newton_irls')
        ->and(pf('Phit\fit_generalized_linear', 'ok ~ x', $bdata, 'binomial', 'probit')['diagnostics']['solver'])->toBe('fisher_scoring')
        ->and(pf('Phit\fit_generalized_linear', 'ok ~ x', $bdata)['diagnostics']['solver'])->toBe('scaled_svd')
        ->and($logical['formula'])->toBe('ok ~ x')
        ->and($logical['response'])->toBe('ok')
        ->and($logical['assign'])->toBe([0, 1])
        ->and($logical['intercept'])->toBeTrue();
});

final class FormulaModelDataSource
{
    public function __construct(private array $c) {}

    public function toModelData(): array
    {
        return $this->c;
    }
}

final class FormulaColumnsSource
{
    public function __construct(private array $c) {}

    public function toColumns(): array
    {
        return $this->c;
    }
}

final class FormulaThrowingSource
{
    public function toModelData(): array
    {
        throw new DomainException('source failed', 77);
    }
}

it('reads objects through toModelData() or toColumns()', function () {
    $base = formulaBase();
    expect(pf('Phit\fit_ols', 'y ~ x + g', new FormulaModelDataSource($base))['coefficients'])->toBe(pf('Phit\fit_ols', 'y ~ x + g', $base)['coefficients'])
        ->and(pf('Phit\model_matrix', 'y ~ x + g', new FormulaColumnsSource($base))['x'])->toBe(pf('Phit\model_matrix', 'y ~ x + g', $base)['x']);
    formulaThrows(fn () => pf('Phit\fit_ols', 'y ~ x', new ArrayObject($base)), 1001, 'toModelData', 'object without methods');
    expect(fn () => pf('Phit\fit_ols', 'y ~ x', new FormulaThrowingSource))->toThrow(DomainException::class, 'source failed');
    expect(thrownCode(fn () => pf('Phit\fit_ols', 'y ~ x', new FormulaThrowingSource)))->toBe(77);
});

it('returns random-effect bars as explicit z and group designs', function () {
    $schools = loadFixture('formula-r.json')['datasets']['schools'];
    $mm = pf('Phit\model_matrix', 'score ~ study_hours + (1 + study_hours | school) + (0 + prior_score || school)', $schools);
    expect($mm['columns'])->toBe(['(Intercept)', 'study_hours'])->and($mm['random'])->toHaveCount(2);
    [$slopes, $independent] = $mm['random'];
    expect($slopes['label'])->toBe('1 + study_hours | school')
        ->and($slopes['group'])->toBe('school')
        ->and($slopes['correlated'])->toBeTrue()
        ->and($slopes['columns'])->toBe(['(Intercept)', 'study_hours'])
        ->and($independent['columns'])->toBe(['prior_score'])
        ->and($independent['correlated'])->toBeFalse()
        ->and($slopes['groups'])->toBe($schools['school'])
        ->and($slopes['z'][3])->toBe([1.0, (float) $schools['study_hours'][3]]);
    $base = formulaBase();
    $nested = pf('Phit\model_matrix', 'y ~ x + (1 | g/f) + (x | g:f)', $base + ['f' => [1, 1, 2, 2, 1, 2]]);
    expect(array_column($nested['random'], 'label'))->toBe(['1 | g', '1 | f:g', 'x | g:f'])
        ->and($nested['random'][1]['groups'][0])->toBe('1:a')
        ->and($nested['random'][2]['groups'][0])->toBe('a:1')
        ->and($nested['random'][0]['groups'][0])->toBe('a')
        ->and(pf('Phit\model_matrix', 'y ~ x + (1 | x)', $base)['random'][0]['groups'][0])->toBe(1);
});

it('follows R term order, naming and coding rules', function () {
    $data = [
        'y' => [1.0, 3.0, 2.0, 5.0, 4.0, null],
        'x' => [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        'a' => ['p', 'q', 'r', 'p', 'q', 'r'],
        'b' => ['u', 'v', 'u', 'v', 'u', 'v'],
    ];
    $columns = function (string $formula) use ($data): array {
        warnings(fn () => pf('Phit\model_matrix', $formula, $data), $result);

        return $result['columns'];
    };
    expect($columns('y ~ a:b'))->toBe(['(Intercept)', 'ap:bu', 'aq:bu', 'ar:bu', 'ap:bv', 'aq:bv', 'ar:bv'])
        ->and($columns('y ~ b:a + a'))->toBe(['(Intercept)', 'aq', 'ar', 'bv:ap', 'bv:aq', 'bv:ar'])
        ->and($columns('y ~ a*b - 1'))->toBe(['ap', 'aq', 'ar', 'bv', 'aq:bv', 'ar:bv'])
        ->and($columns('y ~ a/x'))->toBe(['(Intercept)', 'aq', 'ar', 'ap:x', 'aq:x', 'ar:x'])
        ->and($columns('y ~ x %in% a'))->toBe(['(Intercept)', 'x:ap', 'x:aq', 'x:ar']);
    warnings(fn () => pf('Phit\model_matrix', 'y ~ poly(x, 2, raw = TRUE) + I(x > 3) + log(x, base = 2)', $data), $d);
    expect($d['columns'])->toBe(['(Intercept)', 'poly(x, 2, raw = TRUE)1', 'poly(x, 2, raw = TRUE)2', 'I(x > 3)TRUE', 'log(x, base = 2)'])
        ->and($d['dropped_rows'])->toBe([5])
        ->and($d['assign'])->toBe([0, 1, 1, 2, 3]);
    $ordered = pf('Phit\model_matrix', '~ o', ['o' => ['levels' => ['a', 'b', 'c'], 'codes' => [0, 1, 2], 'ordered' => true]]);
    expect($ordered['columns'])->toBe(['(Intercept)', 'o.L', 'o.Q'])
        ->and($ordered['x'][0][1])->toBeClose(-M_SQRT1_2, 1e-15, 0)
        ->and($ordered['x'][1][2])->toBeClose(-0.816496580927726, 1e-15, 0);
});

it('deparses formulas and numbers as R does', function () {
    $data = ['y' => [1.0, 2.0, 3.5, 4.0], 'x' => [1.0, 2.0, 3.0, 4.5], 'a b' => [2.0, 1.0, 4.0, 3.0]];
    expect(pf('Phit\model_matrix', 'y~log(x,base=2)+I(x^2*3)', $data)['formula'])->toBe('y ~ log(x, base = 2) + I(x^2 * 3)')
        ->and(pf('Phit\model_matrix', 'y ~ I(x * 100000) + I(x * 0.0001) + I(x/3)', $data)['columns'])
        ->toBe(['(Intercept)', 'I(x * 1e+05)', 'I(x * 1e-04)', 'I(x/3)'])
        ->and(pf('Phit\model_matrix', 'y ~ `a b` + -1', $data)['columns'])->toBe(['`a b`'])
        ->and(pf('Phit\model_matrix', 'y ~ factor(x)', $data)['columns'])->toBe(['(Intercept)', 'factor(x)2', 'factor(x)3', 'factor(x)4.5']);
});

// Review finding: a power beyond PHP_INT_MAX was cast to int (0, with a PHP
// warning) and dropped the interactions. The expansion now stops once it no
// longer changes, which matches the extension's terms for any power.
it('expands very large interaction powers to the full interaction', function () {
    $data = ['y' => [1, 2, 3, 4, 5, 6], 'x' => [1, 3, 2, 5, 4, 6], 'z' => [2, 1, 4, 3, 6, 5]];
    $warnings = warnings(fn () => pf('Phit\model_matrix', 'y ~ (x + z)^1e300', $data), $result);
    expect($warnings)->toBe([])
        ->and($result['columns'])->toBe(['(Intercept)', 'x', 'z', 'x:z'])
        ->and(pf('Phit\model_matrix', 'y ~ (x + z)^5', $data)['columns'])->toBe($result['columns']);
});
