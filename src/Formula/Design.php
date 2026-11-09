<?php

declare(strict_types=1);

namespace Phit\Polyfill\Formula;

use Phit\Polyfill\ErrorCode;

/**
 * R's model.frame()/model.matrix() for a formula over named columns,
 * following stats-core formula::model_design: listwise deletion over the variables
 * the formula uses, treatment (unordered) or polynomial (ordered) contrasts,
 * R column names and term order, and random-effect bars as z designs.
 *
 * @internal
 */
final class Design
{
    public string $formula;

    public ?string $responseName;

    /**
     * null, ['numeric', list<float>], ['factor', levels, codes] or
     * ['binomial', successes, failures].
     */
    public ?array $response;

    /** @var list<list<float>> */
    public array $x;

    /** @var list<string> */
    public array $columns;

    /** @var list<int> */
    public array $assign;

    public bool $intercept;

    /** @var list<string> */
    public array $terms;

    /** @var ?list<float> */
    public ?array $offset;

    /** @var list<array{label: string, group: string, columns: list<string>, correlated: bool, z: list<list<float>>, groups: list<int|string>}> */
    public array $random;

    /** @var list<list<float>> */
    public array $extra;

    public int $nInput;

    /** @var list<int> */
    public array $retainedRows;

    /** @var list<int> */
    public array $droppedRows;

    /**
     * What prediction must reuse, as R's xlevels and predvars: factor levels
     * after deletion by variable, and orthogonal poly() coefficients by call.
     *
     * @var array{levels: array<string, list<string>>, poly: array<string, array{alpha: list<float>, norm2: list<float>}>}
     */
    public array $predvars;

    private static function error(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    private static function subset(array $v, array $rows): array
    {
        $out = [];
        foreach ($rows as $i) {
            $out[] = $v[$i];
        }

        return $out;
    }

    /** Values after deletion: ['num', v], ['mat', cols, names] or ['factor', levels, codes, ordered]. */
    private static function code(array $value, array $rows, bool $dropLevels): array
    {
        switch ($value['t']) {
            case 'num':
                return ['num', self::subset($value['v'], $rows)];
            case 'mat':
                return ['mat', array_map(fn ($c) => self::subset($c, $rows), $value['cols']), $value['names']];
            case 'lgl':
                return ['factor', ['FALSE', 'TRUE'], array_map(fn ($b) => (int) $b, self::subset($value['v'], $rows)), false];
            case 'text':
                $kept = self::subset($value['v'], $rows);
                $levels = Evaluator::sortedLevels($kept);
                $index = array_flip($levels);

                return ['factor', $levels, array_map(fn ($s) => $index[$s] ?? 0, $kept), false];
            case 'factor':
                $kept = self::subset($value['codes'], $rows);
                $levels = $value['levels'];
                if ($dropLevels) {
                    [$levels, $kept] = Evaluator::dropUnused($levels, $kept);
                }

                return ['factor', $levels, $kept, $value['ordered']];
            default:
                throw self::error('a string literal cannot be a model variable');
        }
    }

    /**
     * Code a variable against the levels it had when the model was fitted, as
     * predict() does; a level the fit never saw is an error, as in R.
     */
    private static function codeFrozen(string $name, array $value, array $rows, ?array $levels): array
    {
        switch ($value['t']) {
            case 'lgl':
                $labels = array_map(fn ($b) => $b ? 'TRUE' : 'FALSE', self::subset($value['v'], $rows));
                break;
            case 'text':
                $labels = self::subset($value['v'], $rows);
                break;
            case 'factor':
                $labels = array_map(fn ($c) => $value['levels'][$c], self::subset($value['codes'], $rows));
                break;
            default:
                if ($levels !== null) {
                    throw self::error("variable '$name' was a factor when the model was fitted");
                }

                return self::code($value, $rows, false);
        }
        $ordered = $value['t'] === 'factor' && $value['ordered'];
        if ($levels === null) {
            throw self::error("variable '$name' was numeric when the model was fitted");
        }
        $index = [];
        foreach ($levels as $i => $level) {
            $index[$level] ??= $i;
        }
        $codes = [];
        foreach ($labels as $label) {
            $codes[] = $index[$label] ?? throw self::error("factor $name has new level $label");
        }

        return ['factor', array_values($levels), $codes, $ordered];
    }

    /** Contrast matrix (levels x columns) and column suffixes. */
    private static function contrasts(string $name, array $levels, bool $ordered, bool $dummy): array
    {
        $k = count($levels);
        if ($dummy) {
            $m = [];
            for ($i = 0; $i < $k; $i++) {
                for ($j = 0; $j < $k; $j++) {
                    $m[$i][$j] = $i === $j ? 1.0 : 0.0;
                }
            }

            return [$m, $levels];
        }
        if ($k < 2) {
            throw self::error("contrasts can be applied only to factors with 2 or more levels ('$name' has $k after deletion)");
        }
        if ($ordered) {
            $basis = Evaluator::orthonormalPoly(array_map('floatval', range(1, $k)), $k - 1);
            $m = [];
            for ($i = 0; $i < $k; $i++) {
                foreach ($basis as $j => $c) {
                    $m[$i][$j] = $c[$i];
                }
            }
            $names = [];
            for ($d = 1; $d < $k; $d++) {
                $names[] = match ($d) {
                    1 => '.L',
                    2 => '.Q',
                    3 => '.C',
                    default => "^$d",
                };
            }

            return [$m, $names];
        }
        $m = [];
        for ($i = 0; $i < $k; $i++) {
            $m[$i] = [];
            for ($j = 1; $j < $k; $j++) {
                $m[$i][] = $i === $j ? 1.0 : 0.0;
            }
        }

        return [$m, array_slice($levels, 1)];
    }

    private static function isFactor(?array $coded): bool
    {
        return $coded !== null && $coded[0] === 'factor' && count($coded[1]) > 1;
    }

    /**
     * R's model.matrix for one terms object over coded (retained) variables.
     *
     * @return array{columns: list<list<float>>, names: list<string>, assign: list<int>, labels: list<string>}
     */
    private static function buildBlock(Terms $terms, array $coded, int $m): array
    {
        $codes = [];
        foreach ($terms->terms as $t => $term) {
            $codes[$t] = array_map(fn ($v) => $terms->coding($t, $v), $term);
        }
        if (! $terms->intercept) {
            // First factor (by term, then variable order) is coded by indicators.
            foreach ($terms->terms as $t => $term) {
                foreach ($term as $k => $v) {
                    if (self::isFactor($coded[$v])) {
                        $codes[$t][$k] = 2;
                        break 2;
                    }
                }
            }
        }
        $block = [
            'columns' => [],
            'names' => [],
            'assign' => [],
            'labels' => array_map(fn ($t) => $terms->label($t), $terms->terms),
        ];
        if ($terms->intercept) {
            $block['columns'][] = array_fill(0, $m, 1.0);
            $block['names'][] = '(Intercept)';
            $block['assign'][] = 0;
        }
        foreach ($terms->terms as $t => $term) {
            $current = [['', array_fill(0, $m, 1.0)]];
            foreach ($term as $k => $v) {
                $name = $terms->variables[$v][0];
                $value = $coded[$v];
                $parts = [];
                if ($value[0] === 'num') {
                    $parts[] = [$name, $value[1]];
                } elseif ($value[0] === 'mat') {
                    foreach ($value[1] as $c => $column) {
                        $parts[] = [$name.$value[2][$c], $column];
                    }
                } else {
                    [$matrix, $suffixes] = self::contrasts($name, $value[1], $value[3], $codes[$t][$k] === 2);
                    foreach ($suffixes as $j => $suffix) {
                        $parts[] = [$name.$suffix, array_map(fn ($c) => $matrix[$c][$j], $value[2])];
                    }
                }
                $next = [];
                foreach ($parts as [$pname, $pvalues]) {
                    foreach ($current as [$cname, $cvalues]) {
                        $product = [];
                        foreach ($cvalues as $i => $a) {
                            $product[] = $a * $pvalues[$i];
                        }
                        $next[] = [$cname === '' ? $pname : "$cname:$pname", $product];
                    }
                }
                $current = $next;
            }
            foreach ($current as [$name, $values]) {
                $block['names'][] = $name;
                $block['columns'][] = $values;
                $block['assign'][] = $t + 1;
            }
        }

        return $block;
    }

    private static function groupLabel(array $value, int $i): string
    {
        return match ($value['t']) {
            'num' => Parser::formatNumber($value['v'][$i] ?? NAN),
            'lgl' => $value['v'][$i] === true ? 'TRUE' : 'FALSE',
            'text' => $value['v'][$i] ?? '',
            'factor' => $value['codes'][$i] === null ? '' : $value['levels'][$value['codes'][$i]],
            default => '',
        };
    }

    private static function checkFinite(iterable $values, string $what): void
    {
        foreach ($values as $v) {
            if (! is_finite($v)) {
                throw ErrorCode::exception(ErrorCode::NONFINITE, "formula produced infinite values in the $what (for example log(0))");
            }
        }
    }

    /** @param list<list<float>> $rows */
    private static function flatten(array $rows): \Generator
    {
        foreach ($rows as $row) {
            yield from $row;
        }
    }

    /** @param list<list<float>> $columns  @return list<list<float>> */
    private static function rowsOf(array $columns, int $m): array
    {
        $rows = [];
        for ($i = 0; $i < $m; $i++) {
            $row = [];
            foreach ($columns as $c) {
                $row[] = $c[$i];
            }
            $rows[] = $row;
        }

        return $rows;
    }

    /**
     * Parse a formula, evaluate it over data and build the model design.
     *
     * @param  list<list<?float>>  $extra  aligned vectors (weights) that take part in deletion
     */
    public static function build(string $formula, Data $data, bool $dropMissing, bool $dropUnusedLevels, array $extra = [], bool $allowEmpty = false): self
    {
        return self::design($formula, $data, $dropMissing, $dropUnusedLevels, $extra, null, false, $allowEmpty);
    }

    /**
     * The design for new data under a fitted model, as R's predict() builds it:
     * the response is not needed, factors keep the fitted levels, poly() keeps
     * the fitted basis, and rows with missing predictors are reported as
     * dropped rather than raised. Random-effect terms are evaluated only when
     * requested.
     */
    public static function prediction(string $formula, Data $data, array $predvars, bool $random): self
    {
        return self::design($formula, $data, true, false, [], $predvars, $random);
    }

    private static function design(string $formula, Data $data, bool $dropMissing, bool $dropUnusedLevels, array $extra, ?array $predvars, bool $random, bool $allowEmpty = false): self
    {
        $predict = $predvars !== null;
        Evaluator::beginDesign($predict ? $predvars['poly'] : null);
        try {
            return self::designWith($formula, $data, $dropMissing, $dropUnusedLevels, $extra, $predvars, $random, $allowEmpty);
        } finally {
            Evaluator::beginDesign(null);
        }
    }

    private static function designWith(string $formula, Data $data, bool $dropMissing, bool $dropUnusedLevels, array $extra, ?array $predvars, bool $random, bool $allowEmpty = false): self
    {
        $predict = $predvars !== null;
        $levelsOf = fn (string $name): ?array => $predict ? ($predvars['levels'][$name] ?? null) : null;
        $parsed = Parser::parse($formula);
        if ($parsed[0] === 'binary' && $parsed[1] === '~') {
            $responseExpr = $parsed[2];
            $rhs = $parsed[3];
        } elseif ($parsed[0] === 'unary' && $parsed[1] === '~') {
            $responseExpr = null;
            $rhs = $parsed[2];
        } else {
            throw self::error("formula must contain '~', e.g. 'y ~ x'");
        }
        $parts = [];
        $bars = [];
        Terms::splitBars($rhs, $parts, $bars, false);
        $fixed = Terms::build($responseExpr, $parts, $data->names);
        if ($predict && $random === false) {
            $bars = [];
        }
        $specs = [];
        foreach ($bars as $bar) {
            $lhsParts = [];
            $nested = [];
            Terms::splitBars($bar['lhs'], $lhsParts, $nested, false);
            if ($nested !== []) {
                throw self::error('random-effect bars cannot be nested');
            }
            $terms = Terms::build(null, $lhsParts, $data->names);
            if ($terms->offsets !== []) {
                throw self::error('offset() is not allowed inside a random-effect term');
            }
            if ($terms->terms === [] && ! $terms->intercept) {
                throw self::error("random-effect term '".Terms::barLabel($bar)."' has no columns");
            }
            $groups = array_map(fn ($g) => Evaluator::evaluate($g, $data), $bar['groups']);
            foreach ($groups as $g) {
                if ($g['t'] === 'mat' || $g['t'] === 'str') {
                    throw self::error('grouping variables must be single columns');
                }
            }
            $specs[] = ['bar' => $bar, 'terms' => $terms, 'groups' => $groups];
        }
        $n = $data->rows;
        if ($n === 0) {
            throw self::error('data must contain at least one row');
        }
        foreach ($extra as $e) {
            if (count($e) !== $n) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'weights must have one value per original observation');
            }
        }
        $values = [];
        foreach ($fixed->variables as $v => [, $expr]) {
            // predict() does not need the response; it may be absent.
            $values[$v] = $predict && $v === $fixed->response
                ? ['t' => 'num', 'v' => array_fill(0, $n, 0.0)]
                : Evaluator::evaluate($expr, $data);
        }
        $barValues = array_map(
            fn ($s) => array_map(fn ($v) => Evaluator::evaluate($v[1], $data), $s['terms']->variables),
            $specs,
        );
        $retained = [];
        $dropped = [];
        for ($i = 0; $i < $n; $i++) {
            $missing = false;
            foreach ($values as $v) {
                $missing = $missing || Evaluator::missing($v, $i);
            }
            foreach ($extra as $e) {
                $missing = $missing || $e[$i] === null;
            }
            foreach ($barValues as $list) {
                foreach ($list as $v) {
                    $missing = $missing || Evaluator::missing($v, $i);
                }
            }
            foreach ($specs as $s) {
                foreach ($s['groups'] as $g) {
                    $missing = $missing || Evaluator::missing($g, $i);
                }
            }
            if ($missing) {
                $dropped[] = $i;
            } else {
                $retained[] = $i;
            }
        }
        if ($dropped !== [] && ! $dropMissing) {
            throw self::error("missing values found with missing='raise'");
        }
        if ($retained === [] && $predict === false) {
            throw self::error('no complete observations remain after removing missing values');
        }
        $m = count($retained);
        $coded = [];
        foreach ($values as $v => $value) {
            $name = $fixed->variables[$v][0];
            $coded[$v] = match (true) {
                $v === $fixed->response => null,
                $predict => self::codeFrozen($name, $value, $retained, $levelsOf($name)),
                default => self::code($value, $retained, $dropUnusedLevels),
            };
        }
        $block = self::buildBlock($fixed, $coded, $m);
        if ($block['columns'] === [] && !$allowEmpty) {
            throw self::error('the formula has no predictor columns; keep the intercept or add a term');
        }
        $x = self::rowsOf($block['columns'], $m);
        self::checkFinite(self::flatten($x), 'design matrix');
        $offset = null;
        if ($fixed->offsets !== []) {
            $offset = array_fill(0, $m, 0.0);
            foreach ($fixed->offsets as $v) {
                foreach ($coded[$v][1] as $i => $o) {
                    $offset[$i] += $o;
                }
            }
            self::checkFinite($offset, 'offset');
        }
        $response = null;
        if ($fixed->response !== null && $predict === false) {
            $value = $values[$fixed->response];
            if ($value['t'] === 'num') {
                $response = ['numeric', self::subset($value['v'], $retained)];
            } elseif ($value['t'] === 'lgl') {
                $response = ['numeric', array_map(fn ($b) => $b ? 1.0 : 0.0, self::subset($value['v'], $retained))];
            } elseif ($value['t'] === 'mat') {
                if (count($value['cols']) !== 2) {
                    throw self::error('a matrix response must be cbind(successes, failures)');
                }
                $response = ['binomial', self::subset($value['cols'][0], $retained), self::subset($value['cols'][1], $retained)];
            } else {
                $c = self::code($value, $retained, $dropUnusedLevels);
                if ($c[0] !== 'factor') {
                    throw self::error('unsupported response type');
                }
                $response = ['factor', $c[1], $c[2]];
            }
        }
        if ($response !== null && $response[0] === 'numeric') {
            self::checkFinite($response[1], 'response');
        } elseif ($response !== null && $response[0] === 'binomial') {
            self::checkFinite([...$response[1], ...$response[2]], 'response');
        }
        $levels = [];
        $record = function (array $variables, array $coded) use (&$levels): void {
            foreach ($variables as $v => [$name]) {
                $value = $coded[$v] ?? null;
                if ($value !== null && $value[0] === 'factor' && ! array_key_exists($name, $levels)) {
                    $levels[$name] = $value[1];
                }
            }
        };
        $record($fixed->variables, $coded);
        $random = [];
        foreach ($specs as $s => $spec) {
            $barCoded = [];
            foreach ($barValues[$s] as $v => $value) {
                $name = $spec['terms']->variables[$v][0];
                $barCoded[$v] = $predict
                    ? self::codeFrozen($name, $value, $retained, $levelsOf($name))
                    : self::code($value, $retained, $dropUnusedLevels);
            }
            $record($spec['terms']->variables, $barCoded);
            $barBlock = self::buildBlock($spec['terms'], $barCoded, $m);
            $z = self::rowsOf($barBlock['columns'], $m);
            self::checkFinite(self::flatten($z), 'random-effect design');
            $singleInteger = false;
            if (count($spec['groups']) === 1 && $spec['groups'][0]['t'] === 'num') {
                $singleInteger = true;
                foreach ($retained as $i) {
                    $g = $spec['groups'][0]['v'][$i];
                    if ($g === null || fmod($g, 1.0) != 0.0 || is_nan($g) || abs($g) >= 9.0e15) {
                        $singleInteger = false;
                        break;
                    }
                }
            }
            $groups = [];
            foreach ($retained as $i) {
                if ($singleInteger) {
                    $groups[] = (int) $spec['groups'][0]['v'][$i];
                } else {
                    $groups[] = implode(':', array_map(fn ($g) => self::groupLabel($g, $i), $spec['groups']));
                }
            }
            $random[] = [
                'label' => Terms::barLabel($spec['bar']),
                'group' => Terms::barGroupLabel($spec['bar']),
                'columns' => $barBlock['names'],
                'correlated' => $spec['bar']['correlated'],
                'z' => $z,
                'groups' => $groups,
            ];
        }
        $extraKept = array_map(fn ($e) => self::subset($e, $retained), $extra);
        foreach ($extraKept as $e) {
            self::checkFinite($e, 'weights');
        }

        $design = new self;
        $design->formula = Parser::deparse($parsed);
        $design->responseName = $fixed->response === null ? null : $fixed->variables[$fixed->response][0];
        $design->response = $response;
        $design->x = $x;
        $design->columns = $block['names'];
        $design->assign = $block['assign'];
        $design->intercept = $fixed->intercept;
        $design->terms = $block['labels'];
        $design->offset = $offset;
        $design->random = $random;
        $design->extra = $extraKept;
        $design->nInput = $n;
        $design->retainedRows = $retained;
        $design->droppedRows = $dropped;
        $design->predvars = ['levels' => $levels, 'poly' => Evaluator::recordedPoly()];

        return $design;
    }
}
