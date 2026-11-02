<?php

declare(strict_types=1);

namespace Phit\Polyfill\Formula;

/**
 * R terms.formula semantics, following stats-core formula::terms: variables
 * in order of first appearance, terms as sorted sets of variable indices
 * ordered by interaction order, and R's TermCode contrast/dummy heuristic.
 *
 * @internal
 */
final class Terms
{
    /** @var list<array{0: string, 1: array}> variable keys (deparsed) and expressions */
    public array $variables = [];

    public ?int $response = null;

    /** @var list<list<int>> */
    public array $terms = [];

    public bool $intercept = true;

    /** @var list<int> */
    public array $offsets = [];

    /** @var list<string> */
    private array $frameNames = [];

    /** @var list<string> */
    private array $excluded = [];

    /** @param list<int> $term */
    public function label(array $term): string
    {
        return implode(':', array_map(fn (int $v) => $this->variables[$v][0], $term));
    }

    /** 1 = code by contrasts, 2 = code by indicator (dummy) variables (R's TermCode). */
    public function coding(int $index, int $variable): int
    {
        $reduced = array_values(array_filter($this->terms[$index], fn (int $v) => $v !== $variable));
        if ($reduced === []) {
            return 1;
        }
        for ($t = 0; $t < $index; $t++) {
            $all = true;
            foreach ($reduced as $v) {
                if (! in_array($v, $this->terms[$t], true)) {
                    $all = false;
                    break;
                }
            }
            if ($all) {
                return 1;
            }
        }

        return 2;
    }

    /** @param list<int> $a  @param list<int> $b  @return list<int> */
    private static function union(array $a, array $b): array
    {
        $t = array_values(array_unique([...$a, ...$b]));
        sort($t);

        return $t;
    }

    /** Drop empty and repeated terms, keeping first appearance. */
    private static function trim(array $terms): array
    {
        $out = [];
        foreach ($terms as $t) {
            if ($t !== [] && ! in_array($t, $out, true)) {
                $out[] = $t;
            }
        }

        return $out;
    }

    private static function interact(array $left, array $right): array
    {
        $out = [];
        foreach ($left as $l) {
            foreach ($right as $r) {
                $out[] = self::union($l, $r);
            }
        }

        return self::trim($out);
    }

    private function install(array $expr): int
    {
        $key = Parser::deparse($expr);
        foreach ($this->variables as $i => [$k]) {
            if ($k === $key) {
                return $i;
            }
        }
        $this->variables[] = [$key, $expr];

        return count($this->variables) - 1;
    }

    private static function wrapError(array $expr): \Exception
    {
        return Parser::error("invalid model formula: wrap '".Parser::deparse($expr)."' in I(...)");
    }

    /** @return list<list<int>> */
    private function encode(array $expr, bool $negate): array
    {
        switch ($expr[0]) {
            case 'num':
                if ($expr[1] == 1.0) {
                    $this->intercept = ! $negate;

                    return [];
                }
                if ($expr[1] == 0.0) {
                    $this->intercept = $negate;

                    return [];
                }
                throw Parser::error('invalid model formula: numeric terms other than 0 or 1');
            case 'paren':
                return $this->encode($expr[1], $negate);
            case 'unary':
                if ($expr[1] === '+') {
                    return $this->encode($expr[2], $negate);
                }
                if ($expr[1] === '-') {
                    // A leading -term only removes; nothing precedes it.
                    $this->encode($expr[2], ! $negate);

                    return [];
                }
                throw self::wrapError($expr);
            case 'dot':
                if ($this->frameNames === []) {
                    throw Parser::error("'.' in formula requires data columns");
                }
                $out = [];
                foreach ($this->frameNames as $name) {
                    if (! in_array($name, $this->excluded, true)) {
                        $out[] = [$this->install(['name', $name])];
                    }
                }

                return $out;
            case 'call':
                if ($expr[1] === 'offset') {
                    $v = $this->install($expr);
                    if (! in_array($v, $this->offsets, true)) {
                        $this->offsets[] = $v;
                    }

                    return [];
                }

                return [[$this->install($expr)]];
            case 'name':
            case 'bool':
            case 'str':
                return [[$this->install($expr)]];
        }
        // binary
        [, $op, $l, $r] = $expr;
        switch ($op) {
            case '+':
                $left = $this->encode($l, $negate);

                return self::trim([...$left, ...$this->encode($r, $negate)]);
            case '-':
                $left = $this->encode($l, $negate);
                $right = $this->encode($r, ! $negate);

                return array_values(array_filter($left, fn ($t) => ! in_array($t, $right, true)));
            case ':':
                $left = $this->encode($l, $negate);
                $right = $this->encode($r, $negate);

                return self::interact($left, $right);
            case '*':
                $left = $this->encode($l, $negate);
                $right = $this->encode($r, $negate);
                $both = self::interact($left, $right);

                return self::trim([...$left, ...$right, ...$both]);
            case '/':
                $left = $this->encode($l, $negate);
                $right = $this->encode($r, $negate);
                $all = [];
                foreach ($left as $t) {
                    $all = self::union($all, $t);
                }
                $nested = array_map(fn ($t) => self::union($all, $t), $right);

                return self::trim([...$left, ...$nested]);
            case '%in%':
                $left = $this->encode($l, $negate);
                $right = $this->encode($r, $negate);
                $all = [];
                foreach ($right as $t) {
                    $all = self::union($all, $t);
                }

                return self::trim(array_map(fn ($t) => self::union($t, $all), $left));
            case '^':
                if ($r[0] !== 'num' || is_nan($r[1]) || $r[1] < 1.0 || fmod($r[1], 1.0) != 0.0) {
                    throw Parser::error('invalid power in formula: use a positive integer, e.g. (a + b)^2');
                }
                $base = $this->encode($l, $negate);
                $result = $base;
                $acc = $base;
                // Interactions cannot exceed the number of distinct variables, so
                // higher powers add nothing. Capping as the extension does (phit
                // 91ee538) also avoids casting a huge power to int.
                $variables = [];
                foreach ($base as $term) {
                    $variables = self::union($variables, $term);
                }
                $rounds = min($r[1], (float) max(count($variables), 1));
                for ($k = 1.0; $k < $rounds; $k += 1.0) {
                    $acc = self::interact($acc, $base);
                    $result = self::trim([...$result, ...$acc]);
                }

                return $result;
            case '|':
            case '||':
                throw Parser::error('random-effect bars must be parenthesized, additive terms such as (1 | g)');
            default:
                throw self::wrapError($expr);
        }
    }

    /**
     * Encode a right-hand side (bars already removed) plus an optional response.
     *
     * @param  list<array{0: bool, 1: array}>  $parts
     * @param  list<string>  $frameNames
     */
    public static function build(?array $response, array $parts, array $frameNames): self
    {
        $terms = new self;
        $terms->frameNames = $frameNames;
        if ($response !== null) {
            self::collectNames($response, $terms->excluded);
            $terms->response = $terms->install($response);
        }
        $list = [];
        foreach ($parts as [$negate, $part]) {
            $encoded = $terms->encode($part, $negate);
            if ($negate) {
                $list = array_values(array_filter($list, fn ($t) => ! in_array($t, $encoded, true)));
            } else {
                $list = self::trim([...$list, ...$encoded]);
            }
        }
        if ($terms->response !== null) {
            foreach ($list as $t) {
                if (in_array($terms->response, $t, true)) {
                    throw Parser::error('the response cannot also appear on the right-hand side of the formula');
                }
            }
        }
        // R's default keep.order = FALSE: stable sort by interaction order.
        usort($list, fn ($a, $b) => count($a) <=> count($b));
        $terms->terms = $list;

        return $terms;
    }

    private static function collectNames(array $expr, array &$out): void
    {
        switch ($expr[0]) {
            case 'name':
                $out[] = $expr[1];
                break;
            case 'call':
                foreach ($expr[2] as [, $value]) {
                    self::collectNames($value, $out);
                }
                break;
            case 'unary':
                self::collectNames($expr[2], $out);
                break;
            case 'paren':
                self::collectNames($expr[1], $out);
                break;
            case 'binary':
                self::collectNames($expr[2], $out);
                self::collectNames($expr[3], $out);
                break;
        }
    }

    /**
     * Split the right-hand side into fixed additive parts and random-effect bars.
     *
     * @param  list<array{0: bool, 1: array}>  $fixed
     * @param  list<array{lhs: array, groups: list<array>, correlated: bool}>  $bars
     */
    public static function splitBars(array $expr, array &$fixed, array &$bars, bool $negate): void
    {
        if ($expr[0] === 'binary' && $expr[1] === '+') {
            self::splitBars($expr[2], $fixed, $bars, $negate);
            self::splitBars($expr[3], $fixed, $bars, $negate);

            return;
        }
        if ($expr[0] === 'binary' && $expr[1] === '-') {
            self::splitBars($expr[2], $fixed, $bars, $negate);
            self::splitBars($expr[3], $fixed, $bars, ! $negate);

            return;
        }
        if ($expr[0] === 'paren' && self::bar($expr[1], $bars)) {
            return;
        }
        if ($expr[0] === 'binary' && ($expr[1] === '|' || $expr[1] === '||') && self::bar($expr, $bars)) {
            return;
        }
        $fixed[] = [$negate, $expr];
    }

    private static function bar(array $inner, array &$bars): bool
    {
        if ($inner[0] === 'binary' && ($inner[1] === '|' || $inner[1] === '||')) {
            foreach (self::groupFactors($inner[3]) as $groups) {
                $bars[] = ['lhs' => $inner[2], 'groups' => $groups, 'correlated' => $inner[1] === '|'];
            }

            return true;
        }

        return false;
    }

    /** @return list<list<array>> */
    private static function groupFactors(array $expr): array
    {
        if ($expr[0] === 'paren') {
            return self::groupFactors($expr[1]);
        }
        if ($expr[0] === 'binary' && $expr[1] === ':') {
            $left = self::groupFactors($expr[2]);
            $right = self::groupFactors($expr[3]);
            $out = [];
            foreach ($left as $a) {
                foreach ($right as $b) {
                    $out[] = [...$a, ...$b];
                }
            }

            return $out;
        }
        if ($expr[0] === 'binary' && $expr[1] === '/') {
            // lme4: (x | a/b) = (x | a) + (x | b:a), naming the inner factor first.
            $left = self::groupFactors($expr[2]);
            $right = self::groupFactors($expr[3]);
            $all = array_merge(...$left);
            $out = $left;
            foreach ($right as $b) {
                $out[] = [...$b, ...$all];
            }

            return $out;
        }
        if ($expr[0] === 'name' || $expr[0] === 'call') {
            return [[$expr]];
        }
        throw Parser::error("unsupported grouping expression '".Parser::deparse($expr)."' in random-effect term");
    }

    public static function barGroupLabel(array $bar): string
    {
        return implode(':', array_map(Parser::deparse(...), $bar['groups']));
    }

    public static function barLabel(array $bar): string
    {
        return Parser::deparse($bar['lhs']).' '.($bar['correlated'] ? '|' : '||').' '.self::barGroupLabel($bar);
    }
}
