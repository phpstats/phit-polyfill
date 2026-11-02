<?php

declare(strict_types=1);

namespace Phit\Polyfill\Formula;

use Phit\Polyfill\ErrorCode;

/**
 * Tokenizer, Pratt parser, and R-compatible deparser for model formulas,
 * following the extension's stats-core formula::parse.
 *
 * Operator precedence follows R's ?Syntax: ^ > unary -/+ > : > %in% > * / >
 * binary +/- > comparisons > ! > & > | > ~.
 *
 * Expressions are arrays:
 *   ['num', float], ['bool', bool], ['name', string], ['str', string], ['dot'],
 *   ['call', string, list<array{0: ?string, 1: array}>], ['unary', string, array],
 *   ['binary', string, array, array], ['paren', array].
 *
 * @internal
 */
final class Parser
{
    /** Binary operators: [left binding power, right binding power]. */
    private const array BINDING = [
        '~' => [2, 3],
        '|' => [4, 5], '||' => [4, 5],
        '&' => [6, 7],
        '<' => [10, 11], '>' => [10, 11], '<=' => [10, 11], '>=' => [10, 11], '==' => [10, 11], '!=' => [10, 11],
        '+' => [12, 13], '-' => [12, 13],
        '*' => [14, 15], '/' => [14, 15],
        '%in%' => [16, 17],
        ':' => [18, 19],
        '^' => [23, 22],
    ];

    private const array TWO_CHAR = [
        '||' => '||', '&&' => '&', '<=' => '<=', '>=' => '>=', '==' => '==', '!=' => '!=', '**' => '^',
    ];

    private const array RESERVED = [
        'if', 'else', 'repeat', 'while', 'function', 'for', 'next', 'break', 'TRUE', 'FALSE', 'NULL',
        'Inf', 'NaN', 'NA', 'NA_integer_', 'NA_real_', 'NA_character_', 'in', 'T',
    ];

    /** @var list<array> */
    private array $tokens;

    private int $pos = 0;

    public static function error(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    public static function parse(string $text): array
    {
        $parser = new self;
        $parser->tokens = self::tokenize($text);
        if ($parser->tokens === []) {
            throw self::error('formula is empty');
        }
        $expr = $parser->expr(0);
        if ($parser->pos < count($parser->tokens)) {
            throw self::error('unexpected trailing input in formula');
        }

        return $expr;
    }

    private static function isWhitespace(string $c): bool
    {
        return preg_match('/^\p{White_Space}$/u', $c) === 1;
    }

    private static function isAlphabetic(string $c): bool
    {
        return preg_match('/^\p{Alphabetic}$/u', $c) === 1;
    }

    private static function isAlphanumeric(string $c): bool
    {
        return preg_match('/^[\p{Alphabetic}\p{N}]$/u', $c) === 1;
    }

    private static function isDigit(?string $c): bool
    {
        return $c !== null && strlen($c) === 1 && $c >= '0' && $c <= '9';
    }

    /** @return list<array> */
    private static function tokenize(string $text): array
    {
        $chars = preg_split('//u', $text, -1, PREG_SPLIT_NO_EMPTY) ?: [];
        $count = count($chars);
        $out = [];
        $i = 0;
        while ($i < $count) {
            $c = $chars[$i];
            if (self::isWhitespace($c)) {
                $i++;

                continue;
            }
            $next = $chars[$i + 1] ?? null;
            if (self::isDigit($c) || ($c === '.' && self::isDigit($next))) {
                $start = $i;
                while ($i < $count && (self::isDigit($chars[$i]) || $chars[$i] === '.')) {
                    $i++;
                }
                if ($i < $count && ($chars[$i] === 'e' || $chars[$i] === 'E')) {
                    $j = $i + 1;
                    if ($j < $count && ($chars[$j] === '+' || $chars[$j] === '-')) {
                        $j++;
                    }
                    if ($j < $count && self::isDigit($chars[$j])) {
                        $i = $j;
                        while ($i < $count && self::isDigit($chars[$i])) {
                            $i++;
                        }
                    }
                }
                $literal = implode('', array_slice($chars, $start, $i - $start));
                // Rust's f64 parser: digits with at most one '.', optional exponent.
                if (preg_match('/^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/', $literal) !== 1) {
                    throw self::error("invalid number '$literal' in formula");
                }
                $value = (float) $literal;
                // R integer suffix: 2L.
                if ($i < $count && $chars[$i] === 'L') {
                    $i++;
                }
                $out[] = ['num', $value];

                continue;
            }
            if (self::isAlphabetic($c) || $c === '.' || $c === '_') {
                $start = $i;
                while ($i < $count && (self::isAlphanumeric($chars[$i]) || $chars[$i] === '.' || $chars[$i] === '_')) {
                    $i++;
                }
                $out[] = ['name', implode('', array_slice($chars, $start, $i - $start))];

                continue;
            }
            if ($c === '"' || $c === "'") {
                $value = '';
                $i++;
                while (true) {
                    $ch = $chars[$i] ?? null;
                    if ($ch === null) {
                        throw self::error('unterminated string in formula');
                    }
                    if ($ch === $c) {
                        break;
                    }
                    if ($ch === '\\') {
                        $escaped = $chars[$i + 1] ?? throw self::error('unterminated string in formula');
                        $value .= $escaped;
                        $i += 2;
                    } else {
                        $value .= $ch;
                        $i++;
                    }
                }
                $i++;
                $out[] = ['str', $value];

                continue;
            }
            if ($c === '`') {
                $start = $i + 1;
                $i++;
                while ($i < $count && $chars[$i] !== '`') {
                    $i++;
                }
                if ($i >= $count) {
                    throw self::error('unterminated backtick name in formula');
                }
                $name = implode('', array_slice($chars, $start, $i - $start));
                if ($name === '') {
                    throw self::error('empty backtick name in formula');
                }
                $out[] = ['name', $name];
                $i++;

                continue;
            }
            $two = $c.($next ?? '');
            if (isset(self::TWO_CHAR[$two])) {
                $out[] = ['op', self::TWO_CHAR[$two]];
                $i += 2;

                continue;
            }
            $width = 1;
            $token = match ($c) {
                '~', '|', '&', '<', '>', '+', '-', '*', '/', ':', '^' => ['op', $c],
                '!' => ['not'],
                '(' => ['lparen'],
                ')' => ['rparen'],
                ',' => ['comma'],
                '=' => ['assign'],
                '%' => null,
                default => throw self::error("unexpected character '$c' in formula"),
            };
            if ($token === null) {
                if (implode('', array_slice($chars, $i, 4)) === '%in%') {
                    $token = ['op', '%in%'];
                    $width = 4;
                } else {
                    throw self::error('only the %in% operator is supported in formulas');
                }
            }
            $out[] = $token;
            $i += $width;
        }

        return $out;
    }

    private function peek(int $offset = 0): ?array
    {
        return $this->tokens[$this->pos + $offset] ?? null;
    }

    private function next(): ?array
    {
        $token = $this->tokens[$this->pos] ?? null;
        $this->pos++;

        return $token;
    }

    private function prefix(): array
    {
        $token = $this->next();
        if ($token === null) {
            throw self::error('formula ended unexpectedly');
        }
        switch ($token[0]) {
            case 'num':
                return ['num', $token[1]];
            case 'str':
                return ['str', $token[1]];
            case 'name':
                $name = $token[1];
                if (($this->peek()[0] ?? null) === 'lparen') {
                    $this->pos++;
                    $args = [];
                    if (($this->peek()[0] ?? null) === 'rparen') {
                        $this->pos++;

                        return ['call', $name, $args];
                    }
                    while (true) {
                        $argName = null;
                        if (($this->peek()[0] ?? null) === 'name' && ($this->peek(1)[0] ?? null) === 'assign') {
                            $argName = $this->next()[1];
                            $this->pos++;
                        }
                        $value = $this->expr(0);
                        $args[] = [$argName, $value];
                        $sep = $this->next()[0] ?? null;
                        if ($sep === 'comma') {
                            continue;
                        }
                        if ($sep === 'rparen') {
                            break;
                        }
                        throw self::error("expected ',' or ')' in function call");
                    }

                    return ['call', $name, $args];
                }

                return match ($name) {
                    '.' => ['dot'],
                    'TRUE', 'T' => ['bool', true],
                    'FALSE', 'F' => ['bool', false],
                    default => ['name', $name],
                };
            case 'lparen':
                $inner = $this->expr(0);
                if (($this->next()[0] ?? null) !== 'rparen') {
                    throw self::error("expected ')' in formula");
                }

                return ['paren', $inner];
            case 'op':
                return match ($token[1]) {
                    '-' => ['unary', '-', $this->expr(20)],
                    '+' => ['unary', '+', $this->expr(20)],
                    '~' => ['unary', '~', $this->expr(3)],
                    default => throw self::error('unexpected token in formula'),
                };
            case 'not':
                return ['unary', '!', $this->expr(8)];
            default:
                throw self::error('unexpected token in formula');
        }
    }

    private function expr(int $min): array
    {
        $left = $this->prefix();
        while (($token = $this->peek()) !== null && $token[0] === 'op') {
            $op = $token[1];
            [$lbp, $rbp] = self::BINDING[$op];
            if ($lbp < $min) {
                break;
            }
            $this->pos++;
            $right = $this->expr($rbp);
            $left = ['binary', $op, $left, $right];
        }

        return $left;
    }

    /** R's syntactic-name rule; other names are deparsed in backticks. */
    public static function syntactic(string $name): bool
    {
        $chars = preg_split('//u', $name, -1, PREG_SPLIT_NO_EMPTY) ?: [];
        if ($chars === []) {
            return false;
        }
        $first = $chars[0];
        $second = $chars[1] ?? null;
        $starts = self::isAlphabetic($first) || ($first === '.' && ! self::isDigit($second));
        if (! $starts) {
            return false;
        }
        foreach ($chars as $c) {
            if (! (self::isAlphanumeric($c) || $c === '.' || $c === '_')) {
                return false;
            }
        }

        return ! (in_array($name, self::RESERVED, true) && $name !== 'T') && $name !== '...';
    }

    public static function deparseName(string $name): string
    {
        return self::syntactic($name) ? $name : "`$name`";
    }

    /** Rust's `{:.Ne}` scientific formatting: mantissa digits and exponent. */
    private static function scientific(float $value, int $precision): array
    {
        [$mantissa, $exponent] = explode('e', sprintf('%.'.$precision.'e', $value));

        return [$mantissa, (int) $exponent];
    }

    /** Exact decimal digits of an integral double, as Rust prints it with {:.0}. */
    private static function integerDigits(float $value): string
    {
        return sprintf('%.0F', $value);
    }

    /** R's 15-significant-digit number formatting (deparse, as.character, factor labels). */
    public static function formatNumber(float $value): string
    {
        if (is_nan($value)) {
            return 'NaN';
        }
        if (is_infinite($value)) {
            return $value > 0 ? 'Inf' : '-Inf';
        }
        if ($value == 0.0) {
            return '0';
        }
        $rounded = (float) sprintf('%.14e', $value);
        $digits = 15;
        for ($d = 1; $d <= 15; $d++) {
            if ((float) sprintf('%.'.($d - 1).'e', $value) === $rounded) {
                $digits = $d;
                break;
            }
        }
        [$mantissa, $exponent] = self::scientific($value, $digits - 1);
        $sci = $mantissa.'e'.($exponent < 0 ? '-' : '+').str_pad((string) abs($exponent), 2, '0', STR_PAD_LEFT);
        $decimals = max($digits - 1 - $exponent, 0);
        if ($decimals > 40) {
            // The fixed form is longer than any scientific form here.
            return $sci;
        }
        $fixed = $decimals === 0 ? self::integerDigits($value) : sprintf('%.'.$decimals.'F', $value);

        return strlen($fixed) <= strlen($sci) ? $fixed : $sci;
    }

    public static function deparse(array $expr): string
    {
        switch ($expr[0]) {
            case 'num':
                return self::formatNumber($expr[1]);
            case 'bool':
                return $expr[1] ? 'TRUE' : 'FALSE';
            case 'name':
                return self::deparseName($expr[1]);
            case 'str':
                return '"'.str_replace(['\\', '"'], ['\\\\', '\\"'], $expr[1]).'"';
            case 'dot':
                return '.';
            case 'call':
                $args = [];
                foreach ($expr[2] as [$name, $value]) {
                    $args[] = $name !== null ? self::deparseName($name).' = '.self::deparse($value) : self::deparse($value);
                }

                return self::deparseName($expr[1]).'('.implode(', ', $args).')';
            case 'unary':
                return $expr[1].self::deparse($expr[2]);
            case 'paren':
                return '('.self::deparse($expr[1]).')';
            default: // binary
                $op = $expr[1];
                if ($op === '^' || $op === ':' || $op === '/') {
                    return self::deparse($expr[2]).$op.self::deparse($expr[3]);
                }

                return self::deparse($expr[2]).' '.$op.' '.self::deparse($expr[3]);
        }
    }
}
