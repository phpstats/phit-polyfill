<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

use Phit\Polyfill\ErrorCode;

/**
 * The number formatting and table layout of phit's summaries: R's printed number format, p-values, and coefficient
 * tables, and the text tables of summaries and printed results.
 *
 * Written from phit's clean-room specification specs/cleanroom/format.md and the Rust that implements it
 * (stats-core's rformat.rs and summary.rs), operation for operation, so that the output matches the extension byte
 * for byte.
 */
final class Format
{
    /** Significant digits of a coefficient table: R's max(3, getOption("digits") - 3) = 4. */
    public const int DIGITS = 4;

    /** Machine epsilon used by formatPval(). */
    private const float PVAL_EPS = 2.220446049250313e-16;

    /** The largest precision PHP's sprintf() honors; beyond it, decimals come from the exact binary value. */
    private const int SPRINTF_PRECISION = 53;

    /** Unicode White_Space, the characters Rust's str::trim() removes. */
    private const string WHITE_SPACE = '\x{0009}-\x{000D}\x{0020}\x{0085}\x{00A0}\x{1680}\x{2000}-\x{200A}'
        .'\x{2028}\x{2029}\x{202F}\x{205F}\x{3000}';

    /** The digits of one base-10^9 limb. */
    private const int LIMB = 1_000_000_000;

    /**
     * The heading: the title, then the family and link, then the formula, each line ending in a newline.
     *
     * @param  ?array{string, string}  $family
     */
    public static function heading(string $title, ?array $family, ?string $formula): string
    {
        $out = "$title\n";
        if ($family !== null) {
            $out .= "Family: {$family[0]} ({$family[1]})\n";
        }
        if ($formula !== null) {
            $out .= "Formula: $formula\n";
        }

        return $out;
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Decimal conversion: Rust's {:.N} and {:.Ne}.

    /** Whether the sign bit of x is set, so that -0.0 counts as negative. */
    private static function signBit(float $x): bool
    {
        return (unpack('J', pack('E', $x))[1] >> 63) !== 0;
    }

    /**
     * The exact decimal value of a positive finite double: [digits, k] with value = digits × 10^-k.
     *
     * @return array{string, int}
     */
    private static function exactDecimal(float $a): array
    {
        $bits = unpack('J', pack('E', $a))[1];
        $field = ($bits >> 52) & 0x7FF;
        $mantissa = $bits & 0xFFFFFFFFFFFFF;
        if ($field === 0) {
            $e2 = -1074;
        } else {
            $mantissa |= 1 << 52;
            $e2 = $field - 1075;
        }
        $limbs = [$mantissa % self::LIMB, intdiv($mantissa, self::LIMB) % self::LIMB, intdiv($mantissa, self::LIMB * self::LIMB)];
        $k = 0;
        if ($e2 >= 0) {
            for ($left = $e2; $left > 0; $left -= 30) {
                $limbs = self::multiplyLimbs($limbs, 1 << min(30, $left));
            }
        } else {
            $k = -$e2;
            for ($left = $k; $left > 0; $left -= 13) {
                $limbs = self::multiplyLimbs($limbs, 5 ** min(13, $left));
            }
        }
        $top = count($limbs) - 1;
        while ($top > 0 && $limbs[$top] === 0) {
            $top--;
        }
        $digits = (string) $limbs[$top];
        for ($i = $top - 1; $i >= 0; $i--) {
            $digits .= str_pad((string) $limbs[$i], 9, '0', STR_PAD_LEFT);
        }

        return [$digits, $k];
    }

    /**
     * Multiply a little-endian base-10^9 number by a factor below 2^31.
     *
     * @param  list<int>  $limbs
     * @return list<int>
     */
    private static function multiplyLimbs(array $limbs, int $factor): array
    {
        $carry = 0;
        foreach ($limbs as $i => $limb) {
            $product = $limb * $factor + $carry;
            $limbs[$i] = $product % self::LIMB;
            $carry = intdiv($product, self::LIMB);
        }
        while ($carry > 0) {
            $limbs[] = $carry % self::LIMB;
            $carry = intdiv($carry, self::LIMB);
        }

        return $limbs;
    }

    /**
     * The first `keep` digits of a digit string, rounded to nearest on the digits dropped, ties to even. A `keep` of
     * zero or less keeps nothing, leaving "0" or "1".
     */
    private static function roundDigits(string $digits, int $keep): string
    {
        if ($keep >= strlen($digits)) {
            return $digits;
        }
        $kept = $keep > 0 ? substr($digits, 0, $keep) : '0';
        $rest = $keep >= 0 ? substr($digits, $keep) : str_repeat('0', -$keep).$digits;
        $first = $rest[0];
        if ($first > '5') {
            $up = true;
        } elseif ($first < '5') {
            $up = false;
        } elseif (rtrim(substr($rest, 1), '0') !== '') {
            $up = true;
        } else {
            $up = ((int) $kept[strlen($kept) - 1]) % 2 === 1;
        }

        return $up ? self::stepDigits($kept, true) : $kept;
    }

    /** Rust's format!("{:.n}", |x|) for finite x, from the exact binary value. */
    private static function exactFixed(float $a, int $decimals): string
    {
        if ($a == 0.0) {
            $scaled = '0';
        } else {
            [$digits, $k] = self::exactDecimal($a);
            $scaled = $decimals >= $k
                ? $digits.str_repeat('0', $decimals - $k)
                : self::roundDigits($digits, strlen($digits) - ($k - $decimals));
        }
        if (strlen($scaled) < $decimals + 1) {
            $scaled = str_repeat('0', $decimals + 1 - strlen($scaled)).$scaled;
        }
        $integer = substr($scaled, 0, strlen($scaled) - $decimals);

        return $decimals > 0 ? $integer.'.'.substr($scaled, -$decimals) : $integer;
    }

    /**
     * Rust's format!("{:.ne}", |x|) for finite x, from the exact binary value: the mantissa and the exponent.
     *
     * @return array{string, int}
     */
    private static function exactScientific(float $a, int $decimals): array
    {
        if ($a == 0.0) {
            $mantissa = str_repeat('0', $decimals + 1);
            $exp = 0;
        } else {
            [$digits, $k] = self::exactDecimal($a);
            $exp = strlen($digits) - 1 - $k;
            $want = $decimals + 1;
            $mantissa = strlen($digits) <= $want
                ? $digits.str_repeat('0', $want - strlen($digits))
                : self::roundDigits($digits, $want);
            if (strlen($mantissa) > $want) {
                $exp++;
                $mantissa = substr($mantissa, 0, $want);
            }
        }

        return [$decimals > 0 ? $mantissa[0].'.'.substr($mantissa, 1) : $mantissa, $exp];
    }

    /** Rust's format!("{:.n}", x) for finite x; negative zero keeps its sign, as in Rust. */
    private static function fixed(float $x, int $decimals): string
    {
        if ($decimals <= self::SPRINTF_PRECISION) {
            $text = sprintf("%.{$decimals}F", $x);
        } else {
            $text = ($x < 0.0 ? '-' : '').self::exactFixed(abs($x), $decimals);
        }
        if ($x == 0.0 && self::signBit($x) && $text[0] !== '-') {
            $text = '-'.$text;
        }

        return $text;
    }

    /**
     * Rust's format!("{:.ne}", x) for finite nonzero x, split into mantissa text (with its sign) and exponent.
     *
     * @return array{string, int}
     */
    private static function scientific(float $x, int $decimals): array
    {
        if ($decimals <= self::SPRINTF_PRECISION) {
            $text = sprintf("%.{$decimals}e", $x);
            $pos = strpos($text, 'e');

            return [substr($text, 0, $pos), (int) substr($text, $pos + 1)];
        }
        [$mantissa, $exp] = self::exactScientific(abs($x), $decimals);

        return [($x < 0.0 ? '-' : '').$mantissa, $exp];
    }

    /** Rust's format!("{:.n}", x): fixed notation, `NaN`, `inf`, or `-inf`, and `-0.0` for negative zero. */
    public static function rustFixed(float $x, int $decimals): string
    {
        if (is_nan($x)) {
            return 'NaN';
        }
        if (is_infinite($x)) {
            return $x > 0.0 ? 'inf' : '-inf';
        }

        return self::fixed($x, $decimals);
    }

    /** Rust's format!("{v:.n}"), keeping the sign of negative zero; the same as rustFixed(). */
    public static function rustFixedSigned(float $x, int $decimals): string
    {
        return self::rustFixed($x, $decimals);
    }

    /** Rust's f64::powi(), as compiler-rt's __powidf2 evaluates it: repeated squaring, then a reciprocal. */
    private static function powi(float $a, int $b): float
    {
        $reciprocal = $b < 0;
        $r = 1.0;
        while (true) {
            if (($b & 1) !== 0) {
                $r *= $a;
            }
            $b = intdiv($b, 2);
            if ($b === 0) {
                break;
            }
            $a *= $a;
        }

        return $reciprocal ? fdiv(1.0, $r) : $r;
    }

    // ---------------------------------------------------------------------------------------------------------------
    // rformat.rs

    /** Text for a non-finite value. */
    private static function nonFiniteText(float $x): string
    {
        if (is_nan($x)) {
            return 'NaN';
        }

        return $x > 0.0 ? 'Inf' : '-Inf';
    }

    /** Exponent text in C style: a sign and at least two digits. */
    private static function exponentText(int $exp): string
    {
        return ($exp < 0 ? '-' : '+').sprintf('%02d', abs($exp));
    }

    /**
     * Round |x| to `sig` significant digits: the significant digits left after dropping trailing zeros, and the
     * decimal exponent of the rounded value.
     *
     * @return array{int, int}
     */
    private static function sigAndExp(float $x, int $sig): array
    {
        if ($x == 0.0) {
            return [1, 0];
        }
        [$mantissa, $exp] = self::scientific(abs($x), max(1, $sig) - 1);
        $digits = str_replace(['.', '-'], '', $mantissa);
        $kept = strlen(rtrim($digits, '0'));

        return [max(1, $kept === 0 ? 1 : $kept), $exp];
    }

    /** Print x in scientific notation with `decimals` mantissa decimals and a C-style exponent. */
    private static function sciText(float $x, int $decimals): string
    {
        $x = $x == 0.0 ? 0.0 : $x;
        if ($x == 0.0) {
            return ($decimals > 0 ? '0.'.str_repeat('0', $decimals) : '0').'e+00';
        }
        [$mantissa, $exp] = self::scientific($x, $decimals);

        return $mantissa.'e'.self::exponentText($exp);
    }

    /** Print x in fixed notation with `decimals` decimals; negative zero prints as zero. */
    private static function fixedText(float $x, int $decimals): string
    {
        return self::fixed($x == 0.0 ? 0.0 : $x, $decimals);
    }

    /**
     * Pad every cell on the left to the widest cell, or to `minWidth` when that is wider.
     *
     * @param  list<string>  $cells
     * @return list<string>
     */
    private static function padCells(array $cells, int $minWidth): array
    {
        $width = $minWidth;
        foreach ($cells as $c) {
            $width = max($width, self::width($c));
        }

        return array_map(fn (string $c): string => self::padLeft($c, $width), $cells);
    }

    /**
     * R's format(xs, digits = digits): every element in one common notation and padded on the left to one width.
     *
     * @param  list<float>  $xs
     * @return list<string>
     */
    public static function formatReal(array $xs, int $digits): array
    {
        $xs = array_values(array_map(fn ($x): float => (float) $x, $xs));
        $digits = max(1, $digits);

        $anyFinite = false;
        $anyNegative = false;
        $maxSig = 1;
        $intWidth = 1;
        $maxDecimals = 0;
        $wideExponent = false;

        foreach ($xs as $x) {
            if (! is_finite($x)) {
                continue;
            }
            $anyFinite = true;
            $negative = $x < 0.0;
            $anyNegative = $anyNegative || $negative;
            [$sig, $exp] = self::sigAndExp($x, $digits);
            $maxSig = max($maxSig, $sig);
            if ($exp >= 100 || $exp <= -99) {
                $wideExponent = true;
            }

            // Integer digits in fixed notation.
            $intDigits = $exp + 1;
            $probeDecimals = max(0, $digits - 1 - $exp);
            $printed = (float) self::fixedText(abs($x), $probeDecimals);
            if ($printed < self::powi(10.0, $exp)) {
                $intDigits--;
            }
            $width = max(1, $intDigits) + ($negative ? 1 : 0);
            $intWidth = max($intWidth, $width);

            $maxDecimals = max($maxDecimals, $sig - 1 - $exp);
        }

        $notationWidth = 0;
        if (! $anyFinite) {
            $cells = array_map(fn (float $x): string => self::nonFiniteText($x), $xs);
        } else {
            $decimals = max(0, $maxDecimals);
            $fixedWidth = $intWidth + $decimals + ($decimals > 0 ? 1 : 0);
            $mantissaDecimals = $maxSig - 1;
            $sciWidth = ($anyNegative ? 1 : 0) + 1 + ($maxSig > 1 ? 1 + $mantissaDecimals : 0) + ($wideExponent ? 5 : 4);
            $useFixed = $fixedWidth <= $sciWidth;
            $notationWidth = $useFixed ? $fixedWidth : $sciWidth;
            $cells = [];
            foreach ($xs as $x) {
                if (! is_finite($x)) {
                    $cells[] = self::nonFiniteText($x);
                } elseif ($useFixed) {
                    $cells[] = self::fixedText($x, $decimals);
                } else {
                    $cells[] = self::sciText($x, $mantissaDecimals);
                }
            }
        }

        return self::padCells($cells, $notationWidth);
    }

    /** Add one to (`up`) or subtract one from a string of decimal digits. */
    private static function stepDigits(string $digits, bool $up): string
    {
        [$edge, $reset] = $up ? ['9', '0'] : ['0', '9'];
        $i = strlen($digits);
        while (true) {
            if ($i === 0) {
                // Only reachable going up past all nines.
                return '1'.$digits;
            }
            $i--;
            if ($digits[$i] === $edge) {
                $digits[$i] = $reset;
            } else {
                $digits[$i] = (string) ((int) $digits[$i] + ($up ? 1 : -1));

                return $digits;
            }
        }
    }

    /** Round x to `decimals` places: the nearer decimal neighbour, ties to the even last digit. */
    public static function round(float $x, int $decimals): float
    {
        if (! is_finite($x)) {
            return $x;
        }
        $magnitude = abs($x);
        $nearestText = self::fixed($magnitude, $decimals);
        $nearest = (float) $nearestText;
        if ($nearest == $magnitude) {
            // |x| is already a d-decimal value.
            return $x;
        }
        // Scaled digits of the nearest d-decimal value, then its neighbour on the other side of |x|.
        $scaled = str_replace('.', '', $nearestText);
        if ($nearest < $magnitude) {
            [$loDigits, $hiDigits] = [$scaled, self::stepDigits($scaled, true)];
        } else {
            [$loDigits, $hiDigits] = [self::stepDigits($scaled, false), $scaled];
        }
        $lo = (float) "{$loDigits}e-$decimals";
        $hi = (float) "{$hiDigits}e-$decimals";
        $below = $magnitude - $lo;
        $above = $hi - $magnitude;
        $loEven = ((int) $loDigits[strlen($loDigits) - 1]) % 2 === 0;
        $chosen = ($below < $above || ($below == $above && $loEven)) ? $lo : $hi;

        return $x < 0.0 ? -$chosen : $chosen;
    }

    /** Remove trailing zeros, and then a trailing point, from a decimal digit string. */
    private static function stripTrailingZeros(string $text): string
    {
        if (! str_contains($text, '.')) {
            return $text;
        }

        return rtrim(rtrim($text, '0'), '.');
    }

    /** R's formatC(x, digits = digits), which is C's printf("%.{digits}g"). */
    public static function formatG(float $x, int $digits): string
    {
        if (! is_finite($x)) {
            return self::nonFiniteText($x);
        }
        if ($x == 0.0) {
            return '0';
        }
        $precision = max(1, $digits);
        [$mantissa, $exp] = self::scientific($x, $precision - 1);
        if ($exp < -4 || $exp >= $precision) {
            return self::stripTrailingZeros($mantissa).'e'.self::exponentText($exp);
        }

        return self::stripTrailingZeros(self::fixed($x, $precision - 1 - $exp));
    }

    /**
     * R's format.pval(pvalues, digits = digits) with eps = 2.220446049250313e-16.
     *
     * @param  list<float>  $pvalues
     * @return list<string>
     */
    public static function formatPval(array $pvalues, int $digits): array
    {
        $pvalues = array_values(array_map(fn ($p): float => (float) $p, $pvalues));
        $out = array_fill(0, count($pvalues), '');

        $fixedIdx = [];
        $sciIdx = [];
        $tinyIdx = [];
        foreach ($pvalues as $i => $p) {
            if (is_nan($p)) {
                $out[$i] = 'NA';
            } elseif ($p >= self::PVAL_EPS) {
                $floor = floor(log10($p));
                $e = is_finite($floor) ? (int) $floor : PHP_INT_MAX;
                if ($e >= -3 || ($e === -4 && $digits > 1)) {
                    $fixedIdx[] = $i;
                } else {
                    $sciIdx[] = $i;
                }
            } else {
                $tinyIdx[] = $i;
            }
        }

        $widest = 0;
        foreach ([$fixedIdx, $sciIdx] as $group) {
            if ($group === []) {
                continue;
            }
            $values = array_map(fn (int $i): float => $pvalues[$i], $group);
            foreach (self::formatReal($values, $digits) as $j => $cell) {
                $widest = max($widest, self::width($cell));
                $out[$group[$j]] = $cell;
            }
        }

        if ($tinyIdx !== []) {
            $k = max(1, $digits - 2);
            if ($fixedIdx === [] && $sciIdx === []) {
                $separator = $k === 1 ? '' : ' ';
            } else {
                if ($k > 1 && $k + 6 > $widest) {
                    $k = max(1, $widest - 7);
                }
                $separator = ($k === 1 && $widest <= 6) ? '' : ' ';
            }
            $cell = '<'.$separator.self::formatReal([self::PVAL_EPS], $k)[0];
            foreach ($tinyIdx as $i) {
                $out[$i] = $cell;
            }
        }

        return $out;
    }

    /**
     * Format a column with formatReal() at DIGITS, trimming cells and printing NaN as NA.
     *
     * @param  list<float>  $values
     * @return list<string>
     */
    private static function columnCells(array $values): array
    {
        $cells = [];
        foreach (self::formatReal($values, self::DIGITS) as $i => $cell) {
            $cells[] = is_nan($values[$i]) ? 'NA' : self::trim($cell);
        }

        return $cells;
    }

    /**
     * The cells R's printCoefmat(m, digits = 4, signif.stars = FALSE) prints for a table whose columns are an
     * estimate, its standard error, a test statistic, and a p-value. Cells are trimmed of padding.
     *
     * @param  list<float>  $estimates
     * @param  list<float>  $stdErrors
     * @param  list<float>  $statistics
     * @param  list<float>  $pValues
     * @return array{estimates: list<string>, std_errors: list<string>, statistics: list<string>, p_values: list<string>}
     */
    public static function coefficientCells(array $estimates, array $stdErrors, array $statistics, array $pValues): array
    {
        $estimates = array_values(array_map(fn ($v): float => (float) $v, $estimates));
        $stdErrors = array_values(array_map(fn ($v): float => (float) $v, $stdErrors));
        $statistics = array_values(array_map(fn ($v): float => (float) $v, $statistics));
        $joint = array_merge($estimates, $stdErrors);

        $smallest = INF;
        foreach ($joint as $v) {
            if (is_finite($v) && $v != 0.0) {
                $smallest = min($smallest, abs($v));
            }
        }
        $g = is_finite($smallest) ? 1 + (int) floor(log10($smallest)) : 1;
        $decimals = max(1, self::DIGITS - $g);

        $estimateCells = self::columnCells(array_map(fn (float $v): float => self::round($v, $decimals), $joint));
        $stdErrorCells = array_slice($estimateCells, count($estimates));
        $estimateCells = array_slice($estimateCells, 0, count($estimates));

        return [
            'estimates' => $estimateCells,
            'std_errors' => $stdErrorCells,
            'statistics' => self::columnCells(array_map(fn (float $v): float => self::round($v, 3), $statistics)),
            'p_values' => array_map(fn (string $c): string => self::trim($c), self::formatPval($pValues, 3)),
        ];
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Text tables: summary.rs.

    /** The number of characters in a UTF-8 string, as Rust's chars().count(). */
    private static function width(string $s): int
    {
        return strlen($s) - preg_match_all('/[\x80-\xBF]/', $s);
    }

    /** Rust's str::trim(): Unicode white space removed from both ends. */
    private static function trim(string $s): string
    {
        if (! preg_match('/[\x80-\xFF]/', $s)) {
            return trim($s, " \t\n\r\v\f");
        }

        return preg_replace('/^['.self::WHITE_SPACE.']+|['.self::WHITE_SPACE.']+$/u', '', $s) ?? $s;
    }

    /** Rust's str::trim_end(). */
    private static function trimEnd(string $s): string
    {
        if (! preg_match('/[\x80-\xFF]/', $s)) {
            return rtrim($s, " \t\n\r\v\f");
        }

        return preg_replace('/['.self::WHITE_SPACE.']+$/u', '', $s) ?? $s;
    }

    /** Right-align: pad on the left to `width` characters. */
    public static function padLeft(string $s, int $width): string
    {
        return str_repeat(' ', max(0, $width - self::width($s))).$s;
    }

    /** Left-align: pad on the right to `width` characters. */
    private static function padRight(string $s, int $width): string
    {
        return $s.str_repeat(' ', max(0, $width - self::width($s)));
    }

    /**
     * A text table: the first column left-aligned, the rest right-aligned unless listed in `left`, one space apart.
     *
     * @param  list<string>  $header
     * @param  list<list<string>>  $rows
     * @param  list<int>  $left
     */
    public static function table(array $header, array $rows, array $left): string
    {
        $widths = [];
        foreach ($header as $c => $name) {
            $w = self::width($name);
            foreach ($rows as $r) {
                $w = max($w, self::width($r[$c]));
            }
            $widths[$c] = $w;
        }
        $line = function (array $cells) use ($widths, $left): string {
            $parts = [];
            foreach ($cells as $c => $s) {
                $parts[] = ($c === 0 || in_array($c, $left, true)) ? self::padRight($s, $widths[$c]) : self::padLeft($s, $widths[$c]);
            }

            return self::trimEnd(implode(' ', $parts))."\n";
        };
        $out = $line($header);
        foreach ($rows as $row) {
            $out .= $line($row);
        }

        return $out;
    }

    /** R's printCoefmat(digits = 4) without significance stars. */
    public static function coefficientTable(Coefficients $c): string
    {
        $cells = self::coefficientCells($c->estimates, $c->stdErrors, $c->statistics, $c->pValues);
        $rows = [];
        foreach ($c->names as $i => $name) {
            $rows[] = [
                $name,
                $cells['estimates'][$i],
                $cells['std_errors'][$i],
                $cells['statistics'][$i],
                $cells['p_values'][$i],
            ];
        }
        [$stat, $p] = $c->statistic === 't' ? ['t value', 'Pr(>|t|)'] : ['z value', 'Pr(>|z|)'];

        return self::table(['', 'Estimate', 'Std. Error', $stat, $p], $rows, []);
    }

    /**
     * The random-effects table: each term's variances, standard deviations, and correlations, then the residual.
     *
     * @param  list<RandomTerm>  $terms
     */
    public static function randomEffects(array $terms, ?float $residual): string
    {
        $rows = [];
        $variances = [];
        $correlations = [];
        $width = 1;
        if ($terms !== []) {
            $width = max(array_map(fn (RandomTerm $t): int => $t->correlated ? count($t->columns) : 1, $terms));
        }
        $width = max(0, $width - 1);
        foreach ($terms as $term) {
            foreach (array_values($term->columns) as $j => $name) {
                $rows[] = [$j === 0 ? $term->group : '', $name];
                $variances[] = $term->covariance[$j][$j];
                $corr = array_fill(0, $width, '');
                if ($term->correlated) {
                    for ($m = 0; $m < min($j, $width); $m++) {
                        $r = fdiv($term->covariance[$j][$m], sqrt($term->covariance[$j][$j] * $term->covariance[$m][$m]));
                        $corr[$m] = self::rustFixed($r == 0.0 ? 0.0 : $r, 2);
                    }
                }
                $correlations[] = $corr;
            }
        }
        if ($residual !== null) {
            $rows[] = ['Residual', ''];
            $variances[] = $residual;
            $correlations[] = array_fill(0, $width, '');
        }
        $sds = array_map(fn (float $v): float => sqrt($v), $variances);
        $vt = self::formatReal($variances, self::DIGITS);
        $st = self::formatReal($sds, self::DIGITS);
        foreach ($rows as $i => $row) {
            $rows[$i] = [...$row, self::trim($vt[$i]), self::trim($st[$i]), ...$correlations[$i]];
        }
        $header = ['Groups', 'Name', 'Variance', 'Std.Dev.', ...array_fill(0, $width, '')];
        if ($width > 0) {
            $header[4] = 'Corr';
        }

        return self::table($header, $rows, [1]);
    }

    /** A string the extension reads as text: valid UTF-8, or null as Zval::str() gives for any other string. */
    private static function utf8(string $s): ?string
    {
        return preg_match('//u', $s) === 1 ? $s : null;
    }

    /**
     * A table in the layout of R's print.data.frame(): row numbers left-aligned in an unnamed first column, every
     * other column right-aligned under its name, one space apart. A column of strings and nulls is text; any other
     * column is numbers, formatted per column as format(x, digits = 7). A missing number prints NA and a missing
     * string <NA>. Columns are split into blocks whose lines hold at most 79 characters, each block with its own
     * header and row numbers; a column wider than that is a block of its own.
     *
     * @param  array<string|int, mixed>  $columns
     */
    public static function dataFrame(array $columns): string
    {
        $names = [];
        $texts = [];
        $rows = null;
        foreach ($columns as $key => $column) {
            $names[] = (string) $key;
            if (! is_array($column)) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'table columns must be arrays');
            }
            $rows ??= count($column);
            $isText = true;
            foreach ($column as $v) {
                if ($v !== null && ! is_string($v)) {
                    $isText = false;
                    break;
                }
            }
            if ($isText) {
                $texts[] = array_map(fn (?string $v): string => ($v === null ? null : self::utf8($v)) ?? '<NA>', array_values($column));

                continue;
            }
            $values = array_map(
                fn (mixed $v): ?float => is_float($v) ? $v : (is_int($v) ? (float) $v : null),
                array_values($column),
            );
            $present = array_values(array_filter($values, fn (?float $v): bool => $v !== null));
            $formatted = self::formatReal($present, 7);
            $next = 0;
            $cells = [];
            foreach ($values as $v) {
                $cells[] = $v === null ? 'NA' : ($formatted[$next++] ?? '');
            }
            $texts[] = $cells;
        }
        $rows ??= 0;
        if ($rows === 0) {
            return '<0 rows> '.implode(', ', $names)."\n";
        }
        $labels = array_map('strval', range(1, $rows));
        $labelWidth = strlen($labels[$rows - 1]);
        $widths = [];
        foreach ($names as $c => $name) {
            $w = self::width($name);
            foreach ($texts[$c] as $cell) {
                $w = max($w, self::width($cell));
            }
            $widths[$c] = $w;
        }

        // Blocks of columns, filled greedily while a line stays within 79 characters.
        $blocks = [];
        $start = 0;
        $count = count($names);
        while ($start < $count) {
            $end = $start + 1;
            $line = $labelWidth + 1 + $widths[$start];
            while ($end < $count && $line + 1 + $widths[$end] <= 79) {
                $line += 1 + $widths[$end];
                $end++;
            }
            $blocks[] = [$start, $end];
            $start = $end;
        }

        $out = '';
        foreach ($blocks as [$from, $to]) {
            $out .= str_repeat(' ', $labelWidth);
            for ($c = $from; $c < $to; $c++) {
                $out .= ' '.self::padLeft($names[$c], $widths[$c]);
            }
            $out .= "\n";
            foreach ($labels as $r => $label) {
                $out .= str_pad($label, $labelWidth);
                for ($c = $from; $c < $to; $c++) {
                    $out .= ' '.self::padLeft(self::trim($texts[$c][$r] ?? ''), $widths[$c]);
                }
                $out .= "\n";
            }
        }

        return $out;
    }
}
