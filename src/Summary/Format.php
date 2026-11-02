<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/**
 * The formatting every printer shares, from stats-core's summary.rs: R's
 * formatReal() (through format()), format.pval(), and C's %g (through
 * formatC()), built on correctly rounded decimal conversion, which PHP's
 * sprintf('%.*F') and sprintf('%.*e') give; and the coefficient and
 * random-effects tables.
 */
final class Format
{
    /** Significant digits for coefficient tables, as R's max(3, getOption("digits") - 3). */
    public const int DIGITS = 4;

    /** sprintf() truncates a longer precision to this many digits. */
    private const int MAX_PRECISION = 53;

    /**
     * The heading lines: the title, then the family and the formula when known.
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

    /** sprintf('%.{decimals}F'), with decimals beyond sprintf's limit handled as Rust's {:.*} does. */
    private static function fixed(float $x, int $decimals): string
    {
        if ($decimals <= self::MAX_PRECISION) {
            return sprintf("%.{$decimals}F", $x);
        }
        // Only round() asks for this many decimals; see there.
        throw new \LogicException('fixed() precision out of range');
    }

    /** Decimal digits of |x| rounded to `digits` significant digits: [mantissa digits, decimal exponent]. */
    private static function sigDigits(float $x, int $digits): array
    {
        [$mantissa, $exponent] = explode('e', sprintf('%.'.($digits - 1).'e', abs($x)));

        return [str_replace('.', '', $mantissa), (int) $exponent];
    }

    /** R's scientific(): [negative, kpower, nsig, widens]. */
    private static function scientific(float $x, int $digits): array
    {
        if ($x == 0.0) {
            return [false, 0, 1, false];
        }
        [$mantissa, $kpower] = self::sigDigits($x, $digits);
        $nsig = $digits;
        while ($nsig > 1 && $mantissa[$nsig - 1] === '0') {
            $nsig--;
        }
        // Fixed notation needs one digit fewer when |x|, shown with the decimals it would get there, stays below
        // 10^kpower (99996 at 4 digits prints as 99996, where 99999.6 rounds to 100000 and prints as 1e+05).
        $widens = false;
        if ($kpower > 0 && $kpower <= 22) {
            $shown = (float) self::fixed(abs($x), max(0, $digits - 1 - $kpower));
            $widens = $shown < 10.0 ** $kpower;
        }

        return [$x < 0.0, $kpower, $nsig, $widens];
    }

    /** An exponent in R's style: a sign and at least two digits. */
    private static function exponent(int $e): string
    {
        return 'e'.($e < 0 ? '-' : '+').sprintf('%02d', abs($e));
    }

    private static function scientificText(float $x, int $decimals): string
    {
        [$mantissa, $e] = explode('e', sprintf("%.{$decimals}e", $x));

        return $mantissa.self::exponent((int) $e);
    }

    private static function special(float $x): ?string
    {
        return match (true) {
            is_nan($x) => 'NaN',
            $x === INF => 'Inf',
            $x === -INF => '-Inf',
            default => null,
        };
    }

    /** Rust's {:.d} for any double: NaN and the infinities as Rust prints them. */
    public static function rustFixed(float $x, int $decimals): string
    {
        return match (true) {
            is_nan($x) => 'NaN',
            $x === INF => 'inf',
            $x === -INF => '-inf',
            default => self::fixed($x == 0.0 ? 0.0 : $x, $decimals),
        };
    }

    /**
     * R's format(x, digits = digits) for a vector: one common layout, fixed or scientific, padded to one width.
     *
     * @param  list<float>  $xs
     * @return list<string>
     */
    public static function formatReal(array $xs, int $digits): array
    {
        [$neg, $mxsl, $rgt, $mxns, $mxe, $mne] = [false, 1, 0, 0, PHP_INT_MIN, PHP_INT_MAX];
        $finite = false;
        foreach ($xs as $x) {
            if (! is_finite($x)) {
                continue;
            }
            $finite = true;
            [$negative, $kpower, $nsig, $widens] = self::scientific($x, $digits);
            $left = $kpower + 1 - (int) $widens;
            $sleft = (int) $negative + ($left <= 0 ? 1 : $left);
            $rgt = max($rgt, $nsig - $left);
            $neg = $neg || $negative;
            $mxsl = max($mxsl, $sleft);
            $mxns = max($mxns, $nsig);
            $mxe = max($mxe, $kpower);
            $mne = min($mne, $kpower);
        }
        $rgt = max($rgt, 0);
        $fixedWidth = max($mxsl, 1) + $rgt + (int) ($rgt > 0);
        $scientificLayout = null;
        if ($finite) {
            $e = $mxe >= 100 || $mne <= -99 ? 2 : 1;
            $decimals = $mxns - 1;
            $width = (int) $neg + (int) ($decimals > 0) + $decimals + 4 + $e;
            if ($fixedWidth > $width) {
                $scientificLayout = $decimals;
            }
        }
        $texts = [];
        foreach ($xs as $x) {
            $special = self::special($x);
            if ($special !== null) {
                $texts[] = $special;

                continue;
            }
            $x = $x == 0.0 ? 0.0 : $x;
            $texts[] = $scientificLayout !== null ? self::scientificText($x, $scientificLayout) : self::fixed($x, $rgt);
        }
        $width = max(array_map(self::width(...), $texts) ?: [0]);

        return array_map(fn (string $s): string => self::padLeft($s, $width), $texts);
    }

    /**
     * R's round(x, decimals), as the correctly rounded decimal.
     *
     * Beyond sprintf's 53 decimals: a value with 17 or more significant digits
     * there is unchanged, because its correctly rounded decimal parses back to
     * itself. A smaller value keeps the significant digits it has.
     */
    public static function round(float $x, int $decimals): float
    {
        if (! is_finite($x)) {
            return $x;
        }
        if ($decimals <= self::MAX_PRECISION) {
            return (float) self::fixed($x, $decimals);
        }
        if ($x == 0.0) {
            return $x;
        }
        [, $e] = self::sigDigits($x, 17);
        $significant = $decimals + $e + 1;
        if ($significant >= 17) {
            return $x;
        }
        if ($significant >= 1) {
            [$mantissa, $e] = self::sigDigits($x, $significant);

            return ($x < 0.0 ? -1.0 : 1.0) * (float) ($mantissa.'e'.($e - $significant + 1));
        }
        // Below one unit in the last place: zero, or one unit when |x| passes half of it.
        $half = $significant === 0 && strcmp(substr(sprintf('%.40e', abs($x)), 0, 42), '5.'.str_repeat('0', 40)) > 0;

        return ($x < 0.0 ? -1.0 : 1.0) * ($half ? (float) "1e-$decimals" : 0.0);
    }

    /** C's %.{digits}g, as R's formatC(x, digits = digits). */
    public static function formatG(float $x, int $digits): string
    {
        $special = self::special($x);
        if ($special !== null) {
            return $special;
        }
        if ($x == 0.0) {
            return '0';
        }
        [, $e] = self::sigDigits($x, $digits);
        $trim = fn (string $s): string => str_contains($s, '.') ? rtrim(rtrim($s, '0'), '.') : $s;
        if ($e < -4 || $e >= $digits) {
            [$mantissa] = explode('e', sprintf('%.'.($digits - 1).'e', $x));

            return $trim($mantissa).self::exponent($e);
        }

        return $trim(self::fixed($x, max(0, $digits - 1 - $e)));
    }

    /**
     * R's format.pval(pv, digits = digits) with eps = machine epsilon.
     *
     * @param  list<float>  $pvalues
     * @return list<string>
     */
    public static function formatPval(array $pvalues, int $digits): array
    {
        $eps = PHP_FLOAT_EPSILON;
        $out = array_fill(0, count($pvalues), '');
        $normal = array_keys(array_filter($pvalues, fn (float $p): bool => ! is_nan($p) && $p >= $eps));
        $fixed = array_values(array_filter($normal, function (int $i) use ($pvalues, $digits): bool {
            $expo = floor(log10($pvalues[$i]));

            return $expo >= -3.0 || ($expo == -4.0 && $digits > 1);
        }));
        $other = array_values(array_diff($normal, $fixed));
        foreach ([$fixed, $other] as $group) {
            $texts = self::formatReal(array_map(fn (int $i): float => $pvalues[$i], $group), $digits);
            foreach ($group as $k => $i) {
                $out[$i] = $texts[$k];
            }
        }
        $tiny = array_keys(array_filter($pvalues, fn (float $p): bool => ! is_nan($p) && $p < $eps));
        if ($tiny !== []) {
            $epsDigits = max($digits - 2, 1);
            if ($normal === []) {
                $sep = $epsDigits === 1 ? '' : ' ';
            } else {
                $nc = max(array_map(fn (int $i): int => self::width($out[$i]), $normal));
                if ($epsDigits > 1 && $epsDigits + 6 > $nc) {
                    $epsDigits = max($nc - 7, 1);
                }
                $sep = $epsDigits === 1 && $nc <= 6 ? '' : ' ';
            }
            $text = '<'.$sep.self::formatReal([$eps], $epsDigits)[0];
            foreach ($tiny as $i) {
                $out[$i] = $text;
            }
        }
        foreach ($pvalues as $i => $p) {
            if (is_nan($p)) {
                $out[$i] = 'NA';
            }
        }

        return $out;
    }

    // Tables.

    /** Characters, as Rust counts them. */
    private static function width(string $s): int
    {
        return preg_match_all('/./su', $s);
    }

    public static function padLeft(string $s, int $width): string
    {
        return str_repeat(' ', max(0, $width - self::width($s))).$s;
    }

    private static function padRight(string $s, int $width): string
    {
        return $s.str_repeat(' ', max(0, $width - self::width($s)));
    }

    /**
     * A text table: the first column left-aligned, the rest right-aligned unless listed in $left, one space apart.
     *
     * @param  list<string>  $header
     * @param  list<list<string>>  $rows
     * @param  list<int>  $left
     */
    public static function table(array $header, array $rows, array $left): string
    {
        $widths = [];
        foreach ($header as $c => $name) {
            $widths[$c] = max([self::width($name), ...array_map(fn (array $r): int => self::width($r[$c]), $rows)]);
        }
        $line = function (array $cells) use ($widths, $left): string {
            $parts = [];
            foreach ($cells as $c => $s) {
                $parts[] = $c === 0 || in_array($c, $left, true) ? self::padRight($s, $widths[$c]) : self::padLeft($s, $widths[$c]);
            }

            return rtrim(implode(' ', $parts), " \t\n\r\0\x0B\x0C")."\n";
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
        $digits = self::DIGITS;
        $testDigits = $digits - 1;
        $both = [...$c->estimates, ...$c->stdErrors];
        $sizes = array_filter(array_map(abs(...), $both), fn (float $x): bool => is_finite($x) && $x != 0.0);
        $digmin = $sizes === [] ? 1 : 1 + (int) floor(log10(min($sizes)));
        $decimals = max(1, $digits - $digmin);
        $estimateText = self::formatReal(array_map(fn (float $x): float => self::round($x, $decimals), $both), $digits);
        $p = count($c->estimates);
        $statistics = self::formatReal(array_map(fn (float $x): float => self::round($x, $testDigits), $c->statistics), $digits);
        $pValues = self::formatPval($c->pValues, $testDigits);
        $na = fn (string $s, float $x): string => is_nan($x) ? 'NA' : trim($s);
        $rows = [];
        foreach ($c->names as $i => $name) {
            $rows[] = [
                $name,
                $na($estimateText[$i], $c->estimates[$i]),
                $na($estimateText[$p + $i], $c->stdErrors[$i]),
                $na($statistics[$i], $c->statistics[$i]),
                trim($pValues[$i]),
            ];
        }
        [$stat, $pr] = $c->statistic === 't' ? ['t value', 'Pr(>|t|)'] : ['z value', 'Pr(>|z|)'];

        return self::table(['', 'Estimate', 'Std. Error', $stat, $pr], $rows, []);
    }

    /** @param  list<RandomTerm>  $terms */
    public static function randomEffects(array $terms, ?float $residual): string
    {
        $rows = [];
        $variances = [];
        $correlations = [];
        $width = max(array_map(fn (RandomTerm $t): int => $t->correlated ? count($t->columns) : 1, $terms) ?: [1]) - 1;
        foreach ($terms as $term) {
            foreach ($term->columns as $j => $name) {
                $rows[] = [$j === 0 ? $term->group : '', $name];
                $cov = $term->covariance;
                $variances[] = $cov[$j][$j];
                $corr = array_fill(0, max(0, $width), '');
                if ($term->correlated) {
                    for ($m = 0; $m < min($j, $width); $m++) {
                        $corr[$m] = self::rustFixed(fdiv($cov[$j][$m], sqrt($cov[$j][$j] * $cov[$m][$m])), 2);
                    }
                }
                $correlations[] = $corr;
            }
        }
        if ($residual !== null) {
            $rows[] = ['Residual', ''];
            $variances[] = $residual;
            $correlations[] = array_fill(0, max(0, $width), '');
        }
        $vt = self::formatReal($variances, self::DIGITS);
        $st = self::formatReal(array_map(sqrt(...), $variances), self::DIGITS);
        foreach ($rows as $i => &$row) {
            $row[] = trim($vt[$i]);
            $row[] = trim($st[$i]);
            array_push($row, ...$correlations[$i]);
        }
        unset($row);
        $header = ['Groups', 'Name', 'Variance', 'Std.Dev.', ...array_fill(0, max(0, $width), '')];
        if ($width > 0) {
            $header[4] = 'Corr';
        }

        return self::table($header, $rows, [1]);
    }

    /** Rust's {:.d} keeping the sign of -0.0, as format!("{v:.1}") prints it. */
    public static function rustFixedSigned(float $x, int $decimals): string
    {
        if ($x == 0.0 && ! is_nan($x)) {
            return (fdiv(1.0, $x) < 0 ? '-' : '').self::fixed(0.0, $decimals);
        }

        return self::rustFixed($x, $decimals);
    }

    /**
     * A table of named columns in the layout of R's print.data.frame(): row numbers left-aligned in an unnamed first
     * column, every other column right-aligned under its name, one space apart. A column is text when every value is
     * a string or null, and numeric otherwise, formatted per column as format(x, digits = 7). A missing number prints
     * NA and a missing string <NA>. As R does at its default width of 80, columns are split into blocks whose lines
     * hold at most 79 characters, each block with its own header and row numbers.
     *
     * @param  array<string|int, mixed>  $columns
     */
    public static function dataFrame(array $columns): string
    {
        $names = [];
        $texts = [];
        foreach ($columns as $key => $column) {
            $names[] = (string) $key;
            if (! is_array($column)) {
                throw \Phit\Polyfill\ErrorCode::exception(\Phit\Polyfill\ErrorCode::INVALID_INPUT, 'table columns must be arrays');
            }
            $values = array_values($column);
            if (array_filter($values, fn (mixed $v): bool => $v !== null && ! is_string($v)) === []) {
                // A string that is not UTF-8 is not text to the extension, so it prints as missing.
                $texts[] = array_map(fn (mixed $v): string => is_string($v) && preg_match('//u', $v) === 1 ? $v : '<NA>', $values);

                continue;
            }
            $numbers = array_map(fn (mixed $v): ?float => is_float($v) || is_int($v) ? (float) $v : null, $values);
            $formatted = self::formatReal(array_values(array_filter($numbers, fn (?float $v): bool => $v !== null)), 7);
            $cells = [];
            foreach ($numbers as $v) {
                $cells[] = $v === null ? 'NA' : array_shift($formatted) ?? '';
            }
            $texts[] = $cells;
        }
        $rows = count($texts[0] ?? []);
        if ($rows === 0) {
            return '<0 rows> '.implode(', ', $names)."\n";
        }
        $labelWidth = strlen((string) $rows);
        $widths = [];
        foreach ($names as $c => $name) {
            $widths[$c] = max([self::width($name), ...array_map(self::width(...), $texts[$c])]);
        }
        // Blocks of columns, filled greedily while a line stays within 79 characters.
        $blocks = [];
        $start = 0;
        while ($start < count($names)) {
            $end = $start + 1;
            $line = $labelWidth + 1 + $widths[$start];
            while ($end < count($names) && $line + 1 + $widths[$end] <= 79) {
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
            for ($r = 0; $r < $rows; $r++) {
                $out .= str_pad((string) ($r + 1), $labelWidth);
                for ($c = $from; $c < $to; $c++) {
                    $out .= ' '.self::padLeft(trim($texts[$c][$r]), $widths[$c]);
                }
                $out .= "\n";
            }
        }

        return $out;
    }
}
