<?php

declare(strict_types=1);

use Phit\Polyfill\Summary\Format;

// The tables of phit's specs/cleanroom/data (copied to tests/fixtures/cleanroom)
// for R's printed number format (format.md). Every case must match exactly;
// each check reports its mismatches, or for round its largest difference, on
// STDERR.

if (! function_exists('formatCell')) {
    /** A %.17g cell of the R tables. */
    function formatCell(string|int|float $cell): float
    {
        if (! is_string($cell)) {
            return (float) $cell;
        }

        return match ($cell) {
            'Inf' => INF,
            '-Inf' => -INF,
            'NaN', 'NA' => NAN,
            default => (float) $cell,
        };
    }

    /**
     * A vector of the R tables, which write a vector of length one as a scalar.
     *
     * @return list<string|int|float>
     */
    function formatList(mixed $cells): array
    {
        return is_array($cells) ? $cells : [$cells];
    }

    function formatReport(string $table, int $rows, string $result): void
    {
        fwrite(STDERR, sprintf("%s: %d rows, %s\n", $table, $rows, $result));
    }

    /**
     * The trimmed cells of R's printed coefficient table, by column. Columns are
     * right-aligned under their headers, so each ends where its header ends.
     *
     * @param  list<string>  $printed
     * @return list<list<string>>
     */
    function formatCoefmatColumns(array $printed): array
    {
        $header = $printed[0];
        $ends = [];
        $from = 0;
        foreach (['Estimate', 'Std. Error', ' value', 'Pr(>|'] as $label) {
            $at = strpos($header, $label, $from);
            expect($at)->not->toBeFalse($header);
            $end = $label === 'Pr(>|' ? strlen($header) : $at + strlen($label);
            $ends[] = $end;
            $from = $end;
        }
        $columns = [[], [], [], []];
        foreach (array_slice($printed, 1) as $line) {
            $start = strpos($line, ' '); // after the row name
            foreach ($ends as $j => $end) {
                $columns[$j][] = trim(substr($line, $start, $end - $start));
                $start = $end;
            }
        }

        return $columns;
    }
}

it('matches R format()', function () {
    $cases = loadFixture('cleanroom/format.json');
    expect($cases)->toHaveCount(605);
    foreach ($cases as $case) {
        $got = Format::formatReal(array_map(formatCell(...), formatList($case['x'])), $case['digits']);
        expect($got)->toBe(formatList($case['out']), json_encode($case));
    }
    formatReport('format', count($cases), '0 mismatches');
});

it('matches R round() exactly, sign of zero included', function () {
    $cases = loadFixture('cleanroom/round.json');
    expect($cases)->toHaveCount(300);
    $bits = fn (float $x): string => bin2hex(pack('E', $x));
    $worst = 0.0;
    foreach ($cases as $case) {
        $got = Format::round(formatCell($case['x']), $case['digits']);
        $expected = formatCell($case['out']);
        expect($bits($got))->toBe($bits($expected), json_encode($case).' got '.var_export($got, true));
        $worst = max($worst, $expected == 0.0 ? abs($got) : abs($got - $expected) / abs($expected));
    }
    formatReport('round', count($cases), sprintf('largest relative difference %.3e', $worst));
});

it('matches R formatC() once its padding is trimmed', function () {
    $cases = loadFixture('cleanroom/format_g.json');
    expect($cases)->toHaveCount(90);
    foreach ($cases as $case) {
        expect(Format::formatG(formatCell($case['x']), $case['digits']))->toBe(trim($case['out']), json_encode($case));
    }
    formatReport('format_g', count($cases), '0 mismatches');
});

it('matches R format.pval(), padding included', function () {
    $cases = loadFixture('cleanroom/format_pval.json');
    expect($cases)->toHaveCount(200);
    foreach ($cases as $case) {
        $got = Format::formatPval(array_map(formatCell(...), formatList($case['p'])), $case['digits']);
        expect($got)->toBe(formatList($case['out']), json_encode($case));
    }
    formatReport('format_pval', count($cases), '0 mismatches');
});

it('matches the cells R printCoefmat() prints', function () {
    $cases = loadFixture('cleanroom/coefmat.json');
    expect($cases)->toHaveCount(120);
    foreach ($cases as $case) {
        $cells = Format::coefficientCells(
            array_map(formatCell(...), formatList($case['estimate'])),
            array_map(formatCell(...), formatList($case['std_error'])),
            array_map(formatCell(...), formatList($case['statistic'])),
            array_map(formatCell(...), formatList($case['p_value'])),
        );
        $expected = formatCoefmatColumns($case['printed']);
        $label = implode("\n", $case['printed']);
        expect($cells['estimates'])->toBe($expected[0], $label)
            ->and($cells['std_errors'])->toBe($expected[1], $label)
            ->and($cells['statistics'])->toBe($expected[2], $label)
            ->and($cells['p_values'])->toBe($expected[3], $label);
    }
    formatReport('coefmat', count($cases), '0 mismatches');
});
