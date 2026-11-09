<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Distribution\Normal;
use Phit\Polyfill\Distribution\StudentT;
use Phit\Polyfill\Summary\FitReader;
use Phit\Polyfill\Summary\Format;
use Phit\Polyfill\Summary\Printer;

/**
 * Phit\summary() prints a fit in the layout of R's print(summary(fit)), and
 * Phit\tidy() returns its coefficient table as columns, as broom's tidy()
 * does. Both pick the fit's printer from its `model` key through PRINTERS, as
 * the extension's printer.rs does; the layouts are the Summary\Printer
 * implementations, and the shared formatting is Summary\Format.
 */
final class Summary
{
    /**
     * One printer per value of a fit's `model` key, which every model function
     * sets. A new model kind adds one entry here and one Printer.
     *
     * @var array<string, class-string<Printer>>
     */
    public const array PRINTERS = [
        'linear' => Summary\Linear::class,
        'iv' => Summary\Linear::class,
        'ridge' => Summary\Ridge::class,
        'generalized_linear' => Summary\GeneralizedLinear::class,
        'mixed_gaussian' => Summary\MixedGaussian::class,
        'mixed_logistic' => Summary\MixedLogistic::class,
    ];

    /** Phit\summary */
    public static function summary(mixed $fit): string
    {
        $data = self::fitArgument($fit);
        $text = self::render(self::printer($data));
        if (($data['model'] ?? null) !== 'iv') return $text;
        $text = str_replace('Linear model fit by least squares', 'Instrumental-variable model fit by two-stage least squares', $text);
        $text = preg_replace('/^F-statistic:.*\n/m', '', $text);
        $text .= "\nDiagnostic tests (iid covariance):\nTest  df1  df2  statistic  p-value\n";
        $fmt = fn ($v) => $v === null ? 'NA' : Format::formatG((float)$v, 4);
        $row = function ($label, $a, $b, $stat, $p) use (&$text, $fmt) {
            $text .= "$label  ".$fmt($a).'  '.$fmt($b).'  '.$fmt($stat).'  '.$fmt($p)."\n";
        };
        foreach ($data['first_stage'] as $name => $stage) {
            $row("Weak instruments ($name)", $stage['df1'], $stage['df2'], $stage['f_statistic'], $stage['p_value']);
            $row("Sanderson-Windmeijer ($name)", $stage['sw_df1'], $stage['sw_df2'], $stage['sw_f_statistic'], $stage['sw_p_value']);
        }
        $wu = $data['diagnostics']['wu_hausman'];
        $row('Wu-Hausman', $wu['df1'], $wu['df2'], $wu['statistic'], $wu['p_value']);
        if (($sargan = $data['diagnostics']['sargan']) !== null) $row('Sargan', $sargan['df'], null, $sargan['statistic'], $sargan['p_value']);
        return $text;
    }

    /** Phit\tidy */
    public static function tidy(mixed $fit, bool $confInt = false, float $level = 0.95, string $effects = 'fixed'): array
    {
        // An explicit NaN test: opcache rewrites a negated float comparison (docs/notes/opcache-nan.md).
        if (is_nan($level) || $level <= 0.0 || $level >= 1.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'level must be in (0, 1)');
        }
        if ($effects !== 'fixed' && $effects !== 'ran_vals') {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "effects must be 'fixed' or 'ran_vals'");
        }
        $printer = self::printer(self::fitArgument($fit));

        $out = $effects === 'fixed' ? self::tidyFixed($printer, $confInt, $level) : self::tidyRandom($printer, $confInt, $level);
        if ($printer instanceof Summary\Ridge) foreach (['std_error','statistic','p_value','conf_low','conf_high'] as $key) {
            if (isset($out[$key])) $out[$key] = array_map(fn ($v) => is_nan($v) ? null : $v, $out[$key]);
        }
        return $out;
    }

    /** The fixed effects, as broom's tidy(). */
    private static function tidyFixed(Printer $printer, bool $confInt, float $level): array
    {
        $c = $printer->coefficients();
        $out = [
            'term' => $c->names,
            'estimate' => $c->estimates,
            'std_error' => $c->stdErrors,
            'statistic' => $c->statistics,
            'p_value' => $c->pValues,
        ];
        if ($confInt) {
            [$out['conf_low'], $out['conf_high']] = self::confidenceIntervals($c->estimates, $c->stdErrors, $printer->intervalDf(), $level);
        }

        return $out;
    }

    /**
     * Each group's conditional modes and standard errors, as broom.mixed's tidy(effects = "ran_vals"): one row per
     * term, column, and group, in that order. Intervals use the normal distribution.
     */
    private static function tidyRandom(Printer $printer, bool $confInt, float $level): array
    {
        $terms = $printer->randomTerms()
            ?? throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "effects 'ran_vals' needs a mixed model fit");
        $out = ['group' => [], 'level' => [], 'term' => [], 'estimate' => [], 'std_error' => []];
        foreach ($terms as $t) {
            $values = $t->values ?? throw FitReader::notAFit('missing random effects');
            foreach ($t->columns as $j => $name) {
                foreach ($values->levels as $g => $label) {
                    $out['group'][] = $t->group;
                    $out['level'][] = $label;
                    $out['term'][] = $name;
                    $out['estimate'][] = $values->estimates[$g][$j];
                    $out['std_error'][] = $values->stdErrors[$g][$j];
                }
            }
        }
        if ($confInt) {
            [$out['conf_low'], $out['conf_high']] = self::confidenceIntervals($out['estimate'], $out['std_error'], null, $level);
        }

        return $out;
    }

    /** The text of print(summary(fit)) for any printer, ending with one newline. */
    public static function render(Printer $p): string
    {
        return $p->heading().$p->beforeCoefficients()."\n{$p->coefficientsTitle()}:\n"
            .Format::coefficientTable($p->coefficients()).$p->afterCoefficients();
    }

    /** The result array behind a Fit|array argument. */
    private static function fitArgument(mixed $fit): array
    {
        return Objects\Wrap::data($fit) ?? throw FitReader::notAFit('not a fit object or array');
    }

    private static function printer(array $fit): Printer
    {
        $kind = FitReader::text($fit, 'model') ?? throw FitReader::notAFit('missing model');
        $class = self::PRINTERS[$kind] ?? throw FitReader::notAFit("unknown model '$kind'");

        return $class::fromFit($fit);
    }

    /**
     * Wald confidence intervals: t with $df degrees of freedom when given, else the standard normal.
     *
     * @param  list<float>  $estimates
     * @param  list<float>  $stdErrors
     * @return array{list<float>, list<float>}
     */
    private static function confidenceIntervals(array $estimates, array $stdErrors, ?float $df, float $level): array
    {
        $p = 0.5 + $level / 2.0;
        try {
            $q = $df === null ? Normal::create(0.0, 1.0)->quantileAt($p) : StudentT::create($df)->quantileAt($p);
        } catch (\Exception) {
            $q = NAN;
        }
        $low = [];
        $high = [];
        foreach ($estimates as $i => $b) {
            $low[] = $b - $q * $stdErrors[$i];
            $high[] = $b + $q * $stdErrors[$i];
        }

        return [$low, $high];
    }
}
