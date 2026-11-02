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
        'generalized_linear' => Summary\GeneralizedLinear::class,
        'mixed_gaussian' => Summary\MixedGaussian::class,
        'mixed_logistic' => Summary\MixedLogistic::class,
    ];

    /** Phit\summary */
    public static function summary(array $fit): string
    {
        return self::render(self::printer($fit));
    }

    /** Phit\tidy */
    public static function tidy(array $fit, bool $confInt = false, float $level = 0.95): array
    {
        // An explicit NaN test: opcache rewrites a negated float comparison (docs/notes/opcache-nan.md).
        if (is_nan($level) || $level <= 0.0 || $level >= 1.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'level must be in (0, 1)');
        }
        $printer = self::printer($fit);
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

    /** The text of print(summary(fit)) for any printer, ending with one newline. */
    public static function render(Printer $p): string
    {
        return $p->heading().$p->beforeCoefficients()."\n{$p->coefficientsTitle()}:\n"
            .Format::coefficientTable($p->coefficients()).$p->afterCoefficients();
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
