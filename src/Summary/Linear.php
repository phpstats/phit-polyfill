<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

use Phit\Polyfill\Distribution\FisherF;

/** A linear model fitted by least squares, as R's summary.lm. The `model` key is "linear". */
final readonly class Linear implements Printer
{
    /** @param  ?int  $fDf  numerator degrees of freedom of the overall F test, when the intercept is known */
    public function __construct(
        private ?string $formula,
        private float $sigma,
        private int $dfResidual,
        private float $rSquared,
        private float $adjRSquared,
        private ?int $fDf,
        private Coefficients $coefficients,
    ) {}

    public static function fromFit(array $fit): self
    {
        $coefficients = FitReader::coefficients($fit);
        $intercept = FitReader::intercept($fit);

        return new self(
            FitReader::text($fit, 'formula'),
            FitReader::requiredFloat($fit, 'sigma'),
            FitReader::dfResidual($fit),
            FitReader::requiredFloat($fit, 'r_squared'),
            FitReader::requiredFloat($fit, 'adj_r_squared'),
            $intercept === null ? null : count($coefficients->estimates) - (int) $intercept,
            $coefficients,
        );
    }

    public function heading(): string
    {
        return Format::heading('Linear model fit by least squares', null, $this->formula);
    }

    public function beforeCoefficients(): string
    {
        return '';
    }

    public function coefficientsTitle(): string
    {
        return 'Coefficients';
    }

    public function coefficients(): Coefficients
    {
        return $this->coefficients;
    }

    public function afterCoefficients(): string
    {
        $df = $this->dfResidual;
        $out = "\nResidual standard error: ".Format::formatReal([$this->sigma], Format::DIGITS)[0]." on $df degrees of freedom\n";
        $out .= 'Multiple R-squared: '.Format::formatG($this->rSquared, Format::DIGITS)
            .', Adjusted R-squared: '.Format::formatG($this->adjRSquared, Format::DIGITS)."\n";
        $numerator = $this->fDf;
        if ($numerator !== null && $numerator > 0) {
            $r2 = $this->rSquared;
            $f = fdiv(fdiv($r2, (float) $numerator), fdiv(1.0 - $r2, (float) $df));
            try {
                $p = FisherF::create((float) $numerator, (float) $df)->sfAt($f);
            } catch (\Exception) {
                $p = NAN;
            }
            $out .= 'F-statistic: '.Format::formatG($f, Format::DIGITS)." on $numerator and $df DF, p-value: "
                .Format::formatPval([$p], Format::DIGITS)[0]."\n";
        }

        return $out;
    }

    public function intervalDf(): ?float
    {
        return (float) $this->dfResidual;
    }

    public function randomTerms(): ?array
    {
        return null;
    }
}
