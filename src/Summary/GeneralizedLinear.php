<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/**
 * A generalized linear model fitted by maximum likelihood, as R's summary.glm,
 * for any family and link. The `model` key is "generalized_linear".
 */
final readonly class GeneralizedLinear implements Printer
{
    /** @param  ?array{string, string}  $family  the family and link, when the fit reports them */
    public function __construct(
        private ?array $family,
        private ?string $formula,
        private float $dispersion,
        private float $nullDeviance,
        private float $deviance,
        private ?int $dfNull,
        private int $dfResidual,
        private float $aic,
        private Coefficients $coefficients,
    ) {}

    public static function fromFit(array $fit): self
    {
        $coefficients = FitReader::coefficients($fit);
        $n = FitReader::nObservations($fit);
        $family = FitReader::text($fit, 'family');
        $link = FitReader::text($fit, 'link');
        $intercept = FitReader::intercept($fit);

        return new self(
            $family !== null && $link !== null ? [$family, $link] : null,
            FitReader::text($fit, 'formula'),
            // A Gaussian GLM estimates its dispersion and reports t statistics; the others fix it at one.
            $coefficients->statistic === 't' ? FitReader::requiredFloat($fit, 'residual_variance') : 1.0,
            FitReader::requiredFloat($fit, 'null_deviance'),
            FitReader::requiredFloat($fit, 'deviance'),
            $intercept === null ? null : $n - (int) $intercept,
            FitReader::dfResidual($fit),
            FitReader::requiredFloat($fit, 'aic'),
            $coefficients,
        );
    }

    public function heading(): string
    {
        return Format::heading('Generalized linear model fit by maximum likelihood', $this->family, $this->formula);
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
        $family = $this->family === null ? '' : " for {$this->family[0]} family";
        $out = "\n(Dispersion parameter$family taken to be ".Format::formatReal([$this->dispersion], 7)[0].")\n\n";
        $deviances = Format::formatReal([$this->nullDeviance, $this->deviance], Format::DIGITS + 1);
        if ($this->dfNull !== null) {
            $dfs = Format::formatReal([(float) $this->dfNull, (float) $this->dfResidual], 15);
            $out .= "    Null deviance: {$deviances[0]} on {$dfs[0]} degrees of freedom\n";
            $out .= "Residual deviance: {$deviances[1]} on {$dfs[1]} degrees of freedom\n";
        } else {
            $out .= "    Null deviance: {$deviances[0]}\n";
            $out .= "Residual deviance: {$deviances[1]} on {$this->dfResidual} degrees of freedom\n";
        }

        return $out.'AIC: '.Format::formatReal([$this->aic], Format::DIGITS + 1)[0]."\n";
    }

    public function intervalDf(): ?float
    {
        return $this->coefficients->statistic === 't' ? (float) $this->dfResidual : null;
    }

    public function randomTerms(): ?array
    {
        return null;
    }
}
