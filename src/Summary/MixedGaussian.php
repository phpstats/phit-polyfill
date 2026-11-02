<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/** A Gaussian linear mixed model, as lme4's summary.merMod for lmer(). The `model` key is "mixed_gaussian". */
final readonly class MixedGaussian implements Printer
{
    public function __construct(
        private ?string $formula,
        private bool $reml,
        private MixedFit $fit,
        private float $residualVariance,
        private Coefficients $coefficients,
    ) {}

    public static function fromFit(array $fit): self
    {
        return new self(
            FitReader::text($fit, 'formula'),
            FitReader::text($fit, 'method') === 'REML',
            MixedFit::fromFit($fit),
            FitReader::requiredFloat($fit, 'residual_variance'),
            FitReader::coefficients($fit),
        );
    }

    public function heading(): string
    {
        return Format::heading('Linear mixed model fit by '.($this->reml ? 'REML' : 'maximum likelihood'), null, $this->formula);
    }

    public function beforeCoefficients(): string
    {
        $criteria = $this->reml
            ? "\nREML criterion at convergence: ".Format::formatReal([-2.0 * $this->fit->logLikelihood()], Format::DIGITS + 1)[0]."\n"
            : $this->fit->criteria();

        return $criteria.$this->fit->random($this->residualVariance);
    }

    public function coefficientsTitle(): string
    {
        return 'Fixed effects';
    }

    public function coefficients(): Coefficients
    {
        return $this->coefficients;
    }

    public function afterCoefficients(): string
    {
        return $this->fit->singularNote();
    }

    public function intervalDf(): ?float
    {
        return null;
    }
}
