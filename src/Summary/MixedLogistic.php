<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/** A binomial logit mixed model, as lme4's summary.merMod for glmer(). The `model` key is "mixed_logistic". */
final readonly class MixedLogistic implements Printer
{
    /** @param  string  $integration  the integration method, such as "Laplace approximation" */
    public function __construct(
        private ?string $formula,
        private string $integration,
        private MixedFit $fit,
        private Coefficients $coefficients,
    ) {}

    public static function fromFit(array $fit): self
    {
        $integration = match (FitReader::text($fit, 'integration')) {
            'laplace' => 'Laplace approximation',
            null => throw FitReader::notAFit('missing integration'),
            default => 'adaptive Gauss-Hermite quadrature, '.(FitReader::long($fit, 'quadrature_points') ?? 0).' points',
        };

        return new self(FitReader::text($fit, 'formula'), $integration, MixedFit::fromFit($fit), FitReader::coefficients($fit));
    }

    public function heading(): string
    {
        return Format::heading(
            "Generalized linear mixed model fit by maximum likelihood ({$this->integration})",
            ['binomial', 'logit'],
            $this->formula,
        );
    }

    public function beforeCoefficients(): string
    {
        return $this->fit->criteria().$this->fit->random(null);
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

    public function randomTerms(): ?array
    {
        return $this->fit->terms();
    }
}
