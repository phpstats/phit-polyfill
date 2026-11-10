<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/** Coefficient inference and labeling for sandwich covariance. */
final readonly class Robust implements Printer
{
    public function __construct(private Printer $base,private string $label,private float $df) {}
    public static function fromFit(array $fit): self {return new self(Linear::fromFit($fit),$fit['vcov_label'],(float)$fit['inference_df']);}
    public function heading(): string {return $this->base->heading();}
    public function beforeCoefficients(): string {return $this->base->beforeCoefficients()."Covariance: {$this->label}\nInference: t with {$this->df} degrees of freedom\n";}
    public function coefficientsTitle(): string {return $this->base->coefficientsTitle();}
    public function coefficients(): Coefficients {return $this->base->coefficients();}
    public function afterCoefficients(): string {return preg_replace('/^F-statistic:.*\n/m','',$this->base->afterCoefficients());}
    public function intervalDf(): ?float {return $this->df;}
    public function randomTerms(): ?array {return null;}
}
