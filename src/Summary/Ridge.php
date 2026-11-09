<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/** Conditional ridge uncertainty; coefficients at bounds are unavailable. */
final readonly class Ridge implements Printer
{
    private function __construct(private array $fit, private Coefficients $table) {}
    public static function fromFit(array $fit): self { return new self($fit, FitReader::coefficients($fit)); }
    public function heading(): string
    {
        return "Gaussian ridge regression\nFormula: ".$this->fit['formula']."\nLambda: ".Format::formatG($this->fit['lambda'],4)."\n";
    }
    public function beforeCoefficients(): string { return ''; }
    public function coefficientsTitle(): string { return 'Coefficients'; }
    public function coefficients(): Coefficients { return $this->table; }
    public function afterCoefficients(): string
    {
        return "\nEffective degrees of freedom: ".Format::formatG($this->fit['df'],4)
            ."\nResidual standard error: ".Format::formatG($this->fit['sigma'],4)
            ."\nR-squared: ".Format::formatG($this->fit['r_squared'],4)
            ."\nConditional SEs and intervals hold lambda fixed. Estimates are biased; intervals are not coverage-valid.\nAt-bound coefficients have unavailable uncertainty.\n";
    }
    public function intervalDf(): ?float { return (float)$this->fit['df_residual']; }
    public function randomTerms(): ?array { return null; }
}
