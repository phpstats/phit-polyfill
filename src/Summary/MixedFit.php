<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/** What every mixed model prints before its fixed effects. */
final readonly class MixedFit
{
    /** @param  list<RandomTerm>  $terms */
    public function __construct(
        private float $logLikelihood,
        private ?float $aic,
        private ?float $bic,
        private int $dfResidual,
        private int $nObservations,
        private array $terms,
        private bool $singular,
    ) {}

    public static function fromFit(array $fit): self
    {
        return new self(
            FitReader::requiredFloat($fit, 'log_likelihood'),
            FitReader::float($fit, 'aic'),
            FitReader::float($fit, 'bic'),
            FitReader::dfResidual($fit),
            FitReader::nObservations($fit),
            FitReader::randomTerms($fit),
            (array_key_exists('singular', $fit) ? $fit['singular'] : ($fit['boundary'] ?? null)) === true,
        );
    }

    public function logLikelihood(): float
    {
        return $this->logLikelihood;
    }

    /** @return list<RandomTerm> */
    public function terms(): array
    {
        return $this->terms;
    }

    /** The ML criteria table, as lme4 prints it: one decimal each. */
    public function criteria(): string
    {
        $names = [];
        $values = [];
        foreach (['AIC' => $this->aic, 'BIC' => $this->bic] as $name => $value) {
            if ($value !== null) {
                $names[] = $name;
                $values[] = $value;
            }
        }
        array_push($names, 'logLik', 'deviance');
        array_push($values, $this->logLikelihood, -2.0 * $this->logLikelihood);
        $texts = array_map(fn (float $v): string => Format::rustFixedSigned($v, 1), $values);
        $names[] = 'df.resid';
        $texts[] = (string) $this->dfResidual;
        $header = [];
        $row = [];
        foreach ($names as $k => $name) {
            $w = max(strlen($name), strlen($texts[$k]));
            $header[] = Format::padLeft($name, $w);
            $row[] = Format::padLeft($texts[$k], $w);
        }

        return "\n".implode(' ', $header)."\n".implode(' ', $row)."\n";
    }

    /** The random-effects table and the observation and group counts. */
    public function random(?float $residual): string
    {
        $seen = [];
        foreach ($this->terms as $term) {
            if (! array_key_exists($term->group, $seen)) {
                $seen[$term->group] = $term->nGroups;
            }
        }
        $groups = [];
        foreach ($seen as $group => $n) {
            $groups[] = "$group, $n";
        }

        return "\nRandom effects:\n".Format::randomEffects($this->terms, $residual)
            ."Number of obs: {$this->nObservations}, groups: ".implode('; ', $groups)."\n";
    }

    public function singularNote(): string
    {
        return $this->singular ? "\nThe fit is singular: a random-effect variance or correlation is at its boundary.\n" : '';
    }
}
