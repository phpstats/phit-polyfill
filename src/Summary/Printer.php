<?php

declare(strict_types=1);

namespace Phit\Polyfill\Summary;

/**
 * How one kind of fitted model prints, in the layout of R's
 * print(summary(fit)): stats_core::summary::Printer. Each model kind
 * implements it, and Phit\Polyfill\Summary::render() lays the parts out the
 * same way for all of them, with the shared formatting in Format.
 */
interface Printer
{
    /** The printer of a fit with this model kind, read from the fit's keys. */
    public static function fromFit(array $fit): self;

    /** The heading: the fitting method, then the family and the formula when known. Each line ends in a newline. */
    public function heading(): string;

    /** Fit statistics printed before the coefficient table. Mixed models print their criteria and random effects here. */
    public function beforeCoefficients(): string;

    /** The coefficient table's title. */
    public function coefficientsTitle(): string;

    /** The coefficient estimates, standard errors, tests, and p-values. tidy() reads the same data. */
    public function coefficients(): Coefficients;

    /** Fit statistics printed after the coefficient table, each line ending in a newline. */
    public function afterCoefficients(): string;

    /** Degrees of freedom for t confidence intervals, or null for normal intervals. */
    public function intervalDf(): ?float;

    /**
     * The random-effect terms of a mixed model, or null for a model without random effects.
     *
     * @return ?list<RandomTerm>
     */
    public function randomTerms(): ?array;
}
