<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Internal\CalledAs;
use Phit\Polyfill\Internal\Deprecation;

/**
 * Phit\fit_ols: the formula fit (phit 0.3's Phit\fit_lm) when the first
 * argument is a string or formula: and data: are named, and otherwise the
 * deprecated matrix fit (Phit\fit_ols_matrix). An explicit null counts as not
 * given, and every error is thrown before the deprecation.
 */
final class FitOls
{
    public static function fit(
        mixed $formula = null,
        mixed $data = null,
        string $missing = 'drop',
        mixed $weights = null,
        ?array $offset = null,
        bool $intercept = true,
        ?array $x = null,
        ?array $y = null,
    ): array {
        // A string that is not UTF-8 is not a formula to the extension either.
        if (is_string($formula) && mb_check_encoding($formula, 'UTF-8')) {
            if ($x !== null || $y !== null) {
                throw self::invalid('Phit\fit_ols takes a formula and data, or a matrix x and y, not both');
            }
            if ($offset !== null || ! $intercept) {
                throw self::invalid('Phit\fit_ols with a formula takes no offset or intercept arguments; write offset() or - 1 in the formula');
            }
            if ($data === null) {
                throw self::invalid('Phit\fit_ols with a formula needs data');
            }

            return Lm::fit($formula, $data, $missing, $weights);
        }
        if ($formula !== null && ! is_array($formula)) {
            throw self::invalid('Phit\fit_ols needs a formula string, or a matrix x as an array');
        }
        if (($formula !== null && $x !== null) || ($data !== null && $y !== null)) {
            throw self::invalid('Phit\fit_ols takes a formula and data, or a matrix x and y, not both');
        }
        $x ??= $formula;
        if ($y === null && $data !== null) {
            if (! is_array($data)) {
                throw self::invalid('Phit\fit_ols with a matrix x needs y as an array');
            }
            $y = $data;
        }
        if ($x === null || $y === null) {
            throw self::invalid('Phit\fit_ols needs a formula and data (or, deprecated, a matrix x and y)');
        }
        if ($weights !== null && ! is_array($weights)) {
            throw self::invalid('Phit\fit_ols with a matrix x needs weights as an array');
        }
        Deprecation::olsMatrix();

        return CalledAs::run('Phit\fit_ols', static fn (): array => Ols::fit($x, $y, $missing, $weights, $offset, $intercept));
    }

    private static function invalid(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }
}
