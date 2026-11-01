<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Distribution\Support;
use Phit\Polyfill\Formula\Design;
use Phit\Polyfill\Formula\Frontend;
use Phit\Polyfill\Internal\Family\LogisticFamily;
use Phit\Polyfill\Internal\Family\PoissonFamily;
use Phit\Polyfill\Internal\Family\ProbitFamily;

/**
 * Generalized linear model from a formula over the existing GLM cores (Phit\fit_generalized_linear).
 */
final class Glm
{
    public static function fit(
        string $formula,
        mixed $data,
        string $family = 'gaussian',
        string $link = '',
        string $missing = 'drop',
        mixed $weights = null,
    ): array {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($family, 'family');
        Frontend::stringArgument($link, 'link');
        Frontend::stringArgument($missing, 'missing');
        $link = match (true) {
            $family === 'gaussian' && ($link === '' || $link === 'identity') => 'identity',
            $family === 'binomial' && ($link === '' || $link === 'logit') => 'logit',
            $family === 'binomial' && $link === 'probit' => 'probit',
            $family === 'poisson' && ($link === '' || $link === 'log') => 'log',
            in_array($family, ['gaussian', 'binomial', 'poisson'], true) => throw ErrorCode::exception(
                ErrorCode::UNSUPPORTED,
                "link '$link' is not supported for family '$family'; supported: gaussian/identity, binomial/logit, binomial/probit, poisson/log",
            ),
            default => throw ErrorCode::exception(
                ErrorCode::UNSUPPORTED,
                "family '$family' is not supported; use gaussian, binomial or poisson",
            ),
        };
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        Frontend::rejectRandom($design, 'Phit\fit_generalized_linear');
        $x = $design->x;
        $n = count($x);
        $offset = $design->offset ?? array_fill(0, $n, 0.0);
        if ($link === 'identity') {
            $y = Frontend::numericResponse($design, 'Phit\fit_generalized_linear');
            $output = Ols::fitWithData($x, $y, $w, $offset);
            $p = count($output['coefficients']);
            // R's glm null deviance: intercept-plus-offset fit, else the offset alone.
            $centered = [];
            foreach ($y as $i => $v) {
                $centered[$i] = $v - $offset[$i];
            }
            $center = 0.0;
            if ($design->intercept) {
                $num = 0.0;
                $den = 0.0;
                foreach ($centered as $i => $v) {
                    $num += $v * $w[$i];
                }
                foreach ($w as $value) {
                    $den += $value;
                }
                $center = $num / $den;
            }
            $null = 0.0;
            foreach ($centered as $i => $v) {
                $d = $v - $center;
                $null += $w[$i] * ($d * $d);
            }
            $output['deviance'] = $output['rss'];
            $output['null_deviance'] = $null;
            $output['aic'] = -2.0 * $output['log_likelihood'] + 2.0 * ($p + 1.0);
        } elseif ($link === 'logit' || $link === 'probit') {
            $y = self::binomialResponse($design, $w);
            $core = $link === 'logit' ? new LogisticFamily : new ProbitFamily;
            $output = $core->fitData($x, $y, $w, $offset);
        } else {
            $y = Frontend::numericResponse($design, 'Phit\fit_generalized_linear');
            foreach ($y as $v) {
                if ($v < 0.0 || fmod($v, 1.0) != 0.0) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'Poisson y must contain nonnegative integer counts');
                }
            }
            $output = (new PoissonFamily())->fitData($x, $y, $w, $offset);
        }
        $output['family'] = $family;
        $output['link'] = $link;
        $output = Frontend::addNames($output, $design);
        $solver = match ($link) {
            'identity' => 'scaled_svd',
            'probit' => 'fisher_scoring',
            default => 'newton_irls',
        };

        return Frontend::finish($output, $design, 'Phit\fit_generalized_linear', $solver);
    }

    /**
     * Binomial outcomes as proportions; cbind(successes, failures) multiplies
     * the prior weights by the trials. A factor's first level is failure.
     *
     * @param  list<float>  $weights
     * @return list<float>
     */
    /** @internal shared with Phit\fit_mixed_logistic */
    public static function binomialResponse(Design $design, array &$weights): array
    {
        $response = $design->response;
        if ($response === null) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "Phit\\fit_generalized_linear needs a response: 'y ~ ...'");
        }
        $y = [];
        if ($response[0] === 'numeric') {
            $y = $response[1];
        } elseif ($response[0] === 'factor') {
            foreach ($response[2] as $c) {
                $y[] = $c !== 0 ? 1.0 : 0.0;
            }
        } else {
            foreach ($response[1] as $i => $s) {
                $f = $response[2][$i];
                if ($s < 0.0 || $f < 0.0 || fmod($s, 1.0) != 0.0 || fmod($f, 1.0) != 0.0 || $s + $f <= 0.0) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'cbind(successes, failures) needs nonnegative integer counts with positive totals');
                }
                $y[] = $s / ($s + $f);
                $weights[$i] *= $s + $f;
            }
        }
        foreach ($y as $v) {
            if (is_nan($v) || $v < 0.0 || $v > 1.0) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binary/grouped-binomial y must be in [0,1]');
            }
        }
        foreach ($weights as $w) {
            if (fmod($w, 1.0) != 0.0) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binomial weights must be integer trial counts');
            }
        }
        foreach ($y as $i => $v) {
            $events = $v * $weights[$i];
            if (abs($events - Support::round($events)) > 1e-8) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'binary/grouped-binomial outcomes require integer trial counts and successes');
            }
        }

        return $y;
    }
}
