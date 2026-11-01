<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Design;
use Phit\Polyfill\Formula\Frontend;
use Phit\Polyfill\Internal\LinearMixed;
use Phit\Polyfill\Internal\LogisticRandomIntercept;
use Phit\Polyfill\Internal\ModelData;

/**
 * lme4-style formula interfaces over the mixed-model cores: Phit\fit_mixed_gaussian and
 * Phit\fit_mixed_logistic (the extension's formula.rs).
 */
final class MixedFormula
{
    /** Phit\fit_mixed_gaussian */
    public static function lmer(
        string $formula,
        mixed $data,
        bool $reml = true,
        string $missing = 'drop',
        mixed $weights = null,
        ?array $start = null,
        int $maxEvaluations = 20000,
        float $tolerance = 1e-7,
        string $optimizer = 'bobyqa',
        float $objectiveTolerance = 1e-12,
    ): array {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($missing, 'missing');
        Frontend::stringArgument($optimizer, 'optimizer');
        if ($optimizer !== 'bobyqa' && $optimizer !== 'nelder_mead') {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, "optimizer must be 'bobyqa' or 'nelder_mead'");
        }
        $startValues = null;
        if ($start !== null) {
            $startValues = [];
            foreach ($start as $value) {
                $value = ModelData::number($value);
                if ($value === null) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'start cannot contain missing values');
                }
                $startValues[] = $value;
            }
        }
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        self::requireRandom($design, 'Phit\fit_mixed_gaussian');
        $y = Frontend::numericResponse($design, 'Phit\fit_mixed_gaussian');
        $terms = [];
        $termLabels = [];
        foreach ($design->random as $term) {
            [$indices, $labels] = ModelData::groupIndices($term['groups']);
            $termLabels[] = $labels;
            $terms[] = ['z' => $term['z'], 'groups' => $indices, 'correlated' => $term['correlated']];
        }
        if ($maxEvaluations < 0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'maxEvaluations must be positive');
        }
        $result = LinearMixed::fit($design->x, $y, $terms, [
            'reml' => $reml,
            'optimizer' => $optimizer,
            'objective_tolerance' => $objectiveTolerance,
            'offset' => $design->offset ?? array_fill(0, count($y), 0.0),
            'weights' => $w,
            'start' => $startValues,
            'max_evaluations' => $maxEvaluations,
            'tolerance' => $tolerance,
        ]);
        $output = MixedModels::lmmOutput($result, $reml, $termLabels);
        foreach ($design->random as $i => $term) {
            $output['random'][$i]['label'] = $term['label'];
            $output['random'][$i]['group'] = $term['group'];
            $output['random'][$i]['columns'] = $term['columns'];
            $output['random'][$i]['correlated'] = $term['correlated'];
        }
        $output = Frontend::addNames($output, $design);

        return Frontend::finish($output, $design, 'Phit\fit_mixed_gaussian', $optimizer);
    }

    /** Phit\fit_glmer, deprecated: Phit\fit_mixed_logistic with a family and link that must be binomial logit. */
    public static function glmer(
        string $formula,
        mixed $data,
        string $family = 'binomial',
        string $link = '',
        string $missing = 'drop',
        mixed $weights = null,
        int $quadraturePoints = 1,
        int $maxEvaluations = 20000,
        float $tolerance = 1e-7,
    ): array {
        Frontend::stringArgument($family, 'family');
        Frontend::stringArgument($link, 'link');
        if ($family !== 'binomial' || ($link !== '' && $link !== 'logit')) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, "Phit\\fit_glmer supports family 'binomial' with the logit link only (got '$family'/'$link')");
        }

        return self::mixedLogistic($formula, $data, $missing, $weights, $quadraturePoints, $maxEvaluations, $tolerance);
    }

    /** Phit\fit_mixed_logistic */
    public static function mixedLogistic(
        string $formula,
        mixed $data,
        string $missing = 'drop',
        mixed $weights = null,
        int $quadraturePoints = 1,
        int $maxEvaluations = 20000,
        float $tolerance = 1e-7,
    ): array {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($missing, 'missing');
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        self::requireRandom($design, 'Phit\fit_mixed_logistic');
        if (count($design->random) !== 1 || $design->random[0]['columns'] !== ['(Intercept)']) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'Phit\fit_mixed_logistic supports exactly one random intercept, (1 | g); random slopes and several terms are not implemented');
        }
        $y = Glm::binomialResponse($design, $w);
        [$indices, $labels] = ModelData::groupIndices($design->random[0]['groups']);
        if ($quadraturePoints < 0) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'quadraturePoints must be positive');
        }
        if ($maxEvaluations < 0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'maxEvaluations must be positive');
        }
        $r = LogisticRandomIntercept::fit($design->x, $y, $indices, [
            'offset' => $design->offset ?? array_fill(0, count($y), 0.0),
            'weights' => $w,
            'quadrature_points' => $quadraturePoints,
            'max_evaluations' => $maxEvaluations,
            'tolerance' => $tolerance,
        ]);
        $output = MixedModels::logisticMixedOutput($r, $labels);
        $output['random'] = [[
            'label' => $design->random[0]['label'],
            'group' => $design->random[0]['group'],
            'columns' => $design->random[0]['columns'],
        ]];
        $output = Frontend::addNames($output, $design);

        return Frontend::finish($output, $design, 'Phit\fit_mixed_logistic', 'bobyqa');
    }

    private static function requireRandom(Design $design, string $function): void
    {
        if ($design->random === []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "$function needs at least one random-effect term such as (1 | g); use Phit\\fit_ols or Phit\\fit_generalized_linear for fixed effects only");
        }
    }
}
