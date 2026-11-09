<?php

declare(strict_types=1);

namespace Phit\Polyfill\Formula;

use Phit\Polyfill\ErrorCode;
use Phit\Polyfill\Internal\ModelData;

/**
 * The PHP glue shared by Phit\model_matrix, Phit\fit_ols and Phit\fit_generalized_linear, in the
 * extension's order: argument checks, missing policy, data, weights, design.
 *
 * @internal
 */
final class Frontend
{
    private static function invalid(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    /** ext-php-rs rejects &str arguments that are not valid UTF-8 before the function runs. */
    public static function stringArgument(string $value, string $name): void
    {
        if (! mb_check_encoding($value, 'UTF-8')) {
            throw new \Exception("Invalid value given for argument `$name`.", 0);
        }
    }

    private static function missingPolicy(string $missing): bool
    {
        return match ($missing) {
            'drop' => true,
            'raise' => false,
            default => throw self::invalid("missing must be 'drop' or 'raise'"),
        };
    }

    /** @return ?list<?float> */
    private static function weights(mixed $weights, Data $data): ?array
    {
        if ($weights === null) {
            return null;
        }
        if (Data::utf8($weights)) {
            $column = $data->column($weights);
            if ($column === null) {
                throw self::invalid("unknown weights column '$weights'");
            }
            if ($column['t'] !== 'num') {
                throw self::invalid("weights column '$weights' must be numeric");
            }
            $values = $column['v'];
        } elseif (is_array($weights)) {
            $values = [];
            foreach ($weights as $w) {
                $values[] = ModelData::number($w);
            }
        } else {
            throw self::invalid('weights must be null, an array or a column name');
        }
        if (count($values) !== $data->rows) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'optional vectors must have one value per original observation');
        }
        foreach ($values as $w) {
            if ($w !== null && $w <= 0.0) {
                throw self::invalid('weights must be positive');
            }
        }

        return $values;
    }

    /** @return array{0: Design, 1: list<float>} the design and the retained weights */
    public static function prepare(string $formula, mixed $data, string $missing, mixed $weights, bool $dropUnusedLevels, bool $allowEmpty = false): array
    {
        $dropMissing = self::missingPolicy($missing);
        $frame = Data::from($data);
        $weights = self::weights($weights, $frame);
        $design = Design::build($formula, $frame, $dropMissing, $dropUnusedLevels, $weights === null ? [] : [$weights], $allowEmpty);

        return [$design, $design->extra[0] ?? array_fill(0, count($design->x), 1.0)];
    }

    public static function rejectRandom(Design $design, string $function): void
    {
        if ($design->random !== []) {
            throw ErrorCode::exception(
                ErrorCode::UNSUPPORTED,
                "$function does not fit random-effect terms such as ({$design->random[0]['label']}); use Phit\\fit_mixed_gaussian, or Phit\\fit_mixed_gaussian_matrix with Phit\\model_matrix()['random']",
            );
        }
    }

    /** @return list<float> */
    public static function numericResponse(Design $design, string $function): array
    {
        if ($design->response === null) {
            throw self::invalid("$function needs a response: 'y ~ ...'");
        }
        if ($design->response[0] !== 'numeric') {
            throw self::invalid("$function needs a numeric or logical response");
        }

        return $design->response[1];
    }

    /** Naming keys added to every formula fit. */
    public static function addNames(array $output, Design $design): array
    {
        $output['formula'] = $design->formula;
        $output['response'] = $design->responseName;
        $output['columns'] = $design->columns;
        $output['terms'] = $design->terms;
        $output['assign'] = $design->assign;
        $output['intercept'] = $design->intercept;
        $named = [];
        foreach ($design->columns as $j => $name) {
            if (! array_key_exists($j, $output['coefficients'] ?? [])) {
                break;
            }
            $named[$name] = $output['coefficients'][$j];
        }
        $output['named'] = $named;
        // Factor levels and poly() coefficients that Phit\predict() reuses.
        $output['predvars'] = $design->predvars;

        return $output;
    }

    public static function finish(array $output, Design $design, string $function, string $solver): array
    {
        return ModelData::accounting($design->nInput, $design->retainedRows, $design->droppedRows)
            ->finish($output, $function, $solver);
    }
}
