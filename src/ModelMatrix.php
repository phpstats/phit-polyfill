<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Design;
use Phit\Polyfill\Formula\Frontend;

/**
 * R's model.matrix() for a formula over named columns (Phit\model_matrix).
 */
final class ModelMatrix
{
    public static function build(string $formula, mixed $data, string $missing = 'drop'): array
    {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($missing, 'missing');
        [$design] = Frontend::prepare($formula, $data, $missing, null, false);
        $out = [];
        $response = $design->response;
        if ($response === null) {
            $out['y'] = null;
        } elseif ($response[0] === 'numeric') {
            $out['y'] = $response[1];
        } elseif ($response[0] === 'factor') {
            $out['y'] = array_map(fn (int $c) => (float) $c, $response[2]);
            $out['response_levels'] = $response[1];
        } else {
            $out['y'] = array_map(fn ($s, $f) => [$s, $f], $response[1], $response[2]);
        }
        $out['x'] = $design->x;
        $out['columns'] = $design->columns;
        $out['intercept'] = $design->intercept;
        $out['terms'] = $design->terms;
        $out['assign'] = $design->assign;
        $out['formula'] = $design->formula;
        $out['response'] = $design->responseName;
        $out['offset'] = $design->offset;
        $out['random'] = $design->random;
        $dropped = count($design->droppedRows);
        $out['n_input'] = $design->nInput;
        $out['n_observations'] = count($design->retainedRows);
        $out['n_dropped'] = $dropped;
        $out['retained_rows'] = $design->retainedRows;
        $out['dropped_rows'] = $design->droppedRows;
        if ($dropped > 0) {
            // E_USER_WARNING with the "fn(): " prefix PHP adds to native warnings.
            trigger_error(
                "Phit\\model_matrix(): dropped $dropped of {$design->nInput} observations with missing values (listwise deletion)",
                E_USER_WARNING,
            );
        }

        return $out;
    }
}
