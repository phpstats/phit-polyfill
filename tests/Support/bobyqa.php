<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Bobyqa;
use Phit\Polyfill\Internal\Bobyqa\Config;
use Phit\Polyfill\Internal\Bobyqa\RestartConfig;

/*
 * Runs the cases of tests/fixtures/bobyqa-reference.json, shared by
 * tests/Unit/BobyqaTest.php and bin/jit-off-experiment.
 */

function bobyqaFloat(string $hex): float
{
    return unpack('E', hex2bin($hex))[1];
}

function bobyqaHex(float $value): string
{
    return bin2hex(pack('E', $value));
}

/** The harness's objective, with its exact operation order. */
function bobyqaObjective(array $p, array $x): float
{
    $n = $p['n'];
    if ($p['barrier'] !== 0) {
        $s = 0.0;
        for ($i = 0; $i < $n; $i++) {
            $s += $x[$i];
        }
        if ($s > $p['thr']) {
            return $p['barrier'] === 1 ? INF : NAN;
        }
    }
    $quad = function (array $x) use ($p, $n): float {
        $f = $p['c'];
        for ($i = 0; $i < $n; $i++) {
            $f += $p['g'][$i] * $x[$i];
        }
        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $n; $j++) {
                $f += 0.5 * $p['h'][$i * $n + $j] * $x[$i] * $x[$j];
            }
        }

        return $f;
    };

    switch ($p['kind']) {
        case 'rosenbrock':
            if ($n === 1) {
                $b = 1.0 - $x[0];

                return $b * $b;
            }
            $f = 0.0;
            for ($i = 0; $i < $n - 1; $i++) {
                $a = $x[$i + 1] - $x[$i] * $x[$i];
                $b = 1.0 - $x[$i];
                $f += 100.0 * $a * $a + $b * $b;
            }

            return $f;
        case 'quad':
            return $quad($x);
        case 'quartic':
            $f = $quad($x);
            for ($i = 0; $i < $n; $i++) {
                $t = $x[$i] - $p['b'][$i];
                $t2 = $t * $t;
                $f += $p['a'][$i] * $t2 * $t2;
            }

            return $f;
        case 'l1':
            $f = $quad($x);
            for ($i = 0; $i < $n; $i++) {
                $f += $p['a'][$i] * abs($x[$i] - $p['b'][$i]);
            }

            return $f;
        case 'step':
            return floor($quad($x) * $p['a'][0]) / $p['a'][0];
        default:
            return INF;
    }
}

function bobyqaProblem(array $row): array
{
    $p = $row;
    foreach (['c', 'thr', 'rho_begin', 'rho_end', 'f_target'] as $key) {
        $p[$key] = bobyqaFloat($row[$key]);
    }
    foreach (['g', 'h', 'a', 'b', 'x0', 'lower', 'upper'] as $key) {
        $p[$key] = array_map(bobyqaFloat(...), $row[$key]);
    }
    $config = new Config($p['n']);
    $config->npt = $p['npt'];
    $config->rhoBegin = $p['rho_begin'];
    $config->rhoEnd = $p['rho_end'];
    $config->maxFun = $p['max_fun'];
    $config->fTarget = $p['f_target'];
    $config->ftolRel = $row['ftol_rel'] === null ? null : bobyqaFloat($row['ftol_rel']);
    $config->ftolAbs = $row['ftol_abs'] === null ? null : bobyqaFloat($row['ftol_abs']);
    $config->primaParity = $row['prima_parity'];
    if ($row['restart'] !== null) {
        $restart = new RestartConfig;
        $restart->maxRestarts = $row['restart']['max_restarts'];
        $restart->improveRelTol = bobyqaFloat($row['restart']['improve_rel_tol']);
        $restart->stallReductions = $row['restart']['stall_reductions'];
        $restart->cycleBudgetFrac = bobyqaFloat($row['restart']['cycle_budget_frac']);
        $config->restart = $restart;
    }

    return [$p, $config];
}

/** One reference case: the outcome, the final point, the trace, and the solver. */
function bobyqaRun(array $row): array
{
    [$p, $config] = bobyqaProblem($row);
    $solver = new Bobyqa($p['n'], $config);
    $trace = '';
    $x = $p['x0'];
    $outcome = $solver->minimize(function (array $point) use (&$trace, $p): float {
        $trace .= implode(',', array_map(bobyqaHex(...), $point)).';';

        return bobyqaObjective($p, $point);
    }, $x, $p['lower'], $p['upper']);

    return [$outcome, $x, $trace, $solver];
}
