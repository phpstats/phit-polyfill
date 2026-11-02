<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\Internal\Bobyqa\Bobyqb;
use Phit\Polyfill\Internal\Bobyqa\BobyqaException;
use Phit\Polyfill\Internal\Bobyqa\Config;
use Phit\Polyfill\Internal\Bobyqa\Linalg;
use Phit\Polyfill\Internal\Bobyqa\Outcome;
use Phit\Polyfill\Internal\Bobyqa\RestartState;
use Phit\Polyfill\Internal\Bobyqa\Status;
use Phit\Polyfill\Internal\Bobyqa\Util;

/**
 * Follows the Rust crate bobyqa 0.3.2 (Powell's BOBYQA, after PRIMA),
 * which the extension's mixed models minimize with. It reproduces the crate's
 * evaluation sequence and result bit for bit: same loops, same summation
 * order, Rust's NaN-losing f64::max/min, and PRIMA's single-precision literals.
 *
 *     $config = new Config($n);          // npt = 2n + 1, rhoBegin 1.0, rhoEnd 1e-6, maxFun 500n
 *     $config->rhoBegin = 0.1;           // the other fields: fTarget, restart, ftolRel,
 *     $config->rhoEnd = $tolerance;      // ftolAbs, primaParity, as in the crate
 *     $config->maxFun = $maxEvaluations;
 *     $solver = new Bobyqa($n, $config); // throws BobyqaException (message: the crate's Display)
 *     $outcome = $solver->minimize($f, $x, $lower, $upper); // $x is overwritten with the best point
 *     $outcome->status;                  // Status enum; ->label() is the crate's Display text
 *     $outcome->nEval;
 *     $outcome->foundFinite();
 *
 * The objective receives a list of floats and may return INF or NAN for
 * infeasible points (they are moderated to FUNCMAX, as in PRIMA).
 *
 * @internal
 */
final class Bobyqa
{
    public const float FUNCMAX = Util::FUNCMAX;

    public const int MAX_RESTARTS_CAP = 10000;

    /** @var list<int> */
    private array $cycleBoundaries = [];

    private readonly Config $config;

    /**
     * @throws BobyqaException when (n, config) is rejected
     */
    public function __construct(private readonly int $n, Config $config)
    {
        $this->config = clone $config;
        if ($config->restart !== null) {
            $this->config->restart = clone $config->restart;
        }
        $c = $this->config;
        if ($n <= 0
            || $c->npt < $n + 2
            || $c->npt > intdiv(($n + 1) * ($n + 2), 2)
            || ! is_finite($c->rhoBegin) || is_nan($c->rhoBegin) || $c->rhoBegin <= 0.0
            || is_nan($c->rhoEnd) || $c->rhoEnd <= 0.0 || $c->rhoEnd > $c->rhoBegin
            || $c->maxFun < $c->npt + 1
            || is_nan($c->fTarget) || $c->fTarget === INF) {
            throw new BobyqaException(Status::InvalidArgs);
        }
        foreach ([$c->ftolRel, $c->ftolAbs] as $tol) {
            if ($tol !== null && (! is_finite($tol) || $tol < 0.0)) {
                throw new BobyqaException(Status::InvalidArgs);
            }
        }
        $restart = $c->restart;
        if ($restart !== null) {
            if ($restart->maxRestarts < 1 || $restart->maxRestarts > self::MAX_RESTARTS_CAP
                || is_nan($restart->improveRelTol) || $restart->improveRelTol < 0.0
                || is_nan($restart->cycleBudgetFrac) || $restart->cycleBudgetFrac < 0.0 || $restart->cycleBudgetFrac > 1.0) {
                throw new BobyqaException(Status::InvalidArgs);
            }
        }
    }

    /**
     * Minimize $f from $x (overwritten with the best point found) subject to
     * $lower <= x <= $upper. Invalid arguments return an InvalidArgs outcome.
     *
     * @param  callable(list<float>): (float|int)  $f
     * @param  list<float>  $x
     * @param  list<float>  $lower
     * @param  list<float>  $upper
     */
    public function minimize(callable $f, array &$x, array $lower, array $upper): Outcome
    {
        $this->cycleBoundaries = [];
        $n = $this->n;
        $c = $this->config;
        $x = array_map(floatval(...), array_values($x));
        $lower = array_map(floatval(...), array_values($lower));
        $upper = array_map(floatval(...), array_values($upper));
        if (! $this->argsAreValid($x, $lower, $upper)) {
            return new Outcome(NAN, 0, Status::InvalidArgs);
        }

        $xl = [];
        $xu = [];
        for ($i = 0; $i < $n; $i++) {
            $xl[] = Linalg::max($lower[$i], -Util::BOUNDMAX);
            $xu[] = Linalg::min($upper[$i], Util::BOUNDMAX);
        }
        for ($i = 0; $i < $n; $i++) {
            $x[$i] = Linalg::max($xl[$i], Linalg::min($xu[$i], Util::moderatex1($x[$i])));
        }
        Util::reviseX0($x, $xl, $xu, $c->rhoBegin);

        $restart = $c->restart !== null ? new RestartState($c->restart) : null;
        $ftol = ($c->ftolRel === null && $c->ftolAbs === null)
            ? null
            : [$c->ftolRel ?? 0.0, $c->ftolAbs ?? 0.0];

        [$fopt, $nf, $info] = Bobyqb::solve(
            $f,
            $c->maxFun,
            $c->npt,
            $c->fTarget,
            $ftol,
            $c->rhoBegin,
            $c->rhoEnd,
            $x,
            $xl,
            $xu,
            $restart,
            ! $c->primaParity && $c->npt === intdiv(($n + 1) * ($n + 2), 2),
        );
        if ($restart !== null) {
            $this->cycleBoundaries = $restart->cycleBoundaries;
        }

        return new Outcome($fopt, $nf, Status::fromInfo($info));
    }

    /** Restarts performed on the last minimize call. */
    public function lastRestartCount(): int
    {
        return count($this->cycleBoundaries);
    }

    /**
     * Cumulative evaluation count at each restart boundary of the last minimize call.
     *
     * @return list<int>
     */
    public function lastCycleBoundaries(): array
    {
        return $this->cycleBoundaries;
    }

    /**
     * Lengths, NaN-free bounds and x, and bounds (after the BOUNDMAX clamp)
     * ordered and at least 2 * rhoBegin apart.
     *
     * @param  list<float>  $x
     * @param  list<float>  $lower
     * @param  list<float>  $upper
     */
    private function argsAreValid(array $x, array $lower, array $upper): bool
    {
        $n = $this->n;
        if (count($x) !== $n || count($lower) !== $n || count($upper) !== $n) {
            return false;
        }
        for ($i = 0; $i < $n; $i++) {
            if (is_nan($lower[$i]) || is_nan($upper[$i]) || is_nan($x[$i])) {
                return false;
            }
        }
        for ($i = 0; $i < $n; $i++) {
            $cl = Linalg::max($lower[$i], -Util::BOUNDMAX);
            $cu = Linalg::min($upper[$i], Util::BOUNDMAX);
            if ($cl > $cu || $cu - $cl < 2.0 * $this->config->rhoBegin) {
                return false;
            }
        }

        return true;
    }
}
