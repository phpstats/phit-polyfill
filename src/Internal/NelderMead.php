<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

use Phit\Polyfill\Distribution\Support;

/**
 * argmin 0.11.0's NelderMead solver driven by its Executor, over Vec<f64>
 * parameters with argmin-math 0.5.1's "vec" operations, reproduced bit for
 * bit (verified against the crate by tests/reference/neldermead).
 *
 * Mapping from the Rust used in stats-core's lmm.rs:
 *
 *   NelderMead::new(simplex)                  NelderMead::new($simplex)
 *     .with_sd_tolerance(tol)?                  ->withSdTolerance($tol)        throws NelderMeadError
 *   Executor::new(problem, solver)
 *     .configure(|s| s.max_iters(m)).run()?   ->run($cost, $m)               throws NelderMeadError
 *   state.get_termination_reason()            $result->termination           NelderMeadTermination
 *   state.get_best_cost() / get_best_param()  $result->bestCost / bestParam
 *   CostFunction::cost -> Err(e)              $cost throws NelderMeadError
 *
 * A NelderMeadError's message is the Rust error's Display text, so
 * "...failed: {e}" becomes "...failed: ".$e->getMessage(). A cost error during
 * the initial evaluation of the simplex is unwrap()ped by the crate and panics;
 * that, and the crate's index/assert panics on a malformed simplex, surface as
 * Distribution\RustPanic with the Rust panic message.
 *
 * Only the default alpha, gamma, rho and sigma are supported, and no observers,
 * checkpoints, timeouts or ctrlc (all off in stats-core). The vertex sort is
 * std's stable sort, which is insertion sort up to 20 vertices. Above that std
 * uses driftsort: without NaN costs any stable sort gives the same order, but
 * with one (NaN compares Equal to everything, which is no total order)
 * driftsort may order differently or panic ("user-provided comparison function
 * does not correctly implement a total order"); this code keeps insertion sort.
 *
 * @internal
 */
final class NelderMead
{
    private const float ALPHA = 1.0;

    private const float GAMMA = 2.0;

    private const float RHO = 0.5;

    private const float SIGMA = 0.5;

    /** @var list<list<float>> */
    private array $params;

    /** @var list<float> */
    private array $costs;

    /** @var \Closure(list<float>): float */
    private \Closure $cost;

    /**
     * @param  list<list<float>>  $simplex
     */
    private function __construct(private readonly array $simplex, private readonly float $sdTolerance) {}

    /**
     * @param  list<list<float>>  $simplex  n + 1 vertices of n coordinates
     */
    public static function new(array $simplex): self
    {
        $vertices = [];
        foreach ($simplex as $vertex) {
            $vertices[] = array_map(static fn (int|float $v): float => (float) $v, array_values($vertex));
        }

        return new self($vertices, PHP_FLOAT_EPSILON);
    }

    public function withSdTolerance(float $tolerance): self
    {
        if ($tolerance < 0.0) {
            throw new NelderMeadError(self::invalidParameter('`Nelder-Mead`: sd_tolerance must be >= 0.'));
        }

        return new self($this->simplex, $tolerance);
    }

    /**
     * Executor::run(). $maxIters stands in for the u64 max_iters (default
     * u64::MAX, for which PHP_INT_MAX is equivalent).
     *
     * @param  callable(list<float>): float  $cost
     *
     * @throws NelderMeadError
     */
    public function run(callable $cost, int $maxIters = PHP_INT_MAX): NelderMeadResult
    {
        if ($maxIters < 0) {
            throw new \InvalidArgumentException('max_iters is a u64');
        }
        $this->cost = static fn (array $x): float => $cost($x);
        $this->params = $this->simplex;
        $this->costs = [];

        // Solver::init evaluates with problem.cost(p).unwrap().
        foreach ($this->params as $i => $p) {
            try {
                $this->costs[$i] = ($this->cost)($p);
            } catch (NelderMeadError $e) {
                Support::panic('called `Result::unwrap()` on an `Err` value: '.$e->getMessage());
            }
        }
        $this->sort();
        if ($this->params === []) {
            Support::panic('index out of bounds: the len is 0 but the index is 0');
        }

        // IterState::new(): best cost +INF, no best param, target cost -INF.
        $bestParam = null;
        $bestCost = INF;
        $iter = 0;
        [$bestParam, $bestCost] = self::update($this->params[0], $this->costs[0], $bestParam, $bestCost);

        while (true) {
            if ($this->converged()) {
                $reason = NelderMeadTermination::SolverConverged;
            } elseif ($iter >= $maxIters) {
                $reason = NelderMeadTermination::MaxItersReached;
            } elseif ($bestCost <= -INF) {
                $reason = NelderMeadTermination::TargetCostReached;
            } else {
                $reason = null;
            }
            if ($reason !== null) {
                return new NelderMeadResult($reason, $bestParam, $bestCost, $iter);
            }
            $this->nextIter();
            [$bestParam, $bestCost] = self::update($this->params[0], $this->costs[0], $bestParam, $bestCost);
            $iter++;
        }
    }

    /**
     * IterState::update(): strictly better, or both infinite with the same sign.
     *
     * @param  list<float>  $param
     * @param  list<float>|null  $bestParam
     * @return array{0: list<float>|null, 1: float}
     */
    private static function update(array $param, float $cost, ?array $bestParam, float $bestCost): array
    {
        if ($cost < $bestCost || (is_infinite($cost) && is_infinite($bestCost) && ($cost > 0.0) === ($bestCost > 0.0))) {
            return [$param, $cost];
        }

        return [$bestParam, $bestCost];
    }

    /** NelderMead::terminate: sample sd of the vertex costs below sd_tolerance. */
    private function converged(): bool
    {
        $n = (float) count($this->costs);
        // f64's Sum folds from -0.0.
        $sum = -0.0;
        foreach ($this->costs as $c) {
            $sum += $c;
        }
        $c0 = $sum / $n;
        $squares = -0.0;
        foreach ($this->costs as $c) {
            $d = $c - $c0;
            $squares += $d * $d;
        }
        $s = sqrt(fdiv(1.0, $n - 1.0) * $squares);

        return $s < $this->sdTolerance;
    }

    private function nextIter(): void
    {
        $len = count($this->params);
        $x0 = $this->centroid();
        if ($len < 2) {
            // params[len - 2] wraps to usize::MAX in a release build.
            Support::panic('index out of bounds: the len is 1 but the index is 18446744073709551615');
        }
        $last = $len - 1;
        $bestCost = $this->costs[0];
        $worst = $this->params[$last];
        $worstCost = $this->costs[$last];
        $secondWorstCost = $this->costs[$len - 2];

        $xr = self::add($x0, self::scale(self::sub($x0, $worst), self::ALPHA));
        $xrCost = ($this->cost)($xr);

        if ($xrCost < $secondWorstCost && $xrCost >= $bestCost) {
            $this->params[$last] = $xr;
            $this->costs[$last] = $xrCost;
        } elseif ($xrCost < $bestCost) {
            $xe = self::add($x0, self::scale(self::sub($xr, $x0), self::GAMMA));
            $xeCost = ($this->cost)($xe);
            if ($xeCost < $xrCost) {
                $this->params[$last] = $xe;
                $this->costs[$last] = $xeCost;
            } else {
                $this->params[$last] = $xr;
                $this->costs[$last] = $xrCost;
            }
        } elseif ($xrCost >= $secondWorstCost) {
            if ($xrCost < $worstCost) {
                $xc = self::add($x0, self::scale(self::sub($xr, $x0), self::RHO));
                $xcCost = ($this->cost)($xc);
                if ($xcCost <= $xrCost) {
                    $this->params[$last] = $xc;
                    $this->costs[$last] = $xcCost;
                } else {
                    $this->shrink();
                }
            } else {
                $xc = self::add($x0, self::scale(self::sub($worst, $x0), self::RHO));
                $xcCost = ($this->cost)($xc);
                if ($xcCost < $worstCost) {
                    $this->params[$last] = $xc;
                    $this->costs[$last] = $xcCost;
                } else {
                    $this->shrink();
                }
            }
        } else {
            // Only a NaN reflected cost gets here.
            throw new NelderMeadError(
                'Potential bug: "`NelderMead`: Reached unreachable point.". This is potentially a bug. '
                .'Please file a report on https://github.com/argmin-rs/argmin/issues',
            );
        }

        $this->sort();
    }

    /**
     * Sum of all vertices but the worst, accumulated into the best in order,
     * then multiplied by the reciprocal of their count.
     *
     * @return list<float>
     */
    private function centroid(): array
    {
        $count = count($this->params) - 1;
        $acc = $this->params[0];
        for ($k = 1; $k < $count; $k++) {
            $acc = self::add($acc, $this->params[$k]);
        }

        return self::scale($acc, fdiv(1.0, (float) $count));
    }

    private function shrink(): void
    {
        $x0 = $this->params[0];
        $len = count($this->params);
        for ($k = 1; $k < $len; $k++) {
            $p = self::add($x0, self::scale(self::sub($this->params[$k], $x0), self::SIGMA));
            $this->params[$k] = $p;
            $this->costs[$k] = ($this->cost)($p);
        }
    }

    /**
     * slice::sort_by(partial_cmp, NaN => Equal): std's stable sort, which for
     * up to 20 elements is insertion sort with is_less = (a < b).
     */
    private function sort(): void
    {
        $len = count($this->costs);
        for ($i = 1; $i < $len; $i++) {
            $cost = $this->costs[$i];
            $param = $this->params[$i];
            $j = $i;
            while ($j > 0 && $cost < $this->costs[$j - 1]) {
                $this->costs[$j] = $this->costs[$j - 1];
                $this->params[$j] = $this->params[$j - 1];
                $j--;
            }
            $this->costs[$j] = $cost;
            $this->params[$j] = $param;
        }
    }

    /**
     * @param  list<float>  $a
     * @param  list<float>  $b
     * @return list<float>
     */
    private static function add(array $a, array $b): array
    {
        self::assertSameLength($a, $b);
        $out = [];
        foreach ($a as $i => $v) {
            $out[] = $v + $b[$i];
        }

        return $out;
    }

    /**
     * @param  list<float>  $a
     * @param  list<float>  $b
     * @return list<float>
     */
    private static function sub(array $a, array $b): array
    {
        self::assertSameLength($a, $b);
        $out = [];
        foreach ($a as $i => $v) {
            $out[] = $v - $b[$i];
        }

        return $out;
    }

    /**
     * @param  list<float>  $a
     * @return list<float>
     */
    private static function scale(array $a, float $factor): array
    {
        $out = [];
        foreach ($a as $v) {
            $out[] = $v * $factor;
        }

        return $out;
    }

    /**
     * argmin-math's Vec ⊕ Vec preconditions, with their panic messages.
     *
     * @param  list<float>  $a
     * @param  list<float>  $b
     */
    private static function assertSameLength(array $a, array $b): void
    {
        $n1 = count($a);
        $n2 = count($b);
        if ($n1 === 0) {
            Support::panic('assertion failed: n1 > 0');
        }
        if ($n2 === 0) {
            Support::panic('assertion failed: n2 > 0');
        }
        if ($n1 !== $n2) {
            Support::panic("assertion `left == right` failed\n  left: $n1\n right: $n2");
        }
    }

    /** Display of ArgminError::InvalidParameter: the text is Debug-quoted. */
    private static function invalidParameter(string $text): string
    {
        return 'Invalid parameter: "'.addcslashes($text, '"\\').'"';
    }
}
