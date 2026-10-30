<?php

declare(strict_types=1);

namespace Phit\Polyfill;

/**
 * One-dimensional minimization and root finding, ported from the phit
 * extension (stats-core optimize.rs plus the PHP glue), which follows R's C
 * Brent_fmin (behind stats::optimize) and R_zeroin2 (behind stats::uniroot)
 * step for step, so the same function is evaluated at the same points.
 *
 * f is any callable taking one float and returning an int or float; an
 * exception thrown by f propagates unchanged. As in R, a nonfinite value of
 * f during the search is replaced (NAN and INF by PHP_FLOAT_MAX, and for
 * uniroot -INF by -PHP_FLOAT_MAX) with one warning per call; the returned
 * objective and f_root are f's own values.
 */
final class Optimize
{
    /** R's default tolerance, .Machine$double.eps^0.25. */
    private const float DEFAULT_TOL = 0.0001220703125;

    private const string REPLACED = 'nonfinite function values were replaced by the largest finite value';

    /**
     * Phit\optimize
     *
     * @return array{minimum: float, objective: float}|array{maximum: float, objective: float}
     */
    public static function optimize(callable $f, float $lower, float $upper, ?float $tol = null, bool $maximum = false): array
    {
        $name = 'Phit\optimize';
        self::interval($lower, $upper);
        $tol = self::tolerance($tol);
        $warned = false;
        $x = self::fmin($lower, $upper, $tol, static function (float $x) use ($f, $name, $maximum, &$warned): float {
            $y = self::call($f, $x, $name);
            $y = $maximum ? -$y : $y;
            if (is_finite($y)) {
                return $y;
            }
            self::replaced($name, $warned);

            return PHP_FLOAT_MAX;
        });
        $objective = self::call($f, $x, $name);

        return [($maximum ? 'maximum' : 'minimum') => $x, 'objective' => $objective];
    }

    /**
     * Phit\uniroot
     *
     * @return array{root: float, f_root: float, iter: int, estim_prec: float}
     */
    public static function uniroot(callable $f, float $lower, float $upper, ?float $tol = null, int $maxiter = 1000): array
    {
        $name = 'Phit\uniroot';
        self::interval($lower, $upper);
        $tol = self::tolerance($tol);
        if ($maxiter < 1) {
            throw self::invalid('maxiter must be at least 1');
        }
        // The extension evaluates both ends before it checks either result,
        // so f(upper) runs even when f(lower) returned a non-number, and an
        // exception thrown there takes precedence. An exception from
        // f(lower) stops before f(upper).
        $raw = $f($lower);
        $fUpper = self::call($f, $upper, $name);
        $fLower = self::number($raw, $name);
        if (is_nan($fLower) || is_nan($fUpper)) {
            throw self::invalid('f is NAN at lower or upper');
        }
        if (self::sign($fLower) * self::sign($fUpper) > 0.0) {
            throw self::invalid('f must have opposite signs at lower and upper');
        }
        $truncate = static fn (float $v): float => $v < -PHP_FLOAT_MAX ? -PHP_FLOAT_MAX : ($v > PHP_FLOAT_MAX ? PHP_FLOAT_MAX : $v);
        $warned = false;
        [$root, $iterations, $precision] = self::zeroin($lower, $upper, $truncate($fLower), $truncate($fUpper), $tol, $maxiter,
            static function (float $x) use ($f, $name, &$warned): float {
                $y = self::call($f, $x, $name);
                if (is_finite($y)) {
                    return $y;
                }
                self::replaced($name, $warned);

                return $y === -INF ? -PHP_FLOAT_MAX : PHP_FLOAT_MAX;
            });
        if ($iterations === null) {
            self::warn($name, "not converged in $maxiter iterations");
            $iterations = $maxiter;
        }
        $fRoot = self::call($f, $root, $name);

        return ['root' => $root, 'f_root' => $fRoot, 'iter' => $iterations, 'estim_prec' => $precision];
    }

    // ---- PHP glue -------------------------------------------------------

    private static function invalid(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    /** f at x, which must be an int or float. */
    private static function call(callable $f, float $x, string $name): float
    {
        return self::number($f($x), $name);
    }

    private static function number(mixed $y, string $name): float
    {
        if (! \is_int($y) && ! \is_float($y)) {
            throw self::invalid("$name(): f must return an int or float");
        }

        return (float) $y;
    }

    /**
     * Userland code cannot raise E_WARNING; E_USER_WARNING with the "fn(): "
     * prefix PHP adds to native warnings is the closest.
     */
    private static function warn(string $name, string $message): void
    {
        trigger_error("$name(): $message", E_USER_WARNING);
    }

    /** Warns once per search that nonfinite values were replaced, as R warns. */
    private static function replaced(string $name, bool &$warned): void
    {
        if (! $warned) {
            $warned = true;
            self::warn($name, self::REPLACED);
        }
    }

    private static function tolerance(?float $tol): float
    {
        $tol ??= self::DEFAULT_TOL;
        if (! is_finite($tol) || $tol <= 0.0) {
            throw self::invalid('tol must be a positive finite number');
        }

        return $tol;
    }

    private static function interval(float $lower, float $upper): void
    {
        if (! is_finite($lower) || ! is_finite($upper)) {
            throw self::invalid('lower and upper must be finite');
        }
        if ($lower >= $upper) {
            throw self::invalid('lower must be less than upper');
        }
    }

    private static function sign(float $v): float
    {
        return $v > 0.0 ? 1.0 : ($v < 0.0 ? -1.0 : 0.0);
    }

    // ---- numerical core -------------------------------------------------

    /**
     * The minimizer of f on [lower, upper] by Brent's golden-section and
     * parabolic-interpolation search, as R's Brent_fmin.
     *
     * @param callable(float): float $f
     */
    private static function fmin(float $lower, float $upper, float $tol, callable $f): float
    {
        // The squared inverse of the golden ratio.
        $c = (3.0 - sqrt(5.0)) * 0.5;
        $eps = sqrt(PHP_FLOAT_EPSILON);
        $a = $lower;
        $b = $upper;
        $v = $a + $c * ($b - $a);
        $w = $v;
        $x = $v;
        $d = 0.0;
        $e = 0.0;
        $fx = $f($x);
        $fv = $fx;
        $fw = $fx;
        $tol3 = $tol / 3.0;
        while (true) {
            $xm = ($a + $b) * 0.5;
            $tol1 = $eps * abs($x) + $tol3;
            $t2 = $tol1 * 2.0;
            if (abs($x - $xm) <= $t2 - ($b - $a) * 0.5) {
                break;
            }
            $p = 0.0;
            $q = 0.0;
            $r = 0.0;
            if (abs($e) > $tol1) {
                // Fit a parabola.
                $r = ($x - $w) * ($fx - $fv);
                $q = ($x - $v) * ($fx - $fw);
                $p = ($x - $v) * $q - ($x - $w) * $r;
                $q = ($q - $r) * 2.0;
                if ($q > 0.0) {
                    $p = -$p;
                } else {
                    $q = -$q;
                }
                $r = $e;
                $e = $d;
            }
            if (abs($p) >= abs($q * 0.5 * $r) || $p <= $q * ($a - $x) || $p >= $q * ($b - $x)) {
                // A golden-section step.
                $e = $x < $xm ? $b - $x : $a - $x;
                $d = $c * $e;
            } else {
                // A parabolic-interpolation step, kept away from the ends.
                $d = $p / $q;
                $u = $x + $d;
                if ($u - $a < $t2 || $b - $u < $t2) {
                    $d = $x >= $xm ? -$tol1 : $tol1;
                }
            }
            // f is never evaluated within tol1 of x.
            if (abs($d) >= $tol1) {
                $u = $x + $d;
            } elseif ($d > 0.0) {
                $u = $x + $tol1;
            } else {
                $u = $x - $tol1;
            }
            $fu = $f($u);
            if ($fu <= $fx) {
                if ($u < $x) {
                    $b = $x;
                } else {
                    $a = $x;
                }
                $v = $w;
                $w = $x;
                $x = $u;
                $fv = $fw;
                $fw = $fx;
                $fx = $fu;
            } else {
                if ($u < $x) {
                    $a = $u;
                } else {
                    $b = $u;
                }
                if ($fu <= $fw || $w == $x) {
                    $v = $w;
                    $fv = $fw;
                    $w = $u;
                    $fw = $fu;
                } elseif ($fu <= $fv || $v == $x || $v == $w) {
                    $v = $u;
                    $fv = $fu;
                }
            }
        }

        return $x;
    }

    /**
     * A root of f in [lower, upper], given f at both ends with opposite
     * signs, by Brent's method as R's R_zeroin2: the root, the iterations
     * used (null when maxiter ran out first) and |c - b| at the end.
     *
     * @param callable(float): float $f
     * @return array{float, ?int, float}
     */
    private static function zeroin(float $lower, float $upper, float $fLower, float $fUpper, float $tol, int $maxiter, callable $f): array
    {
        $a = $lower;
        $b = $upper;
        $fa = $fLower;
        $fb = $fUpper;
        $c = $a;
        $fc = $fa;
        if ($fa == 0.0) {
            return [$a, 0, 0.0];
        }
        if ($fb == 0.0) {
            return [$b, 0, 0.0];
        }
        for ($used = 0; $used <= $maxiter; $used++) {
            $prevStep = $b - $a;
            if (abs($fc) < abs($fb)) {
                // Make b the best approximation.
                $a = $b;
                $b = $c;
                $c = $a;
                $fa = $fb;
                $fb = $fc;
                $fc = $fa;
            }
            $tolAct = 2.0 * PHP_FLOAT_EPSILON * abs($b) + $tol / 2.0;
            $newStep = ($c - $b) / 2.0;
            if (abs($newStep) <= $tolAct || $fb == 0.0) {
                return [$b, $used, abs($c - $b)];
            }
            // Try interpolation when the previous step was large enough and
            // in the right direction.
            if (abs($prevStep) >= $tolAct && abs($fa) > abs($fb)) {
                $cb = $c - $b;
                if ($a == $c) {
                    // Two distinct points: linear interpolation.
                    $t1 = $fb / $fa;
                    $p = $cb * $t1;
                    $q = 1.0 - $t1;
                } else {
                    // Inverse quadratic interpolation.
                    $q = $fa / $fc;
                    $t1 = $fb / $fc;
                    $t2 = $fb / $fa;
                    $p = $t2 * ($cb * $q * ($q - $t1) - ($b - $a) * ($t1 - 1.0));
                    $q = ($q - 1.0) * ($t1 - 1.0) * ($t2 - 1.0);
                }
                if ($p > 0.0) {
                    $q = -$q;
                } else {
                    $p = -$p;
                }
                if ($p < 0.75 * $cb * $q - abs($tolAct * $q) / 2.0 && $p < abs($prevStep * $q / 2.0)) {
                    $newStep = $p / $q;
                }
            }
            if (abs($newStep) < $tolAct) {
                $newStep = $newStep > 0.0 ? $tolAct : -$tolAct;
            }
            $a = $b;
            $fa = $fb;
            $b += $newStep;
            $fb = $f($b);
            if (($fb > 0.0 && $fc > 0.0) || ($fb < 0.0 && $fc < 0.0)) {
                // Keep c on the other side of the root from b.
                $c = $a;
                $fc = $fa;
            }
        }

        return [$b, null, abs($c - $b)];
    }
}
