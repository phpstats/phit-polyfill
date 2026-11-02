<?php

declare(strict_types=1);

namespace Phit\Polyfill;

/**
 * One-dimensional minimization and root finding, as R's stats::optimize and
 * stats::uniroot, following the phit extension: stats-core optimize.rs (Brent's
 * methods, written from the clean-room specification specs/cleanroom/brent.md
 * and netlib's public-domain fmin.f and zeroin.f) and the PHP glue in
 * php-extension optimize.rs.
 *
 * f is any PHP callable taking one float and returning an int or float. An
 * exception thrown by f stops the search and propagates unchanged. A nonfinite
 * value of f during the search is replaced (NAN and INF by PHP_FLOAT_MAX, and
 * for uniroot -INF by -PHP_FLOAT_MAX) with one warning per call; the returned
 * objective and f_root are f's own values.
 */
final class Optimize
{
    /** R's default tolerance for optimize() and uniroot(): .Machine$double.eps^0.25. */
    private const float DEFAULT_TOL = 0.0001220703125;

    /** f64::EPSILON. */
    private const float EPSILON = 2.220446049250313e-16;

    private const string REPLACED = 'nonfinite function values were replaced by the largest finite value';

    /**
     * Minimize (or maximize) f on [lower, upper], as R's optimize().
     *
     * @return array{minimum?: float, maximum?: float, objective: float}
     */
    public static function optimize(callable $f, float $lower, float $upper, ?float $tol = null, bool $maximum = false): array
    {
        $name = 'Phit\optimize';
        self::interval($lower, $upper);
        $tol = self::tolerance($tol);
        $warned = false;
        $x = self::fmin($lower, $upper, $tol, static function (float $x) use ($f, $name, $maximum, &$warned): float {
            $y = self::number($f, $x, $name);
            $y = $maximum ? -$y : $y;
            if (is_finite($y)) {
                return $y;
            }
            self::note($warned, $name);

            return PHP_FLOAT_MAX;
        });
        $objective = self::number($f, $x, $name);

        return [$maximum ? 'maximum' : 'minimum' => $x, 'objective' => $objective];
    }

    /**
     * A root of f in [lower, upper], where f changes sign, as R's uniroot().
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
        // The extension evaluates both ends before looking at either result. An
        // exception from f at lower propagates at once; a non-numeric value there
        // is reported only after f has been called at upper.
        $fLower = self::call($f, $lower);
        $fUpper = self::call($f, $upper);
        if ($fLower === null || $fUpper === null) {
            throw self::badReturn($name);
        }
        if (is_nan($fLower) || is_nan($fUpper)) {
            throw self::invalid('f is NAN at lower or upper');
        }
        if (self::sign($fLower) * self::sign($fUpper) > 0.0) {
            throw self::invalid('f must have opposite signs at lower and upper');
        }
        $warned = false;
        [$root, $iterations, $precision] = self::zeroin(
            $lower,
            $upper,
            self::truncate($fLower),
            self::truncate($fUpper),
            $tol,
            $maxiter,
            static function (float $x) use ($f, $name, &$warned): float {
                $y = self::number($f, $x, $name);
                if (is_finite($y)) {
                    return $y;
                }
                self::note($warned, $name);

                return $y === -INF ? -PHP_FLOAT_MAX : PHP_FLOAT_MAX;
            },
        );
        if ($iterations === null) {
            trigger_error("$name(): not converged in $maxiter iterations", E_USER_WARNING);
            $iterations = $maxiter;
        }
        $fRoot = self::number($f, $root, $name);

        return ['root' => $root, 'f_root' => $fRoot, 'iter' => $iterations, 'estim_prec' => $precision];
    }

    // ---- PHP glue -------------------------------------------------------

    private static function invalid(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    private static function badReturn(string $name): \Exception
    {
        return self::invalid("$name(): f must return an int or float");
    }

    /** f(x) as a float, or null when f returned something other than an int or float. */
    private static function call(callable $f, float $x): ?float
    {
        $value = $f($x);
        if (\is_float($value)) {
            return $value;
        }
        if (\is_int($value)) {
            return (float) $value;
        }

        return null;
    }

    /** f(x) as a float; anything other than an int or float is an error. */
    private static function number(callable $f, float $x, string $name): float
    {
        return self::call($f, $x) ?? throw self::badReturn($name);
    }

    /** Warns once per search that nonfinite values were replaced. */
    private static function note(bool &$warned, string $name): void
    {
        if (! $warned) {
            $warned = true;
            trigger_error("$name(): ".self::REPLACED, E_USER_WARNING);
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
        if ($v > 0.0) {
            return 1.0;
        }
        if ($v < 0.0) {
            return -1.0;
        }

        return 0.0;
    }

    /** v clamped to [-PHP_FLOAT_MAX, PHP_FLOAT_MAX]; v is not NAN. */
    private static function truncate(float $v): float
    {
        if ($v < -PHP_FLOAT_MAX) {
            return -PHP_FLOAT_MAX;
        }
        if ($v > PHP_FLOAT_MAX) {
            return PHP_FLOAT_MAX;
        }

        return $v;
    }

    // ---- numerical core -------------------------------------------------

    /**
     * A minimizer of f on [lower, upper] by Brent's method (golden section with
     * parabolic interpolation). Returns the final x. An exception from f stops
     * the search.
     *
     * @param callable(float): float $f
     */
    private static function fmin(float $lower, float $upper, float $tol, callable $f): float
    {
        $golden = 0.5 * (3.0 - sqrt(5.0));
        $eps = sqrt(self::EPSILON);
        $third = $tol / 3.0;

        $lo = $lower;
        $hi = $upper;
        $x = $lo + $golden * ($hi - $lo);
        $w = $x;
        $v = $x;
        $step = 0.0; // d
        $prevStep = 0.0; // e
        $fx = $f($x);
        $fw = $fx;
        $fv = $fx;

        while (true) {
            $mid = 0.5 * ($lo + $hi);
            $tol1 = $eps * abs($x) + $third;
            $tol2 = 2.0 * $tol1;
            if (abs($x - $mid) <= $tol2 - 0.5 * ($hi - $lo)) {
                return $x;
            }

            $p = 0.0;
            $q = 0.0;
            $r = 0.0;
            if (abs($prevStep) > $tol1) {
                $r = ($x - $w) * ($fx - $fv);
                $q = ($x - $v) * ($fx - $fw);
                $p = ($x - $v) * $q - ($x - $w) * $r;
                $q = 2.0 * ($q - $r);
                if ($q > 0.0) {
                    $p = -$p;
                } else {
                    $q = -$q;
                }
                $r = $prevStep;
                $prevStep = $step;
            }

            // Parabolic interpolation is accepted unless one of the three tests
            // holds; a NAN makes every test false, and so accepts it.
            $parabolic = true;
            if (abs($p) >= abs(0.5 * $q * $r) || $p <= $q * ($lo - $x) || $p >= $q * ($hi - $x)) {
                $parabolic = false;
            }
            if ($parabolic) {
                $step = fdiv($p, $q);
                $u = $x + $step;
                $below = $u - $lo;
                $above = $hi - $u;
                // Not (both at least tol2), with NAN counting as not.
                if ($below < $tol2 || $above < $tol2 || is_nan($below) || is_nan($above)) {
                    $step = $x >= $mid ? -$tol1 : $tol1;
                }
            } else {
                $prevStep = $x >= $mid ? $lo - $x : $hi - $x;
                $step = $golden * $prevStep;
            }

            if (abs($step) >= $tol1) {
                $u = $x + $step;
            } elseif ($step <= 0.0) {
                $u = $x - $tol1;
            } else {
                $u = $x + $tol1;
            }
            $fu = $f($u);

            if ($fu <= $fx) {
                if ($u >= $x) {
                    $lo = $x;
                } else {
                    $hi = $x;
                }
                $v = $w;
                $fv = $fw;
                $w = $x;
                $fw = $fx;
                $x = $u;
                $fx = $fu;
            } else {
                if ($u >= $x) {
                    $hi = $u;
                } else {
                    $lo = $u;
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
    }

    /**
     * A root of f in [lower, upper] by Brent's method, given f(lower) and
     * f(upper), which have opposite signs or one of which is zero. Returns the
     * root, the evaluations made before convergence (null when the search
     * stopped without converging) and the final bracket width |c - b|.
     *
     * @param callable(float): float $f
     * @return array{float, ?int, float}
     */
    private static function zeroin(float $lower, float $upper, float $fLower, float $fUpper, float $tol, int $maxiter, callable $f): array
    {
        if ($fLower == 0.0) {
            return [$lower, 0, 0.0];
        }
        if ($fUpper == 0.0) {
            return [$upper, 0, 0.0];
        }

        $eps = self::EPSILON;
        $maxEvals = $maxiter + 1;
        $evals = 0;

        // b is the current best estimate, a the previous one, c the opposite end of the bracket.
        $a = $lower;
        $b = $upper;
        $fa = $fLower;
        $fb = $fUpper;
        $c = $a;
        $fc = $fa;

        while (true) {
            // Keep b the point with the smaller |f|.
            if (abs($fc) < abs($fb)) {
                $a = $b;
                $b = $c;
                $c = $a;
                $fa = $fb;
                $fb = $fc;
                $fc = $fa;
            }

            $tol1 = 2.0 * $eps * abs($b) + 0.5 * $tol;
            $half = 0.5 * ($c - $b);
            if (abs($half) <= $tol1 || $fb == 0.0) {
                return [$b, $evals, abs($c - $b)];
            }
            if ($evals >= $maxEvals) {
                return [$b, null, abs($c - $b)];
            }

            // The previous step is b - a, the step last taken (after the swap above).
            $last = $b - $a;
            $d = $half;
            if (abs($last) >= $tol1 && abs($fa) > abs($fb)) {
                $s = fdiv($fb, $fa);
                if ($a == $c) {
                    // secant step
                    $p = 2.0 * $half * $s;
                    $q = 1.0 - $s;
                } else {
                    // inverse quadratic interpolation
                    $qa = fdiv($fa, $fc);
                    $rb = fdiv($fb, $fc);
                    $p = $s * (2.0 * $half * $qa * ($qa - $rb) - ($b - $a) * ($rb - 1.0));
                    $q = ($qa - 1.0) * ($rb - 1.0) * ($s - 1.0);
                }
                if ($p > 0.0) {
                    $q = -$q;
                } else {
                    $p = -$p;
                }
                // Accept the interpolated step only if it stays well inside the bracket and shrinks fast enough.
                if (2.0 * $p < 3.0 * $half * $q - abs($tol1 * $q) && $p < abs(0.5 * $last * $q)) {
                    $d = fdiv($p, $q);
                }
            }

            $a = $b;
            $fa = $fb;
            if (abs($d) > $tol1) {
                $b += $d;
            } elseif ($half > 0.0) {
                $b += $tol1;
            } else {
                $b -= $tol1;
            }
            $fb = $f($b);
            $evals++;

            if ($fb * fdiv($fc, abs($fc)) > 0.0) {
                $c = $a;
                $fc = $fa;
            }
        }
    }
}
