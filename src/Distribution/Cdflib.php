<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/**
 * Follows the cdflib 0.4.3 crate routines (TOMS 708 incomplete beta, TOMS 654
 * incomplete gamma, Cody's normal CDF and the CDFLIB helpers) that the
 * phit distribution API uses. Branches, constants and operation order
 * follow the Rust source so results agree with the extension to the last bit
 * wherever the platform libm agrees.
 *
 * @internal
 */
final class Cdflib
{
    private const float NEG_EXPARG = -708.389_334_568_083_540_9;

    private const float POS_EXPARG = 709.775_615_066_259_888_4;

    // ------------------------------------------------------------ helpers

    /** CDFLIB eval_pol: c0 + c1 x + c2 x^2 + ... (ascending coefficients). */
    public static function evalPol(array $c, float $x): float
    {
        $n = count($c);
        $acc = $c[$n - 1];
        for ($i = $n - 2; $i >= 0; $i--) {
            $acc = $acc * $x + $c[$i];
        }

        return $acc;
    }

    /** exp(mu + x) without intermediate overflow. */
    public static function esum(int $mu, float $x): float
    {
        if ($x <= 0.0) {
            if ($mu >= 0) {
                $w = $mu + $x;
                if ($w <= 0.0) {
                    return exp($w);
                }
            }
        } elseif ($mu <= 0) {
            $w = $mu + $x;
            if ($w >= 0.0) {
                return exp($w);
            }
        }

        return exp((float) $mu) * exp($x);
    }

    /** ln(Γ(b) / Γ(a + b)) for b ≥ 8. */
    public static function algdiv(float $a, float $b): float
    {
        $C0 = 0.833333333333333e-1;
        $C1 = -0.277777777760991e-2;
        $C2 = 0.793650666825390e-3;
        $C3 = -0.595202931351870e-3;
        $C4 = 0.837308034031215e-3;
        $C5 = -0.165322962780713e-2;

        if ($b < $a) {
            $h = $b / $a;
            $c = 1.0 / (1.0 + $h);
            $x = $h / (1.0 + $h);
            $d = $a + ($b - 0.5);
        } else {
            $h = $a / $b;
            $c = $h / (1.0 + $h);
            $x = 1.0 / (1.0 + $h);
            $d = $b + ($a - 0.5);
        }

        $x2 = $x * $x;
        $s3 = 1.0 + ($x + $x2);
        $s5 = 1.0 + ($x + $x2 * $s3);
        $s7 = 1.0 + ($x + $x2 * $s5);
        $s9 = 1.0 + ($x + $x2 * $s7);
        $s11 = 1.0 + ($x + $x2 * $s9);

        $t = 1.0 / $b;
        $t *= $t;
        $w = (((($C5 * $s11 * $t + $C4 * $s9) * $t + $C3 * $s7) * $t + $C2 * $s5) * $t + $C1 * $s3) * $t + $C0;
        $w *= ($c / $b);

        $u = $d * self::alnrel($a / $b);
        $v = $a * (log($b) - 1.0);

        return $v < $u ? $w - $v - $u : $w - $u - $v;
    }

    /** ln B(a, b). */
    public static function betaLog(float $a0, float $b0): float
    {
        $E = 0.918938533204673;

        $a = $a0 < $b0 ? $a0 : $b0;
        $b = $a0 < $b0 ? $b0 : $a0;

        if ($a >= 8.0) {
            $w = self::bcorr($a, $b);
            $h = $a / $b;
            $c = $h / (1.0 + $h);
            $u = -(($a - 0.5) * log($c));
            $v = $b * self::alnrel($h);

            return $v < $u
                ? -(0.5 * log($b)) + $E + $w - $v - $u
                : -(0.5 * log($b)) + $E + $w - $u - $v;
        }

        if ($a < 1.0) {
            if ($b < 8.0) {
                return self::gammaLog($a) + (self::gammaLog($b) - self::gammaLog($a + $b));
            }

            return self::gammaLog($a) + self::algdiv($a, $b);
        }

        if ($a <= 2.0) {
            if ($b <= 2.0) {
                return self::gammaLog($a) + self::gammaLog($b) - self::gsumln($a, $b);
            }
            if ($b >= 8.0) {
                return self::gammaLog($a) + self::algdiv($a, $b);
            }

            return self::reduceB($a, $b, 0.0);
        }

        if ($b > 1000.0) {
            $n = (int) ($a - 1.0);
            $w = 1.0;
            for ($i = 0; $i < $n; $i++) {
                $a -= 1.0;
                $w *= $a / (1.0 + $a / $b);
            }

            return log($w) - $n * log($b) + (self::gammaLog($a) + self::algdiv($a, $b));
        }

        $n = (int) ($a - 1.0);
        $w = 1.0;
        for ($i = 0; $i < $n; $i++) {
            $a -= 1.0;
            $h = $a / $b;
            $w *= $h / (1.0 + $h);
        }
        $w = log($w);
        if ($b >= 8.0) {
            return $w + self::gammaLog($a) + self::algdiv($a, $b);
        }

        return self::reduceB($a, $b, $w);
    }

    private static function reduceB(float $a, float $b, float $w): float
    {
        $n = (int) ($b - 1.0);
        $z = 1.0;
        for ($i = 0; $i < $n; $i++) {
            $b -= 1.0;
            $z *= $b / ($a + $b);
        }

        return $w + log($z) + (self::gammaLog($a) + (self::gammaLog($b) - self::gsumln($a, $b)));
    }

    /** Δ(a) + Δ(b) − Δ(a + b) for a, b ≥ 8. */
    public static function bcorr(float $a0, float $b0): float
    {
        $C0 = 0.833333333333333e-1;
        $C1 = -0.277777777760991e-2;
        $C2 = 0.793650666825390e-3;
        $C3 = -0.595202931351870e-3;
        $C4 = 0.837308034031215e-3;
        $C5 = -0.165322962780713e-2;

        $a = $a0 < $b0 ? $a0 : $b0;
        $b = $a0 < $b0 ? $b0 : $a0;
        $h = $a / $b;
        $c = $h / (1.0 + $h);
        $x = 1.0 / (1.0 + $h);
        $x2 = $x * $x;
        $s3 = 1.0 + ($x + $x2);
        $s5 = 1.0 + ($x + $x2 * $s3);
        $s7 = 1.0 + ($x + $x2 * $s5);
        $s9 = 1.0 + ($x + $x2 * $s7);
        $s11 = 1.0 + ($x + $x2 * $s9);
        $t = 1.0 / $b;
        $t *= $t;
        $w = (((($C5 * $s11 * $t + $C4 * $s9) * $t + $C3 * $s7) * $t + $C2 * $s5) * $t + $C1 * $s3) * $t + $C0;
        $w *= ($c / $b);
        $t = 1.0 / $a;
        $t *= $t;

        return ((((($C5 * $t + $C4) * $t + $C3) * $t + $C2) * $t + $C1) * $t + $C0) / $a + $w;
    }

    // ------------------------------------------------------------ incomplete beta pieces

    private static function fpser(float $a, float $b, float $x, float $eps): float
    {
        $result = 1.0;
        if ($a > 1e-3 * $eps) {
            $result = 0.0;
            $t = $a * log($x);
            if ($t < self::NEG_EXPARG) {
                return $result;
            }
            $result = exp($t);
        }
        $result *= $b / $a;
        $tol = $eps / $a;
        $an = $a + 1.0;
        $t = $x;
        $s = $t / $an;
        while (true) {
            $an += 1.0;
            $t *= $x;
            $c = $t / $an;
            $s += $c;
            if (abs($c) <= $tol || is_nan($c)) {
                break;
            }
        }

        return $result * (1.0 + $a * $s);
    }

    private static function apser(float $a, float $b, float $x, float $eps): float
    {
        $G = 0.577215664901533;
        $bx = $b * $x;
        $t = $x - $bx;
        $c = $b * $eps <= 2e-2
            ? log($x) + self::psi($b) + $G + $t
            : log($bx) + $G + $t;
        $tol = 5.0 * $eps * abs($c);
        $j = 1.0;
        $s = 0.0;
        while (true) {
            $j += 1.0;
            $t *= $x - $bx / $j;
            $aj = $t / $j;
            $s += $aj;
            if (abs($aj) <= $tol || is_nan($aj)) {
                break;
            }
        }

        return -($a * ($c + $s));
    }

    private static function betaPser(float $a, float $b, float $x, float $eps): float
    {
        if ($x == 0.0) {
            return 0.0;
        }

        $a0 = $a < $b ? $a : $b;
        if ($a0 >= 1.0) {
            $z = $a * log($x) - self::betaLog($a, $b);
            $result = exp($z) / $a;
        } else {
            $b0 = $a < $b ? $b : $a;
            if ($b0 <= 1.0) {
                $result = $x ** $a;
                if ($result == 0.0) {
                    return 0.0;
                }
                $apb = $a + $b;
                if ($apb <= 1.0) {
                    $z = 1.0 + self::gam1($apb);
                } else {
                    $u = $a + $b - 1.0;
                    $z = (1.0 + self::gam1($u)) / $apb;
                }
                $c = (1.0 + self::gam1($a)) * (1.0 + self::gam1($b)) / $z;
                $result = $result * $c * ($b / $apb);
            } elseif ($b0 < 8.0) {
                $u = self::gammaLn1($a0);
                $m = (int) ($b0 - 1.0);
                $c = 1.0;
                for ($i = 1; $i <= $m; $i++) {
                    $b0 -= 1.0;
                    $c *= $b0 / ($a0 + $b0);
                }
                $u += log($c);
                $z = $a * log($x) - $u;
                $b0 -= 1.0;
                $apb = $a0 + $b0;
                if ($apb <= 1.0) {
                    $t = 1.0 + self::gam1($apb);
                } else {
                    $u = $a0 + $b0 - 1.0;
                    $t = (1.0 + self::gam1($u)) / $apb;
                }
                $result = exp($z) * ($a0 / $a) * (1.0 + self::gam1($b0)) / $t;
            } else {
                $u = self::gammaLn1($a0) + self::algdiv($a0, $b0);
                $z = $a * log($x) - $u;
                $result = ($a0 / $a) * exp($z);
            }
        }

        if ($result == 0.0 || $a <= 0.1 * $eps) {
            return $result;
        }

        $sum = 0.0;
        $n = 0.0;
        $c = 1.0;
        $tol = $eps / $a;
        while (true) {
            $n += 1.0;
            $c = $c * (0.5 + (0.5 - $b / $n)) * $x;
            $w = $c / ($a + $n);
            $sum += $w;
            if (abs($w) <= $tol || is_nan($w)) {
                break;
            }
        }

        return $result * (1.0 + $a * $sum);
    }

    /** x^a y^b / B(a, b). */
    private static function betaRcomp(float $a, float $b, float $x, float $y): float
    {
        $CONST = 0.398942280401433;
        if ($x == 0.0 || $y == 0.0) {
            return 0.0;
        }
        $a0 = $a < $b ? $a : $b;
        if ($a0 < 8.0) {
            if ($x <= 0.375) {
                $lnx = log($x);
                $lny = self::alnrel(-$x);
            } elseif ($y <= 0.375) {
                $lnx = self::alnrel(-$y);
                $lny = log($y);
            } else {
                $lnx = log($x);
                $lny = log($y);
            }
            $z = $a * $lnx + $b * $lny;
            if ($a0 >= 1.0) {
                return exp($z - self::betaLog($a, $b));
            }
            $b0 = $a < $b ? $b : $a;
            if ($b0 <= 1.0) {
                $result = exp($z);
                if ($result == 0.0) {
                    return 0.0;
                }
                $apb = $a + $b;
                if ($apb <= 1.0) {
                    $z = 1.0 + self::gam1($apb);
                } else {
                    $u = $a + $b - 1.0;
                    $z = (1.0 + self::gam1($u)) / $apb;
                }
                $c = (1.0 + self::gam1($a)) * (1.0 + self::gam1($b)) / $z;

                return $result * ($a0 * $c) / (1.0 + $a0 / $b0);
            }
            if ($b0 < 8.0) {
                $u = self::gammaLn1($a0);
                $n = (int) ($b0 - 1.0);
                $c = 1.0;
                for ($i = 1; $i <= $n; $i++) {
                    $b0 -= 1.0;
                    $c *= $b0 / ($a0 + $b0);
                }
                $u += log($c);
                $z -= $u;
                $b0 -= 1.0;
                $apb = $a0 + $b0;
                if ($apb <= 1.0) {
                    $t = 1.0 + self::gam1($apb);
                } else {
                    $u = $a0 + $b0 - 1.0;
                    $t = (1.0 + self::gam1($u)) / $apb;
                }

                return $a0 * exp($z) * (1.0 + self::gam1($b0)) / $t;
            }
            $u = self::gammaLn1($a0) + self::algdiv($a0, $b0);

            return $a0 * exp($z - $u);
        }
        if ($a <= $b) {
            $h = $a / $b;
            $x0 = $h / (1.0 + $h);
            $y0 = 1.0 / (1.0 + $h);
            $lambda = $a - ($a + $b) * $x;
        } else {
            $h = $b / $a;
            $x0 = 1.0 / (1.0 + $h);
            $y0 = $h / (1.0 + $h);
            $lambda = ($a + $b) * $y - $b;
        }
        $e = -($lambda / $a);
        $u = abs($e) <= 0.6 ? self::rlog1($e) : $e - log($x / $x0);
        $e = $lambda / $b;
        $v = abs($e) <= 0.6 ? self::rlog1($e) : $e - log($y / $y0);
        $z = exp(-($a * $u + $b * $v));

        return $CONST * sqrt($b * $x0) * $z * exp(-self::bcorr($a, $b));
    }

    /** exp(mu) x^a y^b / B(a, b). */
    private static function betaRcomp1(int $mu, float $a, float $b, float $x, float $y): float
    {
        $CONST = 0.398942280401433;
        $a0 = $a < $b ? $a : $b;
        if ($a0 >= 8.0) {
            if ($a <= $b) {
                $h = $a / $b;
                $x0 = $h / (1.0 + $h);
                $y0 = 1.0 / (1.0 + $h);
                $lambda = $a - ($a + $b) * $x;
            } else {
                $h = $b / $a;
                $x0 = 1.0 / (1.0 + $h);
                $y0 = $h / (1.0 + $h);
                $lambda = ($a + $b) * $y - $b;
            }
            $e = -($lambda / $a);
            $u = abs($e) <= 0.6 ? self::rlog1($e) : $e - log($x / $x0);
            $e = $lambda / $b;
            $v = abs($e) <= 0.6 ? self::rlog1($e) : $e - log($y / $y0);
            $t4 = -($a * $u + $b * $v);
            $z = self::esum($mu, $t4);

            return $CONST * sqrt($b * $x0) * $z * exp(-self::bcorr($a, $b));
        }

        if ($x <= 0.375) {
            $lnx = log($x);
            $lny = self::alnrel(-$x);
        } elseif ($y <= 0.375) {
            $lnx = self::alnrel(-$y);
            $lny = log($y);
        } else {
            $lnx = log($x);
            $lny = log($y);
        }
        $z = $a * $lnx + $b * $lny;
        if ($a0 >= 1.0) {
            return self::esum($mu, $z - self::betaLog($a, $b));
        }
        $b0 = $a < $b ? $b : $a;
        if ($b0 >= 8.0) {
            $u = self::gammaLn1($a0) + self::algdiv($a0, $b0);

            return $a0 * self::esum($mu, $z - $u);
        }
        if ($b0 > 1.0) {
            $u = self::gammaLn1($a0);
            $n = (int) ($b0 - 1.0);
            $c = 1.0;
            for ($i = 1; $i <= $n; $i++) {
                $b0 -= 1.0;
                $c *= $b0 / ($a0 + $b0);
            }
            $u += log($c);
            $z -= $u;
            $b0 -= 1.0;
            $apb = $a0 + $b0;
            if ($apb <= 1.0) {
                $t = 1.0 + self::gam1($apb);
            } else {
                $u = $a0 + $b0 - 1.0;
                $t = (1.0 + self::gam1($u)) / $apb;
            }

            return $a0 * self::esum($mu, $z) * (1.0 + self::gam1($b0)) / $t;
        }
        $result = self::esum($mu, $z);
        if ($result == 0.0) {
            return 0.0;
        }
        $apb = $a + $b;
        if ($apb <= 1.0) {
            $z = 1.0 + self::gam1($apb);
        } else {
            $u = $a + $b - 1.0;
            $z = (1.0 + self::gam1($u)) / $apb;
        }
        $c = (1.0 + self::gam1($a)) * (1.0 + self::gam1($b)) / $z;

        return $result * ($a0 * $c) / (1.0 + $a0 / $b0);
    }

    /** I_x(a, b) − I_x(a + n, b) for positive integer n. */
    private static function betaUp(float $a, float $b, float $x, float $y, int $n, float $eps): float
    {
        $apb = $a + $b;
        $ap1 = $a + 1.0;
        $mu = 0;
        $d = 1.0;
        if ($n != 1 && $a >= 1.0 && $apb >= 1.1 * $ap1) {
            $mu = (int) abs(self::NEG_EXPARG);
            $k = (int) self::POS_EXPARG;
            if ($k < $mu) {
                $mu = $k;
            }
            $d = exp(-(float) $mu);
        }
        $bup = self::betaRcomp1($mu, $a, $b, $x, $y) / $a;
        if ($n == 1 || $bup == 0.0) {
            return $bup;
        }
        $nm1 = $n - 1;
        $w = $d;
        $k = 0;
        if ($b > 1.0) {
            if ($y <= 1e-4) {
                $k = $nm1;
            } else {
                $r = ($b - 1.0) * $x / $y - $a;
                if ($r >= 1.0) {
                    $t = (float) $nm1;
                    $k = $nm1;
                    if ($r < $t) {
                        $k = (int) $r;
                    }
                }
            }
            for ($i = 1; $i <= $k; $i++) {
                $l = (float) ($i - 1);
                $d *= (($apb + $l) / ($ap1 + $l)) * $x;
                $w += $d;
            }
        }
        $kp1 = $k + 1;
        for ($i = $kp1; $i <= $nm1; $i++) {
            $l = (float) ($i - 1);
            $d *= ($apb + $l) / ($ap1 + $l) * $x;
            $w += $d;
            if ($d <= $eps * $w) {
                break;
            }
        }

        return $bup * $w;
    }

    /** (P(a, x), Q(a, x)) for a ≤ 1, used by beta_grat. */
    private static function gammaRat1(float $a, float $x, float $r, float $eps): array
    {
        if ($a * $x == 0.0) {
            return $x <= $a ? [0.0, 1.0] : [1.0, 0.0];
        }
        if ($a == 0.5) {
            $rtx = sqrt($x);
            if ($x < 0.25) {
                $p = self::errorF($rtx);

                return [$p, 0.5 + (0.5 - $p)];
            }
            $q = self::errorFc($rtx);

            return [0.5 + (0.5 - $q), $q];
        }

        if ($x < 1.1) {
            $an = 3.0;
            $c = $x;
            $sum = $x / ($a + 3.0);
            $tol = 0.1 * $eps / ($a + 1.0);
            while (true) {
                $an += 1.0;
                $c = -($c * ($x / $an));
                $t = $c / ($a + $an);
                $sum += $t;
                if (abs($t) <= $tol || is_nan($t)) {
                    break;
                }
            }
            $j = $a * $x * (($sum / 6.0 - 0.5 / ($a + 2.0)) * $x + 1.0 / ($a + 1.0));
            $z = $a * log($x);
            $h = self::gam1($a);
            $g = 1.0 + $h;
            $useLabel50 = $x < 0.25 ? $z > -0.13394 : $a < $x / 2.59;
            if ($useLabel50) {
                $l = self::rexp($z);
                $w = 0.5 + (0.5 + $l);
                $q = ($w * $j - $l) * $g - $h;
                if ($q < 0.0) {
                    return [1.0, 0.0];
                }

                return [0.5 + (0.5 - $q), $q];
            }
            $w = exp($z);
            $p = $w * $g * (0.5 + (0.5 - $j));

            return [$p, 0.5 + (0.5 - $p)];
        }

        $a2nm1 = 1.0;
        $a2n = 1.0;
        $b2nm1 = $x;
        $b2n = $x + (1.0 - $a);
        $c = 1.0;
        while (true) {
            $a2nm1 = $x * $a2n + $c * $a2nm1;
            $b2nm1 = $x * $b2n + $c * $b2nm1;
            $am0 = $a2nm1 / $b2nm1;
            $c += 1.0;
            $cma = $c - $a;
            $a2n = $a2nm1 + $cma * $a2n;
            $b2n = $b2nm1 + $cma * $b2n;
            $an0 = $a2n / $b2n;
            if (abs($an0 - $am0) < $eps * $an0 || is_nan($an0)) {
                $q = $r * $an0;

                return [0.5 + (0.5 - $q), $q];
            }
        }
    }

    /**
     * Asymptotic expansion of I_x(a, b) for a ≥ 15 and b ≤ 1, added to $w.
     * Returns null on the routine's soft failures (callers fall back to $w).
     */
    private static function betaGrat(float $a, float $b, float $x, float $y, float $wIn, float $eps): ?float
    {
        $bm1 = $b - 0.5 - 0.5;
        $nu = $a + 0.5 * $bm1;
        $lnx = $y > 0.375 ? log($x) : self::alnrel(-$y);
        $z = -($nu * $lnx);
        if ($b * $z == 0.0) {
            return null;
        }

        $r = $b * (1.0 + self::gam1($b)) * exp($b * log($z));
        $r = $r * exp($a * $lnx) * exp(0.5 * $bm1 * $lnx);
        $u = self::algdiv($b, $a) + $b * log($nu);
        $u = $r * exp(-$u);
        if ($u == 0.0) {
            return null;
        }

        [, $q] = self::gammaRat1($b, $z, $r, $eps);
        $v = 1.0 / $nu;
        $v = 0.25 * ($v * $v);
        $t2 = 0.25 * $lnx * $lnx;
        $l = $wIn / $u;
        $j = $q / $r;
        $sum = $j;
        $t = 1.0;
        $cn = 1.0;
        $n2 = 0.0;
        $c = array_fill(0, 30, 0.0);
        $d = array_fill(0, 30, 0.0);

        for ($n = 1; $n <= 30; $n++) {
            $bp2n = $b + $n2;
            $j = ($bp2n * ($bp2n + 1.0) * $j + ($z + $bp2n + 1.0) * $t) * $v;
            $n2 += 2.0;
            $t *= $t2;
            $cn /= $n2 * ($n2 + 1.0);
            $c[$n - 1] = $cn;
            $s = 0.0;
            if ($n != 1) {
                $coef = $b - $n;
                for ($i = 1; $i < $n; $i++) {
                    $s += $coef * $c[$i - 1] * $d[$n - $i - 1];
                    $coef += $b;
                }
            }
            $d[$n - 1] = $bm1 * $cn + $s / $n;
            $dj = $d[$n - 1] * $j;
            $sum += $dj;
            if ($sum <= 0.0) {
                return null;
            }
            if (abs($dj) <= $eps * ($sum + $l)) {
                return $wIn + $u * $sum;
            }
        }

        return $wIn + $u * $sum;
    }

    /** Asymptotic expansion of I_x(a, b) for a, b ≥ 15. */
    private static function betaAsym(float $a, float $b, float $lambda, float $eps): float
    {
        $E0 = 1.12837916709551;
        $E1 = 0.353553390593274;
        $NUM = 20;

        if ($a < $b) {
            $h = $a / $b;
            $r0 = 1.0 / (1.0 + $h);
            $r1 = ($b - $a) / $b;
            $w0 = 1.0 / sqrt($a * (1.0 + $h));
        } else {
            $h = $b / $a;
            $r0 = 1.0 / (1.0 + $h);
            $r1 = ($b - $a) / $a;
            $w0 = 1.0 / sqrt($b * (1.0 + $h));
        }

        $f = $a * self::rlog1(-($lambda / $a)) + $b * self::rlog1($lambda / $b);
        $t = exp(-$f);
        if ($t == 0.0) {
            return 0.0;
        }
        $z0 = sqrt($f);
        $z = 0.5 * ($z0 / $E1);
        $z2 = $f + $f;

        $a0 = array_fill(0, 21, 0.0);
        $b0 = array_fill(0, 21, 0.0);
        $c = array_fill(0, 21, 0.0);
        $d = array_fill(0, 21, 0.0);

        $a0[0] = 2.0 / 3.0 * $r1;
        $c[0] = -(0.5 * $a0[0]);
        $d[0] = -$c[0];
        $j0 = 0.5 / $E0 * self::errorFcScaled($z0);
        $j1 = $E1;
        $sum = $j0 + $d[0] * $w0 * $j1;
        $s = 1.0;
        $h2 = $h * $h;
        $hn = 1.0;
        $w = $w0;
        $znm1 = $z;
        $zn = $z2;

        $n = 2;
        while ($n <= $NUM) {
            $hn *= $h2;
            $a0[$n - 1] = 2.0 * $r0 * (1.0 + $h * $hn) / ($n + 2.0);
            $np1 = $n + 1;
            $s += $hn;
            $a0[$np1 - 1] = 2.0 * $r1 * $s / ($n + 3.0);

            for ($i = $n; $i <= $np1; $i++) {
                $r = -(0.5 * ($i + 1.0));
                $b0[0] = $r * $a0[0];
                for ($m = 2; $m <= $i; $m++) {
                    $bsum = 0.0;
                    for ($jj = 1; $jj < $m; $jj++) {
                        $mmj = $m - $jj;
                        $bsum += ((float) $jj * $r - (float) $mmj) * $a0[$jj - 1] * $b0[$mmj - 1];
                    }
                    $b0[$m - 1] = $r * $a0[$m - 1] + $bsum / $m;
                }
                $c[$i - 1] = $b0[$i - 1] / ($i + 1.0);
                $dsum = 0.0;
                for ($jj = 1; $jj < $i; $jj++) {
                    $imj = $i - $jj;
                    $dsum += $d[$imj - 1] * $c[$jj - 1];
                }
                $d[$i - 1] = -($dsum + $c[$i - 1]);
            }
            $j0 = $E1 * $znm1 + ($n - 1.0) * $j0;
            $j1 = $E1 * $zn + ((float) $n) * $j1;
            $znm1 *= $z2;
            $zn *= $z2;
            $w *= $w0;
            $t0 = $d[$n - 1] * $w * $j0;
            $w *= $w0;
            $t1 = $d[$np1 - 1] * $w * $j1;
            $sum += $t0 + $t1;
            if (abs($t0) + abs($t1) <= $eps * $sum) {
                break;
            }
            $n += 2;
        }

        $u = exp(-self::bcorr($a, $b));

        return $E0 * $t * $u * $sum;
    }

    /** Continued fraction for I_x(a, b) when a, b > 1. */
    private static function betaFrac(float $a, float $b, float $x, float $y, float $lambda, float $eps): float
    {
        $init = self::betaRcomp($a, $b, $x, $y);
        if ($init == 0.0) {
            return 0.0;
        }
        $c = 1.0 + $lambda;
        $c0 = $b / $a;
        $c1 = 1.0 + 1.0 / $a;
        $yp1 = $y + 1.0;
        $n = 0.0;
        $p = 1.0;
        $s = $a + 1.0;
        $an = 0.0;
        $anp1 = 1.0;
        $bn = 1.0;
        $bnp1 = $c / $c1;
        $r = $c1 / $c;

        while (true) {
            $n += 1.0;
            $t = $n / $a;
            $w = $n * ($b - $n) * $x;
            $e1 = $a / $s;
            $alpha = $p * ($p + $c0) * $e1 * $e1 * ($w * $x);
            $e2 = (1.0 + $t) / ($c1 + $t + $t);
            $beta = $n + $w / $s + $e2 * ($c + $n * $yp1);
            $p = 1.0 + $t;
            $s += 2.0;

            $tn = $alpha * $an + $beta * $anp1;
            $an = $anp1;
            $anp1 = $tn;
            $tn = $alpha * $bn + $beta * $bnp1;
            $bn = $bnp1;
            $bnp1 = $tn;
            // Fully qualified calls, and the new value stored before $r0 = $r: PHP
            // 8.5's function JIT adds an entry point after each by-name call and
            // can let the copy $r0 see a call result written into $r. See
            // docs/notes/function-jit.md.
            $next = \fdiv($anp1, $bnp1);
            $r0 = $r;
            $r = $next;

            if (\abs($r - $r0) <= $eps * $r || \is_nan($r)) {
                return $init * $r;
            }

            $an = \fdiv($an, $bnp1);
            $bn = \fdiv($bn, $bnp1);
            $anp1 = $r;
            $bnp1 = 1.0;
        }
    }

    /**
     * Regularized incomplete beta (I_x(a, b), 1 − I_x(a, b)), both tails
     * computed independently. Invalid input panics in the extension.
     *
     * @return array{0: float, 1: float}
     */
    public static function betaInc(float $a, float $b, float $x, float $y): array
    {
        $eps = PHP_FLOAT_EPSILON;
        $fail = static fn (string $why) => Support::panic("beta_inc($a, $b, $x, $y): $why");
        if ($a < 0.0 || $b < 0.0) {
            $fail('a or b is negative');
        }
        if ($a == 0.0 && $b == 0.0) {
            $fail('both a and b are zero');
        }
        if ($x < 0.0 || $x > 1.0) {
            $fail('x must be in [0..1]');
        }
        if ($y < 0.0 || $y > 1.0) {
            $fail('y must be in [0..1]');
        }
        $z = $x + $y - 0.5 - 0.5;
        if (abs($z) > 3.0 * $eps) {
            $fail('x + y must equal 1');
        }
        if ($x == 0.0) {
            if ($a == 0.0) {
                $fail('degenerate: x = 0 and a = 0');
            }

            return [0.0, 1.0];
        }
        if ($y == 0.0) {
            if ($b == 0.0) {
                $fail('degenerate: y = 0 and b = 0');
            }

            return [1.0, 0.0];
        }
        if ($a == 0.0) {
            return [1.0, 0.0];
        }
        if ($b == 0.0) {
            return [0.0, 1.0];
        }
        if (is_nan($a) || is_nan($b) || is_nan($x) || is_nan($y)) {
            return [NAN, NAN];
        }

        $eps = $eps > 1e-15 ? $eps : 1e-15;
        if (($a > $b ? $a : $b) < 1e-3 * $eps) {
            return [$b / ($a + $b), $a / ($a + $b)];
        }

        $ind = false;
        $a0 = $a;
        $b0 = $b;
        $x0 = $x;
        $y0 = $y;

        if (($a0 < $b0 ? $a0 : $b0) <= 1.0) {
            if ($x > 0.5) {
                $ind = true;
                $a0 = $b;
                $b0 = $a;
                $x0 = $y;
                $y0 = $x;
            }
            [$w, $w1] = self::smallBranch($a0, $b0, $x0, $y0, $eps);
        } else {
            $lambda = $a > $b ? ($a + $b) * $y - $b : $a - ($a + $b) * $x;
            if ($lambda < 0.0) {
                $lambda = abs($lambda);
                $ind = true;
                $a0 = $b;
                $b0 = $a;
                $x0 = $y;
                $y0 = $x;
            }
            [$w, $w1] = self::largeBranch($a0, $b0, $x0, $y0, $lambda, $eps);
        }

        return $ind ? [$w1, $w] : [$w, $w1];
    }

    private static function smallBranch(float $a0, float $b0, float $x0, float $y0, float $eps): array
    {
        $epsA = $eps * $a0;
        if ($b0 < ($eps < $epsA ? $eps : $epsA)) {
            $w = self::fpser($a0, $b0, $x0, $eps);

            return [$w, 0.5 + (0.5 - $w)];
        }
        $epsB = $eps * $b0;
        if ($a0 < ($eps < $epsB ? $eps : $epsB) && $b0 * $x0 <= 1.0) {
            $w1 = self::apser($a0, $b0, $x0, $eps);

            return [0.5 + (0.5 - $w1), $w1];
        }
        if (($a0 > $b0 ? $a0 : $b0) > 1.0) {
            if ($b0 <= 1.0) {
                $w = self::betaPser($a0, $b0, $x0, $eps);

                return [$w, 0.5 + (0.5 - $w)];
            }
            if ($x0 >= 0.3) {
                $w1 = self::betaPser($b0, $a0, $y0, $eps);

                return [0.5 + (0.5 - $w1), $w1];
            }
            if ($x0 < 0.1 && ($x0 * $b0) ** $a0 <= 0.7) {
                $w = self::betaPser($a0, $b0, $x0, $eps);

                return [$w, 0.5 + (0.5 - $w)];
            }
            if ($b0 > 15.0) {
                $w1 = self::betaGrat($b0, $a0, $y0, $x0, 0.0, 15.0 * $eps) ?? 0.0;

                return [0.5 + (0.5 - $w1), $w1];
            }
            $n = 20;
            $w1 = self::betaUp($b0, $a0, $y0, $x0, $n, $eps);
            $w1 = self::betaGrat($b0 + $n, $a0, $y0, $x0, $w1, 15.0 * $eps) ?? $w1;

            return [0.5 + (0.5 - $w1), $w1];
        }
        if ($a0 >= (0.2 < $b0 ? 0.2 : $b0)) {
            $w = self::betaPser($a0, $b0, $x0, $eps);

            return [$w, 0.5 + (0.5 - $w)];
        }
        if ($x0 ** $a0 <= 0.9) {
            $w = self::betaPser($a0, $b0, $x0, $eps);

            return [$w, 0.5 + (0.5 - $w)];
        }
        if ($x0 >= 0.3) {
            $w1 = self::betaPser($b0, $a0, $y0, $eps);

            return [0.5 + (0.5 - $w1), $w1];
        }
        $n = 20;
        $w1 = self::betaUp($b0, $a0, $y0, $x0, $n, $eps);
        $w1 = self::betaGrat($b0 + $n, $a0, $y0, $x0, $w1, 15.0 * $eps) ?? $w1;

        return [0.5 + (0.5 - $w1), $w1];
    }

    private static function largeBranch(float $a0, float $b0, float $x0, float $y0, float $lambda, float $eps): array
    {
        if ($b0 < 40.0 && $b0 * $x0 <= 0.7) {
            $w = self::betaPser($a0, $b0, $x0, $eps);

            return [$w, 0.5 + (0.5 - $w)];
        }
        if ($b0 < 40.0) {
            $n = (int) $b0;
            $b0r = $b0 - $n;
            $nUse = $n;
            if ($b0r == 0.0) {
                $nUse--;
                $b0r = 1.0;
            }
            $w = self::betaUp($b0r, $a0, $y0, $x0, $nUse, $eps);
            if ($x0 <= 0.7) {
                $w += self::betaPser($a0, $b0r, $x0, $eps);

                return [$w, 0.5 + (0.5 - $w)];
            }
            $a0r = $a0;
            if ($a0 <= 15.0) {
                $nn = 20;
                $w += self::betaUp($a0r, $b0r, $x0, $y0, $nn, $eps);
                $a0r += $nn;
            }
            $w = self::betaGrat($a0r, $b0r, $x0, $y0, $w, 15.0 * $eps) ?? $w;

            return [$w, 0.5 + (0.5 - $w)];
        }
        if ($a0 <= $b0) {
            $w = ($a0 <= 100.0 || $lambda > 0.03 * $a0)
                ? self::betaFrac($a0, $b0, $x0, $y0, $lambda, 15.0 * $eps)
                : self::betaAsym($a0, $b0, $lambda, 100.0 * $eps);

            return [$w, 0.5 + (0.5 - $w)];
        }
        $w = ($b0 <= 100.0 || $lambda > 0.03 * $b0)
            ? self::betaFrac($a0, $b0, $x0, $y0, $lambda, 15.0 * $eps)
            : self::betaAsym($a0, $b0, $lambda, 100.0 * $eps);

        return [$w, 0.5 + (0.5 - $w)];
    }

    // ------------------------------------------------------------ gamma family

    /** ln(1 + a). */
    public static function alnrel(float $a): float
    {
        if (abs($a) <= 0.375) {
            $t = $a / ($a + 2.0);
            $t2 = $t * $t;
            $num = ((-0.178874546012214e-1 * $t2 + 0.405303492862024) * $t2 + -0.129418923021993e1) * $t2 + 1.0;
            $den = ((-0.845104217945565e-1 * $t2 + 0.747811014037616) * $t2 + -0.162752256355323e1) * $t2 + 1.0;

            return 2.0 * $t * ($num / $den);
        }

        return log(1.0 + $a);
    }

    /** exp(x) − 1. */
    public static function rexp(float $x): float
    {
        if (abs($x) <= 0.15) {
            $num = (0.238082361044469e-1 * $x + 0.914041914819518e-9) * $x + 1.0;
            $den = (((0.595130811860248e-3 * $x + -0.119041179760821e-1) * $x + 0.107141568980644) * $x + -0.499999999085958) * $x + 1.0;

            return $x * ($num / $den);
        }
        $w = exp($x);

        return $x <= 0.0 ? ($w - 0.5) - 0.5 : $w * (0.5 + (0.5 - 1.0 / $w));
    }

    /** x − 1 − ln x. */
    public static function rlog(float $x): float
    {
        $A = 0.566749439387324e-1;
        $B = 0.456512608815524e-1;
        $P0 = 0.333333333333333;
        $P1 = -0.224696413112536;
        $P2 = 0.620886815375787e-2;
        $Q1 = -0.127408923933623e1;
        $Q2 = 0.354508718369557;

        if ($x < 0.61) {
            $r = $x - 0.5 - 0.5;

            return $r - log($x);
        }
        if ($x < 1.57) {
            if ($x < 0.82) {
                $u = ($x - 0.7) / 0.7;
                $w1 = $A - $u * 0.3;
            } elseif ($x < 1.18) {
                $u = $x - 0.5 - 0.5;
                $w1 = 0.0;
            } else {
                $u = 0.75 * $x - 1.0;
                $w1 = $B + $u / 3.0;
            }
            $r = $u / ($u + 2.0);
            $t = $r * $r;
            $w = (($P2 * $t + $P1) * $t + $P0) / (($Q2 * $t + $Q1) * $t + 1.0);

            return 2.0 * $t * (1.0 / (1.0 - $r) - $r * $w) + $w1;
        }
        $r = $x - 0.5 - 0.5;

        return $r - log($x);
    }

    /** x − ln(1 + x). */
    public static function rlog1(float $x): float
    {
        $A = 0.566749439387324e-1;
        $B = 0.456512608815524e-1;
        $P0 = 0.333333333333333;
        $P1 = -0.224696413112536;
        $P2 = 0.620886815375787e-2;
        $Q1 = -0.127408923933623e1;
        $Q2 = 0.354508718369557;

        if (is_nan($x) || $x < -0.39 || $x > 0.57) {
            $w = $x + 0.5 + 0.5;

            return $x - log($w);
        }
        if ($x < -0.18) {
            $h = ($x + 0.3) / 0.7;
            $w1 = $A - $h * 0.3;
        } elseif ($x > 0.18) {
            $h = 0.75 * $x - 0.25;
            $w1 = $B + $h / 3.0;
        } else {
            $h = $x;
            $w1 = 0.0;
        }
        $r = $h / ($h + 2.0);
        $t = $r * $r;
        $w = (($P2 * $t + $P1) * $t + $P0) / (($Q2 * $t + $Q1) * $t + 1.0);

        return 2.0 * $t * (1.0 / (1.0 - $r) - $r * $w) + $w1;
    }

    /** 1/Γ(1 + a) − 1 for −0.5 ≤ a ≤ 1.5. */
    public static function gam1(float $a): float
    {
        $S1 = 0.273076135303957;
        $S2 = 0.559398236957378e-1;
        $P = [0.577215664901533, -0.409078193005776, -0.230975380857675, 0.597275330452234e-1,
            0.766968181649490e-2, -0.514889771323592e-2, 0.589597428611429e-3];
        $Q = [1.0, 0.427569613095214, 0.158451672430138, 0.261132021441447e-1, 0.423244297896961e-2];
        $R = [-0.422784335098468, -0.771330383816272, -0.244757765222226, 0.118378989872749,
            0.930357293360349e-3, -0.118290993445146e-1, 0.223047661158249e-2, 0.266505979058923e-3,
            -0.132674909766242e-3];

        $t = $a;
        $d = $a - 0.5;
        if ($d > 0.0) {
            $t = $d - 0.5;
        }
        if ($t == 0.0) {
            return 0.0;
        }
        if ($t > 0.0) {
            $top = (((((($P[6] * $t + $P[5]) * $t + $P[4]) * $t + $P[3]) * $t + $P[2]) * $t + $P[1]) * $t) + $P[0];
            $bot = ((($Q[4] * $t + $Q[3]) * $t + $Q[2]) * $t + $Q[1]) * $t + 1.0;
            $w = $top / $bot;

            return $d <= 0.0 ? $a * $w : ($t / $a) * (($w - 0.5) - 0.5);
        }
        $top = ((((((($R[8] * $t + $R[7]) * $t + $R[6]) * $t + $R[5]) * $t + $R[4]) * $t + $R[3]) * $t
            + $R[2]) * $t + $R[1]) * $t + $R[0];
        $bot = ($S2 * $t + $S1) * $t + 1.0;
        $w = $top / $bot;

        return $d <= 0.0 ? $a * (($w + 0.5) + 0.5) : $t * $w / $a;
    }

    /** Γ(a); poles and overflow panic as in the crate. */
    public static function gamma(float $a): float
    {
        if (is_nan($a)) {
            return NAN;
        }
        $D = 0.41893853320467274178;
        $PI = 3.1415926535898;
        $P = [0.539637273585445e-3, 0.261939260042690e-2, 0.204493667594920e-1, 0.730981088720487e-1,
            0.279648642639792, 0.553413866010467, 1.0];
        $Q = [-0.832979206704073e-3, 0.470059485860584e-2, 0.225211131035340e-1, -0.170458969313360,
            -0.567902761974940e-1, 0.113062953091122e1, 1.0];
        $R1 = 0.820756370353826e-3;
        $R2 = -0.595156336428591e-3;
        $R3 = 0.793650663183693e-3;
        $R4 = -0.277777777770481e-2;
        $R5 = 0.833333333333333e-1;

        $x = $a;
        if (abs($a) < 15.0) {
            $t = 1.0;
            $m = ((int) $a) - 1;
            if ($m >= 0) {
                for ($i = 1; $i <= $m; $i++) {
                    $x -= 1.0;
                    $t *= $x;
                }
                $x -= 1.0;
            } else {
                $t = $a;
                if ($a <= 0.0) {
                    $mneg = -$m - 1;
                    for ($i = 1; $i <= $mneg; $i++) {
                        $x += 1.0;
                        $t *= $x;
                    }
                    $x = ($x + 0.5) + 0.5;
                    $t *= $x;
                    if ($t == 0.0) {
                        Support::panic("gamma($a): Γ has a pole at $a");
                    }
                }
                if (abs($t) < 1e-30 && abs($t) * PHP_FLOAT_MAX <= 1.0001) {
                    Support::panic("gamma($a): Γ($a) overflows f64");
                }
                if (abs($t) < 1e-30) {
                    return 1.0 / $t;
                }
            }
            $top = $P[0];
            $bot = $Q[0];
            for ($i = 1; $i < 7; $i++) {
                $top = $P[$i] + $x * $top;
                $bot = $Q[$i] + $x * $bot;
            }
            $g = $top / $bot;

            return $a < 1.0 ? $g / $t : $g * $t;
        }

        if ($a >= 1e3) {
            Support::panic("gamma($a): Γ($a) overflows f64");
        }
        if ($a <= -1e3) {
            Support::panic("gamma($a): Γ($a) underflows f64");
        }
        $s = 0.0;
        if ($a <= 0.0) {
            $x = -$a;
            $n = (int) $x;
            $t = $x - $n;
            if ($t > 0.9) {
                $t = 1.0 - $t;
            }
            $s = sin($PI * $t) / $PI;
            if ($n % 2 == 0) {
                $s = -$s;
            }
            if ($s == 0.0) {
                Support::panic("gamma($a): Γ has a pole at $a");
            }
        }
        $t = 1.0 / ($x * $x);
        $g = (((($R1 * $t + $R2) * $t + $R3) * $t + $R4) * $t + $R5) / $x;
        $lnx = log($x);
        $g = $D + $g + ($x - 0.5) * ($lnx - 1.0);
        $w = $g;
        $t = $g - $w;
        if ($w > 0.99999 * self::POS_EXPARG) {
            Support::panic("gamma($a): Γ($a) overflows f64");
        }
        $result = exp($w) * (1.0 + $t);
        if ($a < 0.0) {
            $result = (1.0 / ($result * $s)) / $x;
        }

        return $result;
    }

    /** ln Γ(1 + a) for −0.2 ≤ a ≤ 1.25. */
    public static function gammaLn1(float $a): float
    {
        if ($a < 0.6) {
            $num = (((((( -0.271935708322958e-2 * $a + -0.673562214325671e-1) * $a + -0.402055799310489) * $a
                + -0.780427615533591) * $a + -0.168860593646662) * $a + 0.844203922187225) * $a) + 0.577215664901533;
            $den = ((((((0.667465618796164e-3 * $a + 0.325038868253937e-1) * $a + 0.361951990101499) * $a
                + 0.156875193295039e1) * $a + 0.312755088914843e1) * $a + 0.288743195473681e1) * $a) + 1.0;

            return -($a * ($num / $den));
        }
        $x = $a - 0.5 - 0.5;
        $num = (((((0.497958207639485e-3 * $x + 0.170502484022650e-1) * $x + 0.156513060486551) * $x
            + 0.565221050691933) * $x + 0.848044614534529) * $x) + 0.422784335098467;
        $den = (((((0.116165475989616e-3 * $x + 0.713309612391000e-2) * $x + 0.101552187439830) * $x
            + 0.548042109832463) * $x + 0.124313399877507e1) * $x) + 1.0;

        return $x * ($num / $den);
    }

    /** Digamma ψ(x). */
    public static function psi(float $xx): float
    {
        if (is_nan($xx)) {
            return NAN;
        }
        $DX0 = 1.461632144968362341262659542325721325;
        $PIOV4 = 0.785398163397448;
        $P1 = [0.895385022981970e-2, 0.477762828042627e1, 0.142441585084029e3, 0.118645200713425e4,
            0.363351846806499e4, 0.413810161269013e4, 0.130560269827897e4];
        $P2 = [-0.212940445131011e1, -0.701677227766759e1, -0.448616543918019e1, -0.648157123766197];
        $Q1 = [0.448452573429826e2, 0.520752771467162e3, 0.221000799247830e4, 0.364127349079381e4,
            0.190831076596300e4, 0.691091682714533e-5];
        $Q2 = [0.322703493791143e2, 0.892920700481861e2, 0.546117738103215e2, 0.777788548522962e1];

        $xmax1 = 2_147_483_647.0;
        $xsmall = 1.0e-9;

        $x = $xx;
        $aug = 0.0;
        if ($x == 0.0) {
            Support::panic("psi($xx): ψ has a pole at $xx");
        }
        if ($x < 0.5) {
            if (abs($x) <= $xsmall) {
                $aug = -1.0 / $x;
            } else {
                $w = -$x;
                $sgn = $PIOV4;
                if ($w <= 0.0) {
                    $w = -$w;
                    $sgn = -$sgn;
                }
                if ($w >= $xmax1) {
                    Support::panic("psi($xx): argument is too large in magnitude");
                }
                $nq = (int) $w;
                $w -= $nq;
                $nq = (int) ($w * 4.0);
                $w = 4.0 * ($w - $nq * 0.25);
                $n = intdiv($nq, 2);
                if ($n + $n != $nq) {
                    $w = 1.0 - $w;
                }
                $z = $PIOV4 * $w;
                $m = intdiv($n, 2);
                if ($m + $m != $n) {
                    $sgn = -$sgn;
                }
                $n = intdiv($nq + 1, 2);
                $m2 = intdiv($n, 2);
                if ($m2 + $m2 == $n) {
                    if ($z == 0.0) {
                        Support::panic("psi($xx): ψ has a pole at $xx");
                    }
                    $aug = 4.0 * $sgn * (cos($z) / sin($z));
                } else {
                    $aug = 4.0 * $sgn * (sin($z) / cos($z));
                }
            }
            $x = 1.0 - $x;
        }

        if ($x <= 3.0) {
            $den = $x;
            $upper = $P1[0] * $x;
            for ($i = 1; $i <= 5; $i++) {
                $den = ($den + $Q1[$i - 1]) * $x;
                $upper = ($upper + $P1[$i]) * $x;
            }
            $den = ($upper + $P1[6]) / ($den + $Q1[5]);

            return $den * ($x - $DX0) + $aug;
        }
        if ($x < $xmax1) {
            $w = (1.0 / $x) / $x;
            $den = $w;
            $upper = $P2[0] * $w;
            for ($i = 1; $i <= 3; $i++) {
                $den = ($den + $Q2[$i - 1]) * $w;
                $upper = ($upper + $P2[$i]) * $w;
            }
            $aug = $upper / ($den + $Q2[3]) - 0.5 / $x + $aug;

            return $aug + log($x);
        }

        return $aug + log($x);
    }

    /** ln Γ(a) for a > 0. */
    public static function gammaLog(float $a): float
    {
        if ($a <= 0.8) {
            return self::gammaLn1($a) - log($a);
        }
        if ($a <= 2.25) {
            $t = $a - 0.5 - 0.5;

            return self::gammaLn1($t);
        }
        if ($a < 10.0) {
            $n = (int) ($a - 1.25);
            $t = $a;
            $w = 1.0;
            for ($i = 1; $i <= $n; $i++) {
                $t -= 1.0;
                $w *= $t;
            }

            return self::gammaLn1($t - 1.0) + log($w);
        }
        $t = 1.0 / $a;
        $t *= $t;
        $w = ((((( -0.165322962780713e-2 * $t + 0.837308034031215e-3) * $t + -0.595202931351870e-3) * $t
            + 0.793650666825390e-3) * $t + -0.277777777760991e-2) * $t + 0.833333333333333e-1) / $a;

        return 0.418938533204673 + $w + ($a - 0.5) * (log($a) - 1.0);
    }

    /** ln Γ(a + b) for 1 ≤ a, b ≤ 2. */
    public static function gsumln(float $a, float $b): float
    {
        $x = $a + $b - 2.0;
        if ($x <= 0.25) {
            return self::gammaLn1(1.0 + $x);
        }
        if ($x <= 1.25) {
            return self::gammaLn1($x) + self::alnrel($x);
        }

        return self::gammaLn1($x - 1.0) + log($x * (1.0 + $x));
    }

    // ------------------------------------------------------------ incomplete gamma (TOMS 654, max accuracy)

    private const array TEMME_D0 = [
        0.833333333333333e-1, -0.148148148148148e-1, 0.115740740740741e-2, 0.352733686067019e-3,
        -0.178755144032922e-3, 0.391926317852244e-4, -0.218544851067999e-5, -0.185406221071516e-5,
        0.829671134095309e-6, -0.176659527368261e-6, 0.670785354340150e-8, 0.102618097842403e-7,
        -0.438203601845335e-8,
    ];

    private const array TEMME_D1 = [
        -0.347222222222222e-2, 0.264550264550265e-2, -0.990226337448560e-3, 0.205761316872428e-3,
        -0.401877572016461e-6, -0.180985503344900e-4, 0.764916091608111e-5, -0.161209008945634e-5,
        0.464712780280743e-8, 0.137863344691572e-6, -0.575254560351770e-7, 0.119516285997781e-7,
    ];

    private const array TEMME_D2 = [
        -0.268132716049383e-2, 0.771604938271605e-3, 0.200938786008230e-5, -0.107366532263652e-3,
        0.529234488291201e-4, -0.127606351886187e-4, 0.342357873409614e-7, 0.137219573090629e-5,
        -0.629899213838006e-6, 0.142806142060642e-6,
    ];

    private const array TEMME_D3 = [
        0.229472093621399e-3, -0.469189494395256e-3, 0.267720632062839e-3, -0.756180167188398e-4,
        -0.239650511386730e-6, 0.110826541153473e-4, -0.567495282699160e-5, 0.142309007324359e-5,
    ];

    private const array TEMME_D4 = [
        0.784039221720067e-3, -0.299072480303190e-3, -0.146384525788434e-5, 0.664149821546512e-4,
        -0.396836504717943e-4, 0.113757269706784e-4,
    ];

    private const array TEMME_D5 = [
        -0.697281375836586e-4, 0.277275324495939e-3, -0.199325705161888e-3, 0.679778047793721e-4,
    ];

    private const array TEMME_D6 = [-0.592166437353694e-3, 0.270878209671804e-3];

    private const float D10 = -0.185185185185185e-2;

    private const float D20 = 0.413359788359788e-2;

    private const float D30 = 0.649434156378601e-3;

    private const float D40 = -0.861888290916712e-3;

    private const float D50 = -0.336798553366358e-3;

    private const float D60 = 0.531307936463992e-3;

    private const float D70 = 0.344367606892378e-3;

    private const float RT2PIN = 0.398942280401433;

    /**
     * Regularized incomplete gamma (P(a, x), Q(a, x)), computed independently.
     *
     * @return array{0: float, 1: float}
     */
    public static function gammaInc(float $a, float $x): array
    {
        if (is_nan($a) || is_nan($x)) {
            return [NAN, NAN];
        }
        if ($a < 0.0) {
            Support::panic("gamma_inc($a, $x): parameter a must be non-negative");
        }
        if ($x < 0.0) {
            Support::panic("gamma_inc($a, $x): argument x must be non-negative");
        }
        if ($a == 0.0 && $x == 0.0) {
            Support::panic("gamma_inc($a, $x): both a and x are zero");
        }
        $result = self::gammaIncCore($a, $x);
        if ($result[0] == 2.0) {
            Support::panic("gamma_inc_with_acc($a, $x, Max): indeterminate at a = $a, x = $x (deep asymptotic regime)");
        }

        return $result;
    }

    private static function gammaIncCore(float $a, float $x): array
    {
        $ALOG10 = 2.30258509299405;
        $acc = 5.0e-15 > PHP_FLOAT_EPSILON ? 5.0e-15 : PHP_FLOAT_EPSILON;
        $e0 = 0.25e-3;
        $x0 = 31.0;
        $big = 20.0;
        $e = PHP_FLOAT_EPSILON;

        if ($a * $x == 0.0) {
            return $x <= $a ? [0.0, 1.0] : [1.0, 0.0];
        }

        if ($a < 1.0) {
            if ($a == 0.5) {
                $rtx = sqrt($x);
                if ($x < 0.25) {
                    $ans = self::errorF($rtx);

                    return [$ans, 0.5 + (0.5 - $ans)];
                }
                $qans = self::errorFc($rtx);

                return [0.5 + (0.5 - $qans), $qans];
            }
            if ($x < 1.1) {
                return self::taylorPOverXa($a, $x, $acc);
            }
            $t1 = $a * log($x) - $x;
            $u = $a * exp($t1);
            if ($u == 0.0) {
                return [1.0, 0.0];
            }
            $r = $u * (1.0 + self::gam1($a));

            return self::continuedFraction($a, $x, $r, $acc, $e);
        }

        if ($a < $big) {
            if ($a <= $x && $x < $x0) {
                $twoa = $a + $a;
                $m = (int) $twoa;
                if ($twoa == (float) $m) {
                    $ii = intdiv($m, 2);

                    return self::finiteQHalfInteger($x, $ii, $a == (float) $ii);
                }
            }
            $t1 = $a * log($x) - $x;
            $r = exp($t1) / self::gamma($a);
        } else {
            $l = $x / $a;
            if ($l == 0.0) {
                return [0.0, 1.0];
            }
            $s = 0.5 + (0.5 - $l);
            $z = self::rlog($l);
            if ($z >= 700.0 / $a) {
                if (abs($s) <= 2.0 * $e) {
                    return [2.0, 0.0];
                }

                return $x <= $a ? [0.0, 1.0] : [1.0, 0.0];
            }
            $y = $a * $z;
            $rta = sqrt($a);
            if (abs($s) <= $e0 / $rta) {
                return self::temmeForLEq1($a, $l, $z, $y, $e);
            }
            if (abs($s) <= 0.4) {
                return self::temmeGeneral($a, $l, $z, $y, $rta, $e);
            }
            $t = 1.0 / $a;
            $t *= $t;
            $t1 = (((0.75 * $t - 1.0) * $t + 3.5) * $t - 105.0) / ($a * 1260.0);
            $t1 -= $y;
            $r = self::RT2PIN * $rta * exp($t1);
        }

        if ($r == 0.0) {
            return $x <= $a ? [0.0, 1.0] : [1.0, 0.0];
        }
        if ($x <= ($a > $ALOG10 ? $a : $ALOG10)) {
            return self::taylorPOverR($a, $x, $r, $acc);
        }
        if ($x < $x0) {
            return self::continuedFraction($a, $x, $r, $acc, $e);
        }

        return self::asymptoticExpansionQ($a, $x, $r, $acc);
    }

    private static function taylorPOverXa(float $a, float $x, float $acc): array
    {
        $an = 3.0;
        $c = $x;
        $sum = $x / ($a + 3.0);
        $tol = 3.0 * $acc / ($a + 1.0);
        while (true) {
            $an += 1.0;
            $c = -($c * ($x / $an));
            $t = $c / ($a + $an);
            $sum += $t;
            if (abs($t) <= $tol || is_nan($t)) {
                break;
            }
        }
        $j = $a * $x * (($sum / 6.0 - 0.5 / ($a + 2.0)) * $x + 1.0 / ($a + 1.0));
        $z = $a * log($x);
        $h = self::gam1($a);
        $g = 1.0 + $h;
        $use200 = $x < 0.25 ? $z > -0.13394 : $a < $x / 2.59;
        if ($use200) {
            $l = self::rexp($z);
            $w = 0.5 + (0.5 + $l);
            $qans = ($w * $j - $l) * $g - $h;
            if ($qans < 0.0) {
                return [1.0, 0.0];
            }

            return [0.5 + (0.5 - $qans), $qans];
        }
        $w = exp($z);
        $ans = $w * $g * (0.5 + (0.5 - $j));

        return [$ans, 0.5 + (0.5 - $ans)];
    }

    private static function taylorPOverR(float $a, float $x, float $r, float $acc): array
    {
        $wk = array_fill(0, 20, 0.0);
        $apn = $a + 1.0;
        $t = $x / $apn;
        $wk[0] = $t;
        $filled = 20;
        for ($n = 2; $n <= 20; $n++) {
            $apn += 1.0;
            $t *= $x / $apn;
            if ($t <= 1e-3) {
                $filled = $n;
                break;
            }
            $wk[$n - 1] = $t;
        }
        $sum = $t;
        $tol = 0.5 * $acc;
        while (true) {
            $apn += 1.0;
            $t *= $x / $apn;
            $sum += $t;
            if ($t <= $tol || is_nan($t)) {
                break;
            }
        }
        for ($m = 1; $m <= $filled - 1; $m++) {
            $sum += $wk[$filled - $m - 1];
        }
        $ans = $r / $a * (1.0 + $sum);

        return [$ans, 0.5 + (0.5 - $ans)];
    }

    private static function continuedFraction(float $a, float $x, float $r, float $acc, float $e): array
    {
        $tol = 5.0 * $e > $acc ? 5.0 * $e : $acc;
        $a2n = 1.0;
        $a2nm1 = 1.0;
        $b2nm1 = $x;
        $b2n = $x + (1.0 - $a);
        $c = 1.0;
        while (true) {
            $a2nm1 = $x * $a2n + $c * $a2nm1;
            $b2nm1 = $x * $b2n + $c * $b2nm1;
            $am0 = $a2nm1 / $b2nm1;
            $c += 1.0;
            $cma = $c - $a;
            $a2n = $a2nm1 + $cma * $a2n;
            $b2n = $b2nm1 + $cma * $b2n;
            $an0 = $a2n / $b2n;
            if (abs($an0 - $am0) < $tol * $an0 || is_nan($an0)) {
                $qans = $r * $an0;

                return [0.5 + (0.5 - $qans), $qans];
            }
        }
    }

    private static function asymptoticExpansionQ(float $a, float $x, float $r, float $acc): array
    {
        $wk = array_fill(0, 20, 0.0);
        $amn = $a - 1.0;
        $t = $amn / $x;
        $wk[0] = $t;
        $filled = 20;
        for ($n = 2; $n <= 20; $n++) {
            $amn -= 1.0;
            $t *= $amn / $x;
            if (abs($t) <= 1e-3) {
                $filled = $n;
                break;
            }
            $wk[$n - 1] = $t;
        }
        $sum = $t;
        while (abs($t) > $acc) {
            $amn -= 1.0;
            $t *= $amn / $x;
            $sum += $t;
        }
        for ($m = 1; $m <= $filled - 1; $m++) {
            $sum += $wk[$filled - $m - 1];
        }
        $qans = $r / $x * (1.0 + $sum);

        return [0.5 + (0.5 - $qans), $qans];
    }

    private static function finiteQHalfInteger(float $x, int $i, bool $aEqI): array
    {
        if ($aEqI) {
            $sum = exp(-$x);
            $t = $sum;
            $n = 1;
            $c = 0.0;
        } else {
            $rtx = sqrt($x);
            $sum = self::errorFc($rtx);
            $t = exp(-$x) / (1.77245385090552 * $rtx);
            $n = 0;
            $c = -0.5;
        }
        while ($n != $i) {
            $n++;
            $c += 1.0;
            $t = $x * $t / $c;
            $sum += $t;
        }

        return [0.5 + (0.5 - $sum), $sum];
    }

    /** The c0..c6, d70 polynomial in z and u shared by both Temme paths (short form). */
    private static function temmeShort(float $z, float $u): float
    {
        $D0 = self::TEMME_D0;
        $D1 = self::TEMME_D1;
        $D2 = self::TEMME_D2;
        $D3 = self::TEMME_D3;
        $D4 = self::TEMME_D4;
        $D5 = self::TEMME_D5;
        $D6 = self::TEMME_D6;
        $c0 = (((((($D0[6] * $z + $D0[5]) * $z + $D0[4]) * $z + $D0[3]) * $z + $D0[2]) * $z + $D0[1]) * $z + $D0[0]) * $z - 1.0 / 3.0;
        $c1 = ((((($D1[5] * $z + $D1[4]) * $z + $D1[3]) * $z + $D1[2]) * $z + $D1[1]) * $z + $D1[0]) * $z + self::D10;
        $c2 = (((($D2[4] * $z + $D2[3]) * $z + $D2[2]) * $z + $D2[1]) * $z + $D2[0]) * $z + self::D20;
        $c3 = ((($D3[3] * $z + $D3[2]) * $z + $D3[1]) * $z + $D3[0]) * $z + self::D30;
        $c4 = ($D4[1] * $z + $D4[0]) * $z + self::D40;
        $c5 = ($D5[1] * $z + $D5[0]) * $z + self::D50;
        $c6 = $D6[0] * $z + self::D60;

        return ((((((self::D70 * $u + $c6) * $u + $c5) * $u + $c4) * $u + $c3) * $u + $c2) * $u + $c1) * $u + $c0;
    }

    private static function temmeGeneral(float $a, float $l, float $zIn, float $y, float $rta, float $e): array
    {
        $s = 0.5 + (0.5 - $l);
        if (abs($s) <= 2.0 * $e && 3.28e-3 < $a * $e * $e) {
            return [2.0, 0.0];
        }
        $c = exp(-$y);
        $w = 0.5 * self::errorFcScaled(sqrt($y));
        $u = 1.0 / $a;
        $z = sqrt($zIn + $zIn);
        if ($l < 1.0) {
            $z = -$z;
        }
        if (abs($s) <= 1.0e-3) {
            $t = self::temmeShort($z, $u);
        } else {
            $D0 = self::TEMME_D0;
            $D1 = self::TEMME_D1;
            $D2 = self::TEMME_D2;
            $D3 = self::TEMME_D3;
            $D4 = self::TEMME_D4;
            $D5 = self::TEMME_D5;
            $D6 = self::TEMME_D6;
            $c0 = ((((((((((((($D0[12] * $z + $D0[11]) * $z + $D0[10]) * $z + $D0[9]) * $z + $D0[8]) * $z + $D0[7]) * $z
                + $D0[6]) * $z + $D0[5]) * $z + $D0[4]) * $z + $D0[3]) * $z + $D0[2]) * $z + $D0[1]) * $z + $D0[0]) * $z) - 1.0 / 3.0;
            $c1 = ((((((((((($D1[11] * $z + $D1[10]) * $z + $D1[9]) * $z + $D1[8]) * $z + $D1[7]) * $z + $D1[6]) * $z
                + $D1[5]) * $z + $D1[4]) * $z + $D1[3]) * $z + $D1[2]) * $z + $D1[1]) * $z + $D1[0]) * $z + self::D10;
            $c2 = ((((((((($D2[9] * $z + $D2[8]) * $z + $D2[7]) * $z + $D2[6]) * $z + $D2[5]) * $z + $D2[4]) * $z
                + $D2[3]) * $z + $D2[2]) * $z + $D2[1]) * $z + $D2[0]) * $z + self::D20;
            $c3 = ((((((($D3[7] * $z + $D3[6]) * $z + $D3[5]) * $z + $D3[4]) * $z + $D3[3]) * $z + $D3[2]) * $z
                + $D3[1]) * $z + $D3[0]) * $z + self::D30;
            $c4 = ((((($D4[5] * $z + $D4[4]) * $z + $D4[3]) * $z + $D4[2]) * $z + $D4[1]) * $z + $D4[0]) * $z + self::D40;
            $c5 = ((($D5[3] * $z + $D5[2]) * $z + $D5[1]) * $z + $D5[0]) * $z + self::D50;
            $c6 = ($D6[1] * $z + $D6[0]) * $z + self::D60;
            $t = ((((((self::D70 * $u + $c6) * $u + $c5) * $u + $c4) * $u + $c3) * $u + $c2) * $u + $c1) * $u + $c0;
        }

        if ($l < 1.0) {
            $ans = $c * ($w - self::RT2PIN * $t / $rta);

            return [$ans, 0.5 + (0.5 - $ans)];
        }
        $qans = $c * ($w + self::RT2PIN * $t / $rta);

        return [0.5 + (0.5 - $qans), $qans];
    }

    private static function temmeForLEq1(float $a, float $l, float $zIn, float $y, float $e): array
    {
        if (3.28e-3 < $a * $e * $e) {
            return [2.0, 0.0];
        }
        $RTPI = 1.77245385090552;
        $c = 0.5 + (0.5 - $y);
        $w = (0.5 - sqrt($y) * (0.5 + (0.5 - $y / 3.0)) / $RTPI) / $c;
        $u = 1.0 / $a;
        $z = sqrt($zIn + $zIn);
        if ($l < 1.0) {
            $z = -$z;
        }
        $t = self::temmeShort($z, $u);
        $rta = sqrt($a);
        if ($l < 1.0) {
            $ans = $c * ($w - self::RT2PIN * $t / $rta);

            return [$ans, 0.5 + (0.5 - $ans)];
        }
        $qans = $c * ($w + self::RT2PIN * $t / $rta);

        return [0.5 + (0.5 - $qans), $qans];
    }

    // ------------------------------------------------------------ error function

    private const array ERF_A = [7.71058495001320e-05, -1.33733772997339e-03, 3.23076579225834e-02,
        4.79137145607681e-02, 1.28379167095513e-01];

    private const array ERF_B = [3.01048631703895e-03, 5.38971687740286e-02, 3.75795757275549e-01];

    private const array ERF_P = [-1.36864857382717e-07, 5.64195517478974e-01, 7.21175825088309e+00,
        4.31622272220567e+01, 1.52989285046940e+02, 3.39320816734344e+02, 4.51918953711873e+02,
        3.00459261020162e+02];

    private const array ERF_Q = [1.00000000000000e+00, 1.27827273196294e+01, 7.70001529352295e+01,
        2.77585444743988e+02, 6.38980264465631e+02, 9.31354094850610e+02, 7.90950925327898e+02,
        3.00459260956983e+02];

    private const array ERF_R = [2.10144126479064e+00, 2.62370141675169e+01, 2.13688200555087e+01,
        4.65807828718470e+00, 2.82094791773523e-01];

    private const array ERF_S = [9.41537750555460e+01, 1.87114811799590e+02, 9.90191814623914e+01,
        1.80124575948747e+01];

    private const float ERF_C = 0.564189583547756;

    public static function errorF(float $x): float
    {
        $A = self::ERF_A;
        $B = self::ERF_B;
        $P = self::ERF_P;
        $Q = self::ERF_Q;
        $R = self::ERF_R;
        $S = self::ERF_S;
        $ax = abs($x);
        if ($ax <= 0.5) {
            $t = $x * $x;
            $top = (((($A[0] * $t + $A[1]) * $t + $A[2]) * $t + $A[3]) * $t) + $A[4] + 1.0;
            $bot = (($B[0] * $t + $B[1]) * $t + $B[2]) * $t + 1.0;
            $erf = $ax * ($top / $bot);
        } elseif ($ax <= 4.0) {
            $top = ((((((($P[0] * $ax + $P[1]) * $ax + $P[2]) * $ax + $P[3]) * $ax + $P[4]) * $ax + $P[5]) * $ax + $P[6]) * $ax) + $P[7];
            $bot = ((((((($Q[0] * $ax + $Q[1]) * $ax + $Q[2]) * $ax + $Q[3]) * $ax + $Q[4]) * $ax + $Q[5]) * $ax + $Q[6]) * $ax) + $Q[7];
            $erf = 0.5 + (0.5 - exp(-($x * $x)) * $top / $bot);
        } elseif ($ax < 5.8) {
            $x2 = $x * $x;
            $t = 1.0 / $x2;
            $top = ((($R[0] * $t + $R[1]) * $t + $R[2]) * $t + $R[3]) * $t + $R[4];
            $bot = ((($S[0] * $t + $S[1]) * $t + $S[2]) * $t + $S[3]) * $t + 1.0;
            $erf = (self::ERF_C - $top / ($x2 * $bot)) / $ax;
            $erf = 0.5 + (0.5 - exp(-$x2) * $erf);
        } else {
            $erf = 1.0;
        }
        if ($x < 0.0) {
            $erf = -$erf;
        }

        return $erf;
    }

    public static function errorFc(float $x): float
    {
        return self::errorFcInner($x, false);
    }

    public static function errorFcScaled(float $x): float
    {
        return self::errorFcInner($x, true);
    }

    private static function errorFcInner(float $x, bool $scaled): float
    {
        $A = self::ERF_A;
        $B = self::ERF_B;
        $P = self::ERF_P;
        $Q = self::ERF_Q;
        $R = self::ERF_R;
        $S = self::ERF_S;
        $ax = abs($x);
        if ($ax <= 0.5) {
            $t = $x * $x;
            $top = (((($A[0] * $t + $A[1]) * $t + $A[2]) * $t + $A[3]) * $t) + $A[4] + 1.0;
            $bot = (($B[0] * $t + $B[1]) * $t + $B[2]) * $t + 1.0;
            $erfc = 0.5 + (0.5 - $x * ($top / $bot));
            if ($scaled) {
                $erfc *= exp($t);
            }

            return $erfc;
        }
        if ($ax <= 4.0) {
            $top = ((((((($P[0] * $ax + $P[1]) * $ax + $P[2]) * $ax + $P[3]) * $ax + $P[4]) * $ax + $P[5]) * $ax + $P[6]) * $ax) + $P[7];
            $bot = ((((((($Q[0] * $ax + $Q[1]) * $ax + $Q[2]) * $ax + $Q[3]) * $ax + $Q[4]) * $ax + $Q[5]) * $ax + $Q[6]) * $ax) + $Q[7];
            $erfc = $top / $bot;
        } else {
            if ($x <= -5.6) {
                return $scaled ? 2.0 * exp($x * $x) : 2.0;
            }
            if (! $scaled && ($x > 100.0 || $x * $x > -self::NEG_EXPARG)) {
                return 0.0;
            }
            $t = 1.0 / $x;
            $t *= $t;
            $top = ((($R[0] * $t + $R[1]) * $t + $R[2]) * $t + $R[3]) * $t + $R[4];
            $bot = ((($S[0] * $t + $S[1]) * $t + $S[2]) * $t + $S[3]) * $t + 1.0;
            $erfc = (self::ERF_C - $t * $top / $bot) / $ax;
        }
        if ($scaled) {
            if ($x < 0.0) {
                $erfc = 2.0 * exp($x * $x) - $erfc;
            }
        } else {
            $erfc *= exp(-($x * $x));
            if ($x < 0.0) {
                $erfc = 2.0 - $erfc;
            }
        }

        return $erfc;
    }

    // ------------------------------------------------------------ normal

    /**
     * (Φ(x), 1 − Φ(x)) by Cody's algorithm, both tails computed directly.
     *
     * @return array{0: float, 1: float}
     */
    public static function cumnor(float $x): array
    {
        static $A = [2.2352520354606839287e00, 1.6102823106855587881e02, 1.0676894854603709582e03,
            1.8154981253343561249e04, 6.5682337918207449113e-2];
        static $B = [4.7202581904688241870e01, 9.7609855173777669322e02, 1.0260932208618978205e04,
            4.5507789335026729956e04];
        static $C = [3.9894151208813466764e-1, 8.8831497943883759412e00, 9.3506656132177855979e01,
            5.9727027639480026226e02, 2.4945375852903726711e03, 6.8481904505362823326e03,
            1.1602651437647350124e04, 9.8427148383839780218e03, 1.0765576773720192317e-8];
        static $D = [2.2266688044328115691e01, 2.3538790178262499861e02, 1.5193775994075548050e03,
            6.4855582982667607550e03, 1.8615571640885098091e04, 3.4900952721145977266e04,
            3.8912003286093271411e04, 1.9685429676859990727e04];
        static $P = [2.1589853405795699e-1, 1.274011611602473639e-1, 2.2235277870649807e-2,
            1.421619193227893466e-3, 2.9112874951168792e-5, 2.307344176494017303e-2];
        static $Q = [1.28426009614491121e00, 4.68238212480865118e-1, 6.59881378689285515e-2,
            3.78239633202758244e-3, 7.29751555083966205e-5];
        $SQRPI = 0.39894228040143267794;
        $THRSH = 0.66291;
        $ROOT32 = 5.656854248;

        $eps = 0.5 * PHP_FLOAT_EPSILON;
        $min = PHP_FLOAT_MIN;

        $y = abs($x);
        if ($y <= $THRSH) {
            $xsq = $y > $eps ? $x * $x : 0.0;
            $xnum = $A[4] * $xsq;
            $xden = $xsq;
            for ($i = 0; $i < 3; $i++) {
                $xnum = ($xnum + $A[$i]) * $xsq;
                $xden = ($xden + $B[$i]) * $xsq;
            }
            $r = $x * ($xnum + $A[3]) / ($xden + $B[3]);
            $result = 0.5 + $r;
            $ccum = 0.5 - $r;
        } elseif ($y <= $ROOT32) {
            $xnum = $C[8] * $y;
            $xden = $y;
            for ($i = 0; $i < 7; $i++) {
                $xnum = ($xnum + $C[$i]) * $y;
                $xden = ($xden + $D[$i]) * $y;
            }
            $r = ($xnum + $C[7]) / ($xden + $D[7]);
            $xsq = floor($y * 16.0) / 16.0;
            $del = ($y - $xsq) * ($y + $xsq);
            $result = exp(-0.5 * $xsq * $xsq) * exp(-0.5 * $del) * $r;
            $ccum = 1.0 - $result;
            if ($x > 0.0) {
                [$result, $ccum] = [$ccum, $result];
            }
        } else {
            $xsq = 1.0 / ($x * $x);
            $xnum = $P[5] * $xsq;
            $xden = $xsq;
            for ($i = 0; $i < 4; $i++) {
                $xnum = ($xnum + $P[$i]) * $xsq;
                $xden = ($xden + $Q[$i]) * $xsq;
            }
            $r = $xsq * ($xnum + $P[4]) / ($xden + $Q[4]);
            $r = ($SQRPI - $r) / $y;
            $v = $x * 16.0;
            $xsq = ($v < 0.0 ? ceil($v) : floor($v)) / 16.0;
            $del = ($x - $xsq) * ($x + $xsq);
            $result = exp(-0.5 * $xsq * $xsq) * exp(-0.5 * $del) * $r;
            $ccum = 1.0 - $result;
            if ($x > 0.0) {
                [$result, $ccum] = [$ccum, $result];
            }
        }

        if ($result < $min) {
            $result = 0.0;
        }
        if ($ccum < $min) {
            $ccum = 0.0;
        }

        return [$result, $ccum];
    }

    /** x with Φ(x) = p, from CDFLIB's Newton iteration on cumnor. */
    public static function dinvnr(float $p, float $q): float
    {
        $R2PI = 0.3989422804014326;
        $pp = $p < $q ? $p : $q;
        $strtx = self::stvaln($pp);
        $xcur = $strtx;
        for ($i = 1; $i <= 100; $i++) {
            $cum = self::cumnor($xcur)[0];
            $dx = fdiv($cum - $pp, $R2PI * exp(-0.5 * $xcur * $xcur));
            $xcur -= $dx;
            if (abs(fdiv($dx, $xcur)) < 1.0e-13) {
                return $p <= $q ? $xcur : -$xcur;
            }
        }

        return $p <= $q ? $strtx : -$strtx;
    }

    private static function stvaln(float $p): float
    {
        static $XDEN = [0.993484626060e-1, 0.588581570495e0, 0.531103462366e0, 0.103537752850e0, 0.38560700634e-2];
        static $XNUM = [-0.322232431088e0, -1.000000000000e0, -0.342242088547e0, -0.204231210245e-1, -0.453642210148e-4];
        if ($p <= 0.5) {
            $sign = -1.0;
            $z = $p;
        } else {
            $sign = 1.0;
            $z = 1.0 - $p;
        }
        $y = sqrt(-2.0 * log($z));

        return $sign * ($y + self::evalPol($XNUM, $y) / self::evalPol($XDEN, $y));
    }

    /** Asymptotic Student t quantile used as a Newton start (CDFLIB dt1). */
    public static function dt1(float $p, float $q, float $df): float
    {
        static $COEF = [[1.0, 1.0], [3.0, 16.0, 5.0], [-15.0, 17.0, 19.0, 3.0], [-945.0, -1920.0, 1482.0, 776.0, 79.0]];
        static $DENOM = [4.0, 96.0, 384.0, 92160.0];
        $x = abs(self::dinvnr($p, $q));
        $xx = $x * $x;
        $sum1 = $x;
        $denpow = 1.0;
        for ($i = 0; $i < 4; $i++) {
            $term = self::evalPol($COEF[$i], $xx) * $x;
            $denpow *= $df;
            $sum1 += fdiv($term, $denpow * $DENOM[$i]);
        }

        return $p >= 0.5 ? $sum1 : -$sum1;
    }
}
