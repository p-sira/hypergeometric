//! Real Zhang-Jin hypergeometric kernels translated from SciPy's specfun.
//! See LICENSES/scipy-specfun.txt, xsf-BSD-3-Clause.txt, cephes-iv.txt.
#![allow(
    unused_assignments,
    unused_mut,
    unused_parens,
    non_snake_case,
    clippy::all
)]
use crate::gamma::{digamma as psi, gamma as Gamma, gamma_ratio, lgam_sgn, rgamma};
const HYP2F1_EPS: f64 = 1e-13;
const HYP2F1_ETHRESH: f64 = 1e-12;
const HYP2F1_MAXITER: i32 = 10000;
const MACHEP: f64 = f64::EPSILON * 0.5;
#[inline]
fn neg_one_pow(k: i32) -> f64 {
    if k % 2 == 0 { 1.0 } else { -1.0 }
}

pub(super) fn hyp2f1(mut a: f64, mut b: f64, mut c: f64, mut x: f64) -> f64 {
    let mut d: f64 = 0.0;
    let mut d1: f64 = 0.0;
    let mut d2: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut p: f64 = 0.0;
    let mut q: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut ax: f64 = 0.0;
    let mut ia: f64 = 0.0;
    let mut ib: f64 = 0.0;
    let mut ic: f64 = 0.0;
    let mut id: f64 = 0.0;
    let mut err: f64 = 0.0;
    let mut t1: f64 = 0.0;
    let mut i: i32 = 0;
    let mut aid: i32 = 0;
    let mut neg_int_a: i32 = 0;
    let mut neg_int_b: i32 = 0;
    let mut neg_int_ca_or_cb: i32 = 0;
    err = 0.0f64;
    ax = (x).abs();
    s = (1.0f64 - x);
    ia = (a).round();
    ib = (b).round();
    if (x == 0.0f64) {
        {
            return 1.0f64;
        }
    }
    // Exact identity, evaluated with log1p so tiny |x| and very large
    // negative x do not lose significant digits to 1-x rounding.
    if a == 1.0 && b == 1.0 && c == 2.0 && x < 1.0 {
        return -(-x).ln_1p() / x;
    }
    d = ((c - a) - b);
    id = (d).round();
    if (((a == ((0) as f64)) || (b == ((0) as f64))) && (c != ((0) as f64))) {
        {
            return 1.0f64;
        }
    }
    if a == -1.0 && c != 0.0 {
        return (-b / c).mul_add(x, 1.0);
    }
    if ((a <= ((0) as f64)) && (a == ia)) {
        {
            neg_int_a = 1;
        }
    }
    if ((b <= ((0) as f64)) && (b == ib)) {
        {
            neg_int_b = 1;
        }
    }
    if (((d <= ((-1) as f64)) && !(((d - id).abs() > HYP2F1_EPS) && (s < ((0) as f64))))
        && !((neg_int_a != 0) || (neg_int_b != 0)))
    {
        {
            return ((s).powf(d) * hyp2f1((c - a), (c - b), c, x));
        }
    }
    if (((d <= ((0) as f64)) && (x == ((1) as f64))) && !((neg_int_a != 0) || (neg_int_b != 0))) {
        return f64::INFINITY;
    }
    if ((ax < 1.0f64) || (x == (-1.0f64))) {
        {
            if (b == c) {
                {
                    if (neg_int_b != 0) {
                        {
                            y = hyp2f1_neg_c_equal_bc(a, b, x);
                        }
                    } else {
                        {
                            y = (s).powf((-a));
                        }
                    }
                    return y;
                }
            }
            if (a == c) {
                return if neg_int_a != 0 {
                    hyp2f1_neg_c_equal_bc(b, a, x)
                } else {
                    s.powf(-b)
                };
            }
        }
    }
    if (c <= 0.0f64) {
        {
            ic = (c).round();
            if (c == ic) {
                {
                    if ((neg_int_a != 0) && (ia > ic)) {
                        return hyt2f1(a, b, c, x, std::slice::from_mut(&mut err));
                    }
                    if ((neg_int_b != 0) && (ib > ic)) {
                        return hyt2f1(a, b, c, x, std::slice::from_mut(&mut err));
                    }
                    return f64::INFINITY;
                }
            }
        }
    }
    if ((neg_int_a != 0) || (neg_int_b != 0)) {
        return hyt2f1(a, b, c, x, std::slice::from_mut(&mut err));
    }
    // Pfaff's transformation avoids subtracting the large coefficients
    // of the connection at 1/x. Mixed-sign series use the DD kernel below.
    // Bound x so rounding x/(x-1) cannot erase its distance from one.
    if (-10000.0..-1.0).contains(&x) {
        if a.abs() < b.abs() {
            return s.powf(-a) * hyp2f1(a, c - b, c, x / (x - 1.0));
        }
        return s.powf(-b) * hyp2f1(c - a, b, c, x / (x - 1.0));
    }
    t1 = (b - a).abs();
    if ((x < (-2.0f64)) && ((t1 - (t1).round()).abs() > HYP2F1_EPS)) {
        {
            p = hyp2f1(
                a,
                ((((1) as f64) - c) + a),
                ((((1) as f64) - b) + a),
                (1.0f64 / x),
            );
            q = hyp2f1(
                b,
                ((((1) as f64) - c) + b),
                ((((1) as f64) - a) + b),
                (1.0f64 / x),
            );
            p *= (-x).powf((-a));
            q *= (-x).powf((-b));
            t1 = Gamma(c);
            s = ((t1 * Gamma((b - a))) * (rgamma(b) * rgamma((c - a))));
            y = ((t1 * Gamma((a - b))) * (rgamma(a) * rgamma((c - b))));
            return ((s * p) + (y * q));
        }
    } else {
        if (x < (-1.0f64)) {
            {
                if ((a).abs() < (b).abs()) {
                    {
                        return ((s).powf((-a)) * hyp2f1(a, (c - b), c, (x / (x - ((1) as f64)))));
                    }
                } else {
                    {
                        return ((s).powf((-b)) * hyp2f1(b, (c - a), c, (x / (x - ((1) as f64)))));
                    }
                }
            }
        }
    }
    if (ax > 1.0f64) {
        return f64::INFINITY;
    }
    p = (c - a);
    ia = (p).round();
    if ((ia <= 0.0f64) && ((p - ia).abs() < HYP2F1_EPS)) {
        neg_int_ca_or_cb = 1;
    }
    r = (c - b);
    ib = (r).round();
    if ((ib <= 0.0f64) && ((r - ib).abs() < HYP2F1_EPS)) {
        neg_int_ca_or_cb = 1;
    }
    id = (d).round();
    q = (d - id).abs();
    if (ax == 1.0) {
        {
            if (x > 0.0f64) {
                {
                    if (neg_int_ca_or_cb != 0) {
                        {
                            if (d >= 0.0f64) {
                                return ((s).powf(d)
                                    * hys2f1(
                                        (c - a),
                                        (c - b),
                                        c,
                                        x,
                                        std::slice::from_mut(&mut err),
                                    ));
                            } else {
                                return f64::INFINITY;
                            }
                        }
                    }
                    if (d <= 0.0f64) {
                        return f64::INFINITY;
                    }
                    y = ((Gamma(c) * Gamma(d)) * (rgamma(p) * rgamma(r)));
                    return y;
                }
            }
            if (d <= (-1.0f64)) {
                return f64::INFINITY;
            }
        }
    }
    if (d < 0.0f64) {
        {
            y = hyt2f1(a, b, c, x, std::slice::from_mut(&mut err));
            if (err < HYP2F1_ETHRESH) {
                return y;
            }
            err = 0.0f64;
            aid = ((((2) as f64) - id) as i32);
            e = (c + ((aid) as f64));
            d2 = hyp2f1(a, b, e, x);
            d1 = hyp2f1(a, b, (e + 1.0f64), x);
            q = ((a + b) + 1.0f64);
            i = 0;
            while (i < aid) {
                {
                    r = (e - 1.0f64);
                    y = ((((e * (r - (((2.0f64 * e) - q) * x))) * d2)
                        + ((((e - a) * (e - b)) * x) * d1))
                        / ((e * r) * s));
                    e = r;
                    d1 = d2;
                    d2 = y;
                }
                {
                    i += 1;
                    i
                };
            }
            return y;
        }
    }
    if (neg_int_ca_or_cb != 0) {
        {
            return ((s).powf(d) * hys2f1((c - a), (c - b), c, x, std::slice::from_mut(&mut err)));
        }
    }
    y = hyt2f1(a, b, c, x, std::slice::from_mut(&mut err));
    if (err > HYP2F1_ETHRESH) {
        {}
    }
    return y;
}

pub(super) fn hyp2f1_neg_c_equal_bc(mut a: f64, mut b: f64, mut x: f64) -> f64 {
    let mut k: f64 = 0.0;
    let mut collector: f64 = ((1) as f64);
    let mut sum: f64 = ((1) as f64);
    let mut collector_max: f64 = ((1) as f64);
    if !((b).abs() < 1e5f64) {
        {
            return f64::NAN;
        }
    }
    k = ((1) as f64);
    while (k <= (-b)) {
        {
            collector *= ((((a + k) - ((1) as f64)) * x) / k);
            collector_max = ((collector).abs()).max(collector_max);
            sum += collector;
        }
        {
            k += 1.0;
            k
        };
    }
    if ((1e-16f64 * (((1) as f64) + (collector_max / (sum).abs()))) > 1e-7f64) {
        {
            return f64::NAN;
        }
    }
    return sum;
}

pub(super) fn hyt2f1(mut a: f64, mut b: f64, mut c: f64, mut x: f64, loss: &mut [f64]) -> f64 {
    let mut p: f64 = 0.0;
    let mut q: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut t: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut w: f64 = 0.0;
    let mut d: f64 = 0.0;
    let mut err: f64 = 0.0;
    let mut err1: f64 = 0.0;
    let mut ax: f64 = 0.0;
    let mut id: f64 = 0.0;
    let mut d1: f64 = 0.0;
    let mut d2: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut y1: f64 = 0.0;
    let mut i: i32 = 0;
    let mut aid: i32 = 0;
    let mut sign: i32 = 0;
    let mut ia: i32 = 0;
    let mut ib: i32 = 0;
    let mut neg_int_a: i32 = 0;
    let mut neg_int_b: i32 = 0;
    ia = (((a).round()) as i32);
    ib = (((b).round()) as i32);
    if ((a <= ((0) as f64)) && (a == ia as f64)) {
        {
            neg_int_a = 1;
        }
    }
    if ((b <= ((0) as f64)) && (b == ib as f64)) {
        {
            neg_int_b = 1;
        }
    }
    err = 0.0f64;
    s = (1.0f64 - x);
    if ((x < (-0.5f64)) && !((neg_int_a != 0) || (neg_int_b != 0))) {
        {
            if (b > a) {
                y = ((s).powf((-a))
                    * hys2f1(a, (c - b), c, ((-x) / s), std::slice::from_mut(&mut err)));
            } else {
                y = ((s).powf((-b))
                    * hys2f1((c - a), b, c, ((-x) / s), std::slice::from_mut(&mut err)));
            }
            {
                loss[0] = err;
                return y;
            }
        }
    }
    d = ((c - a) - b);
    id = (d).round();
    if ((x > 0.9f64) && !((neg_int_a != 0) || (neg_int_b != 0))) {
        {
            if ((d - id).abs() > MACHEP) {
                {
                    let mut sgngam: i32 = 0;
                    // Close to x=1, direct summation can require hundreds
                    // or thousands of terms. The continuation at 1-x is
                    // short and avoids that pass when its two coefficients
                    // are well separated from gamma poles. Keep the former
                    // error-driven choice at nearby integer d, or when
                    // parameters cross possible gamma poles.
                    let direct_first = x <= 0.95
                        || a <= 0.0
                        || b <= 0.0
                        || c < a
                        || c < b
                        || c > 32.0
                        || (d - id).abs() < 0.125;
                    if direct_first {
                        y = hys2f1(a, b, c, x, std::slice::from_mut(&mut err));
                        if err < HYP2F1_ETHRESH {
                            loss[0] = err;
                            return y;
                        }
                    }
                    q = hys2f1(a, b, (1.0f64 - d), s, std::slice::from_mut(&mut err));
                    sign = 1;
                    w = lgam_sgn(d, std::slice::from_mut(&mut sgngam));
                    sign *= sgngam;
                    w -= lgam_sgn((c - a), std::slice::from_mut(&mut sgngam));
                    sign *= sgngam;
                    w -= lgam_sgn((c - b), std::slice::from_mut(&mut sgngam));
                    sign *= sgngam;
                    q *= (((sign) as f64) * (w).exp());
                    r = ((s).powf(d)
                        * hys2f1(
                            (c - a),
                            (c - b),
                            (d + 1.0f64),
                            s,
                            std::slice::from_mut(&mut err1),
                        ));
                    sign = 1;
                    w = lgam_sgn((-d), std::slice::from_mut(&mut sgngam));
                    sign *= sgngam;
                    w -= lgam_sgn(a, std::slice::from_mut(&mut sgngam));
                    sign *= sgngam;
                    w -= lgam_sgn(b, std::slice::from_mut(&mut sgngam));
                    sign *= sgngam;
                    r *= (((sign) as f64) * (w).exp());
                    y = (q + r);
                    q = (q).abs();
                    r = (r).abs();
                    if (q > r) {
                        {
                            r = q;
                        }
                    }
                    err += (err1 + ((MACHEP * r) / y));
                    y *= Gamma(c);
                    {
                        loss[0] = err;
                        return y;
                    }
                }
            } else {
                {
                    if (id >= 0.0f64) {
                        {
                            e = d;
                            d1 = d;
                            d2 = 0.0f64;
                            aid = ((id) as i32);
                        }
                    } else {
                        {
                            e = (-d);
                            d1 = 0.0f64;
                            d2 = d;
                            aid = ((-id) as i32);
                        }
                    }
                    ax = (s).ln();
                    let mut psi_one = psi(1.0);
                    let mut psi_e = psi(1.0 + e);
                    let mut psi_a = psi(a + d1);
                    let mut psi_b = psi(b + d1);
                    // For positive arguments, psi(z+1)=psi(z)+1/z removes
                    // four costly digamma evaluations from every term.
                    // Keep direct evaluation across nonpositive arguments
                    // where recurrence through a pole can lose accuracy.
                    let recur_psi = a + d1 > 0.0 && b + d1 > 0.0;
                    y = (((psi_one + psi_e) - psi_a) - psi_b) - ax;
                    y *= rgamma((e + 1.0f64));
                    p = ((((a + d1) * (b + d1)) * s) * rgamma((e + 2.0f64)));
                    t = 1.0f64;
                    loop {
                        {
                            if recur_psi {
                                psi_one += 1.0 / t;
                                psi_e += 1.0 / (t + e);
                                psi_a += 1.0 / ((a + (t - 1.0)) + d1);
                                psi_b += 1.0 / ((b + (t - 1.0)) + d1);
                                r = (((psi_one + psi_e) - psi_a) - psi_b) - ax;
                            } else {
                                r = ((((psi((1.0f64 + t)) + psi(((1.0f64 + t) + e)))
                                    - psi(((a + t) + d1)))
                                    - psi(((b + t) + d1)))
                                    - ax);
                            }
                            q = (p * r);
                            y += q;
                            p *= ((s * ((a + t) + d1)) / (t + 1.0f64));
                            p *= (((b + t) + d1) / ((t + 1.0f64) + e));
                            t += 1.0f64;
                            if (t > ((HYP2F1_MAXITER) as f64)) {
                                {
                                    loss[0] = 1.0f64;
                                    return f64::NAN;
                                }
                            }
                        }
                        if !((y == ((0) as f64)) || ((q / y).abs() > HYP2F1_EPS)) {
                            break;
                        }
                    }
                    if (id == 0.0f64) {
                        {
                            y *= (Gamma(c) / (Gamma(a) * Gamma(b)));
                            {
                                loss[0] = err;
                                return y;
                            }
                        }
                    }
                    y1 = 1.0f64;
                    t = 0.0f64;
                    p = 1.0f64;
                    i = 1;
                    while (i < aid) {
                        {
                            r = ((1.0f64 - e) + t);
                            p *= (((s * ((a + t) + d2)) * ((b + t) + d2)) / r);
                            t += 1.0f64;
                            p /= t;
                            y1 += p;
                        }
                        {
                            i += 1;
                            i
                        };
                    }
                    p = Gamma(c);
                    y1 *= ((Gamma(e) * p) * (rgamma((a + d1)) * rgamma((b + d1))));
                    y *= (p * (rgamma((a + d2)) * rgamma((b + d2))));
                    if ((aid & 1) != 0) {
                        y = (-y);
                    }
                    q = (s).powf(id);
                    if (id > 0.0f64) {
                        y *= q;
                    } else {
                        y1 *= q;
                    }
                    y += y1;
                    {
                        loss[0] = err;
                        return y;
                    }
                }
            }
        }
    }
    y = hys2f1(a, b, c, x, std::slice::from_mut(&mut err));
    loss[0] = err;
    return y;
}

// Binary64 defining series with a first-order rounding-error certificate.
// Each term carries at most about 8(k+1) roundings relative to the exact
// term and Neumaier summation keeps the accumulation error O(eps|sum|), so
// 8 eps sum (k+1)|t_k| bounds the error. Accept only when that certificate
// meets the same 1e-13 target as the DD kernel; otherwise fall back to DD.
fn hys2f1_f64(a: f64, b: f64, c: f64, x: f64) -> Option<f64> {
    let (mut term, mut sum, mut correction) = (1.0_f64, 1.0_f64, 0.0_f64);
    let mut weighted = 1.0_f64;
    let certify = |sum: f64, weighted: f64| {
        (8.0 * f64::EPSILON * weighted <= 1e-13 * sum.abs()).then_some(sum)
    };
    for k in 0..10000 {
        let kf = k as f64;
        let numerator = (a + kf) * (b + kf);
        if numerator == 0.0 {
            return certify(sum + correction, weighted);
        }
        term *= numerator * x / ((c + kf) * (kf + 1.0));
        let next = sum + term;
        correction += if sum.abs() >= term.abs() {
            (sum - next) + term
        } else {
            (term - next) + sum
        };
        sum = next;
        weighted += (kf + 2.0) * term.abs();
        if !weighted.is_finite() || !sum.is_finite() {
            return None;
        }
        // Since 0 <= bound, test the cheap necessary size condition before
        // forming the tail-ratio bound and its divisions.
        if term.abs() > 0.5 * f64::EPSILON * sum.abs() {
            continue;
        }
        let (an, bn, cn) = (a + kf + 1.0, b + kf + 1.0, c + kf + 1.0);
        let bound = x.abs() * (an / cn).max(1.0) * (bn / (kf + 2.0)).max(1.0);
        if an > 0.0
            && bn > 0.0
            && cn > 0.0
            && bound < 1.0
            && term.abs() <= 0.5 * f64::EPSILON * (1.0 - bound) * sum.abs()
        {
            return certify(sum + correction, weighted);
        }
    }
    None
}

// Mixed-sign defining series in DD. In particular, near x=1 a negative
// numerator can cancel many digits before its positive tail converges.
fn hys2f1_dd(a: f64, b: f64, c: f64, x: f64) -> Option<f64> {
    use crate::dd::Dd;
    let (mut term, mut sum) = (Dd::from(1.0), Dd::from(1.0));
    let mut absolute_sum = 1.0;
    for k in 0..10000 {
        let kf = k as f64;
        let numerator = Dd::sum(a, kf) * Dd::sum(b, kf);
        if numerator.hi == 0.0 && numerator.lo == 0.0 {
            return (absolute_sum * 1e-28 <= 1e-13 * sum.hi.abs()).then_some(sum.hi + sum.lo);
        }
        term = (term * numerator)
            .mul_f64(x)
            .div_dd(Dd::sum(c, kf).mul_small(kf + 1.0));
        sum = sum + term;
        absolute_sum += term.hi.abs();
        if !sum.hi.is_finite() || !absolute_sum.is_finite() {
            return None;
        }
        if term.hi.abs() > 1e-25 * sum.hi.abs() {
            continue;
        }
        let (an, bn, cn) = (a + kf + 1.0, b + kf + 1.0, c + kf + 1.0);
        // Each ratio is bounded by its current value or its limit one.
        let bound = x.abs() * (an / cn).max(1.0) * (bn / (kf + 2.0)).max(1.0);
        if an > 0.0
            && bn > 0.0
            && cn > 0.0
            && bound < 1.0
            && term.hi.abs() <= 1e-25 * (1.0 - bound) * sum.hi.abs()
        {
            return (absolute_sum * 1e-28 <= 1e-13 * sum.hi.abs()).then_some(sum.hi + sum.lo);
        }
    }
    None
}

pub(super) fn hys2f1(mut a: f64, mut b: f64, mut c: f64, mut x: f64, loss: &mut [f64]) -> f64 {
    if (a < 0.0 || b < 0.0 || c < 0.0 || x < 0.0)
        && a.abs().max(b.abs()).max(c.abs()) <= 128.0
        && x.abs() < 1.0
    {
        if let Some(value) = hys2f1_f64(a, b, c, x).or_else(|| hys2f1_dd(a, b, c, x)) {
            loss[0] = 1e-13;
            return value;
        }
    }
    let mut f: f64 = 0.0;
    let mut g: f64 = 0.0;
    let mut h: f64 = 0.0;
    let mut k: f64 = 0.0;
    let mut m: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut u: f64 = 0.0;
    let mut umax: f64 = 0.0;
    let mut i: i32 = 0;
    let mut ib: i32 = 0;
    let mut intflag: i32 = 0;
    if ((b).abs() > (a).abs()) {
        {
            f = b;
            b = a;
            a = f;
        }
    }
    ib = (((b).round()) as i32);
    if (((b == ib as f64) && (ib <= 0)) && ((b).abs() < (a).abs())) {
        {
            f = b;
            b = a;
            a = f;
            intflag = 1;
        }
    }
    if (((((a).abs() > ((c).abs() + ((1) as f64))) || (intflag != 0))
        && ((c - a).abs() > ((2) as f64)))
        && ((a).abs() > ((2) as f64)))
    {
        {
            return hyp2f1ra(a, b, c, x, loss);
        }
    }
    i = 0;
    umax = 0.0f64;
    f = a;
    g = b;
    h = c;
    s = 1.0f64;
    u = 1.0f64;
    k = 0.0f64;
    loop {
        {
            if (h == 0.0) {
                {
                    loss[0] = 1.0f64;
                    return f64::INFINITY;
                }
            }
            m = (k + 1.0f64);
            u = (u * ((((f + k) * (g + k)) * x) / ((h + k) * m)));
            s += u;
            k = (u).abs();
            if (k > umax) {
                umax = k;
            }
            k = m;
            if ({
                i += 1;
                i
            } > HYP2F1_MAXITER)
            {
                {
                    loss[0] = 1.0f64;
                    return s;
                }
            }
        }
        if !((s == ((0) as f64)) || ((u / s).abs() > MACHEP)) {
            break;
        }
    }
    loss[0] = (((MACHEP * umax) / (s).abs()) + (MACHEP * ((i) as f64)));
    return s;
}

pub(super) fn hyp2f1ra(mut a: f64, mut b: f64, mut c: f64, mut x: f64, loss: &mut [f64]) -> f64 {
    let mut f2: f64 = 0.0;
    let mut f1: f64 = 0.0;
    let mut f0: f64 = 0.0;
    let mut n: i32 = 0;
    let mut t: f64 = 0.0;
    let mut err: f64 = 0.0;
    let mut da: f64 = 0.0;
    if (((c < ((0) as f64)) && (a <= c)) || ((c >= ((0) as f64)) && (a >= c))) {
        {
            da = (a - c).round();
        }
    } else {
        {
            da = (a).round();
        }
    }
    t = (a - da);
    loss[0] = ((0) as f64);
    debug_assert!((da != ((0) as f64)));
    if ((da).abs() > ((HYP2F1_MAXITER) as f64)) {
        {
            loss[0] = 1.0f64;
            return f64::NAN;
        }
    }
    if (da < ((0) as f64)) {
        {
            f2 = ((0) as f64);
            f1 = hys2f1(t, b, c, x, std::slice::from_mut(&mut err));
            loss[0] += err;
            f0 = hys2f1((t - ((1) as f64)), b, c, x, std::slice::from_mut(&mut err));
            loss[0] += err;
            t -= ((1) as f64);
            n = 1;
            while (((n) as f64) < (-da)) {
                {
                    f2 = f1;
                    f1 = f0;
                    f0 = ((((-((((((2) as f64) * t) - c) - (t * x)) + (b * x))) / (c - t)) * f1)
                        - (((t * (x - ((1) as f64))) / (c - t)) * f2));
                    t -= ((1) as f64);
                }
                {
                    n += 1;
                    n
                };
            }
        }
    } else {
        {
            f2 = ((0) as f64);
            f1 = hys2f1(t, b, c, x, std::slice::from_mut(&mut err));
            loss[0] += err;
            f0 = hys2f1((t + ((1) as f64)), b, c, x, std::slice::from_mut(&mut err));
            loss[0] += err;
            t += ((1) as f64);
            n = 1;
            while (((n) as f64) < da) {
                {
                    f2 = f1;
                    f1 = f0;
                    f0 = ((-((((((((2) as f64) * t) - c) - (t * x)) + (b * x)) * f1)
                        + ((c - t) * f2)))
                        / (t * (x - ((1) as f64))));
                    t += ((1) as f64);
                }
                {
                    n += 1;
                    n
                };
            }
        }
    }
    return f0;
}

pub(super) fn chgm(mut x: f64, mut a: f64, mut b: f64) -> f64 {
    let mut i: i32 = 0;
    let mut j: i32 = 0;
    let mut la: i32 = 0;
    let mut n: i32 = 0;
    let mut nl: i32 = 0;
    let mut a0: f64 = a;
    let mut a1: f64 = a;
    let mut x0: f64 = x;
    let mut y0: f64 = 0.0;
    let mut y1: f64 = 0.0;
    let mut hg1: f64 = 0.0;
    let mut hg2: f64 = 0.0;
    let mut r1: f64 = 0.0;
    let mut r2: f64 = 0.0;
    let mut rg: f64 = 0.0;
    let mut xg: f64 = 0.0;
    let mut sum1: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut pi: f64 = 3.141592653589793f64;
    let mut hg: f64 = 0.0f64;
    if (x < 0.0f64) {
        {
            a = (b - a);
            a0 = a;
            x = (x).abs();
        }
    }
    nl = 0;
    la = 0;
    if (a >= 2.0f64) {
        {
            nl = 1;
            la = ((a) as i32);
            a -= ((la + 1) as f64);
        }
    }
    y0 = 0.0f64;
    y1 = 0.0f64;
    n = 0;
    while (n < (nl + 1)) {
        {
            if (a0 >= 2.0f64) {
                {
                    a += 1.0f64;
                }
            }
            if ((x <= (30.0f64 + (b).abs())) || (a < 0.0f64)) {
                {
                    hg = 1.0f64;
                    rg = 1.0f64;
                    j = 1;
                    while (j < 501) {
                        {
                            rg = (((rg * ((a + ((j) as f64)) - 1.0f64))
                                / (((j) as f64) * ((b + ((j) as f64)) - 1.0f64)))
                                * x);
                            hg += rg;
                            if ((hg != 0.0f64) && ((rg / hg).abs() < 1e-15f64)) {
                                {
                                    if (x0 < 0.0f64) {
                                        {
                                            hg *= (x0).exp();
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                        {
                            j += 1;
                            j
                        };
                    }
                }
            } else {
                {
                    xg = (b - a);
                    sum1 = 1.0f64;
                    sum2 = 1.0f64;
                    r1 = 1.0f64;
                    r2 = 1.0f64;
                    i = 1;
                    while (i < 17) {
                        {
                            r1 = ((((-r1) * ((a + ((i) as f64)) - 1.0f64))
                                * ((a - b) + ((i) as f64)))
                                / (x * ((i) as f64)));
                            r2 = ((((-r2) * (((b - a) + ((i) as f64)) - 1.0f64))
                                * (a - ((i) as f64)))
                                / (x * ((i) as f64)));
                            sum1 += r1;
                            sum2 += r2;
                        }
                        {
                            i += 1;
                            i
                        };
                    }
                    if (x0 >= 0.0f64) {
                        {
                            hg1 = (((gamma_ratio(b, xg, 0.0f64) * (x).powf((-a)))
                                * (pi * a).cos())
                                * sum1);
                            hg2 = ((gamma_ratio(b, a, x) * (x).powf((a - b))) * sum2);
                        }
                    } else {
                        {
                            hg1 = (((gamma_ratio(b, xg, x0) * (x).powf((-a))) * (pi * a).cos())
                                * sum1);
                            hg2 = ((gamma_ratio(b, a, 0.0f64) * (x).powf((a - b))) * sum2);
                        }
                    }
                    hg = (hg1 + hg2);
                }
            }
            if (n == 0) {
                {
                    y0 = hg;
                }
            }
            if (n == 1) {
                {
                    y1 = hg;
                }
            }
        }
        {
            n += 1;
            n
        };
    }
    if (a0 >= 2.0f64) {
        {
            i = 1;
            while (i < la) {
                {
                    hg = ((((((2.0f64 * a) - b) + x) * y1) + ((b - a) * y0)) / a);
                    y0 = y1;
                    y1 = hg;
                    a += 1.0f64;
                }
                {
                    i += 1;
                    i
                };
            }
        }
    }
    a = a1;
    x = x0;
    return hg;
}

// U(-n,b,x) = (-1)^n n! L_n^(b-1)(x), where n is a nonnegative integer.
// The Laguerre three-term recurrence avoids cancellation between the two
// hypergeometric series in the legacy chgu dispatch. Keep its use to a
// checked moderate-parameter region: higher orders/arguments may overflow
// intermediates or accumulate substantial phase error near polynomial zeros.
fn hyperu_neg_integer_polynomial(n: usize, b: f64, x: f64) -> Option<f64> {
    if n == 0 {
        return Some(1.0);
    }
    let mut prev = 1.0;
    let mut curr = x - b;
    for k in 2..=n {
        let k = k as f64;
        let next = (x - (2.0 * k + b - 2.0)).mul_add(curr, -(k - 1.0) * (k + b - 2.0) * prev);
        if !next.is_finite() {
            return None;
        }
        prev = curr;
        curr = next;
    }
    Some(curr)
}

pub(super) fn chgu(mut x: f64, mut a: f64, mut b: f64, md: &mut [i32], isfer: &mut [i32]) -> f64 {
    if a < 0.0 && a == a.trunc() && a >= -128.0 && (0.0..=16.0).contains(&b) && x > 0.0 && x <= 10.0
    {
        if let Some(value) = hyperu_neg_integer_polynomial((-a) as usize, b, x) {
            md[0] = 2;
            isfer[0] = 0;
            return value;
        }
    }
    let mut il1: i32 = 0;
    let mut il2: i32 = 0;
    let mut il3: i32 = 0;
    let mut bl1: i32 = 0;
    let mut bl2: i32 = 0;
    let mut bl3: i32 = 0;
    let mut bn: i32 = 0;
    let mut id1: i32 = 0;
    let mut id: i32 = 0;
    let mut aa: f64 = 0.0;
    let mut hu: f64 = 0.0f64;
    let mut hu1: f64 = 0.0;
    let mut b00: f64 = 0.0;
    aa = ((a - b) + 1.0f64);
    isfer[0] = 0;
    il1 = ((a == (((a) as i32) as f64)) && (a <= 0.0f64)) as i32;
    il2 = ((aa == (((aa) as i32) as f64)) && (aa <= 0.0f64)) as i32;
    il3 = (((a * ((a - b) + 1.0f64)).abs() / x) <= 2.0f64) as i32;
    bl1 = ((x <= 5.0f64) || ((x <= 10.0f64) && (a <= 2.0f64))) as i32;
    bl2 = (((x > 5.0f64) && (x <= 12.5f64)) && ((a >= 1.0f64) && (b >= (a + 4.0f64)))) as i32;
    bl3 = (((x > 12.5f64) && (a >= 5.0f64)) && (b >= (a + 5.0f64))) as i32;
    bn = ((b == (((b) as i32) as f64)) && (b != 0.0f64)) as i32;
    id = (-100);
    hu1 = 0.0f64;
    if (b != (((b) as i32) as f64)) {
        {
            hu = chgus(x, a, b, std::slice::from_mut(&mut id1));
            md[0] = 1;
            if (id1 >= 9) {
                {
                    return hu;
                }
            }
            hu1 = hu;
        }
    }
    if (((il1 != 0) || (il2 != 0)) || (il3 != 0)) {
        {
            hu = chgul(x, a, b, std::slice::from_mut(&mut id));
            md[0] = 2;
            if (id >= 9) {
                {
                    return hu;
                }
            }
            if (id1 > id) {
                {
                    md[0] = 1;
                    id = id1;
                    hu = hu1;
                }
            }
        }
    }
    if (a >= 1.0f64) {
        {
            if ((bn != 0) && (((bl1 != 0) || (bl2 != 0)) || (bl3 != 0))) {
                {
                    hu = chgubi(x, a, b, std::slice::from_mut(&mut id));
                    md[0] = 3;
                }
            } else {
                {
                    hu = chguit(x, a, b, std::slice::from_mut(&mut id));
                    md[0] = 4;
                }
            }
        }
    } else {
        {
            if (b <= a) {
                {
                    b00 = b;
                    a -= (b - 1.0f64);
                    b = (2.0f64 - b);
                    hu = chguit(x, a, b, std::slice::from_mut(&mut id));
                    hu *= (x).powf((1.0f64 - b00));
                    md[0] = 4;
                }
            } else {
                if ((bn != 0) && (il1 == 0)) {
                    {
                        hu = chgubi(x, a, b, std::slice::from_mut(&mut id));
                        md[0] = 3;
                    }
                }
            }
        }
    }
    if (id < 6) {
        {
            isfer[0] = 6;
        }
    }
    return hu;
}

pub(super) fn chgus(mut x: f64, mut a: f64, mut b: f64, id: &mut [i32]) -> f64 {
    let mut j: i32 = 0;
    let mut d1: f64 = 0.0;
    let mut d2: f64 = 0.0;
    let mut ga: f64 = 0.0;
    let mut gb: f64 = 0.0;
    let mut gab: f64 = 0.0;
    let mut gb2: f64 = 0.0;
    let mut h0: f64 = 0.0;
    let mut hmax: f64 = 0.0;
    let mut hmin: f64 = 0.0;
    let mut hu: f64 = 0.0;
    let mut hu0: f64 = 0.0;
    let mut hua: f64 = 0.0;
    let mut r1: f64 = 0.0;
    let mut r2: f64 = 0.0;
    let mut pi: f64 = 3.141592653589793f64;
    id[0] = 100;
    ga = gamma2(a);
    gb = gamma2(b);
    gab = gamma2(((1.0f64 + a) - b));
    gb2 = gamma2((2.0f64 - b));
    hu0 = (pi / (pi * b).sin());
    r1 = (hu0 / (gab * gb));
    r2 = ((hu0 * (x).powf((1.0f64 - b))) / (ga * gb2));
    hu = (r1 - r2);
    hmax = 0.0f64;
    hmin = 1e300f64;
    h0 = 0.0f64;
    j = 1;
    while (j < 151) {
        {
            r1 = (((r1 * ((a + ((j) as f64)) - 1.0f64))
                / (((j) as f64) * ((b + ((j) as f64)) - 1.0f64)))
                * x);
            r2 = (((r2 * ((a - b) + ((j) as f64)))
                / (((j) as f64) * ((1.0f64 - b) + ((j) as f64))))
                * x);
            hu += (r1 - r2);
            hua = (hu).abs();
            if (hua > hmax) {
                {
                    hmax = hua;
                }
            }
            if (hua < hmin) {
                {
                    hmin = hua;
                }
            }
            if ((hu - h0).abs() < ((hu).abs() * 1e-15f64)) {
                {
                    break;
                }
            }
            h0 = hu;
        }
        {
            j += 1;
            j
        };
    }
    d1 = (hmax).log10();
    d2 = 0.0f64;
    if (hmin != 0.0f64) {
        {
            d2 = (hmin).log10();
        }
    }
    id[0] = ((((15) as f64)
        - if ((d1 - d2) < ((0) as f64)) {
            (d2 - d1)
        } else {
            (d1 - d2)
        }) as i32);
    return hu;
}

pub(super) fn gamma2(mut x: f64) -> f64 {
    let mut ga: f64 = 0.0;
    let mut gr: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut z: f64 = 0.0;
    let mut k: i32 = 0;
    let mut m: i32 = 0;
    let mut pi: f64 = 3.141592653589793f64;
    let mut g = [
        1.0000000000000000e+00f64,
        0.5772156649015329e+00f64,
        (-0.6558780715202538e+00f64),
        (-0.4200263503409520e-01f64),
        0.1665386113822915e+00f64,
        (-0.4219773455554430e-01f64),
        (-0.9621971527877000e-02f64),
        0.7218943246663000e-02f64,
        (-0.1165167591859100e-02f64),
        (-0.2152416741149000e-03f64),
        0.1280502823882000e-03f64,
        (-0.2013485478070000e-04f64),
        (-0.1250493482100000e-05f64),
        0.1133027232000000e-05f64,
        (-0.2056338417000000e-06f64),
        0.6116095000000000e-08f64,
        0.5002007500000000e-08f64,
        (-0.1181274600000000e-08f64),
        0.1043427000000000e-09f64,
        0.7782300000000000e-11f64,
        (-0.3696800000000000e-11f64),
        0.5100000000000000e-12f64,
        (-0.2060000000000000e-13f64),
        (-0.5400000000000000e-14f64),
        0.1400000000000000e-14f64,
        0.1000000000000000e-15f64,
    ];
    if (x == (((x) as i32) as f64)) {
        {
            if (x > 0.0f64) {
                {
                    ga = 1.0f64;
                    m = ((x - ((1) as f64)) as i32);
                    k = 2;
                    while (k < (m + 1)) {
                        {
                            ga *= ((k) as f64);
                        }
                        {
                            k += 1;
                            k
                        };
                    }
                }
            } else {
                {
                    ga = 1e300f64;
                }
            }
        }
    } else {
        {
            r = 1.0f64;
            if ((x).abs() > 1.0f64) {
                {
                    z = (x).abs();
                    m = ((z) as i32);
                    k = 1;
                    while (k < (m + 1)) {
                        {
                            r *= (z - ((k) as f64));
                        }
                        {
                            k += 1;
                            k
                        };
                    }
                    z -= ((m) as f64);
                }
            } else {
                {
                    z = x;
                }
            }
            gr = g[(25) as usize];
            k = 25;
            while (k > 0) {
                {
                    gr *= z;
                    gr += g[(k - 1) as usize];
                }
                {
                    k -= 1;
                    k
                };
            }
            ga = (1.0f64 / (gr * z));
            if ((x).abs() > 1.0f64) {
                {
                    ga *= r;
                    if (x < 0.0f64) {
                        {
                            ga = ((-pi) / ((x * ga) * (pi * x).sin()));
                        }
                    }
                }
            }
        }
    }
    return ga;
}

pub(super) fn chgul(mut x: f64, mut a: f64, mut b: f64, id: &mut [i32]) -> f64 {
    let mut il1: i32 = 0;
    let mut il2: i32 = 0;
    let mut k: i32 = 0;
    let mut nm: i32 = 0;
    let mut aa: f64 = 0.0;
    let mut hu: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut r0: f64 = 0.0f64;
    let mut ra: f64 = 0.0f64;
    id[0] = (-100);
    aa = ((a - b) + 1.0f64);
    il1 = ((a == (((a) as i32) as f64)) && (a <= 0.0f64)) as i32;
    il2 = ((aa == (((aa) as i32) as f64)) && (aa <= 0.0f64)) as i32;
    nm = 0;
    if (il1 != 0) {
        {
            nm = (((a).abs()) as i32);
        }
    }
    if (il2 != 0) {
        {
            nm = (((aa).abs()) as i32);
        }
    }
    if ((il1 != 0) || (il2 != 0)) {
        {
            hu = 1.0f64;
            r = 1.0f64;
            k = 1;
            while (k <= nm) {
                {
                    r = ((((-r) * ((a + ((k) as f64)) - 1.0f64)) * ((a - b) + ((k) as f64)))
                        / (((k) as f64) * x));
                    hu += r;
                }
                {
                    k += 1;
                    k
                };
            }
            hu *= (x).powf((-a));
            id[0] = 10;
        }
    } else {
        {
            hu = 1.0f64;
            r = 1.0f64;
            k = 1;
            while (k <= 25) {
                {
                    r = ((((-r) * ((a + ((k) as f64)) - 1.0f64)) * ((a - b) + ((k) as f64)))
                        / (((k) as f64) * x));
                    ra = (r).abs();
                    if (((k > 5) && (ra >= r0)) || (ra < 1e-15f64)) {
                        {
                            break;
                        }
                    }
                    r0 = ra;
                    hu += r;
                }
                {
                    k += 1;
                    k
                };
            }
            id[0] = ((((ra).log10()).abs()) as i32);
            hu *= (x).powf((-a));
        }
    }
    return hu;
}

pub(super) fn chgubi(mut x: f64, mut a: f64, mut b: f64, id: &mut [i32]) -> f64 {
    let mut id1: i32 = 0;
    let mut id2: i32 = 0;
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut m: i32 = 0;
    let mut n: i32 = 0;
    let mut a0: f64 = 0.0;
    let mut a1: f64 = 0.0;
    let mut a2: f64 = 0.0;
    let mut da1: f64 = 0.0;
    let mut da2: f64 = 0.0;
    let mut db1: f64 = 0.0;
    let mut db2: f64 = 0.0;
    let mut ga: f64 = 0.0;
    let mut ga1: f64 = 0.0;
    let mut h0: f64 = 0.0;
    let mut hm1: f64 = 0.0;
    let mut hm2: f64 = 0.0;
    let mut hm3: f64 = 0.0;
    let mut hmax: f64 = 0.0;
    let mut hmin: f64 = 0.0;
    let mut hu: f64 = 0.0;
    let mut hu1: f64 = 0.0;
    let mut hu2: f64 = 0.0;
    let mut hw: f64 = 0.0;
    let mut ps: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut rn: f64 = 0.0;
    let mut rn1: f64 = 0.0;
    let mut s0: f64 = 0.0;
    let mut s1: f64 = 0.0;
    let mut s2: f64 = 0.0;
    let mut sa: f64 = 0.0;
    let mut sb: f64 = 0.0;
    let mut ua: f64 = 0.0;
    let mut ub: f64 = 0.0;
    let mut el: f64 = 0.5772156649015329f64;
    id[0] = (-100);
    n = (((b - ((1) as f64)).abs()) as i32);
    rn1 = 1.0f64;
    rn = 1.0f64;
    j = 1;
    while (j <= n) {
        {
            rn *= ((j) as f64);
            if (j == (n - 1)) {
                {
                    rn1 = rn;
                }
            }
        }
        {
            j += 1;
            j
        };
    }
    ps = psi_spec(a);
    ga = gamma2(a);
    if (b > 0.0f64) {
        {
            a0 = a;
            a1 = (a - ((n) as f64));
            a2 = a1;
            ga1 = gamma2(a1);
            ua = (neg_one_pow(n - 1) / (rn * ga1));
            ub = ((rn1 / ga) * (x).powf(((-n) as f64)));
        }
    } else {
        {
            a0 = (a + ((n) as f64));
            a1 = a0;
            a2 = a;
            ga1 = gamma2(a1);
            ua = ((neg_one_pow(n - 1) / (rn * ga)) * (x).powf(((n) as f64)));
            ub = (rn1 / ga1);
        }
    }
    hm1 = 1.0f64;
    r = 1.0f64;
    hmax = 0.0f64;
    hmin = 1e300f64;
    h0 = 0.0f64;
    k = 1;
    while (k <= 150) {
        {
            r = (((r * ((a0 + ((k) as f64)) - ((1) as f64))) * x) / (((n + k) * k) as f64));
            hm1 += r;
            hu1 = (hm1).abs();
            if (hu1 > hmax) {
                {
                    hmax = hu1;
                }
            }
            if (hu1 < hmin) {
                {
                    hmin = hu1;
                }
            }
            if ((hm1 - h0).abs() < ((hm1).abs() * 1.0e-15f64)) {
                {
                    break;
                }
            }
            h0 = hm1;
        }
        {
            k += 1;
            k
        };
    }
    da1 = (hmax).log10();
    da2 = ((0) as f64);
    if (hmin != ((0) as f64)) {
        {
            da2 = (hmin).log10();
        }
    }
    id[0] = (15 - (((da1 - da2).abs()) as i32));
    hm1 *= (x).ln();
    s0 = ((0) as f64);
    m = 1;
    while (m <= n) {
        {
            if (b >= ((0) as f64)) {
                {
                    s0 -= (1.0f64 / ((m) as f64));
                }
            }
            if (b < ((0) as f64)) {
                {
                    s0 += ((1.0f64 - a) / (((m) as f64) * ((a + ((m) as f64)) - ((1) as f64))));
                }
            }
        }
        {
            m += 1;
            m
        };
    }
    hm2 = ((ps + (((2) as f64) * el)) + s0);
    r = ((1) as f64);
    hmax = ((0) as f64);
    hmin = 1.0e+300f64;
    // Extend sums across terms instead of restarting from m=1 at each k.
    // For b<0, s0 is precisely the starting prefix of s1.
    let mut running_s1 = if b > 0.0 { 0.0 } else { s0 };
    let mut running_s2 = 0.0;
    if b > 0.0 && n > 8 {
        for m in 1..=n {
            running_s2 += 1.0 / (m + 1) as f64;
        }
    }
    k = 1;
    while (k <= 150) {
        {
            if (b > ((0) as f64)) {
                {
                    let kf = k as f64;
                    running_s1 -= ((kf + 2.0 * a) - 2.0) / (kf * ((kf + a) - 1.0));
                    s1 = running_s1;
                    s2 = if n > 8 {
                        running_s2
                    } else {
                        let mut sum = 0.0;
                        for m in 1..=n {
                            sum += 1.0 / (k + m) as f64;
                        }
                        sum
                    };
                }
            } else {
                {
                    let mf = (n + k) as f64;
                    running_s1 += (1.0 - a) / (mf * ((mf + a) - 1.0));
                    running_s2 += 1.0 / k as f64;
                    s1 = running_s1;
                    s2 = running_s2;
                }
            }
            hw = ((((((2) as f64) * el) + ps) + s1) - s2);
            r = (((r * ((a0 + ((k) as f64)) - ((1) as f64))) * x) / (((n + k) * k) as f64));
            hm2 += (r * hw);
            hu2 = (hm2).abs();
            if (hu2 > hmax) {
                {
                    hmax = hu2;
                }
            }
            if (hu2 < hmin) {
                {
                    hmin = hu2;
                }
            }
            if (((hm2 - h0) / hm2).abs() < 1.0e-15f64) {
                {
                    break;
                }
            }
            h0 = hm2;
            if b > 0.0 && n > 8 {
                running_s2 += 1.0 / (k + n + 1) as f64 - 1.0 / (k + 1) as f64;
            }
        }
        {
            k += 1;
            k
        };
    }
    db1 = (hmax).log10();
    db2 = 0.0f64;
    if (hmin != 0.0f64) {
        {
            db2 = (hmin).log10();
        }
    }
    id1 = (15 - (((db1 - db2).abs()) as i32));
    if (id1 < id[0]) {
        {
            id[0] = id1;
        }
    }
    hm3 = 1.0f64;
    if (n == 0) {
        {
            hm3 = 0.0f64;
        }
    }
    r = 1.0f64;
    k = 1;
    while (k < n) {
        {
            r = (((r * ((a2 + ((k) as f64)) - 1.0f64)) / (((k - n) * k) as f64)) * x);
            hm3 += r;
        }
        {
            k += 1;
            k
        };
    }
    sa = (ua * (hm1 + hm2));
    sb = (ub * hm3);
    hu = (sa + sb);
    id2 = 0;
    if (sa != 0.0f64) {
        {
            id1 = ((((sa).abs()).log10()) as i32);
        }
    }
    if (hu != 0.0f64) {
        {
            id2 = ((((hu).abs()).log10()) as i32);
        }
    }
    if ((sa * sb) < 0.0f64) {
        {
            id[0] -= ((((id1 - id2) as f64).abs()) as i32);
        }
    }
    return hu;
}

pub(super) fn psi_spec(mut x: f64) -> f64 {
    let mut k: i32 = 0;
    let mut n: i32 = 0;
    let mut ps: f64 = 0.0;
    let mut s: f64 = 0.0f64;
    let mut x2: f64 = 0.0;
    let mut xa: f64 = (x).abs();
    let mut pi: f64 = 3.141592653589793f64;
    let mut el: f64 = 0.5772156649015329f64;
    let mut a = [
        (-0.8333333333333e-01f64),
        0.83333333333333333e-02f64,
        (-0.39682539682539683e-02f64),
        0.41666666666666667e-02f64,
        (-0.75757575757575758e-02f64),
        0.21092796092796093e-01f64,
        (-0.83333333333333333e-01f64),
        0.4432598039215686f64,
    ];
    if ((x == (((x) as i32) as f64)) && (x <= 0.0f64)) {
        {
            return 1e300f64;
        }
    } else {
        if (xa == (((xa) as i32) as f64)) {
            {
                n = ((xa) as i32);
                k = 1;
                while (k < n) {
                    {
                        s += (1.0f64 / ((k) as f64));
                    }
                    {
                        k += 1;
                        k
                    };
                }
                ps = ((-el) + s);
            }
        } else {
            if ((xa + 0.5f64) == (((xa + 0.5f64) as i32) as f64)) {
                {
                    n = ((xa - 0.5f64) as i32);
                    k = 1;
                    while (k < (n + 1)) {
                        {
                            s += (1.0f64 / ((2.0f64 * ((k) as f64)) - 1.0f64));
                        }
                        {
                            k += 1;
                            k
                        };
                    }
                    ps = (((-el) + (2.0f64 * s)) - 1.386294361119891f64);
                }
            } else {
                {
                    if (xa < 10.0f64) {
                        {
                            n = ((10.0f64 - (((xa) as i32) as f64)) as i32);
                            k = 0;
                            while (k < n) {
                                {
                                    s += (1.0f64 / (xa + ((k) as f64)));
                                }
                                {
                                    k += 1;
                                    k
                                };
                            }
                            xa += ((n) as f64);
                        }
                    }
                    x2 = (1.0f64 / (xa * xa));
                    ps = ((xa).ln() - (0.5f64 / xa));
                    ps += (x2
                        * ((((((((((((((a[(7) as usize] * x2) + a[(6) as usize]) * x2)
                            + a[(5) as usize])
                            * x2)
                            + a[(4) as usize])
                            * x2)
                            + a[(3) as usize])
                            * x2)
                            + a[(2) as usize])
                            * x2)
                            + a[(1) as usize])
                            * x2)
                            + a[(0) as usize]));
                    ps -= s;
                }
            }
        }
    }
    if (x < 0.0f64) {
        {
            ps -= (((pi * (pi * x).cos()) / (pi * x).sin()) + (1.0f64 / x));
        }
    }
    return ps;
}

pub(super) fn chguit(mut x: f64, mut a: f64, mut b: f64, id: &mut [i32]) -> f64 {
    let mut k: i32 = 0;
    let mut j: i32 = 0;
    let mut m: i32 = 0;
    let mut a1: f64 = 0.0;
    let mut b1: f64 = 0.0;
    let mut c: f64 = 0.0;
    let mut d: f64 = 0.0;
    let mut f1: f64 = 0.0;
    let mut f2: f64 = 0.0;
    let mut g: f64 = 0.0;
    let mut ga: f64 = 0.0;
    let mut hu: f64 = 0.0;
    let mut hu0: f64 = 0.0;
    let mut hu1: f64 = 0.0;
    let mut hu2: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut t1: f64 = 0.0;
    let mut t2: f64 = 0.0;
    let mut t3: f64 = 0.0;
    let mut t4: f64 = 0.0;
    let mut t = [
        0.259597723012478e-01f64,
        0.778093339495366e-01f64,
        0.129449135396945e+00f64,
        0.180739964873425e+00f64,
        0.231543551376029e+00f64,
        0.281722937423262e+00f64,
        0.331142848268448e+00f64,
        0.379670056576798e+00f64,
        0.427173741583078e+00f64,
        0.473525841761707e+00f64,
        0.518601400058570e+00f64,
        0.562278900753945e+00f64,
        0.604440597048510e+00f64,
        0.644972828489477e+00f64,
        0.683766327381356e+00f64,
        0.720716513355730e+00f64,
        0.755723775306586e+00f64,
        0.788693739932264e+00f64,
        0.819537526162146e+00f64,
        0.848171984785930e+00f64,
        0.874519922646898e+00f64,
        0.898510310810046e+00f64,
        0.920078476177628e+00f64,
        0.939166276116423e+00f64,
        0.955722255839996e+00f64,
        0.969701788765053e+00f64,
        0.981067201752598e+00f64,
        0.989787895222222e+00f64,
        0.995840525118838e+00f64,
        0.999210123227436e+00f64,
    ];
    let mut w = [
        0.519078776312206e-01f64,
        0.517679431749102e-01f64,
        0.514884515009810e-01f64,
        0.510701560698557e-01f64,
        0.505141845325094e-01f64,
        0.498220356905502e-01f64,
        0.489955754557568e-01f64,
        0.480370318199712e-01f64,
        0.469489888489122e-01f64,
        0.457343797161145e-01f64,
        0.443964787957872e-01f64,
        0.429388928359356e-01f64,
        0.413655512355848e-01f64,
        0.396806954523808e-01f64,
        0.378888675692434e-01f64,
        0.359948980510845e-01f64,
        0.340038927249464e-01f64,
        0.319212190192963e-01f64,
        0.297524915007890e-01f64,
        0.275035567499248e-01f64,
        0.251804776215213e-01f64,
        0.227895169439978e-01f64,
        0.203371207294572e-01f64,
        0.178299010142074e-01f64,
        0.152746185967848e-01f64,
        0.126781664768159e-01f64,
        0.100475571822880e-01f64,
        0.738993116334531e-02f64,
        0.471272992695363e-02f64,
        0.202681196887362e-02f64,
    ];
    id[0] = 9;
    a1 = (a - 1.0f64);
    b1 = ((b - a) - 1.0f64);
    c = (12.0f64 / x);
    hu0 = 0.0f64;
    m = 10;
    while (m <= 100) {
        {
            hu1 = 0.0f64;
            g = ((0.5f64 * c) / ((m) as f64));
            d = g;
            j = 1;
            while (j < (m + 1)) {
                {
                    s = 0.0f64;
                    k = 1;
                    while (k <= 30) {
                        {
                            t1 = (d + (g * t[(k - 1) as usize]));
                            t2 = (d - (g * t[(k - 1) as usize]));
                            f1 = ((((-x) * t1).exp() * (t1).powf(a1)) * (1.0f64 + t1).powf(b1));
                            f2 = ((((-x) * t2).exp() * (t2).powf(a1)) * (1.0f64 + t2).powf(b1));
                            s += (w[(k - 1) as usize] * (f1 + f2));
                        }
                        {
                            k += 1;
                            k
                        };
                    }
                    hu1 += (s * g);
                    d += (2.0f64 * g);
                }
                {
                    j += 1;
                    j
                };
            }
            if ((1.0f64 - (hu0 / hu1)).abs() < 1.0e-9f64) {
                {
                    break;
                }
            }
            hu0 = hu1;
        }
        m += 5;
    }
    ga = gamma2(a);
    hu1 /= ga;
    m = 2;
    while (m <= 10) {
        {
            hu2 = 0.0f64;
            g = (0.5f64 / ((m) as f64));
            d = g;
            j = 1;
            while (j <= m) {
                {
                    s = 0.0f64;
                    k = 1;
                    while (k <= 30) {
                        {
                            t1 = (d + (g * t[(k - 1) as usize]));
                            t2 = (d - (g * t[(k - 1) as usize]));
                            t3 = (c / (1.0f64 - t1));
                            t4 = (c / (1.0f64 - t2));
                            f1 = (((((t3 * t3) / c) * ((-x) * t3).exp()) * (t3).powf(a1))
                                * (1.0f64 + t3).powf(b1));
                            f2 = (((((t4 * t4) / c) * ((-x) * t4).exp()) * (t4).powf(a1))
                                * (1.0f64 + t4).powf(b1));
                            s += (w[(k - 1) as usize] * (f1 + f2));
                        }
                        {
                            k += 1;
                            k
                        };
                    }
                    hu2 += (s * g);
                    d += (2.0f64 * g);
                }
                {
                    j += 1;
                    j
                };
            }
            if ((1.0f64 - (hu0 / hu2)).abs() < 1.0e-9f64) {
                {
                    break;
                }
            }
            hu0 = hu2;
        }
        m += 2;
    }
    ga = gamma2(a);
    hu2 /= ga;
    hu = (hu1 + hu2);
    return hu;
}

#[cfg(test)]
mod optimization_tests {
    use super::*;

    // Original integer-d continuation of hyt2f1, using four independently
    // evaluated digamma values in each term.
    fn reference_hyt2f1_integer_d(a: f64, b: f64, c: f64, x: f64) -> f64 {
        let d = (c - a) - b;
        let id = d.round();
        assert_eq!(d, id);
        let s = 1.0 - x;
        let (e, d1, d2, aid) = if id >= 0.0 {
            (d, d, 0.0, id as i32)
        } else {
            (-d, 0.0, d, (-id) as i32)
        };
        let ax = s.ln();
        let mut y = (((psi(1.0) + psi(1.0 + e)) - psi(a + d1)) - psi(b + d1)) - ax;
        y *= rgamma(e + 1.0);
        let mut p = ((a + d1) * (b + d1) * s) * rgamma(e + 2.0);
        let mut t = 1.0;
        loop {
            let r = (((psi(1.0 + t) + psi(1.0 + t + e)) - psi(a + t + d1)) - psi(b + t + d1)) - ax;
            let q = p * r;
            y += q;
            p *= (s * (a + t + d1)) / (t + 1.0);
            p *= (b + t + d1) / (t + 1.0 + e);
            t += 1.0;
            if t > HYP2F1_MAXITER as f64 {
                return f64::NAN;
            }
            if y != 0.0 && (q / y).abs() <= HYP2F1_EPS {
                break;
            }
        }
        if id == 0.0 {
            return y * (Gamma(c) / (Gamma(a) * Gamma(b)));
        }
        let mut y1 = 1.0;
        let mut t = 0.0;
        let mut p = 1.0;
        for _ in 1..aid {
            let r = (1.0 - e) + t;
            p *= (s * (a + t + d2) * (b + t + d2)) / r;
            t += 1.0;
            p /= t;
            y1 += p;
        }
        let g = Gamma(c);
        y1 *= (Gamma(e) * g) * (rgamma(a + d1) * rgamma(b + d1));
        y *= g * (rgamma(a + d2) * rgamma(b + d2));
        if (aid & 1) != 0 {
            y = -y;
        }
        let q = s.powf(id);
        if id > 0.0 {
            y *= q;
        } else {
            y1 *= q;
        }
        y + y1
    }

    #[test]
    fn hyt2f1_digamma_recurrence_matches_original_continuation() {
        for &(a, b, c) in &[
            (0.5, 0.5, 1.0),
            (0.5, 0.75, 1.25),
            (1.0, 1.0, 2.0),
            (0.5, 1.25, 2.75),
            (1.25, 2.0, 4.25),
            (0.5, 1.25, 9.75),
            (0.5, 0.75, 0.25),
            (1.5, 1.25, 0.75),
            (-0.25, 0.75, 0.5),
            (-0.75, 1.25, 0.5),
        ] {
            for &x in &[0.91, 0.98, 0.999, 0.99999999999] {
                let mut loss = [0.0];
                let got = hyt2f1(a, b, c, x, &mut loss);
                let reference = reference_hyt2f1_integer_d(a, b, c, x);
                assert!(
                    (got - reference).abs() <= 3e-12 * reference.abs().max(1.0),
                    "a={a} b={b} c={c} x={x}: {got} vs {reference}"
                );
            }
        }
    }

    #[test]
    fn integer_d_continuation_matches_logarithmic_identity() {
        // 2F1(1,1;2;x) = -ln(1-x)/x, from integrating the geometric series.
        for &x in &[0.91, 0.99, 0.999999, 0.99999999999] {
            let mut loss = [0.0];
            let result = hyt2f1(1.0, 1.0, 2.0, x, &mut loss);
            let expected = -(-x).ln_1p() / x;
            assert!((result - expected).abs() <= 2e-12 * expected.abs(), "x={x}");
        }
    }

    // Pre-optimization chgubi with the original nested sums. Keep this
    // independent of the streaming implementation to catch changes in its
    // convergence decisions, final value, and reported number of digits.
    fn reference_chgubi(x: f64, a: f64, b: f64) -> (f64, i32) {
        let n = (b - 1.0).abs() as i32;
        let (mut rn, mut rn1) = (1.0, 1.0);
        for j in 1..=n {
            rn *= j as f64;
            if j == n - 1 {
                rn1 = rn;
            }
        }
        let ps = psi_spec(a);
        let ga = gamma2(a);
        let (a0, a2, ua, ub) = if b > 0.0 {
            let a1 = a - n as f64;
            (
                a,
                a1,
                neg_one_pow(n - 1) / (rn * gamma2(a1)),
                (rn1 / ga) * x.powf(-(n as f64)),
            )
        } else {
            let a1 = a + n as f64;
            (
                a1,
                a,
                (neg_one_pow(n - 1) / (rn * ga)) * x.powf(n as f64),
                rn1 / gamma2(a1),
            )
        };
        let (mut hm1, mut r, mut hmax, mut hmin, mut h0) =
            (1.0_f64, 1.0_f64, 0.0_f64, 1e300_f64, 0.0_f64);
        for k in 1..=150 {
            r = ((r * (a0 + k as f64 - 1.0)) * x) / ((n + k) * k) as f64;
            hm1 += r;
            hmax = hmax.max(hm1.abs());
            hmin = hmin.min(hm1.abs());
            if (hm1 - h0).abs() < hm1.abs() * 1e-15 {
                break;
            }
            h0 = hm1;
        }
        let da1 = hmax.log10();
        let da2 = if hmin != 0.0 { hmin.log10() } else { 0.0 };
        let mut digits = 15 - ((da1 - da2).abs() as i32);
        hm1 *= x.ln();

        let mut s0 = 0.0;
        for m in 1..=n {
            let mf = m as f64;
            if b >= 0.0 {
                s0 -= 1.0 / mf;
            } else {
                s0 += (1.0 - a) / (mf * ((a + mf) - 1.0));
            }
        }
        let mut hm2 = (ps + 2.0 * 0.5772156649015329) + s0;
        r = 1.0;
        hmax = 0.0;
        hmin = 1e300;
        for k in 1..=150 {
            let (mut s1, mut s2) = (0.0, 0.0);
            if b > 0.0 {
                for m in 1..=k {
                    let mf = m as f64;
                    s1 -= ((mf + 2.0 * a) - 2.0) / (mf * ((mf + a) - 1.0));
                }
                for m in 1..=n {
                    s2 += 1.0 / (k + m) as f64;
                }
            } else {
                for m in 1..=(k + n) {
                    let mf = m as f64;
                    s1 += (1.0 - a) / (mf * ((mf + a) - 1.0));
                }
                for m in 1..=k {
                    s2 += 1.0 / m as f64;
                }
            }
            let hw = ((2.0 * 0.5772156649015329 + ps) + s1) - s2;
            r = ((r * (a0 + k as f64 - 1.0)) * x) / ((n + k) * k) as f64;
            hm2 += r * hw;
            hmax = hmax.max(hm2.abs());
            hmin = hmin.min(hm2.abs());
            if ((hm2 - h0) / hm2).abs() < 1e-15 {
                break;
            }
            h0 = hm2;
        }
        let db1 = hmax.log10();
        let db2 = if hmin != 0.0 { hmin.log10() } else { 0.0 };
        digits = digits.min(15 - ((db1 - db2).abs() as i32));

        let mut hm3 = if n == 0 { 0.0 } else { 1.0 };
        r = 1.0;
        for k in 1..n {
            r = ((r * (a2 + k as f64 - 1.0)) / ((k - n) * k) as f64) * x;
            hm3 += r;
        }
        let sa = ua * (hm1 + hm2);
        let sb = ub * hm3;
        let hu = sa + sb;
        let mut id1 = 0;
        let mut id2 = 0;
        if sa != 0.0 {
            id1 = sa.abs().log10() as i32;
        }
        if hu != 0.0 {
            id2 = hu.abs().log10() as i32;
        }
        if sa * sb < 0.0 {
            digits -= (id1 - id2).abs();
        }
        (hu, digits)
    }

    #[test]
    fn chgubi_streaming_matches_original_nested_sums() {
        for &b in &[
            -16.0, -9.0, -3.0, -1.0, 1.0, 2.0, 3.0, 9.0, 10.0, 16.0, 32.0,
        ] {
            for &a in &[0.5, 1.25, 2.0, 2.75, 5.5] {
                for &x in &[0.125, 0.75, 2.0, 6.0] {
                    let (reference, expected_digits) = reference_chgubi(x, a, b);
                    let mut digits = [0];
                    let result = chgubi(x, a, b, &mut digits);
                    if reference.is_nan() {
                        assert!(result.is_nan(), "a={a} b={b} x={x}");
                    } else if reference.is_infinite() {
                        assert_eq!(result, reference, "a={a} b={b} x={x}");
                    } else {
                        assert!(
                            (result - reference).abs() <= 3e-12 * reference.abs().max(1.0),
                            "a={a} b={b} x={x}: {result} vs {reference}"
                        );
                    }
                    assert!(
                        (digits[0] - expected_digits).abs() <= 1,
                        "precision estimate a={a} b={b} x={x}: {} vs {expected_digits}",
                        digits[0]
                    );
                }
            }
        }
    }

    #[test]
    fn chgubi_positive_integer_b_matches_finite_laplace_sum() {
        // U(1,b,x) = integral exp(-xt)(1+t)^(b-2) dt for integer b>=2.
        // Expand the polynomial and integrate its powers independently.
        for b in [2_i32, 3, 9, 10] {
            for &x in &[1.0, 2.0, 5.0, 10.0] {
                let mut term = 1.0 / x;
                let mut expected = term;
                for k in 1..=(b - 2) {
                    term *= (b - 1 - k) as f64 / x;
                    expected += term;
                }
                let mut digits = [0];
                let result = chgubi(x, 1.0, b as f64, &mut digits);
                assert!(
                    (result - expected).abs() <= 2e-9 * expected,
                    "b={b} x={x}: {result} vs {expected}"
                );
            }
        }
    }
}
