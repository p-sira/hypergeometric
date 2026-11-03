//! Signed Lanczos log-Gamma, reciprocal gamma and real digamma.
//! Lanczos coefficients follow the public-domain numerical approximation
//! tabulated by Paul Godfrey; see DLMF §5 for Gamma reflection/recurrence.
use core::f64::consts::PI;

const COEFFS: [f64; 9] = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
];
fn positive_lgamma(x: f64) -> f64 {
    debug_assert!(x > 0.0);
    if x < 0.5 {
        return positive_lgamma(x + 1.0) - x.ln();
    }
    let z = x - 1.0;
    let mut sum = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        sum += c / (z + i as f64);
    }
    let t = z + 7.5;
    0.91893853320467274178032973640562 + (z + 0.5) * t.ln() - t + sum.ln()
}
pub(crate) fn lgamma_sign(x: f64) -> (f64, f64) {
    if x.is_nan() || x == f64::NEG_INFINITY || (x <= 0.0 && x == x.floor()) {
        return (f64::INFINITY, f64::NAN);
    }
    if x == f64::INFINITY {
        return (f64::INFINITY, 1.0);
    }
    if x > 0.0 {
        return (positive_lgamma(x), 1.0);
    }
    // Reduce the argument before multiplying by pi. Direct sin(pi*x) loses
    // nearly all relative precision for x next to a negative integer pole:
    // at x = -100 - 1 ulp, the rounding error in pi*x can exceed the phase.
    // A noninteger binary64 x has |x| < 2^52, so the nearest integer and its
    // parity are represented exactly.
    let nearest = x.round();
    let phase = (PI * (x - nearest)).sin();
    let s = if nearest % 2.0 == 0.0 { phase } else { -phase };
    if s == 0.0 {
        return (f64::INFINITY, f64::NAN);
    }
    (
        PI.ln() - s.abs().ln() - positive_lgamma(1.0 - x),
        s.signum(),
    )
}
#[inline]
pub(crate) fn gamma(x: f64) -> f64 {
    let (l, s) = lgamma_sign(x);
    s * l.exp()
}
#[inline]
pub(crate) fn rgamma(x: f64) -> f64 {
    if x <= 0.0 && x == x.floor() {
        return 0.0;
    }
    let (l, s) = lgamma_sign(x);
    s * (-l).exp()
}
pub(crate) fn gamma_ratio(a: f64, b: f64, extra: f64) -> f64 {
    let (la, sa) = lgamma_sign(a);
    let (lb, sb) = lgamma_sign(b);
    (la - lb + extra).exp() * sa * sb
}
pub(crate) fn lgam_sgn(x: f64, sign: &mut [i32]) -> f64 {
    let (l, s) = lgamma_sign(x);
    sign[0] = s as i32;
    l
}
/// Real digamma from reflection, recurrence, Bernoulli asymptotics.
pub(crate) fn digamma(x: f64) -> f64 {
    if !x.is_finite() || (x <= 0.0 && x == x.floor()) {
        return f64::NAN;
    }
    if x < 0.0 {
        // cot(pi*x) has period one. Reduce first to retain the binary64
        // distance to a negative integer pole (as with log-Gamma above).
        return digamma(1.0 - x) - PI / (PI * (x - x.round())).tan();
    }
    let (mut y, mut x) = (0.0, x);
    while x < 12.0 {
        y -= 1.0 / x;
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    y + x.ln()
        - 0.5 * inv
        - inv2
            * (1.0 / 12.0
                - inv2
                    * (1.0 / 120.0
                        - inv2 * (1.0 / 252.0 - inv2 * (1.0 / 240.0 - inv2 * (5.0 / 660.0)))))
}

#[cfg(test)]
mod tests {
    use super::{digamma, lgamma_sign};

    #[test]
    fn gamma_reflection_near_negative_integer_poles() {
        // 90-digit mpmath 1.3.0, evaluated at the exact binary64 inputs.
        for &(x, want_log_abs, want_sign, want_digamma) in &[
            (
                -1.0000000000000002,
                36.04365338911715,
                1.0,
                4503599627370496.0,
            ),
            (
                -0.9999999999999999,
                36.7368005696771,
                -1.0,
                -9007199254740992.0,
            ),
            (
                -100.00000000000001,
                -331.85460524980607,
                -1.0,
                70368744177668.61,
            ),
            (
                -99.99999999999999,
                -331.85460524980596,
                1.0,
                -70368744177659.39,
            ),
            (
                -1000.0000000000001,
                -5882.322849724086,
                -1.0,
                8796093022214.908,
            ),
            (
                -999.9999999999999,
                -5882.322849724085,
                1.0,
                -8796093022201.092,
            ),
        ] {
            let (got_log, got_sign) = lgamma_sign(x);
            assert_eq!(got_sign, want_sign, "log-gamma sign at {x}");
            assert!(
                (got_log - want_log_abs).abs() < 1e-12_f64.max(2e-14 * want_log_abs.abs()),
                "log-gamma at {x}: {got_log} vs {want_log_abs}"
            );
            let got_psi = digamma(x);
            assert!(
                (got_psi - want_digamma).abs() < 4e-15 * want_digamma.abs(),
                "digamma at {x}: {got_psi} vs {want_digamma}"
            );
        }
    }
}
