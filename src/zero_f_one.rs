//! Confluent limit function, compensated and double-double power series.
use crate::{dd::Dd, gamma};
type StrErr = &'static str;
const ERROR: StrErr = "0F1 arguments must be finite; parameter cannot be a nonpositive integer";
fn validate(v: f64) -> Result<(), StrErr> {
    if !v.is_finite() || (v <= 0.0 && v == v.floor()) {
        Err(ERROR)
    } else {
        Ok(())
    }
}
/// Double-double power series for real `z`. Returns `None` when the series
/// is outside its range or loses too much precision to cancellation.
/// Full-precision anchors retain terms to 1e-32 and require a 1e-26 relative
/// roundoff estimate; ordinary scalar evaluation keeps its existing thresholds.
pub(super) fn series_real<const FULL_PRECISION: bool>(v: f64, z: f64) -> Option<Dd> {
    let norm = z.abs();
    if norm > 100.0 && !(v > 0.0 && norm / v < 32.0) && !(v < 0.0 && norm <= 10000.0 && v > -2500.0)
    {
        return None;
    }
    let mut real = Dd::from(1.0);
    // Positive and negative terms are accumulated separately, so every
    // accumulation is the cheaper same-sign double-double addition; one full
    // addition combines them at the end.
    let (mut positive, mut negative) = (real, Dd::from(0.0));
    let mut peak: f64 = 1.0;
    for k in 1..=4096 {
        // term_k = term_(k-1) z / ((v + k - 1) k); v + k - 1 is exact in
        // double-double and the integer factor k needs no promotion.
        let kf = k as f64;
        let denominator = Dd::sum(v, kf - 1.0).mul_f64(kf);
        real = real.mul_f64(z).div_dd(denominator);
        if real.hi >= 0.0 {
            positive = positive.add_same_sign(real);
        } else {
            negative = negative.add_same_sign(real);
        }
        let term = real.hi.abs();
        let magnitude = (positive.hi + negative.hi).abs();
        peak = peak.max(term);
        if !magnitude.is_finite() {
            return None;
        }
        // Once the next ratio |z| / ((v + k) (k + 1)) is below one half, the
        // remaining tail is at most twice the current term, so terms below
        // 1e-20 of the sum are beyond binary64 resolution.
        let tail_bounded = (v + kf) * (kf + 1.0) > 2.0 * norm;
        let tolerance = if FULL_PRECISION { 1e-32 } else { 1e-20 };
        if (tail_bounded && term <= tolerance * magnitude)
            || term <= 1e-32 * magnitude
            || term == 0.0
        {
            let sum = positive + negative;
            let accuracy = if FULL_PRECISION { 1e-26 } else { 1e-13 };
            if peak * 1e-30 > accuracy * sum.hi.abs() {
                return None;
            }
            return Some(sum);
        }
    }
    None
}

// Certify compensated binary64 summation using an absolute-term forward-error
// bound. The same bound works for alternating terms, provided the final result
// is not too small relative to the term magnitudes; cancellation-sensitive
// cases fall through to the existing double-double/Bessel kernels.
fn compensated_series(v: f64, z: f64) -> Option<f64> {
    if v <= 0.0 || !(-10.0..=100.0).contains(&z) {
        return None;
    }
    let (mut term, mut sum, mut correction, mut weighted) = (1.0, 1.0, 0.0, 2.0);
    for k in 1..=256 {
        let kf = k as f64;
        term *= z / ((v + (kf - 1.0)) * kf);
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
        // Four rounded operations per term, and two eps for compensated sum.
        weighted += (4.0 * kf + 2.0) * term.abs();
        if !sum.is_finite() {
            return None;
        }
        if (v + kf) * (kf + 1.0) > 2.0 * z.abs() && term.abs() <= 0.125 * f64::EPSILON * sum.abs() {
            // gamma_(4k) <= (4k)*eps/(1-1024*eps), since k <= 256.
            let error = f64::EPSILON * weighted / (1.0 - 1024.0 * f64::EPSILON) + 2.0 * term.abs();
            return (error <= 16.0 * f64::EPSILON * sum.abs()).then_some(sum);
        }
    }
    None
}

// The defining 0F1 series has only positive terms when b,z>0. Summing in
// logarithmic form avoids both premature term overflow and underflow of the
// normalized large-order modified Bessel I. This path is used only when the
// faster compensated/double-double sums cannot certify their result.
fn positive_log_series(b: f64, z: f64) -> Option<f64> {
    if b <= 0.0 || z <= 0.0 {
        return None;
    }
    let mut log_term: f64 = 0.0;
    let mut log_sum: f64 = 0.0;
    let log_z = z.ln();
    for k in 1..=16384 {
        let kf = k as f64;
        let next_b = b + (kf - 1.0);
        if !next_b.is_finite() {
            return Some(log_sum.exp());
        }
        log_term += log_z - next_b.ln() - kf.ln();
        if log_term == f64::NEG_INFINITY {
            return Some(log_sum.exp());
        }
        let largest = log_sum.max(log_term);
        log_sum = largest + ((log_sum - largest).exp() + (log_term - largest).exp()).ln();
        if log_sum > f64::MAX.ln() {
            return Some(f64::INFINITY);
        }
        if z / ((b + kf) * (kf + 1.0)) < 0.5 && log_term < log_sum - 38.0 {
            return Some(log_sum.exp());
        }
    }
    None
}

// Uniform Debye expansion of normalized J_(b-1) for large positive b
// and 2*sqrt(-z) < (b-1)/sqrt(2). The regular Bessel J underflows in
// this regime even when 0F1 is ordinary-sized (e.g. b=100000,z=-1e7).
// Stirling cancellation is performed analytically, using log1p for the
// small exponent difference. DLMF 10.19.3, 10.41.7 and 5.11.1.
fn large_order_negative(b: f64, z: f64) -> Option<f64> {
    if b < 500.0 || z >= 0.0 {
        return None;
    }
    let nu = b - 1.0;
    let q = (2.0 * (-z).sqrt() / nu).powi(2);
    let limit = if b >= 1000.0 { 0.5 } else { 0.25 };
    if q >= limit {
        return None;
    }
    let s = (1.0 - q).sqrt();
    let d = -q / (1.0 + s);
    let p = 1.0 / s;
    let p2 = p * p;
    let p3 = p2 * p;
    let p4 = p2 * p2;
    let p5 = p4 * p;
    let p6 = p3 * p3;
    let p7 = p6 * p;
    let p9 = p7 * p2;
    let u1 = (3.0 * p - 5.0 * p3) / 24.0;
    let u2 = (81.0 * p2 - 462.0 * p4 + 385.0 * p6) / 1152.0;
    let u3 = (30375.0 * p3 - 369603.0 * p5 + 765765.0 * p7 - 425425.0 * p9) / 414720.0;
    let correction = 1.0 / (12.0 * nu) - 1.0 / (360.0 * nu.powi(3));
    let exponent = nu * (d - (d / 2.0).ln_1p()) - 0.5 * s.ln()
        + correction
        + (u1 / nu + u2 / nu.powi(2) + u3 / nu.powi(3)).ln_1p();
    Some(exponent.exp())
}

// DLMF 10.39.9 / 10.27.6:
//   0F1(;b; z) = Gamma(b) * |z|^((1-b)/2) * Bessel_(b-1)(2 sqrt|z|),
// with Bessel=J for z<0, I for z>0. Evaluate the exponential factors in
// logarithmic form and use exp(-t) I_v(t) for positive z.
fn bessel_continuation(b: f64, z: f64) -> Option<f64> {
    let order = b - 1.0;
    if order == b {
        return None;
    }
    let root = z.abs().sqrt();
    let t = 2.0 * root;
    let basis = if z < 0.0 {
        bessel::jv(order, t)
    } else {
        bessel::ive(order, t)
    };
    if !basis.is_finite() || basis == 0.0 {
        return None;
    }
    let (gamma_log, gamma_sign) = gamma::lgamma_sign(b);
    let exponent =
        gamma_log + (1.0 - b) * root.ln() + basis.abs().ln() + if z > 0.0 { t } else { 0.0 };
    if exponent.is_nan() {
        return None;
    }
    Some(gamma_sign * basis.signum() * exponent.exp())
}

fn real_kernel(v: f64, z: f64) -> Result<f64, StrErr> {
    validate(v)?;
    if !z.is_finite() {
        return Err(ERROR);
    }
    if z == 0.0 {
        return Ok(1.0);
    }
    // Half-integer parameters reduce exactly to elementary trigonometric
    // functions. These identities remain stable at zeros where the general
    // alternating series cannot certify any relative accuracy.
    if v == 0.5 {
        let root = z.abs().sqrt();
        let t = 2.0 * root;
        return Ok(if z > 0.0 {
            t.cosh()
        } else {
            // Near a cosine zero, rounding sqrt(-z) before evaluating cos
            // loses the entire result. Recover the low part of the root via
            // (-z) - root², computed with fused multiply-add, and correct
            // cos to first order. The omitted term is O(ulp(root)²).
            let root_correction = (-root).mul_add(root, -z) / root;
            // At enormous |z| the spacing between successive f64 roots is
            // larger than a radian, so a linear phase correction is invalid.
            // Return the bounded trigonometric value in that regime.
            if root_correction.abs() <= 0.125 {
                t.cos() - t.sin() * root_correction
            } else {
                t.cos()
            }
        });
    }
    if v == 1.5 {
        let t = 2.0 * z.abs().sqrt();
        return Ok(if z > 0.0 {
            if t > 700.0 {
                // sinh(t) overflows before sinh(t)/t does.
                (t - (2.0 * t).ln()).exp()
            } else {
                t.sinh() / t
            }
        } else {
            t.sin() / t
        });
    }
    if v == -0.5 {
        let t = 2.0 * z.abs().sqrt();
        return Ok(if z > 0.0 {
            t.cosh() - t * t.sinh()
        } else {
            t.cos() + t * t.sin()
        });
    }
    // Large negative arguments with |z| > v have strongly alternating
    // terms. The strict binary64 error bound usually rejects this region;
    // go directly to the higher-precision series to avoid a redundant pass.
    let cancellation_prone = v > 0.0 && z < -3.0 && v < -z;
    if !cancellation_prone && let Some(value) = compensated_series(v, z) {
        return Ok(value);
    }
    if let Some(value) = series_real::<false>(v, z) {
        return Ok(value.hi);
    }
    if let Some(value) = positive_log_series(v, z) {
        return Ok(value);
    }
    if let Some(value) = large_order_negative(v, z) {
        return Ok(value);
    }
    if let Some(value) = bessel_continuation(v, z) {
        return Ok(value);
    }
    Err("0F1 cannot be resolved accurately by the bounded real series or Bessel continuation")
}

/// Real 0F1(;v;z), NaN for poles and unresolved cancellation.
pub(crate) fn hyp0f1(v: f64, z: f64) -> f64 {
    real_kernel(v, z).unwrap_or(f64::NAN)
}
