//! Confluent limit function, compensated and double-double power series.
use crate::dd::Dd;
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
    if norm > 100.0 && !(v > 0.0 && norm / v < 16.0) {
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
        return Ok(if z > 0.0 { t.sinh() / t } else { t.sin() / t });
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
    // The Bessel fallback in sci-special was deliberately not transferred:
    // this self-contained crate reports unsupported cancellation/argument
    // regions explicitly rather than importing a cyclic backend dependency.
    Err("0F1 cannot be resolved accurately by the bounded real series")
}

/// Real 0F1(;v;z), NaN for poles and unresolved cancellation.
pub(crate) fn hyp0f1(v: f64, z: f64) -> f64 {
    real_kernel(v, z).unwrap_or(f64::NAN)
}
