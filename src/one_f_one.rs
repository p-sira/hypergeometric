//! Sci legacy confluent hypergeometric dispatcher and Kummer continuation.
//! Adapted from sci-special's pure Rust implementation. See crate LICENSE.
use crate::dd::Dd;

pub(crate) fn hyp1f1(a: f64, b: f64, x: f64) -> f64 {
    if !a.is_finite() || !b.is_finite() || !x.is_finite() || a.abs() > 100_000.0 {
        return f64::NAN;
    }
    // A nonpositive integer denominator is a pole unless a negative integer
    // numerator terminates the series before the vanishing denominator term.
    // This includes a=b=-n: M(-n,-n,x) is a polynomial, not exp(x).
    if b <= 0.0 && b == b.floor() && !(a < 0.0 && a == a.floor() && a >= b) {
        return f64::INFINITY;
    }
    if a == 0.0 || x == 0.0 {
        return 1.0;
    }
    if a == b && (a > 0.0 || a != a.floor()) {
        return x.exp();
    }
    if b > 0.0 && a == b + 1.0 && a - 1.0 == b {
        // M(b+1,b,x) = exp(x) (1 + x/b). This avoids a power series or
        // outward continuation. Check both directions so an inexact b+1
        // (e.g. 0.1 + 1.0) cannot spuriously satisfy the identity for the
        // exact binary64 input parameters. Near the zero x=-b, add before
        // dividing to avoid cancellation in 1 + x/b.
        let factor = (b + x) / b;
        if factor == 0.0 {
            return 0.0;
        }
        if factor.is_finite() && (-700.0..=700.0).contains(&x) {
            return x.exp() * factor;
        }
        // Keep tiny results representable when exp(x) underflows, and
        // avoid overflow in x/b when b is subnormal.
        let (sign, log_factor) = if factor.is_finite() {
            (factor.signum(), factor.abs().ln())
        } else {
            (x.signum(), x.abs().ln() - b.ln())
        };
        return sign * (x + log_factor).exp();
    }
    // M(-n,b,x) is a normalized Laguerre polynomial. Its defining sum
    // can cancel tens of digits; the three-term recurrence avoids those
    // enormous alternating terms. DD also preserves accuracy near zeros.
    if a < 0.0 && a == a.floor() {
        return laguerre((-a) as usize, b, x);
    }
    // Forward recurrence in -a avoids the exponentially large alternating
    // terms of the defining series. Limit the seed and recurrence work;
    // uncertain results still use the existing series fallback.
    if (-128.0..-16.0).contains(&a)
        && (2.0..=64.0).contains(&b)
        && (8.0..=32.0).contains(&x)
        && let Some(value) = if [a, b, x].iter().all(|v| v * 16.0 == (v * 16.0).trunc()) {
            negative_numerator::<true>(a, b, x)
        } else {
            negative_numerator::<false>(a, b, x)
        }
    {
        return value;
    }
    // The Kummer-transformed series still has catastrophic cancellation when
    // b-a is large and negative. Continue the analytic solution outwards instead
    // of summing exponentially larger alternating terms at the target.
    if a > 0.0 && b < a && (-1000.0..0.0).contains(&x) && (-64.0..=64.0).contains(&b) {
        if b > 0.0 && a <= 64.0 && x >= -12.0 {
            // The transformed series can still be cheap at small |x|.
            // Certify its cancellation before accepting it: near zeros, or
            // when large terms cancel, keep the outward continuation.
            let (mut term, mut sum, mut correction, mut absolute_sum) = (1.0, 1.0, 0.0, 1.0);
            for k in 1..=256 {
                let kf = k as f64;
                term *= (b - a + kf - 1.0) * -x / ((b + kf - 1.0) * kf);
                let adjusted = term - correction;
                let next = sum + adjusted;
                correction = (next - sum) - adjusted;
                sum = next;
                absolute_sum += term.abs();
                if !sum.is_finite() || !absolute_sum.is_finite() {
                    break;
                }
                let next_ratio = (b - a + kf).abs() * -x / ((b + kf) * (kf + 1.0));
                // The consecutive-term ratio has passed its maximum by
                // k > 2(a-b)+b+1, so the remaining tail is geometric.
                if kf > 2.0 * (a - b) + b + 1.0
                    && next_ratio <= 0.5
                    && term.abs() <= 0.125 * f64::EPSILON * sum.abs()
                {
                    let error_bound = (8.0 * kf + 2.0) * f64::EPSILON * absolute_sum;
                    if error_bound <= 1e-12 * sum.abs() {
                        return x.exp() * sum;
                    }
                    break;
                }
            }
        }
        return crate::kummer::negative_argument(a, b, x);
    }
    // In the moderate-argument region the defining series is more accurate
    // than the truncated large-x expansion. Kummer's transformation makes
    // the argument positive and avoids alternating exponential cancellation.
    if x.abs() <= 250.0 && a.abs().max(b.abs()) <= 128.0 {
        let (alpha, argument, scale) = if x < 0.0 {
            (b - a, -x, x.exp())
        } else {
            (a, x, 1.0)
        };
        // When both rising-factorial parameters are positive, every term in
        // the transformed series is positive. Ordinary accumulation avoids
        // the compensated-add dependency chain without cancellation loss.
        if alpha > 0.0 && b > 0.0 {
            let (mut term, mut sum) = (1.0, 1.0);
            for k in 0..1000 {
                let k = k as f64;
                term *= (alpha + k) * argument / ((b + k) * (k + 1.0));
                sum += term;
                if !sum.is_finite() {
                    break;
                }
                if term <= f64::EPSILON * sum {
                    if alpha <= b {
                        // Retain the proven inexpensive stopping rule for
                        // the original parameter regime.
                        if k > argument + alpha {
                            return scale * sum;
                        }
                    } else {
                        // For alpha > b the ratios decrease monotonically.
                        // Bound the geometric tail to stop safely before
                        // the old conservative k > argument + alpha limit.
                        let ratio = (alpha + k + 1.0) * argument / ((b + k + 1.0) * (k + 2.0));
                        if ratio < 1.0 && term * ratio <= f64::EPSILON * sum * (1.0 - ratio) {
                            return scale * sum;
                        }
                    }
                }
            }
        }
        let alpha = if x < 0.0 {
            Dd::from(b) - Dd::from(a)
        } else {
            Dd::from(a)
        };
        if let Some(value) =
            signed_series_f64(alpha.hi, b, argument).or_else(|| signed_series(alpha, b, argument))
        {
            return scale * value;
        }
        // An uncertified alternating series must not fall back to the same
        // binary64 sum in the legacy kernel.
        return f64::NAN;
    }
    crate::specfun::chgm(x, a, b)
}

// DLMF 13.3.1. Seed M(s,b,x), M(s+1,b,x), 0<s<1, from one
// positive DD series: M(s+1)=sum + sum(k*term_k)/s. Keeping the low
// parts of the seeds is essential; binary64 seeds lose many digits here.
fn negative_numerator<const DYADIC: bool>(a: f64, b: f64, x: f64) -> Option<f64> {
    numerator_recurrence::<DYADIC, true>(a, b, x)
        .or_else(|| numerator_recurrence::<DYADIC, false>(a, b, x))
}

fn numerator_recurrence<const DYADIC: bool, const HYBRID: bool>(
    a: f64,
    b: f64,
    x: f64,
) -> Option<f64> {
    let steps = (-a).ceil() as usize;
    let s = a + steps as f64;
    let (mut term, mut sum, mut weighted) = (Dd::from(1.0), Dd::from(1.0), Dd::from(0.0));
    let seed_tolerance = if HYBRID { 1e-26 } else { 1e-32 };
    let mut converged = false;
    for k in 0..256 {
        let kf = k as f64;
        term = if DYADIC {
            // Multiples of 1/16 have exact small integer coefficients.
            // Combine the numerator factors when x is an integer as well.
            let denominator = (16.0 * b + 16.0 * kf) * (kf + 1.0);
            if x == x.trunc() {
                term.mul_div_small((16.0 * s + 16.0 * kf) * x, denominator)
            } else {
                term.mul_f64(x)
                    .mul_div_small(16.0 * s + 16.0 * kf, denominator)
            }
        } else {
            (term * Dd::sum(s, kf))
                .mul_f64(x)
                .div_dd(Dd::sum(b, kf).mul_small(kf + 1.0))
        };
        sum = sum.add_same_sign(term);
        let weight = term.mul_small(kf + 1.0);
        weighted = weighted.add_same_sign(weight);
        if kf > 2.0 * x && weight.hi <= seed_tolerance * weighted.hi {
            converged = true;
            break;
        }
    }
    if !converged {
        return None;
    }
    let mut previous = sum + weighted.div_dd(Dd::from(s));
    let mut current = sum;
    // For k>2x the seed tail ratios are below 1/2. Ten times the
    // stopping tolerance covers both tails; 1e-28 covers DD roundoff.
    let seed_error = 10.0 * seed_tolerance + 1e-28;
    let seed_previous_error = seed_error * previous.hi.abs();
    let seed_current_error = seed_error * current.hi.abs();
    let mut coefficients = [(0.0, 0.0, 0.0); 128];
    for (j, coefficients) in coefficients.iter_mut().enumerate().take(steps) {
        // After the initial large values have decayed, try inexpensive
        // binary64 steps. The adjoint check includes their larger rounding
        // errors; if it fails, repeat the recurrence entirely in DD.
        let ordinary = HYBRID && j as f64 > x * 0.5 + 4.0;
        let (next, ac, bc) = if ordinary {
            let aj = s - j as f64;
            let denominator = b - aj;
            let left = b - 2.0 * aj - x;
            (
                Dd::from((left * current.hi + aj * previous.hi) / denominator),
                left / denominator,
                aj / denominator,
            )
        } else if DYADIC {
            let aj = s - j as f64;
            let denominator = b - aj;
            let left = b - 2.0 * aj - x;
            (
                (current.mul_small(16.0 * left) + previous.mul_small(16.0 * aj))
                    .div_small(16.0 * denominator),
                left / denominator,
                aj / denominator,
            )
        } else {
            let aj = Dd::sum(s, -(j as f64));
            let denominator = Dd::from(b) - aj;
            let left = Dd::from(b) - aj.mul_small(2.0) - Dd::from(x);
            (
                (left * current + aj * previous).div_dd(denominator),
                left.hi / denominator.hi,
                aj.hi / denominator.hi,
            )
        };
        *coefficients = (
            ac,
            bc,
            (if ordinary { 16.0 * f64::EPSILON } else { 1e-29 })
                * ((ac * current.hi).abs() + (bc * previous.hi).abs()),
        );
        previous = current;
        current = next;
    }
    // Propagate sensitivities backwards, retaining their signs. Forward
    // absolute bounds would spuriously grow exponentially in the oscillatory
    // region. The adjoint weights bound each local DD rounding contribution.
    let (mut wc, mut wp, mut error) = (1.0_f64, 0.0_f64, 0.0_f64);
    for &(ac, bc, local_error) in coefficients[..steps].iter().rev() {
        error += wc.abs() * local_error;
        (wc, wp) = (wc * ac + wp, wc * bc);
    }
    error += wc.abs() * seed_current_error + wp.abs() * seed_previous_error;
    let value = current.hi + current.lo;
    (value.is_finite() && error <= 1e-12 * value.abs()).then_some(value)
}

// Normalized Laguerre recurrence, with common binary rescaling so a large
// intermediate cannot overflow before division by b+n.
fn laguerre(degree: usize, b: f64, x: f64) -> f64 {
    let ratio = x / b;
    if !ratio.is_finite() {
        return if degree == 1 { 1.0 - ratio } else { f64::NAN };
    }
    let (mut previous, mut current) = (Dd::from(1.0), Dd::from(1.0) - Dd::from(x) / Dd::from(b));
    let mut exponent = 0;
    for n in 1..degree {
        let nf = n as f64;
        let denominator = Dd::sum(b, nf);
        let factor = Dd::from(1.0) + Dd::sum(nf, -x).div_dd(denominator);
        let previous_factor = Dd::from(nf).div_dd(denominator);
        let coefficient = factor.hi.abs().max(previous_factor.hi.abs()).max(1.0);
        if !coefficient.is_finite() {
            return f64::NAN;
        }
        while current.hi.abs().max(previous.hi.abs()) > 1e280 / coefficient {
            let scale = 2.0_f64.powi(-512);
            previous = previous.mul_f64(scale);
            current = current.mul_f64(scale);
            exponent += 512;
        }
        let next = factor * current - previous_factor * previous;
        previous = current;
        current = next;
    }
    let value = current.hi + current.lo;
    if exponent == 0 {
        value
    } else {
        value.signum() * (value.abs().ln() + f64::from(exponent) * std::f64::consts::LN_2).exp()
    }
}

// Binary64 pre-pass for `signed_series` with the same stopping rule. Each
// term carries at most about 8(k+1) roundings (including rounding of `a`)
// and Neumaier summation keeps accumulation error O(eps|sum|), so accept
// only when 8 eps sum (k+1)|t_k| meets the DD kernel's 1e-13 target.
fn signed_series_f64(a: f64, b: f64, x: f64) -> Option<f64> {
    let (mut term, mut sum, mut correction) = (1.0_f64, 1.0_f64, 0.0_f64);
    let mut weighted = 1.0_f64;
    for k in 0..1000 {
        let kf = k as f64;
        term *= (a + kf) * x / ((b + kf) * (kf + 1.0));
        let next = sum + term;
        correction += if sum.abs() >= term.abs() {
            (sum - next) + term
        } else {
            (term - next) + sum
        };
        sum = next;
        weighted += (kf + 2.0) * term.abs();
        if !sum.is_finite() || !weighted.is_finite() {
            return None;
        }
        if term.abs() > 0.0625 * f64::EPSILON * sum.abs() {
            continue;
        }
        let next_a = a + kf + 1.0;
        let next_b = b + kf + 1.0;
        let ratio_bound = x / (kf + 2.0) * (next_a / next_b).max(1.0);
        if next_a > 0.0
            && next_b > 0.0
            && ratio_bound < 0.5
            && term.abs() <= 0.0625 * f64::EPSILON * sum.abs()
        {
            let sum = sum + correction;
            return (8.0 * f64::EPSILON * weighted <= 1e-13 * sum.abs()).then_some(sum);
        }
    }
    None
}

// Double-double defining series for mixed-sign rising factorials. Wait until
// both parameters are positive and the remaining ratios are bounded by 1/2;
// a tiny term before a nearby numerator zero does not bound the later tail.
fn signed_series(a: Dd, b: f64, x: f64) -> Option<f64> {
    let (mut term, mut sum) = (Dd::from(1.0), Dd::from(1.0));
    let mut absolute_sum = 1.0;
    for k in 0..1000 {
        let kf = k as f64;
        term = (term * a.add_f64(kf))
            .mul_f64(x)
            .div_dd(Dd::sum(b, kf).mul_small(kf + 1.0));
        sum = sum + term;
        absolute_sum += term.hi.abs();
        if !sum.hi.is_finite() || !absolute_sum.is_finite() {
            return None;
        }
        if term.hi.abs() > 1e-25 * sum.hi.abs() {
            continue;
        }
        let next_a = a.add_f64(kf + 1.0).hi;
        let next_b = b + kf + 1.0;
        let ratio_bound = x / (kf + 2.0) * (next_a / next_b).max(1.0);
        if next_a > 0.0
            && next_b > 0.0
            && ratio_bound < 0.5
            && term.hi.abs() <= 1e-25 * sum.hi.abs()
        {
            return (absolute_sum * 1e-29 <= 1e-13 * sum.hi.abs()).then_some(sum.hi + sum.lo);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::hyp1f1;

    #[test]
    fn negative_numerator_cancellation() {
        let value =
            super::negative_numerator::<true>(-50.5, 2.0, 30.0).expect("certified recurrence");
        assert!((value / 3610.292112991606 - 1.0).abs() < 1e-13, "{value}");
    }

    #[test]
    fn recurrence_rejects_uncertainty_near_a_zero() {
        // 100-digit mpmath root, rounded to binary64. The function value is
        // about 4.10e-12, too small to certify relatively from these seeds.
        assert!(super::negative_numerator::<false>(-50.5, 2.0, 29.60332694942178).is_none());
    }

    #[test]
    fn positive_series_stops_after_the_peak() {
        // Kummer's identity turns M(64,2,z) into the finite polynomial
        // exp(z) M(-62,2,-z); all polynomial terms are positive for z > 0.
        for z in [0.75, 2.0, 10.0] {
            let (mut term, mut polynomial) = (1.0, 1.0);
            for k in 0..62 {
                let k = f64::from(k);
                term *= (k - 62.0) * -z / ((k + 2.0) * (k + 1.0));
                polynomial += term;
            }
            let expected = z.exp() * polynomial;
            let actual = hyp1f1(64.0, 2.0, z);
            assert!(
                (actual / expected - 1.0).abs() <= 3e-14,
                "z={z}: {actual} vs {expected}"
            );
        }
    }

    #[test]
    fn positive_series_small_numerator_relative_to_denominator() {
        // M(1,2,z) = expm1(z)/z, including arguments where many positive
        // terms must be added before the series enters its decreasing tail.
        for z in [2.0_f64, 100.0, 200.0] {
            let expected = z.exp_m1() / z;
            let actual = hyp1f1(1.0, 2.0, z);
            assert!(
                (actual / expected - 1.0).abs() <= 3e-14,
                "z={z}: {actual} vs {expected}"
            );
        }
    }
}
