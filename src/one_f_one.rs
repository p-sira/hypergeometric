//! Sci legacy confluent hypergeometric dispatcher and Kummer continuation.
//! Adapted from sci-special's pure Rust implementation. See crate LICENSE.
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
    // Nonpositive integer a terminates the defining hypergeometric series.
    // Evaluating that polynomial directly avoids the gamma poles introduced
    // by the general asymptotic formulas for negative x.
    if a < 0.0 && a == a.floor() {
        let (mut term, mut sum) = (1.0, 1.0);
        for k in 0..(-a) as usize {
            let k = k as f64;
            term *= (a + k) * x / ((b + k) * (k + 1.0));
            sum += term;
        }
        return sum;
    }
    // The Kummer-transformed series still has catastrophic cancellation when
    // b-a is large and negative. Continue the analytic solution outwards instead
    // of summing exponentially larger alternating terms at the target.
    if a > 0.0 && b < a && (-1000.0..0.0).contains(&x) && (-8.0..=64.0).contains(&b) {
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
        let (mut term, mut sum, mut correction) = (1.0, 1.0, 0.0);
        for k in 0..1000 {
            let k = k as f64;
            term *= (alpha + k) * argument / ((b + k) * (k + 1.0));
            let adjusted = term - correction;
            let next = sum + adjusted;
            correction = (next - sum) - adjusted;
            sum = next;
            if term.abs() <= f64::EPSILON * sum.abs() {
                return scale * sum;
            }
        }
    }
    crate::specfun::chgm(x, a, b)
}

#[cfg(test)]
mod tests {
    use super::hyp1f1;

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
