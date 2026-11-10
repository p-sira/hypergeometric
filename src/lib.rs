//! Real hypergeometric and Wright functions in pure Rust.
//!
//! Functions return NaN for invalid arguments or numerical regions where
//! a bounded calculation cannot reliably resolve the real value. See individual
//! APIs and README for documented work limits. `hyp0f1` and `wrightomega`
//! accept `f32` and `f64` through `num_traits::Float`, evaluate in `f64`, and
//! convert the result back to the input type. Other functions use `f64`.
//!
//! # Example
//! ```
//! use hypergeometric::{hyp0f1, hyp1f1, hyp2f1, hyperu, wright_bessel, wrightomega};
//! assert!((hyp0f1(2.5_f64, 0.0)-1.0).abs() < 1e-14);
//! assert!((hyp1f1(1.0, 1.0, 0.5)-0.5_f64.exp()).abs() < 1e-14);
//! assert!((hyp2f1(1.0, 1.0, 2.0, 0.5)-2.0*2.0_f64.ln()).abs() < 1e-12);
//! assert!((hyperu(0.0, 2.5, 1.0)-1.0).abs() < 1e-14);
//! assert!(wright_bessel(0.5, 1.5, 1.0).is_finite());
//! let w: f64 = wrightomega(1.0);
//! assert!((w+w.ln()-1.0).abs()<1e-13);
//! ```
#![allow(clippy::excessive_precision)]
mod dd;
mod gamma;
mod kummer;
mod one_f_one;
mod specfun;
mod zero_f_one;

/// Real confluent limit function 0F1(;b;z) = Σ_{k≥0} z^k / ((b)_k k!).
///
/// # Parameters and domain
/// `b` is the dimensionless lower parameter and `z` the real argument.
/// Both must be finite. For `b = 0,-1,-2,…`, the function has a pole and
/// this API returns NaN. In particular `hyp0f1(b,0) = 1` for valid `b`.
///
/// # Numerical behavior
/// Compensated and double-double power series cover the moderate range and
/// near-parameter-pole regions. Larger arguments use the real Bessel
/// identities `Γ(b)|z|^((1-b)/2) J_(b-1)(2√|z|)` for negative `z`, and
/// the corresponding modified Bessel `I` for positive `z` with exponential
/// scaling. A positive-term logarithmic series prevents overflow/underflow of
/// intermediate values for `b,z>0`. Large positive `b` with negative `z`
/// also uses a uniform large-order expansion when `J` underflows.
///
/// Can return signed infinities when the mathematical value exceeds binary64
/// range; a bounded method can still return NaN in unresolved extreme or
/// highly oscillatory regions, particularly close to Bessel zeros.
/// Near zeros, relative error can grow even when absolute error is small.
/// Work is bounded to 16,384 logarithmic-series terms and 4,096 double-double
/// terms. No complex analytic continuation is provided. Accepts
/// `num_traits::Float` (including `f32` and `f64`), computes in `f64`, and
/// converts the result to `T`; this does not provide additional precision.
///
/// # Example
/// ```
/// use hypergeometric::hyp0f1;
/// assert_eq!(hyp0f1(0.5, 0.0), 1.0);
/// assert!((hyp0f1(0.5, -1.0) - 2.0_f64.cos()).abs() < 1e-14);
/// assert!(hyp0f1(0.0_f64, -1.0).is_nan());
/// ```
///
/// References: DLMF [16.2](https://dlmf.nist.gov/16.2),
/// [10.39](https://dlmf.nist.gov/10.39),
/// [10.19](https://dlmf.nist.gov/10.19).
pub fn hyp0f1<T: num_traits::Float>(b: T, z: T) -> T {
    let b = b.to_f64().unwrap_or(f64::NAN);
    let z = z.to_f64().unwrap_or(f64::NAN);
    T::from(zero_f_one::hyp0f1(b, z)).unwrap_or(T::nan())
}

/// Kummer confluent function M(a,b,z) = Σ_{k≥0} (a)_k z^k/((b)_k k!).
///
/// Finite real parameters |a|≤100000 and finite z, with poles at b=0,-1,…
/// except where a terminating series removes the pole. Negative arguments
/// use Kummer transformation and outward Taylor continuation.
/// Reference: DLMF [13.2](https://dlmf.nist.gov/13.2).
pub fn hyp1f1(a: f64, b: f64, z: f64) -> f64 {
    one_f_one::hyp1f1(a, b, z)
}

/// Gauss function 2F1(a,b;c;z)=Σ_{k≥0}(a)_k(b)_k z^k/((c)_k k!).
///
/// Real continuation follows the Zhang–Jin/Cephes specfun algorithm;
/// parameter magnitudes above 100000 and non-finite inputs return NaN.
/// Real continuation beyond the branch cut z>1 is only supported where the
/// series terminates. At other singularities numerical overflow may result.
/// Reference: DLMF [15.2](https://dlmf.nist.gov/15.2).
pub fn hyp2f1(a: f64, b: f64, c: f64, z: f64) -> f64 {
    if [a, b, c, z].iter().any(|v| !v.is_finite()) || a.abs().max(b.abs()).max(c.abs()) > 100000.0 {
        return f64::NAN;
    }
    specfun::hyp2f1(a, b, c, z)
}

/// Tricomi's function U(a,b,z), the recessive Kummer solution.
///
/// For a>0 and z>0, U(a,b,z)=Γ(a)⁻¹∫₀^∞e^(-zt)t^(a-1)(1+t)^(b-a-1)dt.
/// Finite real |a|,|b|≤100000, z≥0. Uses Zhang–Jin series, recurrence,
/// asymptotic, and quadrature dispatch with bounded work; some regions have
/// only about nine valid digits. Real z<0 is outside the supported domain.
/// Reference: DLMF [13.4](https://dlmf.nist.gov/13.4).
pub fn hyperu(a: f64, b: f64, z: f64) -> f64 {
    if [a, b, z].iter().any(|v| !v.is_finite()) || a.abs().max(b.abs()) > 100000.0 || z < 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return 1.0;
    }
    if z == 0.0 {
        // U(-n,b,z)=(-1)^n n! L_n^(b-1)(z), so its endpoint value
        // is the terminating polynomial (-1)^n (b)_n, including b>=1.
        if a < 0.0 && a == a.floor() && -a <= 4096.0 {
            let mut value = 1.0;
            for k in 0..(-a) as usize {
                value *= -(b + k as f64);
            }
            return value;
        }
        if b >= 1.0 {
            return f64::INFINITY;
        }
        // Gamma(1-b)/Gamma(1+a-b); reciprocal gamma preserves removable
        // denominator poles as exactly zero for a negative integer.
        return gamma::gamma(1.0 - b) * gamma::rgamma(1.0 + a - b);
    }
    if b == 0.0
        && a > 0.0
        && a <= 4.0
        && z <= 20.0
        && let Some(value) = hyperu_zero_b_integral(a, z)
    {
        return value;
    }
    if b == 0.0 {
        let transformed = z * hyperu(a + 1.0, 2.0, z);
        if transformed.is_finite() {
            return transformed;
        }
    }
    if (b == 1.0 && z < 1.0 && a > -0.25 && a < 0.3) || (b == 2.0 && z < 1.0 && a > 1.0 && a < 2.0)
    {
        return (z + 2.0 + 2.0 * a - b) * hyperu(a + 1.0, b, z)
            - (a + 1.0) * (a + 2.0 - b) * hyperu(a + 2.0, b, z);
    }
    let mut md = [0];
    let mut error = [0];
    specfun::chgu(z, a, b, &mut md, &mut error)
}

// DLMF 13.4.4, for b=0 and a>0, substitute t=u/(1-u), then v=u^a:
// U(a,0,z)= 1/Gamma(a+1) integral_0^1 exp(-z*v^(1/a)/(1-v^(1/a)))dv.
// The v substitution removes the integrable endpoint singularity for small a.
// Adaptive Simpson has fixed depth and a global evaluation budget.
fn hyperu_zero_b_integral(a: f64, z: f64) -> Option<f64> {
    fn f(v: f64, a: f64, z: f64) -> f64 {
        if v >= 1.0 {
            return 0.0;
        }
        if a <= 1.0 {
            let u = v.powf(1.0 / a);
            (-z * u / (1.0 - u)).exp()
        } else {
            let p = (3.0 / a).ceil().max(1.0);
            let u = v.powf(p);
            p * v.powf(p * a - 1.0) * (-z * u / (1.0 - u)).exp()
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn recurse(
        a: f64,
        z: f64,
        l: f64,
        r: f64,
        fl: f64,
        fm: f64,
        fr: f64,
        whole: f64,
        tol: f64,
        depth: usize,
        budget: &mut usize,
    ) -> Option<f64> {
        if *budget < 2 {
            return None;
        }
        *budget -= 2;
        let m = (l + r) * 0.5;
        let lm = (l + m) * 0.5;
        let mr = (m + r) * 0.5;
        let f_lm = f(lm, a, z);
        let f_mr = f(mr, a, z);
        let left = (m - l) * (fl + 4.0 * f_lm + fm) / 6.0;
        let right = (r - m) * (fm + 4.0 * f_mr + fr) / 6.0;
        let delta = left + right - whole;
        if delta.abs() <= 15.0 * tol {
            return Some(left + right + delta / 15.0);
        }
        if depth == 0 {
            return None;
        }
        let left = recurse(a, z, l, m, fl, f_lm, fm, left, tol * 0.5, depth - 1, budget)?;
        let right = recurse(
            a,
            z,
            m,
            r,
            fm,
            f_mr,
            fr,
            right,
            tol * 0.5,
            depth - 1,
            budget,
        )?;
        Some(left + right)
    }
    let fl = f(0.0, a, z);
    let fr = f(1.0, a, z);
    let fm = f(0.5, a, z);
    let mut budget = 65536;
    let integral = recurse(
        a,
        z,
        0.0,
        1.0,
        fl,
        fm,
        fr,
        (fl + 4.0 * fm + fr) / 6.0,
        1e-13,
        26,
        &mut budget,
    )?;
    Some(
        integral
            / if a <= 1.0 {
                gamma::gamma(a + 1.0)
            } else {
                gamma::gamma(a)
            },
    )
}

/// Wright generalized Bessel function Φ(a,b;z)=Σ_{k≥0}z^k/(k!Γ(ak+b)).
///
/// The real positive-sum domain is a,b,z≥0 with finite arguments. At most
/// 100000 series terms are evaluated, using log-sum-exp for dynamic range.
/// Reference: [Wright function](https://dlmf.nist.gov/10.46).
pub fn wright_bessel(a: f64, b: f64, z: f64) -> f64 {
    log_wright_bessel(a, b, z).exp()
}

// 90-digit mpmath ln(k!) references; avoid two Lanczos evaluations per Wright term.
const WRIGHT_LOG_FACTORIAL: [f64; 65] = [
    0.0,
    0.0,
    std::f64::consts::LN_2,
    1.79175946922805496e+00,
    3.17805383034794575e+00,
    4.78749174278204581e+00,
    6.57925121201010121e+00,
    8.52516136106541467e+00,
    1.06046029027452509e+01,
    1.28018274800814691e+01,
    1.51044125730755159e+01,
    1.75023078458738865e+01,
    1.99872144956618847e+01,
    2.25521638531234245e+01,
    2.51912211827386798e+01,
    2.78992713838408903e+01,
    3.06718601060806719e+01,
    3.35050734501368908e+01,
    3.63954452080330526e+01,
    3.93398841871994946e+01,
    4.23356164607534851e+01,
    4.53801388984769076e+01,
    4.84711813518352272e+01,
    5.16066755677643769e+01,
    5.47847293981123187e+01,
    5.80036052229805179e+01,
    6.12617017610020014e+01,
    6.45575386270063376e+01,
    6.78897431371815401e+01,
    7.12570389671680147e+01,
    7.46582363488301581e+01,
    7.80922235533153071e+01,
    8.15579594561150429e+01,
    8.50544670175815156e+01,
    8.85808275421976816e+01,
    9.21361756036870929e+01,
    9.57196945421432019e+01,
    9.93306124547874276e+01,
    1.02968198614513810e+02,
    1.06631760260643460e+02,
    1.10320639714757391e+02,
    1.14034211781461707e+02,
    1.17771881399745070e+02,
    1.21533081515438639e+02,
    1.25317271149356898e+02,
    1.29123933639127216e+02,
    1.32952575035616320e+02,
    1.36802722637326355e+02,
    1.40673923648234251e+02,
    1.44565743946344895e+02,
    1.48477766951773020e+02,
    1.52409592584497346e+02,
    1.56360836303078798e+02,
    1.60331128216630901e+02,
    1.64320112263195170e+02,
    1.68327445448427653e+02,
    1.72352797139162789e+02,
    1.76395848406997345e+02,
    1.80456291417543781e+02,
    1.84533828861449479e+02,
    1.88628173423671598e+02,
    1.92739047287844897e+02,
    1.96866181672890008e+02,
    2.01009316399281516e+02,
    2.05168199482641199e+02,
];

/// Logarithm of the generalized Wright Bessel function on a,b,z≥0.
/// Returns -∞ when the positive-domain series evaluates to exactly zero.
pub fn log_wright_bessel(a: f64, b: f64, z: f64) -> f64 {
    if [a, b, z].iter().any(|v| !v.is_finite()) || a < 0.0 || b < 0.0 || z < 0.0 {
        return f64::NAN;
    }
    let lgam = |p: f64| gamma::lgamma_sign(p).0;
    if a == 0.0 {
        return z - lgam(b);
    }
    if z == 0.0 {
        return -lgam(b);
    }
    // For a=1, Gamma(b+k)=Gamma(b+k-1)*(b+k-1). Reusing this
    // recurrence removes a Lanczos evaluation for every series term.
    // Restrict the fast path to moderate arguments, where incrementally
    // updating the logarithms does not accumulate significant roundoff.
    if a == 1.0 && b <= 1000.0 && z <= 1000.0 {
        let lz = z.ln();
        // At b=0 the k=0 reciprocal-Gamma term vanishes; k=1 is z.
        let (mut exponent, start) = if b == 0.0 { (lz, 1) } else { (-lgam(b), 0) };
        let (mut scale, mut sum) = (exponent, 1.0_f64);
        for k in (start + 1)..100000 {
            let kf = k as f64;
            let previous = exponent;
            exponent += lz - kf.ln() - (b + kf - 1.0).ln();
            let term = if exponent > scale {
                sum *= (scale - exponent).exp();
                scale = exponent;
                1.0
            } else {
                (exponent - scale).exp()
            };
            sum += term;
            if exponent < previous && term <= f64::EPSILON * sum {
                return scale + sum.ln();
            }
        }
        return f64::NAN;
    }
    let lz = z.ln();
    let (mut scale, mut sum, mut prev) = (f64::NEG_INFINITY, 0.0_f64, f64::NEG_INFINITY);
    for k in 0..100000 {
        let log_factorial = WRIGHT_LOG_FACTORIAL
            .get(k)
            .copied()
            .unwrap_or_else(|| lgam(k as f64 + 1.0));
        let k = k as f64;
        let exponent = k * lz - log_factorial - lgam(a * k + b);
        if exponent == f64::NEG_INFINITY {
            // k=0 may vanish at b=0, while later terms can be finite.
            // For k>=1, -infinity means log-Gamma has overflowed at a
            // positive, increasing argument. All subsequent terms vanish,
            // so the accumulated sum is already converged.
            if k > 0.0 {
                return scale + sum.ln();
            }
            continue;
        }
        let term = if exponent > scale {
            sum *= (scale - exponent).exp();
            scale = exponent;
            1.0
        } else {
            (exponent - scale).exp()
        };
        sum += term;
        if exponent < prev && term <= f64::EPSILON * sum {
            return scale + sum.ln();
        }
        prev = exponent;
    }
    f64::NAN
}

/// Real Wright omega w solving w + ln(w) = x with w≥0.
///
/// Evaluates principal Lambert W of exp(x) without overflow for large x.
/// For negative infinity returns +0; for positive infinity returns +infinity.
/// Accepts `num_traits::Float` (including `f32` and `f64`), computes in
/// `f64`, and converts the result to `T`.
/// Reference: Lawrence et al., ACM TOMS Algorithm 917 (BSD-licensed).
pub fn wrightomega<T: num_traits::Float>(x: T) -> T {
    let x = x.to_f64().unwrap_or(f64::NAN);
    T::from(wrightomega_f64(x)).unwrap_or(T::nan())
}

fn wrightomega_f64(x: f64) -> f64 {
    if x == f64::NEG_INFINITY {
        return 0.0;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x.is_nan() {
        return x;
    }
    if x < -50.0 {
        return x.exp();
    }
    if x > 1e20 {
        return x;
    }
    let mut w = if x < -2.0 {
        x.exp()
    } else if x < 1.0 {
        (2.0 * (x - 1.0) / 3.0).exp()
    } else {
        let l = x.ln();
        x - l + l / x
    };
    for _ in 0..2 {
        let r = (x - w) - w.ln();
        let p = w + 1.0;
        let t = 2.0 * p * (p + 2.0 * r / 3.0);
        w *= 1.0 + (r / p) * (t - r) / (t - 2.0 * r);
    }
    w
}
