#![allow(dead_code)]
//! Shared double-double arithmetic for cancellation-sensitive series.
use core::ops::{Add, Div, Mul, Neg, Sub};
#[allow(unused_imports)]
#[derive(Clone, Copy)]
pub(crate) struct Dd {
    pub(crate) hi: f64,
    pub(crate) lo: f64,
}
impl From<f64> for Dd {
    fn from(hi: f64) -> Self {
        Self { hi, lo: 0.0 }
    }
}
#[inline(always)]
fn two_sum(a: f64, b: f64) -> Dd {
    let hi = a + b;
    let v = hi - a;
    Dd {
        hi,
        lo: (a - (hi - v)) + (b - v),
    }
}
/// Error-free product: `a*b = hi + lo` exactly, barring overflow or underflow.
#[inline(always)]
fn two_prod(a: f64, b: f64) -> Dd {
    let hi = a * b;
    Dd {
        hi,
        lo: a.mul_add(b, -hi),
    }
}
/// Dekker split for bounded double-double series. Callers keep |a| below
/// 1e280, so multiplying by 2^27+1 cannot overflow.
#[inline(always)]
pub(crate) fn split(a: f64) -> (f64, f64) {
    let c = 134_217_729.0 * a;
    let hi = c - (c - a);
    (hi, a - hi)
}

/// A fixed multiplier for bounded double-double series. Splitting its leading
/// part once avoids repeated FMA dispatch. Operands and intermediate products
/// must stay within the normal binary64 range and below 1e280 in magnitude.
#[derive(Clone, Copy)]
pub(crate) struct DdFactor {
    value: Dd,
    hi: f64,
    lo: f64,
}
impl DdFactor {
    #[inline]
    pub(crate) fn new(value: Dd) -> Self {
        let (hi, lo) = split(value.hi);
        Self { value, hi, lo }
    }
    #[inline(always)]
    fn product_error(self, a: f64, p: f64) -> f64 {
        let (hi, lo) = split(a);
        ((hi * self.hi - p) + hi * self.lo + lo * self.hi) + lo * self.lo
    }
    #[inline(always)]
    pub(crate) fn mul(self, a: Dd) -> Dd {
        let p = a.hi * self.value.hi;
        let error = self.product_error(a.hi, p);
        let lo = error + (a.hi * self.value.lo + a.lo * self.value.hi);
        let sum = p + lo;
        Dd {
            hi: sum,
            lo: lo - (sum - p),
        }
    }
}

/// Reuse a split binary64 multiplier in a series when both leading parts
/// are safely bounded. Outside this range retain the FMA-based product.
#[derive(Clone, Copy)]
pub(crate) enum SeriesFactor {
    Split(DdFactor),
    Fma(f64),
}
impl SeriesFactor {
    #[inline]
    pub(crate) fn new(value: f64) -> Self {
        if (1e-100..=1e100).contains(&value.abs()) {
            Self::Split(DdFactor::new(Dd::from(value)))
        } else {
            Self::Fma(value)
        }
    }
    #[inline(always)]
    pub(crate) fn product_error(self, value: f64, product: f64) -> f64 {
        match self {
            Self::Split(factor) if (1e-100..=1e100).contains(&value.abs()) => {
                factor.product_error(value, product)
            }
            Self::Split(factor) => value.mul_add(factor.value.hi, -product),
            Self::Fma(factor) => value.mul_add(factor, -product),
        }
    }
    #[inline(always)]
    pub(crate) fn mul(self, term: Dd) -> Dd {
        match self {
            Self::Split(factor) if (1e-100..=1e100).contains(&term.hi.abs()) => factor.mul(term),
            Self::Split(factor) => term.mul_f64(factor.value.hi),
            Self::Fma(value) => term.mul_f64(value),
        }
    }
}

/// Exact product residual for an integer factor of magnitude < 2^26 and a
/// bounded normal leading part, with no FMA dispatch.
#[inline(always)]
pub(crate) fn small_product_error(value: f64, n: f64, product: f64) -> f64 {
    let (hi, lo) = split(value);
    (hi * n - product) + lo * n
}

impl Dd {
    /// Exact sum of two binary64 values.
    #[inline(always)]
    pub(crate) fn sum(a: f64, b: f64) -> Self {
        two_sum(a, b)
    }
    /// Quotient by a double-double divisor, to about 2^-104 relative
    /// accuracy. One correction step suffices because `hi - q*divisor.hi`
    /// is exact; cheaper than `/`, which forms the residual with full
    /// double-double subtraction and addition.
    #[inline(always)]
    pub(crate) fn div_dd(self, rhs: Self) -> Self {
        let q = self.hi / rhs.hi;
        let product = rhs.mul_f64(q);
        let residual = (self.hi - product.hi) + (self.lo - product.lo);
        let correction = residual / rhs.hi;
        let hi = q + correction;
        Self {
            hi,
            lo: correction - (hi - q),
        }
    }
    /// Sum with a binary64 value; cheaper than promoting it to `Dd`.
    #[inline(always)]
    pub(crate) fn add_f64(self, rhs: f64) -> Self {
        let leading = two_sum(self.hi, rhs);
        two_sum(leading.hi, leading.lo + self.lo)
    }
    /// Product with a binary64 factor; cheaper than promoting it to `Dd`.
    #[inline(always)]
    pub(crate) fn mul_f64(self, rhs: f64) -> Self {
        let product = two_prod(self.hi, rhs);
        two_sum(product.hi, product.lo + self.lo * rhs)
    }
    /// Sum of two values of the same sign. Without cancellation, one
    /// error-free addition of the leading parts suffices; cheaper than `+`.
    #[inline(always)]
    pub(crate) fn add_same_sign(self, rhs: Self) -> Self {
        let leading = two_sum(self.hi, rhs.hi);
        let hi = leading.hi;
        let lo = leading.lo + (self.lo + rhs.lo);
        let sum = hi + lo;
        Self {
            hi: sum,
            lo: lo - (sum - hi),
        }
    }
    /// Product by an exact integer of magnitude < 2^26. The split products
    /// are exact, avoiding runtime FMA dispatch. Requires |self.hi| < 1e280.
    #[inline(always)]
    pub(crate) fn mul_small(self, n: f64) -> Self {
        let p = self.hi * n;
        let lo = small_product_error(self.hi, n, p) + self.lo * n;
        let sum = p + lo;
        Self {
            hi: sum,
            lo: lo - (sum - p),
        }
    }
    /// Product and quotient by exact small integers, with one normalization.
    /// Requires the same bounds as `mul_small` and `div_small`.
    #[inline(always)]
    pub(crate) fn mul_div_small(self, numerator: f64, denominator: f64) -> Self {
        let product = self.hi * numerator;
        let error = small_product_error(self.hi, numerator, product) + self.lo * numerator;
        let q = product / denominator;
        let p = q * denominator;
        let remainder = (product - p) - small_product_error(q, denominator, p);
        let correction = (remainder + error) / denominator;
        let hi = q + correction;
        Self {
            hi,
            lo: correction - (hi - q),
        }
    }

    /// Quotient by an exact nonzero integer of magnitude < 2^26.
    /// Requires |self.hi / n| < 1e280.
    #[inline(always)]
    pub(crate) fn div_small(self, n: f64) -> Self {
        let q = self.hi / n;
        let p = q * n;
        let error = small_product_error(q, n, p);
        let correction = (((self.hi - p) - error) + self.lo) / n;
        let sum = q + correction;
        Self {
            hi: sum,
            lo: correction - (sum - q),
        }
    }
}
impl Add for Dd {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let leading = two_sum(self.hi, rhs.hi);
        let trailing = two_sum(self.lo, rhs.lo);
        let partial = two_sum(leading.hi, leading.lo + trailing.hi);
        two_sum(partial.hi, partial.lo + trailing.lo)
    }
}
impl Neg for Dd {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}
impl Sub for Dd {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}
impl Mul for Dd {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let product = two_prod(self.hi, rhs.hi);
        two_sum(product.hi, product.lo + self.hi * rhs.lo + self.lo * rhs.hi)
    }
}
impl Div for Dd {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let q = Dd::from(self.hi / rhs.hi);
        let residual = self - rhs * q;
        q + Dd::from((residual.hi + residual.lo) / rhs.hi)
    }
}

#[cfg(test)]
mod tests {
    use super::{Dd, SeriesFactor};

    #[test]
    fn fixed_series_factor_preserves_product_residual() {
        // Compare against the hardware/library FMA residual, including
        // subnormal operands and values outside the split-safe range.
        for x in [
            -1e100,
            1e101,
            1e-100,
            1e-101,
            0.0,
            f64::from_bits(1),
            f64::MIN_POSITIVE,
            -0.743356278913,
        ] {
            let factor = SeriesFactor::new(x);
            for a in [
                -1e100,
                1e101,
                1e-100,
                1e-101,
                f64::MIN_POSITIVE,
                f64::from_bits(1),
                1.8672581934,
            ] {
                let term = Dd::sum(a, a * 1e-17);
                let product = term.hi * x;
                assert_eq!(
                    factor.product_error(term.hi, product),
                    term.hi.mul_add(x, -product)
                );
                let actual = factor.mul(term);
                let reference = term.mul_f64(x);
                assert_eq!(actual.hi, reference.hi);
                assert_eq!(actual.lo, reference.lo);
            }
        }
    }
}
