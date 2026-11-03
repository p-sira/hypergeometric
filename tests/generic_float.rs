//! Match sci-special's generic scalar contracts without changing the f64 kernels.
use hypergeometric::{hyp0f1, wrightomega};
use num_traits::Float;

fn generic_call<T: Float>(b: T, z: T, x: T) -> (T, T) {
    (hyp0f1(b, z), wrightomega(x))
}

#[test]
fn accepts_f32_and_f64_through_float() {
    let (h32, w32): (f32, f32) = generic_call(0.5_f32, -1.0, 1.0);
    assert!((h32 - 2.0_f32.cos()).abs() < 2e-7);
    assert!((w32 + w32.ln() - 1.0).abs() < 3e-7);

    let (h64, w64): (f64, f64) = generic_call(0.5_f64, -1.0, 1.0);
    assert!((h64 - 2.0_f64.cos()).abs() < 1e-14);
    assert!((w64 + w64.ln() - 1.0).abs() < 1e-13);

    // Conversion follows the established f64 kernels, rounded to the
    // requested output type, rather than computing a separate f32 series.
    assert_eq!(h32, h64 as f32);
    assert_eq!(w32, w64 as f32);
}

#[test]
fn f32_domain_overflow_and_underflow_match_f64_conversion() {
    assert!(hyp0f1(0.0_f32, -1.0).is_nan());
    assert!(hyp0f1(f32::INFINITY, 1.0).is_nan());
    assert!(hyp0f1(1.0_f32, f32::NAN).is_nan());
    assert!(wrightomega(f32::NAN).is_nan());
    assert_eq!(wrightomega(f32::NEG_INFINITY), 0.0_f32);
    assert_eq!(wrightomega(f32::INFINITY), f32::INFINITY);

    // This result is representable in binary64 but overflows binary32.
    assert!(hyp0f1(0.5_f64, 2500.0).is_finite());
    assert_eq!(hyp0f1(0.5_f32, 2500.0), f32::INFINITY);

    // Wright omega remains nonzero at a subnormal f32 magnitude.
    let small = wrightomega(-100.0_f32);
    assert!(small > 0.0 && small.is_subnormal());
    assert_eq!(small, (-100.0_f64).exp() as f32);
    assert_eq!(wrightomega(-200.0_f32), 0.0_f32);
}
