//! Small nonzero values must be validated relative to their magnitude.
use hypergeometric::{
    hyp0f1, hyp1f1, hyp2f1, hyperu, log_wright_bessel, wright_bessel, wrightomega,
};

fn relative(got: f64, want: f64, rtol: f64) {
    assert!(want.is_finite() && want != 0.0 && got.is_finite());
    assert!(
        (got - want).abs() <= rtol * want.abs(),
        "got {got:e}, expected {want:e}, relative error {:e}",
        (got - want).abs() / want.abs()
    );
}

#[test]
fn kummer_equal_parameters_preserve_exponentially_small_values() {
    // M(a,a,z) = exp(z); a zero result would have passed old absolute floors.
    for &a in &[0.1, 1.0, 10.0] {
        for &z in &[-200.0, -120.0, -40.0] {
            relative(hyp1f1(a, a, z), z.exp(), 3e-14);
        }
    }
    for &z in &[-0.9, -0.25, 0.25, 0.9] {
        relative(hyp2f1(1.0, 1.0, 2.0, z), -(1.0 - z).ln() / z, 3e-12);
    }
}

#[test]
fn wright_functions_preserve_small_positive_outputs() {
    relative(wrightomega(-100.0), (-100.0_f64).exp(), 2e-14);
    // Φ(0,b;z)=exp(z)/Γ(b). Here Γ(51)=50!.
    let factorial = (1..=50).fold(1.0_f64, |acc, n| acc * f64::from(n));
    let want = 0.5_f64.exp() / factorial;
    assert!(want > 0.0 && want < 1e-60);
    relative(wright_bessel(0.0, 51.0, 0.5), want, 3e-13);
    relative(log_wright_bessel(0.0, 51.0, 0.5).exp(), want, 3e-13);
}

#[test]
fn hyperu_contiguous_identity_at_small_nonzero_values() {
    // DLMF 13.2.40: U(a,b,z)=z^(1-b) U(1+a-b,2-b,z).
    for &z in &[2.0, 10.0, 40.0] {
        let lhs = hyperu(12.0, -3.2, z);
        let rhs = z.powf(4.2) * hyperu(16.2, 5.2, z);
        assert!(lhs > 0.0 && lhs < 1e-8, "unexpected U scale {lhs:e}");
        relative(lhs, rhs, 5e-8);
    }
}

#[test]
fn highly_oscillatory_half_integer_0f1_stays_bounded() {
    // Root-precision correction is local: on huge arguments an ulp of sqrt
    // spans many radians and a linear phase correction must not diverge.
    for &z in &[-1e8, -1e100, -1e308] {
        let value = hyp0f1(0.5, z);
        assert!(
            value.is_finite() && value.abs() <= 1.0,
            "0F1(.5,{z})={value}"
        );
    }
}
