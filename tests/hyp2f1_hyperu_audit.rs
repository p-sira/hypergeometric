//! Independent 100-digit mpmath reference values for the Gauss and Tricomi
//! functions. U(-n,b,z) values come from mpmath's generalized Laguerre
//! polynomial; the 2F1 references use high-precision log1p.
#![allow(clippy::excessive_precision)]
use hypergeometric::{hyp2f1, hyperu};

fn close(got: f64, want: f64, rtol: f64) {
    assert!(got.is_finite(), "expected finite {want:e}, got {got:e}");
    assert!(
        (got - want).abs() <= rtol * want.abs(),
        "got {got:.17e}, expected {want:.17e}; relative error {:.3e}",
        (got - want).abs() / want.abs()
    );
}

#[test]
fn two_f_one_logarithmic_identity_across_scales() {
    // 2F1(1,1;2;z) = -log(1-z)/z.  The references were evaluated from
    // the exact binary64 z inputs at 100 digits of mpmath precision.
    let reference = [
        (-1e300, 6.9077552789821366899e-298),
        (-1e10, 2.3025850930040456840e-9),
        (-10.0, 0.2397895272798370544),
        (-1e-12, 0.9999999999995),
        (-1e-300, 1.0),
        (1e-300, 1.0),
        (1e-12, 1.0000000000005),
        (0.1, 1.0536051565782630155),
        (0.5, 1.3862943611198906188),
        (0.999999999999, 27.631043237920989003),
        (1.0 - f64::EPSILON, 36.043653389117164093),
    ];
    for (z, want) in reference {
        close(hyp2f1(1.0, 1.0, 2.0, z), want, 5e-14);
    }
    assert_eq!(hyp2f1(1.0, 1.0, 2.0, 0.0), 1.0);
    assert_eq!(hyp2f1(1.0, 1.0, 2.0, 1.0), f64::INFINITY);
}

#[test]
fn negative_integer_u_matches_high_precision_laguerre() {
    // Independent 100-digit mpmath: (-1)^n * n! * laguerre(n,b-1,z).
    let reference = [
        (1, 3.5, 2.0, -1.5),
        (2, 2.0, 0.05, 5.7025),
        (10, 1.5, 0.25, -622594.9998083115),
        (20, 2.5, 0.1, 6.996386106082293505e19),
        (30, 4.5, 1.0, 3.560678909690723248e34),
        (50, 3.5, 2.0, -8.578345147783045073e65),
        (100, 3.5, 2.0, 1.712131469552552809e159),
        (100, 10.0, 2.0, -2.045331566021620668e165),
        (100, 4.5, 5.0, 7.664912929648721514e159),
        (128, 2.5, 0.1, -1.317608246558791942e217),
    ];
    for (n, b, z, want) in reference {
        close(hyperu(-(n as f64), b, z), want, 6e-12);
    }
}

#[test]
fn negative_integer_u_polynomial_zeros_and_zero_lower_parameter() {
    // U(-1,b,z) = z-b and U(-2,b,z) = z^2-2(b+1)z+b(b+1).
    assert_eq!(hyperu(-1.0, 2.0, 2.0), 0.0);
    close(hyperu(-1.0, 0.0, 0.2), 0.2, 5e-15);
    let near_root = 2.0 - 2.0_f64.sqrt();
    assert!(hyperu(-2.0, 1.0, near_root).abs() < 1e-14);
    close(hyperu(-2.0, 0.0, 0.05), -0.0975, 5e-14);
}
