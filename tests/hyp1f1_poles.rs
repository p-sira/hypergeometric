//! Independently checked against SciPy 1.17 and mpmath 1.3.
use hypergeometric::hyp1f1;

#[test]
fn equal_negative_integer_parameters_terminate() {
    // M(-n,-n,z) = sum_(k=0)^n z^k/k!, rather than exp(z).
    for &(n, z, expected) in &[
        (1.0, 1.0, 2.0),
        (2.0, 1.0, 2.5),
        (3.0, 2.0, 6.333_333_333_333_333),
        (2.0, -2.0, 1.0),
    ] {
        let got = hyp1f1(-n, -n, z);
        assert!(
            (got - expected).abs() <= 2e-15 * expected.abs(),
            "n={n}, z={z}: {got} vs {expected}"
        );
    }
}

#[test]
fn equal_noninteger_parameters_remain_exponential() {
    // The equal-parameter identity applies on either side of zero provided
    // the denominator is not an integer pole.
    for &(a, z) in &[(-0.5, -700.0), (-2.5, 50.0), (0.5, -50.0)] {
        assert_eq!(hyp1f1(a, a, z), z.exp());
    }
}

#[test]
fn nonpositive_integer_denominator_poles() {
    for &(a, b) in &[(0.0, 0.0), (0.0, -2.0), (-2.0, -1.0), (2.0, 0.0)] {
        for &z in &[0.0, 1.0] {
            assert!(hyp1f1(a, b, z).is_infinite(), "a={a}, b={b}, z={z}");
        }
    }
    assert_eq!(hyp1f1(-1.0, -2.0, 1.0), 1.5);
    assert_eq!(hyp1f1(-2.0, -2.0, 0.0), 1.0);
}

#[test]
fn contiguous_exponential_identity_preserves_small_values() {
    // M(b+1,b,z) = exp(z)*(1+z/b). Independent mpmath references are rounded
    // to binary64, including values for which exp(z) itself underflows.
    let cases: &[(f64, f64, f64)] = &[
        (1.0, 1.0, 5.43656365691809),
        (1.0, -745.0, -2.1e-321),
        (1.0, -750.0, -1.5e-323),
        (2.0, -2.0, 0.0),
        (2.0, -2.0000000000000004, -3.005046949869086e-17),
        (1.5, -1.5000000000000002, -3.302989883801137e-17),
    ];
    for &(b, z, expected) in cases {
        let actual = hyp1f1(b + 1.0, b, z);
        let tolerance = (2e-14 * expected.abs()).max(f64::from_bits(1));
        assert!(
            (actual - expected).abs() <= tolerance,
            "b={b}, z={z}: got {actual:e} vs {expected:e}"
        );
    }
    // Even if (b+z)/b overflows, the final mathematical result can vanish.
    assert_eq!(hyp1f1(1.5, 0.5, -1e308), 0.0);
    // Rounding b+1 can change the mathematical parameter difference. For
    // these exact binary64 inputs, M(1.1,0.1,-0.1) is not exactly zero.
    let near_zero = hyp1f1(1.1, 0.1, -0.1);
    let expected = -7.186292542843042e-17;
    assert!((near_zero / expected - 1.0).abs() < 2e-11, "{near_zero:e}");
}
