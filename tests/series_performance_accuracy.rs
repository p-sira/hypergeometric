use hypergeometric::{hyp0f1, hyp1f1, hyp2f1};

fn check(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance * expected.abs(),
        "got {actual:.17e}, expected {expected:.17e}"
    );
}

// Independent 145-digit Decimal defining series evaluated at exact binary64
// inputs. Includes cancellation, positive series, and both sides of the
// bounded small-integer denominator dispatch.
#[test]
fn zero_f_one_bounded_series() {
    for (b, z, expected) in [
        (0.75, -2.0, -0.4669282148173254),
        (1.25, -4.0, -0.2856390277004589),
        (2.75, -10.0, -0.06866084682669293),
        (0.25, 9.0, 527.6762638934789),
        (0.24999999999999997, 9.0, 527.676263893479),
        (128.0, -100.0, 0.4567427874589717),
        (128.00000000000003, -100.0, 0.4567427874589718),
    ] {
        check(hyp0f1(b, z), expected, 1e-13);
    }
}

#[test]
fn one_f_one_cancellation_series() {
    for (a, b, z, expected) in [
        (
            -1.3624479232633842,
            2.779414522533803,
            2.333742,
            -1.0386817049695845e-8,
        ),
        (
            -2.2967123841077215,
            4.798076079147022,
            2.930414325277413,
            0.03587642839603475,
        ),
        (-1.7, -0.4, 2.25, -4.442896063964281),
    ] {
        check(hyp1f1(a, b, z), expected, 1e-13);
    }
}

#[test]
fn two_f_one_mixed_sign_series() {
    for (a, b, c, z, expected) in [
        (-1.7, 1.9, 0.75, 0.88, -0.39257433418382487),
        (-0.75, 1.8, 0.5, 0.8, -2.2088281286066076),
        (-1.5, -0.4, 1.2, -0.8, 0.5794223742266561),
        (
            -1.033198585223257,
            -1.274591987378873,
            1.0445164367736997,
            -0.7908204562605475,
            0.004613859851103431,
        ),
    ] {
        check(hyp2f1(a, b, c, z), expected, 1e-12);
    }
}
