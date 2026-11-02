use hypergeometric::hyp0f1;

#[test]
fn large_argument_reference() {
    // Arbitrary precision mpmath hypergeometric values, supplied independently.
    let data: &[(f64, f64, f64)] = &[
        (0.25, -1000.0, 3.3485985350492914),
        (0.25, -10000.0, 5.076412235732962),
        (1.25, -1000.0, 0.027640104044790954),
        (1.25, -10000.0, -0.010038357459017138),
        (4.0, -1000.0, 5.674432067916388e-6),
        (100.0, -1000.0, 2.554798380348188e-5),
        (-0.3, -1000.0, 4.283222028923697),
        (-0.3, 1000.0, -5.607614125355963e28),
        (0.25, 200.0, 3.7902228868978423e12),
        (0.25, 1000.0, 7.093970113691315e27),
        (50.0, 10000.0, 3.1136144040041786e47),
    ];
    for &(b, z, expected) in data {
        let got = hyp0f1(b, z);
        assert!(
            (got / expected - 1.0).abs() < 2e-10,
            "0F1(;{b};{z})={got:?} expected {expected:?}"
        );
    }
}

#[test]
fn wide_range_independent_high_precision() {
    // 180-digit mpmath 1.3.0 references for the exact f64 inputs, including
    // near b-poles, large positive/negative arguments, and large positive b.
    let cases = include_str!("scipy/hyp0f1_large_inputs.csv");
    let values = include_str!("fixtures/hyp0f1_large_expected.csv");
    let input_lines = cases.lines().skip(1);
    let reference_lines = values.lines().skip(1);
    assert_eq!(input_lines.clone().count(), reference_lines.clone().count());
    assert!(input_lines.clone().count() >= 190);
    for (row, expected) in input_lines.zip(reference_lines) {
        let (b, z) = row.split_once(',').unwrap();
        let b = b.parse::<f64>().unwrap();
        let z = z.parse::<f64>().unwrap();
        let expected = expected.parse::<f64>().unwrap();
        let actual = hyp0f1(b, z);
        let rel = (actual / expected - 1.0).abs();
        assert!(
            rel <= 2e-10,
            "0F1(;{b};{z}): got {actual}, reference {expected}, relative {rel}"
        );
    }
}

#[test]
fn large_b_negative_argument_and_overflow() {
    // The normalized Bessel value can be finite even when J_{b-1} underflows.
    let b = 100000.0;
    let got = hyp0f1(b, -1e7);
    let expected = 3.5384111609128018e-44;
    assert!(
        (got / expected - 1.0).abs() < 2e-10,
        "got {got}, expected {expected}"
    );
    assert!(hyp0f1(0.25, 1e6).is_infinite());
    assert!(hyp0f1(0.25, -1e4).is_finite());
    // Exact sinh identity stays finite even when the intermediate sinh
    // overflows: at z = 126025, 2 sqrt(z) = 710.
    let t = 710.0_f64;
    let expected = (t - (2.0 * t).ln()).exp();
    assert!((hyp0f1(1.5, 126025.0) / expected - 1.0).abs() < 1e-13);
}

#[test]
fn large_argument_contiguous_relation() {
    // 0F1(;b-1;z) = 0F1(;b;z) + z/[b(b-1)] * 0F1(;b+1;z).
    // Scale by the participating terms so cancellation does not create
    // a vacuous fixed absolute tolerance for small nonzero values.
    for b in [1.25_f64, 4.0, 10.0, 100.0, 500.0, 1000.0, 100000.0] {
        for z in [-10000.0_f64, -1000.0, -200.0, 200.0, 1000.0, 10000.0] {
            let lhs = hyp0f1(b - 1.0, z);
            let first = hyp0f1(b, z);
            let second = z / (b * (b - 1.0)) * hyp0f1(b + 1.0, z);
            let scale = lhs.abs() + first.abs() + second.abs();
            if scale.is_finite() {
                assert!(lhs.is_finite() && first.is_finite() && second.is_finite());
                assert!(
                    (lhs - first - second).abs() <= 2e-9 * scale,
                    "contiguous recurrence b={b}, z={z}: {lhs} vs {first}+{second}"
                );
            }
        }
    }
}
