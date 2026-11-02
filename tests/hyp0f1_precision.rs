// Exact binary64 0F1 inputs and independently generated 105-digit references.
use hypergeometric::hyp0f1;

#[test]
fn negative_argument_shortcut_matches_independent_high_precision() {
    // 105-digit Decimal defining series evaluated at exactly rounded binary64
    // inputs, including the two adjacent floats around a cosine zero.
    // Ordinary absolute floors can hide sign and relative errors here.
    let input = include_str!("scipy/hyp0f1_precision_inputs.csv");
    let reference = include_str!("fixtures/hyp0f1_precision_expected.csv");
    let input_lines = input.lines().skip(1);
    let reference_lines = reference.lines().skip(1);
    assert_eq!(input_lines.clone().count(), reference_lines.clone().count());
    for (row, value) in input_lines.zip(reference_lines) {
        let (v, z) = row.split_once(',').unwrap();
        let v = v.parse::<f64>().unwrap();
        let z = z.parse::<f64>().unwrap();
        let expected = value.parse::<f64>().unwrap();
        let actual = hyp0f1(v, z);
        assert!(
            (actual - expected).abs() <= 3e-12 * expected.abs(),
            "0F1(;{v};{z}): got {actual}, reference {expected}"
        );
    }
}

#[test]
fn negative_argument_shortcut_respects_recurrence_and_poles() {
    // Contiguous relation from the defining Pochhammer series:
    // F(;b-1;z) = F(;b;z) + z F(;b+1;z)/(b(b-1)).
    for &v in &[0.5_f64, 1.125, 2.5, 5.0, 10.0, 25.0] {
        for &z in &[-100.0, -25.0, -10.0, -3.0001, -3.0, -2.9999, 2.0, 20.0] {
            let lhs = hyp0f1(v - 1.0, z);
            let first = hyp0f1(v, z);
            let second = z / (v * (v - 1.0)) * hyp0f1(v + 1.0, z);
            let rhs = first + second;
            // Contiguous relation may cancel nearly to zero. Its natural
            // error scale is the sum of magnitudes of its three terms.
            let scale = lhs.abs() + first.abs() + second.abs();
            assert!(
                (lhs - rhs).abs() <= 2e-10 * scale,
                "0F1 recurrence v={v}, z={z}: {lhs} vs {rhs}"
            );
        }
    }
    for &v in &[0.0_f64, -1.0, -2.0, -3.0] {
        assert!(hyp0f1(v, -10.0).is_nan());
    }
    assert!(hyp0f1(1.0, f64::NAN).is_nan());
}
