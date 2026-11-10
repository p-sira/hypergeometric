//! Independent references generated at 90 and 180 decimal digits from the
//! exact binary64 inputs. See reference/generate_accuracy_audit.py.
use hypergeometric::{
    hyp0f1, hyp1f1, hyp2f1, hyperu, log_wright_bessel, wright_bessel, wrightomega,
};

#[test]
fn high_precision_grid_and_boundary_regressions() {
    let mut count = 0;
    for line in include_str!("fixtures/numerical_audit.csv").lines() {
        let (name, values) = line.split_once(',').unwrap();
        let values: Vec<f64> = values.split(',').map(|s| s.parse().unwrap()).collect();
        let want = *values.last().unwrap();
        let got = match name {
            "hyp0f1" => hyp0f1(values[0], values[1]),
            "hyp1f1" => hyp1f1(values[0], values[1], values[2]),
            "hyp2f1" => hyp2f1(values[0], values[1], values[2], values[3]),
            "hyperu" => hyperu(values[0], values[1], values[2]),
            "log_wright_bessel" => log_wright_bessel(values[0], values[1], values[2]),
            "wrightomega" => wrightomega(values[0]),
            _ => panic!("unknown API: {name}"),
        };
        // The legacy U kernel documents roughly nine valid digits on some
        // paths. Other sampled regions should retain at least eleven.
        let rtol = if name == "hyperu" { 3e-8 } else { 1e-11 };
        assert_reference(got, want, rtol, line);
        if name == "log_wright_bessel" {
            assert_reference(
                wright_bessel(values[0], values[1], values[2]),
                want.exp(),
                1e-11,
                line,
            );
        }
        count += 1;
    }
    assert!(count >= 1900, "unexpected fixture count: {count}");
}

fn assert_reference(got: f64, want: f64, rtol: f64, case: &str) {
    if want == 0.0 || want.is_infinite() {
        assert_eq!(got, want, "{case}");
    } else {
        assert!(
            got.is_finite() && (got - want).abs() <= rtol * want.abs(),
            "{case}: got {got:.17e}, relative error {:.3e}",
            (got - want).abs() / want.abs()
        );
    }
}
