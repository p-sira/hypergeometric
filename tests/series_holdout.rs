//! Independent 140-digit Decimal references outside the Sci benchmark grids.
use hypergeometric::{hyp1f1, hyp2f1};

#[test]
fn independent_mixed_sign_series_holdout() {
    let mut count = 0;
    for line in include_str!("fixtures/series_holdout.csv").lines() {
        let values: Vec<_> = line.split(',').collect();
        let (got, expected) = match values[0] {
            "hyp1f1" => (
                hyp1f1(
                    values[1].parse().unwrap(),
                    values[2].parse().unwrap(),
                    values[3].parse().unwrap(),
                ),
                values[4].parse::<f64>().unwrap(),
            ),
            "hyp2f1" => (
                hyp2f1(
                    values[1].parse().unwrap(),
                    values[2].parse().unwrap(),
                    values[3].parse().unwrap(),
                    values[4].parse().unwrap(),
                ),
                values[5].parse::<f64>().unwrap(),
            ),
            _ => panic!("unexpected function"),
        };
        assert!(
            got.is_finite() && (got - expected).abs() <= 1e-11 * expected.abs(),
            "{line}: got {got:.17e}, reference {expected:.17e}"
        );
        count += 1;
    }
    assert_eq!(count, 320);
}
