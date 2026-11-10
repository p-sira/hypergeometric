//! Independent Decimal references for the exact 256-input Sci benchmark grid.
use hypergeometric::hyp2f1;

#[test]
fn scipy_workload_matches_exact_binary64_references() {
    let mut count = 0;
    for (index, line) in include_str!("fixtures/hyp2f1_bench_exact.csv")
        .lines()
        .skip(1)
        .enumerate()
    {
        let values: Vec<f64> = line
            .split(',')
            .map(|value| value.parse().unwrap())
            .collect();
        let (a, b, c, z, expected) = (values[0], values[1], values[2], values[3], values[4]);
        let actual = hyp2f1(a, b, c, z);
        assert!(
            actual.is_finite() && (actual - expected).abs() <= 1e-11 * expected.abs(),
            "case {index}, ({a},{b},{c},{z}): got {actual:.17e}, reference {expected:.17e}"
        );
        count += 1;
    }
    assert_eq!(count, 256);
}
