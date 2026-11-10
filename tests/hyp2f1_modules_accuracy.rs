//! High-precision Decimal references for Sci's exact 1,000-input 2F1 workload.
use hypergeometric::hyp2f1;

#[test]
fn sci_modules_workload() {
    let mut count = 0;
    for (i, line) in include_str!("fixtures/hyp2f1_modules_exact.csv")
        .lines()
        .skip(1)
        .enumerate()
    {
        let (z, expected) = line.split_once(',').unwrap();
        let z: f64 = z.parse().unwrap();
        let expected: f64 = expected.parse().unwrap();
        let actual = hyp2f1(0.5, 1.5, 2.5, z);
        assert!(
            actual.is_finite() && (actual - expected).abs() <= 1e-12 * expected.abs(),
            "case {i}: z={z}, actual={actual:.17e}, expected={expected:.17e}"
        );
        count += 1;
    }
    assert_eq!(count, 1000);
}
