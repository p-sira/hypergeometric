//! Decimal-based independent high precision regression references.
use hypergeometric::hyp1f1;

#[test]
fn binary64_benchmark_120_digit_series() {
    let mut count = 0;
    for (i, row) in include_str!("fixtures/hyp1f1_bench_exact.csv")
        .lines()
        .skip(1)
        .enumerate()
    {
        let row: Vec<_> = row.split(',').map(|v| v.parse::<f64>().unwrap()).collect();
        let (a, b, z, want) = (row[0], row[1], row[2], row[3]);
        let got = hyp1f1(a, b, z);
        assert!(
            (got - want).abs() <= 3e-13 * want.abs().max(1e-15),
            "case {i}, ({a},{b},{z}): got {got}, reference {want}"
        );
        count += 1;
    }
    assert_eq!(count, 256);
}

#[test]
fn cancellation_100_digit_series_and_scipy_comparison() {
    let mut count = 0;
    let mut misses = Vec::new();
    for (i, row) in include_str!("fixtures/hyp1f1_cancellation.csv")
        .lines()
        .skip(1)
        .enumerate()
    {
        let row: Vec<_> = row.split(',').map(|v| v.parse::<f64>().unwrap()).collect();
        let (a, b, z, want) = (row[0], row[1], row[2], row[3]);
        let got = hyp1f1(a, b, z);
        if (got - want).abs() > 2e-10 * want.abs().max(1e-15) {
            misses.push(format!(
                "{i}, ({a},{b},{z}): got {got}, reference {want}, ratio {}",
                got / want
            ));
        }
        count += 1;
    }
    assert_eq!(count, 239);
    assert!(
        misses.is_empty(),
        "{} mismatches:\n{}",
        misses.len(),
        misses
            .iter()
            .take(25)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}
