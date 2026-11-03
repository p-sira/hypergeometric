//! Independently evaluated Wright functions at Sci's 64 CodSpeed inputs.
use hypergeometric::{log_wright_bessel, wright_bessel};

#[test]
fn wright_bessel_matches_independent_benchmark_references() {
    for row in include_str!("fixtures/wright_bench.csv").lines().skip(1) {
        let values: Vec<f64> = row.split(',').map(|s| s.parse().unwrap()).collect();
        let [a, b, z, expected_log] = <[f64; 4]>::try_from(values).unwrap();
        let actual_log = log_wright_bessel(a, b, z);
        assert!(
            (actual_log - expected_log).abs() <= 4e-13 * expected_log.abs().max(1.0),
            "log_wright_bessel({a},{b},{z})={actual_log}, ref={expected_log}"
        );
        let expected = expected_log.exp();
        let actual = wright_bessel(a, b, z);
        assert!(
            (actual / expected - 1.0).abs() <= 4e-13,
            "wright_bessel({a},{b},{z})={actual}, ref={expected}"
        );
    }
}
