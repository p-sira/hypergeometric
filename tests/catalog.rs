//! Independent SciPy 1.18.1 fixture regression tests: input and expected are separate.
use hypergeometric::{
    hyp0f1, hyp1f1, hyp2f1, hyperu, log_wright_bessel, wright_bessel, wrightomega,
};
use std::collections::HashMap;

#[test]
fn original_sci_inventory_against_independent_scipy() {
    let mut expected = HashMap::new();
    for line in include_str!("fixtures/catalog_expected.csv")
        .lines()
        .skip(1)
    {
        let (id, value) = line.split_once(',').unwrap();
        expected.insert(id, value.parse::<f64>().unwrap());
    }
    let mut total = 0;
    let mut failures = Vec::new();
    for line in include_str!("fixtures/catalog_inputs.csv").lines().skip(1) {
        let fields: Vec<_> = line.split(',').collect();
        let (id, api) = (fields[0], fields[1]);
        let parse = |i: usize| fields[i].parse::<f64>().unwrap_or(f64::NAN);
        let (a, b, c, z) = (parse(2), parse(3), parse(4), parse(5));
        let got = match api {
            "hyp0f1" => hyp0f1(b, z),
            "hyp1f1" => hyp1f1(a, b, z),
            "hyp2f1" => hyp2f1(a, b, c, z),
            "hyperu" => hyperu(a, b, z),
            "wright_bessel" => wright_bessel(a, b, z),
            "log_wright_bessel" => log_wright_bessel(a, b, z),
            "wrightomega" => wrightomega(z),
            _ => panic!("unknown function {api}"),
        };
        let want = expected[id];
        // Relaxed tolerance for the decades-old Zhang-Jin U kernel; documented
        // nine-digit paths can lose precision and singularities need special checks.
        let rel = match api {
            "hyperu" => 3e-8,
            "hyp1f1" => 5e-9,
            "hyp2f1" => 5e-9,
            _ => 1e-10,
        };
        let abs = match api {
            "hyperu" => 1e-11,
            "hyp1f1" => 2e-10,
            _ => 1e-12,
        };
        let valid = if want.is_nan() {
            got.is_nan()
        } else if want.is_infinite() {
            got == want
        } else {
            (got - want).abs() <= abs + rel * want.abs()
        };
        if !valid {
            failures.push(format!(
                "row {id} {api} ({a},{b},{c},{z}) -> got {got:e}, expected {want:e}, relative error {}",
                (got - want).abs() / want.abs().max(1e-300)
            ));
        }
        total += 1;
    }
    assert!(total >= 1000, "insufficient coverage: {total}");
    assert_eq!(total, expected.len());
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures
            .iter()
            .take(65)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}
