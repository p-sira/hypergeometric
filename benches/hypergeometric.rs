//! CodSpeed/divan benchmarks for the public hypergeometric and Wright APIs.
//!
//! Each function is benchmarked on a few representative numerical regimes
//! (series, transformations, asymptotics), plus a sweep over the SciPy
//! reference input grid in `tests/fixtures/inputs.csv`.

use divan::black_box;
use hypergeometric::{
    hyp0f1, hyp1f1, hyp2f1, hyperu, log_wright_bessel, wright_bessel, wrightomega,
};

fn main() {
    divan::main();
}

mod hyp0f1_bench {
    use super::*;

    /// (b, z) cases: small series, near-pole b, large negative/positive z (Bessel paths).
    const CASES: &[(&str, f64, f64)] = &[
        ("series", 2.5, 1.5),
        ("near_pole_b", -2.999999, 3.0),
        ("bessel_j", 0.5, -500.0),
        ("bessel_i", 1.5, 500.0),
        ("large_b_neg_z", 200.0, -1.0e5),
    ];

    #[divan::bench(args = CASES.iter().map(|c| c.0))]
    fn hyp0f1_case(name: &str) -> f64 {
        let &(_, b, z) = CASES.iter().find(|c| c.0 == name).unwrap();
        hyp0f1(black_box(b), black_box(z))
    }

    #[divan::bench]
    fn f32_series() -> f32 {
        hyp0f1(black_box(2.5_f32), black_box(1.5_f32))
    }
}

mod hyp1f1_bench {
    use super::*;

    /// (a, b, z) cases: direct series, Kummer transformation, large argument, cancellation.
    const CASES: &[(&str, f64, f64, f64)] = &[
        ("series", 1.5, 2.5, 0.75),
        ("negative_z", 2.059, 4.488, -4.29),
        ("large_z", 0.5, 1.5, 200.0),
        ("cancellation", -50.5, 2.0, 30.0),
        ("large_a", 1000.0, 10.0, 5.0),
    ];

    #[divan::bench(args = CASES.iter().map(|c| c.0))]
    fn hyp1f1_case(name: &str) -> f64 {
        let &(_, a, b, z) = CASES.iter().find(|c| c.0 == name).unwrap();
        hyp1f1(black_box(a), black_box(b), black_box(z))
    }
}

mod hyp2f1_bench {
    use super::*;

    /// (a, b, c, z) cases across the real continuation regions.
    const CASES: &[(&str, f64, f64, f64, f64)] = &[
        ("small_z", 1.0, 1.0, 2.0, 0.5),
        ("near_one", 0.5, 1.5, 2.25, 0.95),
        ("negative_z", 1.25, 2.5, 3.5, -5.0),
        ("integer_c_minus_ab", 1.0, 2.0, 3.0, 0.9),
        ("integer_d_nontrivial", 1.25, 0.75, 2.0, 0.999),
        ("terminating", -10.0, 3.5, 2.0, 2.5),
    ];

    #[divan::bench(args = CASES.iter().map(|c| c.0))]
    fn hyp2f1_case(name: &str) -> f64 {
        let &(_, a, b, c, z) = CASES.iter().find(|c| c.0 == name).unwrap();
        hyp2f1(black_box(a), black_box(b), black_box(c), black_box(z))
    }
}

mod hyperu_bench {
    use super::*;

    /// (a, b, z) cases: series, b=0 quadrature, recurrence, asymptotic.
    const CASES: &[(&str, f64, f64, f64)] = &[
        ("small_z", 1.5, 2.5, 0.5),
        ("b_zero_integral", 1.5, 0.0, 2.0),
        ("recurrence", 0.1, 1.0, 0.5),
        ("large_z", 2.0, 3.0, 50.0),
        ("negative_a", -3.5, 1.25, 4.0),
        ("integer_b_positive", 1.25, 16.0, 0.75),
        ("integer_b_negative", 1.25, -16.0, 0.75),
    ];

    #[divan::bench(args = CASES.iter().map(|c| c.0))]
    fn hyperu_case(name: &str) -> f64 {
        let &(_, a, b, z) = CASES.iter().find(|c| c.0 == name).unwrap();
        hyperu(black_box(a), black_box(b), black_box(z))
    }
}

mod wright_bench {
    use super::*;

    /// (a, b, z) cases for the Wright generalized Bessel function.
    const CASES: &[(&str, f64, f64, f64)] = &[
        ("small_z", 0.5, 1.5, 1.0),
        ("moderate_z", 1.0, 2.0, 50.0),
        ("large_z_small_a", 0.1, 1.0, 100.0),
    ];

    #[divan::bench(args = CASES.iter().map(|c| c.0))]
    fn wright_bessel_case(name: &str) -> f64 {
        let &(_, a, b, z) = CASES.iter().find(|c| c.0 == name).unwrap();
        wright_bessel(black_box(a), black_box(b), black_box(z))
    }

    #[divan::bench(args = CASES.iter().map(|c| c.0))]
    fn log_wright_bessel_case(name: &str) -> f64 {
        let &(_, a, b, z) = CASES.iter().find(|c| c.0 == name).unwrap();
        log_wright_bessel(black_box(a), black_box(b), black_box(z))
    }

    #[divan::bench(args = [-60.0, -1.0, 1.0, 1.0e3])]
    fn wrightomega_f64(x: f64) -> f64 {
        wrightomega(black_box(x))
    }
}

/// Batch evaluation over the SciPy reference input grid, grouped by function.
mod reference_grid {
    use super::*;

    const INPUTS: &str = include_str!("../tests/fixtures/inputs.csv");

    struct Row {
        a: f64,
        b: f64,
        c: f64,
        z: f64,
    }

    fn parse(v: &str) -> f64 {
        if v.is_empty() {
            f64::NAN
        } else {
            v.parse().unwrap()
        }
    }

    fn rows(func: &str) -> Vec<Row> {
        INPUTS
            .lines()
            .skip(1)
            .filter_map(|line| {
                let f: Vec<&str> = line.split(',').collect();
                (f[1] == func).then(|| Row {
                    a: parse(f[2]),
                    b: parse(f[3]),
                    c: parse(f[4]),
                    z: parse(f[5]),
                })
            })
            .collect()
    }

    fn sweep(bencher: divan::Bencher, func: &str, eval: fn(&Row) -> f64) {
        let data = rows(func);
        assert!(!data.is_empty());
        bencher.bench_local(|| {
            let mut acc = 0.0;
            for row in black_box(&data) {
                acc += black_box(eval(row));
            }
            acc
        });
    }

    #[divan::bench]
    fn hyp0f1_grid(bencher: divan::Bencher) {
        sweep(bencher, "hyp0f1", |r| hyp0f1(r.b, r.z));
    }

    #[divan::bench]
    fn hyp1f1_grid(bencher: divan::Bencher) {
        sweep(bencher, "hyp1f1", |r| hyp1f1(r.a, r.b, r.z));
    }

    #[divan::bench]
    fn hyp2f1_grid(bencher: divan::Bencher) {
        sweep(bencher, "hyp2f1", |r| hyp2f1(r.a, r.b, r.c, r.z));
    }

    #[divan::bench]
    fn hyperu_grid(bencher: divan::Bencher) {
        sweep(bencher, "hyperu", |r| hyperu(r.a, r.b, r.z));
    }

    #[divan::bench]
    fn wright_bessel_grid(bencher: divan::Bencher) {
        sweep(bencher, "wright_bessel", |r| wright_bessel(r.a, r.b, r.z));
    }

    #[divan::bench]
    fn wrightomega_grid(bencher: divan::Bencher) {
        sweep(bencher, "wrightomega", |r| wrightomega(r.z));
    }
}
