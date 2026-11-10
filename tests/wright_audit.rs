//! Wright reference values independently evaluated with mpmath (90 digits)
//! and SciPy's `wrightomega`.
#![allow(clippy::excessive_precision)]
use hypergeometric::{log_wright_bessel, wright_bessel, wrightomega};

#[test]
fn wright_bessel_matches_independent_high_precision_values() {
    // (a, b, z, Phi(a,b;z), ln(Phi(a,b;z))) from a 90-digit mpmath sum.
    let cases = [
        (0.5, 1.5, 20.0, 75886.172737090194, 11.236989769410466),
        (0.1, 0.0, 1.0, 0.2943338753714321, -1.223040525468938),
        (1.0, 1.0, 10.0, 90.47595439632761, 4.505084118123957),
        (2.0, 0.0, 1.0, 1.0847305124172486, 0.08133158048977667),
        (5.0, 0.0, 1.0, 0.04166804453453966, -3.1780207620657555),
    ];
    for (a, b, z, want, want_log) in cases {
        let got = wright_bessel(a, b, z);
        let got_log = log_wright_bessel(a, b, z);
        assert!(
            (got - want).abs() <= 2e-13 * want,
            "Phi({a},{b};{z})={got:e}"
        );
        assert!(
            (got_log - want_log).abs() <= 2e-13,
            "ln Phi({a},{b};{z})={got_log:e}"
        );
    }
}

#[test]
fn wright_unit_step_gamma_recurrence_matches_bessel_identity() {
    // 90-digit mpmath 1.3.0, using the independent identity
    // Phi(1,b;z) = z^((1-b)/2) I_(b-1)(2 sqrt(z)).
    let reference_log = [
        (0.0, 0.01, -4.6001743457241728),
        (0.0, 1.0, 0.46413447354615972),
        (0.0, 100.0, 19.86653971551339),
        (0.1, 50.0, 13.638479295234283),
        (0.5, 100.0, 18.734487876515356),
        (1.0, 10.0, 4.505084118123957),
        (1.0, 1000.0, 60.255094572585392),
        (2.0, 50.0, 9.915088174547206),
        (3.5, 100.0, 11.672998497591175),
        (15.5, 1000.0, 8.50584263881101),
        (0.001, 0.05, -2.950963014992705),
        (100.0, 0.5, -359.12920549332966),
    ];
    for (b, z, want_log) in reference_log {
        let got = log_wright_bessel(1.0, b, z);
        assert!(
            (got - want_log).abs() <= 3e-12 * want_log.abs().max(1.0),
            "ln Phi(1,{b};{z})={got:.17e} vs {want_log:.17e}"
        );
    }
    // The public domain check rejects b<0 before selecting the recurrence.
    for b in [-0.5, -1.0, -100.0] {
        assert!(log_wright_bessel(1.0, b, 2.0).is_nan());
    }
}

#[test]
fn gamma_overflow_stops_a_converged_wright_series() {
    // At this scale the first term is the only term that can affect f64.
    // Every later log-Gamma overflows, although the mathematical sum is finite.
    for &a in &[1e307, 1e308, f64::MAX] {
        assert_eq!(log_wright_bessel(a, 1.0, 1.0), 0.0);
        assert_eq!(wright_bessel(a, 1.0, 1.0), 1.0);
        assert!((log_wright_bessel(a, 3.0, 1.0) + 2.0_f64.ln()).abs() < 2e-14);
        assert!((wright_bessel(a, 3.0, 1.0) - 0.5).abs() < 2e-14);
        // b=0 removes the k=0 term; all remaining terms underflow, too.
        assert_eq!(log_wright_bessel(a, 0.0, 1.0), f64::NEG_INFINITY);
        assert_eq!(wright_bessel(a, 0.0, 1.0), 0.0);
    }
}

#[test]
fn wrightomega_matches_independent_reference() {
    let cases: [(f64, f64); 9] = [
        (-50.0, 1.9287498479639178e-22),
        (-20.0, 2.061153618190205e-9),
        (-2.0, 0.1200282389876412),
        (-1.0, 0.27846454276107374),
        (0.0, 0.5671432904097838),
        (0.5, 0.7662486081617502),
        (2.0, 1.5571455989976113),
        (10.0, 7.9294200950196965),
        (100.0, 95.44148664557584),
    ];
    for (x, want) in cases {
        let got = wrightomega(x);
        assert!((got - want).abs() <= 8e-15 * want, "omega({x})={got:e}");
    }
}

#[test]
#[ignore = "manual performance comparison: cargo test --release --test wright_audit -- --ignored --nocapture"]
fn targeted_wright_release_performance() {
    use std::hint::black_box;
    use std::time::Instant;

    for &(name, a, b, z, iterations) in &[
        ("moderate", 0.5, 1.5, 20.0, 20_000),
        ("long_series", 0.01, 0.1, 1000.0, 2_000),
        ("gamma_overflow", 1e308, 1.0, 1.0, 200),
    ] {
        for _ in 0..16 {
            black_box(log_wright_bessel(black_box(a), black_box(b), black_box(z)));
        }
        let start = Instant::now();
        let mut last = 0.0;
        for _ in 0..iterations {
            last = black_box(log_wright_bessel(black_box(a), black_box(b), black_box(z)));
        }
        let ns = start.elapsed().as_secs_f64() * 1e9 / iterations as f64;
        println!("{name}: {ns:.1} ns/call (n={iterations}, result={last})");
    }
}
