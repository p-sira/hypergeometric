# hypergeometric

<p align="center">
    <a href="LICENSE" style="text-decoration:none">
        <img src="https://img.shields.io/badge/License-BSD--3--Clause-brightgreen.svg" alt="License">
    </a>
    <a href="https://crates.io/crates/hypergeometric" style="text-decoration:none">
        <img src="https://img.shields.io/crates/v/hypergeometric" alt="Crate">
    </a>
    <a href="https://docs.rs/hypergeometric" style="text-decoration:none">
        <img src="https://img.shields.io/badge/Docs-docs.rs-blue" alt="Documentation">
    </a>
    <a href="https://app.codspeed.io/p-sira/hypergeometric?utm_source=badge" style="text-decoration:none">
        <img src="https://img.shields.io/endpoint?url=https://codspeed.io/badge.json" alt="CodSpeed"/>
    </a>
</p>

Pure Rust library for real hypergeometric and Wright special functions, extracted from [Sci](https://github.com/p-sira/sci).

- Evaluate the confluent limit function 0F1, Kummer's 1F1, and Gauss's 2F1.
- Compute Tricomi's confluent hypergeometric function U.
- Evaluate the Wright generalized Bessel function, its logarithm, and Wright omega.

```rust
use hypergeometric::{hyp0f1, hyp1f1, hyp2f1, hyperu, wright_bessel, log_wright_bessel, wrightomega};

assert_eq!(hyp0f1(2.5_f64, 0.0), 1.0);
assert!((hyp1f1(1.0, 1.0, 0.5) - 0.5_f64.exp()).abs() < 1e-14);
assert!((hyp2f1(1.0, 1.0, 2.0, 0.5) - 4.0_f64.ln()).abs() < 1e-12);
assert_eq!(hyperu(0.0, 2.5, 1.0), 1.0);
assert!((wright_bessel(0.0, 1.0, 1.0) - 1.0_f64.exp()).abs() < 1e-14);
assert!((log_wright_bessel(0.0, 1.0, 1.0) - 1.0).abs() < 1e-14);
let w: f64 = wrightomega(1.0);
assert!((w + w.ln() - 1.0).abs() < 1e-13);
```

`hyp0f1` and `wrightomega` accept both `f32` and `f64` through `num_traits::Float`; they evaluate in `f64` and return the input type. The remaining functions use `f64`. All functions operate on real arguments; complex continuation is not supported. Invalid inputs and unresolved numerical regions may return `NaN`, while poles or overflow may produce infinities. The crate uses `bessel` for large-argument 0F1 continuation. Numerical accuracy is checked against independent SciPy, Decimal, and mpmath reference fixtures in `tests/`.

Adapted Sci and SciPy kernels and Wright omega's Algorithm 917 retain their license notices in [`LICENSE`](LICENSE) and [`LICENSES/`](LICENSES/).
