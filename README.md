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

let a = hyp0f1(2.5_f64, 1.0); // 1.4615741153700916
let b = hyp1f1(1.0, 2.0, 0.5); // 1.2974425414002562
let c = hyp2f1(1.0, 2.0, 3.0, 0.5); // 1.545177444479562
let d = hyperu(1.0, 1.0, 1.0); // 0.5963473623231941
let e = wright_bessel(1.0, 2.0, 1.0); // 1.5906368546373284
let f = log_wright_bessel(1.0, 2.0, 1.0); // 0.46413447354615933
let w: f64 = wrightomega(0.0); // 0.5671432904097838
```

`hyp0f1` and `wrightomega` accept both `f32` and `f64` through `num_traits::Float`; they evaluate in `f64` and return the input type. The remaining functions use `f64`. All functions operate on real arguments; complex continuation is not supported. Invalid inputs and unresolved numerical regions may return `NaN`, while poles or overflow may produce infinities. The crate uses `bessel` for large-argument 0F1 continuation. 

For the full accuracy report, see [sci-testing](https://github.com/p-sira/sci-testing/blob/main/docs/accuracy/special.md#scispecialhypergeometric), since [sci-special](https://github.com/p-sira/sci/tree/main/crates/sci-special) re-exports related implementations from `hypergeometric`.
