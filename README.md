# hypergeometric

A standalone, pure Rust implementation of real hypergeometric and Wright functions extracted from [Sci](https://github.com/p-sira/sci). The crate has no runtime dependencies, FFI, or dependency on Sci's facade or other extracted crates.

```rust
use hypergeometric::{hyp0f1, hyp1f1, hyp2f1, hyperu, wright_bessel, log_wright_bessel, wrightomega};
let m = hyp1f1(1.5, 2.5, 0.5);
let u = hyperu(1.5, 2.5, 0.5);
assert!(m.is_finite() && u.is_finite());
```

## API and mathematical conventions

| Function | Mathematical definition | Main real domain |
| --- | --- | --- |
| `hyp0f1(b, z)` | $\sum_{k\ge 0} z^k/((b)_k k!)$ | Finite b,z, b ≠ 0,-1,…; bounded compensated series |
| `hyp1f1(a,b,z)` | $\sum_{k\ge 0}(a)_k z^k/((b)_k k!)$ | Finite a,b,z and bounded |a|, with terminating-polynomial exceptions at b poles |
| `hyp2f1(a,b,c,z)` | $\sum_{k\ge 0}(a)_k(b)_k z^k/((c)_k k!)$ | Real analytic continuation for z≤1, with terminating exceptions |
| `hyperu(a,b,z)` | $\Gamma(a)^{-1}\int_0^\infty e^{-zt}t^{a-1}(1+t)^{b-a-1}dt$ (a,z>0) | z≥0, real a,b |
| `wright_bessel(a,b,z)` | $\sum_{k\ge0} z^k/[k!\Gamma(ak+b)]$ | a,b,z≥0 |
| `log_wright_bessel(a,b,z)` | $\log\Phi(a,b;z)$ | a,b,z≥0 |
| `wrightomega(x)` | Real w≥0 satisfying $w+\log w=x$ | All real x, including infinities |

All public APIs return `f64`; invalid arguments yield NaN, poles may yield +∞, and overflowing values may yield infinities. Their implementations bound iteration, and inputs with large cancellation can return NaN instead of an unreliable finite value. `hyp0f1` uses exact half-integer identities with a compensated square-root phase near trigonometric zeros, and a double-double series on moderate arguments; its previous Sci-specific Bessel fallback has not been transferred (for |z|>100, support depends on the relation between b and z). `hyp1f1`, `hyp2f1`, and `hyperu` use established historical series and continuation kernels with known cancellation limitations; `hyperu(a,0,z)` uses independent adaptive quadrature for 0<a≤4 and 0<z≤20; in some `hyperu` integral regions about nine valid decimal digits are expected. Complex branches are not supported. For z>1, the Gauss function is real-valued only for selected configurations; the principal complex continuation is not implemented. Maximum parameter magnitude for general functions is 100000.

## Development and verification

`tests/fixtures/inputs.csv` and `tests/fixtures/expected.csv` are independent: expected values are produced by **SciPy 1.18.1 / NumPy 2.5.4** (with independent `scipy.integrate.quad` and DLMF contiguous recurrences for U parameters where SciPy's direct function is inaccurate or undefined) and never computed from Rust. To rebuild:

```sh
python3 tests/scipy/make_inputs.py
uv run --script tests/scipy/generate_hypergeometric.py
python3 tests/scipy/make_catalog_inputs.py
uv run --script tests/scipy/generate_hypergeometric.py catalog
python3 tests/scipy/generate_hyp0f1_precision.py
uv run --script tests/scipy/generate_hyp1f1_cancellation.py
python3 tests/scipy/generate_hyp1f1_benchmark_decimal.py
cargo test
cargo clippy --all-targets -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps
```

Independent identity tests cover recurrence/normalization and boundary behavior. The fixture coverage is a finite sampled domain, not a proof of accuracy across every mathematically valid parameter range. The original 1,118-case Sci hypergeometric/Wright catalog is copied into `tests/scipy/sci_special_cases_hypergeometric.jsonl`; its separate inputs/expected CSVs are generated with `make_catalog_inputs.py` and `generate_hypergeometric.py catalog`. The checked-in reference suites also include 1,721 grid cases, 239 independently precision-checked 1F1 cancellation cases, 256 independent 120-digit Decimal benchmark cases, and 40 independently generated 105-digit Decimal 0F1 references (including adjacent binary64 inputs around a cosine zero). The 0F1 reference generator uses `Decimal.from_float` to preserve exact input rounding. Relative tests have no fixed absolute floor for nonzero values; the 0F1 contiguous-recurrence test scales its error by the magnitudes of the terms that may cancel. The test data and generator scripts remain in Git but are omitted from the published crate package.

SciPy 1.18.1 `special.hyperu(0.1,0,0.05)` differs from independent quadrature by about 8.6e-7; the generator takes the smooth transformed quadrature as its reference on that affected parameter slice. Regression references for known `hyp1f1` cancellation may additionally be generated using Sci's original high-precision Decimal scripts.

## Provenance and licenses

- `src/specfun.rs`: isolated real Zhang–Jin kernels, translated and adapted from SciPy `xsf/cpu/specfun` (BSD attribution and original Zhang–Jin notice in `LICENSES/scipy-specfun.txt`, `LICENSES/xsf-BSD-3-Clause.txt`).
- `src/kummer.rs`, `src/one_f_one.rs`, `src/zero_f_one.rs`, `src/dd.rs`: original Sci Rust kernels; Sci's BSD license in `LICENSE`.
- `src/lib.rs`: Wright Bessel positive-series and Wright omega real algorithms from Sci; Wright omega descends from BSD Algorithm 917 by Piers Lawrence (`LICENSES/wright.txt`).
- `src/gamma.rs`: self-contained Lanczos/reflection and digamma helper, using standard identities from [DLMF 5](https://dlmf.nist.gov/5).

BSD notices and attribution are supplied for redistributors. The extracted crate can be built and tested independently from Sci. It has not been published to crates.io.
