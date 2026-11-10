# Numerical accuracy audit

The audit compares the crate against mpmath 1.3.0 at the exact binary64
input values. Each reference rounds identically to binary64 when evaluated
at both 90 and 180 decimal digits. The deterministic grid uses seed
20261110; additional cases exercise zeros, near-pole parameters, terminating
series, endpoints, and overflow. The baseline is `main` at `48b630d`.

The regression suite contains 1,942 reference cases and additionally checks
`wright_bessel` against the exponential of the 100 logarithmic references.
Relative tolerance is `1e-11`, except for `hyperu`, whose legacy kernel has
roughly nine valid digits on some paths and is checked at `3e-8`. Exact zeros
and infinities must match exactly; there is no absolute-error floor.

| Function | Reference cases | Baseline failures | Fixed failures | Largest finite relative error after fixes |
| --- | ---: | ---: | ---: | ---: |
| `hyp0f1` | 315 | 6 | 0 | 2.50e-14 |
| `hyp1f1` | 607 | 146 | 0 | 2.59e-14 |
| `hyp2f1` | 613 | 100 | 0 | 3.23e-13 |
| `hyperu` | 207 | 5 | 0 | 5.58e-10 |
| `wrightomega` | 100 | 0 | 0 | 2.17e-16 |
| `log_wright_bessel` | 100 | 0 | 0 | 2.42e-14 |

Failures include wrong zeros, NaNs, and infinities as well as excessive
relative error. The largest finite error column excludes nonfinite results.

## Findings and fixes

- **Terminating `1F1`:** the defining sum can cancel over 40 decimal digits.
  For example, `M(-97,58.58075296254387,178.2752193528607)` previously
  returned about `-9.26e29` instead of `7.15e-11`. A normalized Laguerre
  recurrence in the existing double-double arithmetic avoids those terms;
  common binary scaling limits intermediate overflow. Mixed-sign
  nonterminating series now use double-double arithmetic, a geometric tail
  bound, and a cancellation check. Negative-argument continuation covers
  a wider range of negative lower parameters.
- **`2F1` near singularities:** tolerance-based equality incorrectly treated
  inputs immediately below `z=1` as the endpoint, and parameters adjacent
  to nonpositive integers as exact poles or terminating numerators. Exact
  comparisons now select those mathematical cases. The linear polynomial
  preserves tiny nonzero results, and equal negative numerator/denominator
  parameters have symmetric terminating behavior.
- **`2F1` continuation:** the connection at `1/z` can subtract large terms
  for ordinary negative arguments. Bounded negative arguments now use
  Pfaff's transformation. Mixed-sign and alternating series use
  double-double summation with a tail and cancellation check before the
  legacy recurrence is considered.
- **Elementary half-integer `0F1`:** rounding the square root corrupts
  phases near zeros, including the `b=3/2` sine identity. All three
  elementary identities share a residual-based phase correction and use
  angle addition. The `b=-1/2` positive-argument identity avoids infinity
  minus infinity on overflow.
- **`U` at zero:** divergence follows the sign of reciprocal gamma.
  Finite gamma ratios are computed in the logarithmic domain. Terminating
  endpoint polynomials retain exact zeros even when earlier factors
  overflow, and cover the documented numerator range.

The identities follow [DLMF §13.3](https://dlmf.nist.gov/13.3),
[Kummer's transformations](https://dlmf.nist.gov/13.2#vii),
[Gauss series](https://dlmf.nist.gov/15.2), and
[Pfaff's transformations](https://dlmf.nist.gov/15.8#E1).

## Reproduction and limits

```sh
python -m venv /tmp/hypergeometric-reference
/tmp/hypergeometric-reference/bin/pip install mpmath==1.3.0
/tmp/hypergeometric-reference/bin/python tests/reference/generate_accuracy_audit.py
cargo test --test numerical_audit
cargo test --release --test numerical_audit
```

The committed references require no Python dependency when running Rust tests.

This is sampled coverage, not an accuracy guarantee over all finite inputs.
Existing work limits, unsupported complex continuation, and unresolved extreme
regions still apply. Relative accuracy close to zeros remains sensitive to
conditioning. Large oscillatory `0F1` arguments whose square-root spacing
exceeds the phase-correction range still use a bounded trigonometric fallback.
Some `U` paths retain their documented lower precision. Double-double series
and recurrence paths may take more time than the former binary64 sums.
