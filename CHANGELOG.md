# Changelog

## 0.1

### 0.1.1

**Performance Improvements**
- Reduce integer-difference 2F1 continuation costs with digamma recurrence updates; stream integer-parameter Tricomi-U harmonic sums instead of rebuilding nested prefixes. Simplify large-positive 0F1 log-sum-exp accumulation while retaining overflow checks and high-precision reference tests.
- Skip the costly preliminary direct series for well-conditioned noninteger 2F1 continuation near 1; use the short connection series at 1-z. Accumulate `log_wright_bessel` with a log-gamma recurrence when a=1, stop positive-term 1F1 series with a bounded geometric tail, and streamline the double-double negative-argument Kummer recurrence.
- `wright_bessel`, `log_wright_bessel`: Cache 65 independently computed logarithmic factorials for the common low-order Wright series and avoid a redundant exponential in the log-sum-exp update; retain the bounded fallback for higher terms ([#4](https://github.com/p-sira/hypergeometric/pull/4)).

**Tests**
- Add 90-digit independent mpmath references for the Wright series arguments in Sci's Divan benchmark, including values and logarithms.
- Add high-precision 2F1 continuation, Wright recurrence, and negative-argument Kummer checks, analytic 1F1 identities, and targeted CodSpeed benchmarks for the newly optimized paths.
- Retain the existing domain, convergence and numerical-accuracy regressions.

### 0.1.0

**New Features**
- Initial standalone extraction of real 0F1, 1F1, 2F1, Tricomi U, Wright generalized Bessel, and Wright omega functions from Sci.
- Add independent numerical reference datasets and accuracy tests.
