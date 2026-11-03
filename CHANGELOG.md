# Changelog

## 0.1

### 0.1.1 (unreleased)

**Performance Improvements**
- `wright_bessel`, `log_wright_bessel`: Cache 65 independently computed logarithmic factorials for the common low-order Wright series and avoid a redundant exponential in the log-sum-exp update; retain the bounded fallback for higher terms ([#4](https://github.com/p-sira/hypergeometric/pull/4)).

**Tests**
- Add 90-digit independent mpmath references for the Wright series arguments in Sci's Divan benchmark, including values and logarithms.
- Retain the existing domain, convergence and numerical-accuracy regressions.

### 0.1.0

**New Features**
- Initial standalone extraction of real 0F1, 1F1, 2F1, Tricomi U, Wright generalized Bessel, and Wright omega functions from Sci.
- Add independent numerical reference datasets and accuracy tests.
