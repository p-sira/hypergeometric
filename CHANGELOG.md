# Changelog

## 0.1.2

### Performance

- Speed up `hyp2f1` continuation for integer parameter differences and noninteger cases near 1 ([#5](https://github.com/p-sira/hypergeometric/pull/5)).
- Reduce series evaluation costs in `hyp0f1` for large positive arguments, `hyp1f1` for positive and negative arguments, and `hyperu` for integer parameters ([#5](https://github.com/p-sira/hypergeometric/pull/5)).
- Speed up `log_wright_bessel` when `a = 1` using a log-gamma recurrence ([#5](https://github.com/p-sira/hypergeometric/pull/5)).

## 0.1.1

### Performance

- Speed up low-order series in `wright_bessel` and `log_wright_bessel` by caching logarithmic factorials and reducing exponential evaluations ([#4](https://github.com/p-sira/hypergeometric/pull/4)).

## 0.1.0

### Added

- Initial release of real 0F1, 1F1, 2F1, Tricomi U, Wright generalized Bessel, and Wright omega functions, extracted from Sci.
