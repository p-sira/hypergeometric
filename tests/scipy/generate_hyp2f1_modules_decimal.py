"""Generate exact binary64 references for Sci's 1,000-call 2F1 benchmark.

Run from the repository root:
    python3 tests/scipy/generate_hyp2f1_modules_decimal.py

Only Python's standard library is required. References use the defining
series at 120 and 180 decimal digits, with matching binary64 rounding.
"""

import csv
from pathlib import Path

from generate_hyp2f1_benchmark_decimal import A, C, MASK, hyp2f1


def benchmark_arguments():
    state = (9 * A + C) & MASK
    for _ in range(1000):
        state = (state * A + C) & MASK
        v = ((state >> 11) / float(1 << 53)) * 2.0 - 1.0
        yield 0.9 * v


def main() -> None:
    destination = Path(__file__).parents[1] / "fixtures/hyp2f1_modules_exact.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["z", "expected"])
        for z in benchmark_arguments():
            low = hyp2f1(0.5, 1.5, 2.5, z, 120)
            high = hyp2f1(0.5, 1.5, 2.5, z, 180)
            assert float(low) == float(high), z
            writer.writerow([repr(z), format(high, ".55g")])


if __name__ == "__main__":
    main()
