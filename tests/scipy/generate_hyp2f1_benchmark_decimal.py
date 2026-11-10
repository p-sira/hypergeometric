"""Generate independent high-precision 2F1 references for Sci's benchmark grid.

Only the Python standard library is required. Decimal.from_float preserves
each exact binary64 parameter, and the defining 2F1 series converges for
the benchmark's |z| < 0.9. Values are checked at 120 and 180 digits.

Run: python tests/scipy/generate_hyp2f1_benchmark_decimal.py
"""

import csv
from decimal import Decimal, localcontext
from pathlib import Path


MASK = (1 << 64) - 1
A = 6364136223846793005
C = 1442695040888963407


def values(seed: int, lo: float, hi: float) -> list[float]:
    state = (seed * A + C) & MASK
    result = []
    for _ in range(256):
        state = (state * A + C) & MASK
        v = ((state >> 11) / float(1 << 53)) * 2.0 - 1.0
        u = 0.5 * max(1e-6, min(v + 1.0, 2.0 - 1e-6))
        result.append(lo + (hi - lo) * u)
    return result


def hyp2f1(a: float, b: float, c: float, z: float, precision: int) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = precision
        aa, bb, cc, zz = map(Decimal.from_float, (a, b, c, z))
        term = total = Decimal(1)
        tolerance = Decimal(10) ** (-(precision - 20))
        for k in range(1, 6000):
            term = term * (aa + k - 1) * (bb + k - 1) * zz / (
                (cc + k - 1) * k
            )
            total += term
            if k > 50 and abs(term) < tolerance * max(
                abs(total), Decimal("1e-80")
            ):
                return total
    raise ArithmeticError(f"2F1 series did not converge for {(a, b, c, z)}")


def main() -> None:
    destination = Path(__file__).parents[1] / "fixtures/hyp2f1_bench_exact.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["a", "b", "c", "z", "expected"])
        for a, b, c, z in zip(
            values(55, -2.0, 2.0),
            values(56, -2.0, 2.0),
            values(57, 0.5, 5.0),
            values(58, -0.9, 0.9),
        ):
            low = hyp2f1(a, b, c, z, 120)
            high = hyp2f1(a, b, c, z, 180)
            assert float(low) == float(high), (a, b, c, z)
            writer.writerow(
                [repr(a), repr(b), repr(c), repr(z), format(high, ".55g")]
            )


if __name__ == "__main__":
    main()
