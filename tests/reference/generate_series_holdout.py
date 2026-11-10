"""Reproduce independent Decimal references outside the Sci benchmark grids.

The test inputs span both signs of the numerator parameters and arguments,
including cancellation-sensitive and slowly converging 2F1 series.
All arithmetic starts from the exact binary64 input values.
"""

import csv
from decimal import Decimal, localcontext
from pathlib import Path
from random import Random


def confluent(a: float, b: float, z: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = 140
        aa, bb, zz = map(Decimal.from_float, (a, b, z))
        term = total = Decimal(1)
        for k in range(1, 1200):
            term = term * (aa + k - 1) * zz / ((bb + k - 1) * k)
            total += term
            if k > 50 and abs(term) < Decimal("1e-115") * max(
                abs(total), Decimal("1e-70")
            ):
                return total
    raise ArithmeticError(f"1F1 series did not converge for {(a, b, z)}")


def gauss(a: float, b: float, c: float, z: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = 140
        aa, bb, cc, zz = map(Decimal.from_float, (a, b, c, z))
        term = total = Decimal(1)
        for k in range(1, 9000):
            term = term * (aa + k - 1) * (bb + k - 1) * zz / (
                (cc + k - 1) * k
            )
            total += term
            if k > 70 and abs(term) < Decimal("1e-110") * max(
                abs(total), Decimal("1e-70")
            ):
                return total
    raise ArithmeticError(f"2F1 series did not converge for {(a, b, c, z)}")


def inputs():
    rng = Random(20261112)
    for _ in range(160):
        a = rng.uniform(-8.0, 8.0)
        b = rng.uniform(0.6, 9.0)
        z = rng.uniform(-8.0, 8.0)
        yield "hyp1f1", (a, b, z)
    for _ in range(160):
        a = rng.uniform(-5.0, 5.0)
        b = rng.uniform(-5.0, 5.0)
        c = rng.uniform(0.6, 9.0)
        z = rng.uniform(-0.96, 0.96)
        yield "hyp2f1", (a, b, c, z)


def main():
    path = Path(__file__).parents[1] / "fixtures/series_holdout.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        for name, args in inputs():
            value = confluent(*args) if name == "hyp1f1" else gauss(*args)
            writer.writerow([name, *map(repr, args), format(value, ".55g")])


if __name__ == "__main__":
    main()
