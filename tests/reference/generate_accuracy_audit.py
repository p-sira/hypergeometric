"""Regenerate numerical_audit.csv with mpmath 1.3.0 (test-only dependency).

Run: python -m pip install mpmath==1.3.0
     python tests/reference/generate_accuracy_audit.py

References use exact binary64 inputs at both 90 and 180 decimal digits.
Complex values are excluded from the real API audit. No Rust output is used.
"""
import csv
import math
from pathlib import Path
import random

import mpmath as mp


def cases():
    rng = random.Random(20261110)
    for _ in range(300):
        b = rng.choice([-0.5, 0.5, 1.5, rng.uniform(-20, 100)])
        yield "hyp0f1", [b, rng.uniform(-500, 500)]
    for _ in range(600):
        a = rng.choice([-float(rng.randrange(1, 101)), rng.uniform(-20, 50)])
        yield "hyp1f1", [a, rng.uniform(-10, 60), rng.uniform(-200, 200)]
    for _ in range(600):
        a, b, c = rng.uniform(-10, 20), rng.uniform(-10, 20), rng.uniform(-10, 30)
        z = rng.choice([rng.uniform(-10, 0.9), 1 - 10 ** rng.uniform(-15, -2)])
        yield "hyp2f1", [a, b, c, z]
    for _ in range(200):
        yield "hyperu", [rng.uniform(-5, 10), rng.uniform(-5, 10), 10 ** rng.uniform(-4, 2)]
    for _ in range(100):
        yield "wrightomega", [rng.uniform(-745, 10000)]
    for _ in range(100):
        yield "log_wright_bessel", [10 ** rng.uniform(-2, 1), rng.uniform(0, 20), 10 ** rng.uniform(-3, 3)]
    for b in [-0.5, 0.5, 1.5]:
        for z in [-2.4674011002723395, -2.46740110027234, -1e28, 126025.0, 1e6]:
            yield "hyp0f1", [b, z]
    for a, b, z in [(-100.0, 40.0, 150.0), (-100.0, 40.0, -150.0),
                    (-2.0, 1e300, 1.0), (-100.0, 2.0, -10000.0),
                    (-1.0, 2.0, math.nextafter(2.0, 0.0)),
                    (-2.0, -2.0, -2.0), (-50.0, -50.0, -20.0)]:
        yield "hyp1f1", [a, b, z]
    for a, b, c, z in [
        (-1.0, 1.0, 1.0, math.nextafter(1.0, 0.0)),
        (-2.0, 1.0, -2.0, 0.5),
        (1.0, -2.0, -2.0, 0.5),
        (-2.0, -3.0, -2.0, 0.5),
        (math.nextafter(-1.0, 0.0), 2.0, 3.0, 0.99),
        (math.nextafter(-2.0, 0.0), 2.0, 3.0, 0.99),
        (0.5, 2.0, math.nextafter(-1.0, 0.0), 0.5),
        (0.5, 2.0, math.nextafter(-1.0, -math.inf), 0.5),
        (0.5, 2.0, 1e-14, 0.5),
        (0.5, 1.0, math.nextafter(1.0, 0.0), 0.9),
        (0.5, 0.5, 2.0, math.nextafter(1.0, 0.0)),
        (0.5, 0.5, 0.75, math.nextafter(1.0, 0.0)),
        (2.0, 3.0, 4.0, math.nextafter(1.0, 0.0)),
    ]:
        yield "hyp2f1", [a, b, c, z]
    for a, b in [(-0.5, 1.0), (-0.5, 2.0), (-1.5, 2.0),
                 (1000.0, -1000.0), (-2.0, -1000.0),
                 (-1000.0, -999.0), (-5000.0, -4999.0)]:
        yield "hyperu", [a, b, 0.0]


def wright_log(a, b, z):
    total = mp.mpf(0)
    for k in range(5000):
        term = z**k * mp.rgamma(a * k + b) / mp.factorial(k)
        total += term
        if k > 2 and term < mp.eps * total:
            return mp.log(total)
    raise ValueError("Wright reference did not converge")


def reference(name, args):
    args = list(map(mp.mpf, args))
    if name == "wrightomega":
        return mp.lambertw(mp.exp(args[0]))
    if name == "log_wright_bessel":
        return wright_log(*args)
    if name == "hyperu" and args[2] == 0:
        a, b, _ = args
        if a < 0 and mp.isint(a):
            return (-1)**int(-a) * mp.rf(b, int(-a))
        if b >= 1:
            return mp.sign(mp.rgamma(a)) * mp.inf
        return mp.gamma(1 - b) * mp.rgamma(1 + a - b)
    if name == "hyp2f1":
        a, b, c, z = args
        degrees = [int(-v) for v in [a, b] if v <= 0 and mp.isint(v)]
        if degrees and c <= 0 and mp.isint(c) and min(degrees) <= -c:
            return mp.fsum(mp.rf(a, k) * mp.rf(b, k) * z**k
                           / (mp.rf(c, k) * mp.factorial(k))
                           for k in range(min(degrees) + 1))
    return getattr(mp, name)(*args)


def main():
    rows = []
    excluded = 0
    for name, args in cases():
        try:
            with mp.workdps(90):
                low = float(reference(name, args))
            with mp.workdps(180):
                high = float(reference(name, args))
        except TypeError:  # mpc result outside the supported real domain
            excluded += 1
            continue
        assert low == high, (name, args, low, high)
        rows.append([name, *args, high])
    path = Path(__file__).resolve().parents[1] / "fixtures" / "numerical_audit.csv"
    with path.open("w", newline="") as output:
        csv.writer(output, lineterminator="\n").writerows(rows)
    print(f"Wrote {len(rows)} cases; excluded {excluded} complex results")


if __name__ == "__main__":
    main()
