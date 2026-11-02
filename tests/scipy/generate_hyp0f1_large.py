#!/usr/bin/env python3
"""Independent mpmath 1.3.0 references for large, negative and near-pole 0F1.

mpmath calculates the hypergeometric function at 180 digits. Inputs are first
rounded to binary64 (then reconstructed exactly) so near-pole comparisons
measure the same mathematical function as the Rust API. Requires mpmath 1.3.0.
"""
import csv
import math
from pathlib import Path
import mpmath as mp

assert mp.__version__ == "1.3.0", "Fixture generator requires mpmath 1.3.0"
mp.mp.dps = 180
root = Path(__file__).resolve().parents[1]
regular_b = [0.125, 0.25, 0.75, 1.25, 2.5, 4.0, 10.0, 50.0, 100.0, 250.0,
             500.0, 1000.0, 100000.0, -0.125, -0.3, -2.3, -20.125,
             -1 + 2**-40, -1 - 2**-40, -3 + 2**-40, -3 - 2**-40,
             1e-8, -1e-8]
regular_z = [-100000.0, -10000.0, -1000.0, -200.0,
             200.0, 1000.0, 10000.0, 100000.0]
pairs = [(b, z) for b in regular_b for z in regular_z]
pairs += [(math.nextafter(pole, side), z)
          for pole in (-1.0, -3.0, -100.0)
          for side in (-math.inf, math.inf)
          for z in (-100000.0, -1000.0, 1000.0, 100000.0)]
pairs += [(100000.0, -1e7), (1000.0, -1e5), (1000.0, -5e4),
          (500.0, -1e4), (100000.0, -1e6), (500.0, -1.5e4),
          (1000.0, -1.2e5), (250.0, -10000.0)]
with (root / 'scipy/hyp0f1_large_inputs.csv').open('w', newline='') as inputs, \
        (root / 'fixtures/hyp0f1_large_expected.csv').open('w', newline='') as refs:
    iw = csv.writer(inputs, lineterminator='\n')
    rw = csv.writer(refs, lineterminator='\n')
    iw.writerow(['b', 'z'])
    rw.writerow(['value'])
    for b, z in pairs:
        value = mp.hyper([], [mp.mpf(b)], mp.mpf(z))
        if abs(value) > mp.mpf('1.0e308') or (value and abs(value) < mp.mpf('1e-300')):
            continue
        assert mp.isfinite(value)
        iw.writerow([repr(b), repr(z)])
        rw.writerow([repr(float(value))])
