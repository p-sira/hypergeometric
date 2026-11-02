#!/usr/bin/env python3
"""Independent 105-decimal-digit 0F1 series for exact binary64 arguments.

Each CSV decimal string is first rounded to its f64 input; Decimal.from_float
reconstructs that exact value. Particularly near zeros, evaluating the printed
string as an exact real gives a different function value. No Sci/SciPy or Rust
implementation is used in reference construction.
"""
import csv
from decimal import Decimal, localcontext
from pathlib import Path

root=Path(__file__).resolve().parents[1]
with (root/'scipy/hyp0f1_precision_inputs.csv').open(newline='') as inp, \
        (root/'fixtures/hyp0f1_precision_expected.csv').open('w',newline='') as out:
    reader=csv.DictReader(inp)
    writer=csv.writer(out,lineterminator="\n")
    writer.writerow(['value'])
    for row in reader:
        with localcontext() as ctx:
            ctx.prec=105
            b=Decimal.from_float(float(row['parameter']))
            z=Decimal.from_float(float(row['z']))
            term=total=Decimal(1)
            for k in range(1,10000):
                term*=z/((b+k-1)*k)
                total+=term
                if k>abs(z)+abs(b)+10 and abs(term)<Decimal('1e-85')*max(abs(total),Decimal('1e-200')):
                    break
            else:
                raise RuntimeError(('series did not converge',row))
            writer.writerow([repr(float(total))])
