#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1"]
# ///
"""Rebuild independent SciPy expected results from the committed input CSV.

Requires Python 3.13, scipy==1.18.1, numpy==2.5.4. The computation imports
only SciPy and NumPy; it neither imports nor executes the Rust crate.
"""
import csv
from decimal import Decimal,localcontext
from pathlib import Path
import scipy
import numpy
from scipy import special, integrate
import math

assert scipy.__version__ == "1.18.1", scipy.__version__
assert numpy.__version__ == "2.5.4", numpy.__version__
root = Path(__file__).resolve().parents[1] / "fixtures"
def hyperu_quadrature(a, b, z):
    """Independent DLMF 13.4.4 integral, smoothed with elementary substitutions."""
    if a <= 1.0:
        # t=u/(1-u), v=u**a removes t**(a-1) at the lower endpoint.
        def f(v):
            if v >= 1.0:
                return 0.0
            u = v**(1.0/a)
            return math.exp(-z*u/(1.0-u)) * ((1.0-u)**(-b))
        return integrate.quad(f,0.0,1.0,epsabs=1e-12,epsrel=1e-12,limit=500)[0]/special.gamma(a+1.0)
    # u=v**p makes the lower endpoint smooth when a>1.
    p=max(1,math.ceil(3.0/a))
    def f(v):
        if v >= 1.0:
            return 0.0
        u=v**p
        return p*v**(p*a-1.0)*math.exp(-z*u/(1.0-u))*(1.0-u)**(-b)
    return integrate.quad(f,0.0,1.0,epsabs=1e-12,epsrel=1e-12,limit=500)[0]/special.gamma(a)

def hyperu_reference(row):
    a,b,z = (float(row[k]) for k in ("a", "b", "z"))
    # The original SciPy U kernel returns NaN/wildly wrong values near b=0;
    # quadrature of its defining positive-a integral supplies an independent oracle.
    if b == 0.0 and a < 0.0 and z>0.0:
        # Recurrence DLMF 13.3.7 shifts negative a to the integral domain.
        # U(a,b,z)=(z+2a+2-b) U(a+1,b,z)
        #          -(a+1)(a+2-b) U(a+2,b,z).
        n=math.floor(-a)+1
        f1=hyperu_quadrature(a+n,b,z)
        f2=hyperu_quadrature(a+n+1,b,z)
        for k in range(n-1,-1,-1):
            t=a+k
            f0=(z+2.0*t+2.0-b)*f1-(t+1.0)*(t+2.0-b)*f2
            f2,f1=f1,f0
        return f1
    if a>0 and z>0 and (
        (b == 0.0 and a<=4.0 and z<=20.0)
        or (a==12.0 and b==10.0 and z==0.01)
        or (a==1.01 and b==2.0 and z==0.01)
    ):
        return hyperu_quadrature(a,b,z)
    return special.hyperu(a,b,z)

def hyp0f1_reference(b,z):
    if b >= 50.0 and abs(z) <= 100.0:
        with localcontext() as context:
            context.prec=100
            bb,zz=Decimal.from_float(b), Decimal.from_float(z)
            term=total=Decimal(1)
            for k in range(1,10000):
                term*=zz/((bb+k-1)*k)
                total+=term
                if abs(term)<=abs(total)*Decimal('1e-90'):
                    return float(total)
            raise RuntimeError('0f1 reference work limit')
    return special.hyp0f1(b,z)

functions = {
    "hyp0f1": lambda r: hyp0f1_reference(float(r["b"]), float(r["z"])),
    "hyp1f1": lambda r: special.hyp1f1(float(r["a"]), float(r["b"]), float(r["z"])),
    "hyp2f1": lambda r: special.hyp2f1(float(r["a"]), float(r["b"]), float(r["c"]), float(r["z"])),
    "hyperu": hyperu_reference,
    "wright_bessel": lambda r: special.wright_bessel(float(r["a"]), float(r["b"]), float(r["z"])),
    "log_wright_bessel": lambda r: special.log_wright_bessel(float(r["a"]), float(r["b"]), float(r["z"])),
    "wrightomega": lambda r: special.wrightomega(float(r["z"])),
}
import sys
suffix='catalog_' if len(sys.argv)>1 and sys.argv[1]=='catalog' else ''
with (root / (suffix+"inputs.csv")).open(newline="") as fin, (root / (suffix+"expected.csv")).open("w", newline="") as fout:
    reader = csv.DictReader(fin)
    writer = csv.writer(fout)
    writer.writerow(["id", "value"])
    for row in reader:
        value = functions[row["func"]](row)
        writer.writerow([row["id"], repr(float(value))])
