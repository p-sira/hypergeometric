"""Generate 90-digit independent Wright references for Sci's CodSpeed inputs."""
from pathlib import Path
import mpmath as mp
mp.mp.dps = 90
mask = (1 << 64) - 1

def params(n, seed, lo, hi):
    state = (seed * 6364136223846793005 + 1442695040888963407) & mask
    for _ in range(n):
        state = (state * 6364136223846793005 + 1442695040888963407) & mask
        v = (((state >> 11) / 2**53) * 2 - 1)
        u = 0.5 * min(max(v + 1, 1e-6), 2 - 1e-6)
        yield lo + (hi - lo) * u

def wright(a, b, z):
    a, b, z = map(mp.mpf, (a,b,z))
    if a == 0: return mp.exp(z) * mp.rgamma(b)
    if z == 0: return mp.rgamma(b)
    total = mp.mpf(0)
    for k in range(2000):
        term = z**k * mp.rgamma(a*k+b) / mp.factorial(k)
        total += term
        if k > 30 and term <= mp.mpf('1e-87') * total:
            return total
    raise RuntimeError('reference failed to converge')

out = Path(__file__).resolve().parents[1] / 'fixtures' / 'wright_bench.csv'
with out.open('w') as f:
    f.write('a,b,z,log_value\n')
    for seeds in [(82,83,84),(85,86,87)]:
        arrs = [list(params(32,seeds[0],0.0,1.0)),
                list(params(32,seeds[1],0.0,2.0)),
                list(params(32,seeds[2],0.0,10.0))]
        for a,b,z in zip(*arrs):
            log_ref = mp.log(wright(a,b,z))
            f.write(f'{a:.17g},{b:.17g},{z:.17g},{mp.nstr(log_ref,31)}\n')
print(out, '64 inputs')
