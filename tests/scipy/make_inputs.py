#!/usr/bin/env python3
"""Create deterministic domain, boundary, and stress inputs independently of Rust."""
import csv
from pathlib import Path
from itertools import product
root = Path(__file__).resolve().parents[1] / "fixtures"
rows=[]
def add(fn, a='', b='', c='', z=''):
    rows.append([f'{len(rows):05d}', fn, a, b, c, z])
for b,z in product([.125,.5,1.,2.5,5.,10.,25.,-.3,-2.3], [-100.,-50.,-25.,-10.,-3.,-.5,0.,.5,2.,10.,50.,100.]):
    add('hyp0f1',b=b,z=z)
for a,b,z in product([-3.,-.5,0.,.5,1.5,3.,10.,65.1], [-3.2,.1,1.,2.5,4.,10.], [-40.,-10.,-2.,0.,2.,10.,40.]):
    add('hyp1f1',a=a,b=b,z=z)
for a,b,c,z in product([-2.,.5,1.,2.3,4.], [-1.,.5,1.,3.], [.5,1.,2.,3.5,5.], [-5.,-.9,-.25,0.,.25,.75,.99,1.]):
    add('hyp2f1',a=a,b=b,c=c,z=z)
for a,b,z in product([-2.,0.,.1,.5,1.,2.,4.], [0.,.5,1.,2.,3.,4.], [.05,.5,2.,10.,100.]):
    add('hyperu',a=a,b=b,z=z)
for a,b,z in product([0.,.1,.5,1.,2.], [0.,.1,.5,1.,3.], [0.,.5,2.,10.,20.]):
    add('wright_bessel',a=a,b=b,z=z)
    add('log_wright_bessel',a=a,b=b,z=z)
for z in [-100.,-50.,-20.,-5.,-1.,0.,1.,10.,100.,1000.]:
    add('wrightomega',z=z)
for fn,args in [
    ('hyp0f1',('',-1.,'',.5)),('hyp0f1',('',0.,'',1.)),
    ('hyp1f1',(1.,-1.,'',.5)),('hyp2f1',(1.,1.,0.,.5)),
    ('hyperu',(1.,1.,'',-1.)),('wright_bessel',(-1.,1.,'',2.)),
    ('log_wright_bessel',(1.,-1.,'',2.))]:
    add(fn,*args)
with (root/'inputs.csv').open('w',newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['id','func','a','b','c','z'])
    writer.writerows(rows)
print(len(rows),'input rows')
