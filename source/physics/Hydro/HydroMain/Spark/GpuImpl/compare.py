#!/usr/bin/env python3
"""Compare a RIKYU batch CSV to NewImpl SSPRK3 output at identical resolution/time."""
import argparse
import csv
import math
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('reference'); p.add_argument('candidate')
p.add_argument('--tolerance',type=float,default=1e-9)
a=p.parse_args()
def read(path):
    with open(path) as f:
        return [[float(v) for v in row.values()] for row in csv.DictReader(f)]
r=read(a.reference); c=read(a.candidate)
if not r or len(r)!=len(c):
    raise SystemExit('Different or empty profiles')
if any(len(row)!=6 or not all(math.isfinite(v) for v in row) for row in r+c):
    raise SystemExit('Malformed/nonfinite profile')
if any(row[1]<=0 or row[3]<=0 or row[4]<=0 for row in c):
    raise SystemExit('Nonpositive state')
err=max(abs(x-y) for rr,cc in zip(r,c) for x,y in zip(rr,cc))
print('maximum absolute profile difference:',err)
if err>a.tolerance:
    raise SystemExit('Profile comparison failed')
