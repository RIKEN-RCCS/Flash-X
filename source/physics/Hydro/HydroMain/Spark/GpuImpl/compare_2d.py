#!/usr/bin/env python3
"""Compare 2D profiles across decompositions, rotated axes, or a planar 1D reference."""
import argparse
import csv
import math
FIELDS=['x','y','density','velocity_x','velocity_y','pressure','internal_energy','total_energy']

def read(path,fields=FIELDS):
    with open(path) as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=fields:
            raise ValueError('Unexpected CSV columns: '+str(reader.fieldnames))
        rows=[[float(row[k]) for k in fields] for row in reader]
    if not rows or any(not all(math.isfinite(v) for v in row) for row in rows):
        raise ValueError('Empty or nonfinite profile')
    return rows

def compare(reference,candidate,tolerance=1e-9,transpose=False,axis=None):
    rows=read(candidate)
    if any(r[2]<=0 or r[5]<=0 or r[6]<=0 for r in rows):
        raise ValueError('Nonpositive candidate state')
    error=0.
    if axis:
        one=read(reference,['x','density','velocity','pressure','internal_energy','total_energy'])
        profile={round(r[0],12):r for r in one}
        for r in rows:
            position=r[0 if axis=='x' else 1]
            if round(position,12) not in profile:
                raise ValueError('Reference resolution does not match planar direction')
            q=profile[round(position,12)]
            values=[r[2],r[3 if axis=='x' else 4],r[5],r[6],r[7]]
            expected=q[1:]
            error=max(error,abs(r[4 if axis=='x' else 3]),*(abs(a-b) for a,b in zip(values,expected)))
    else:
        ref=read(reference)
        if transpose:
            rows=[[r[1],r[0],r[2],r[4],r[3],*r[5:]] for r in rows]
        ref.sort(key=lambda r:(r[1],r[0])); rows.sort(key=lambda r:(r[1],r[0]))
        if len(ref)!=len(rows):
            raise ValueError('Different profile sizes')
        error=max(abs(a-b) for rr,cc in zip(ref,rows) for a,b in zip(rr,cc))
    if error>tolerance:
        raise ValueError(f'Profile difference {error} exceeds {tolerance}')
    return error

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference'); p.add_argument('candidate')
    p.add_argument('--tolerance',type=float,default=1e-9)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--transpose',action='store_true')
    group.add_argument('--reference-1d',choices=['x','y'])
    a=p.parse_args()
    print('maximum absolute profile difference:',compare(a.reference,a.candidate,a.tolerance,a.transpose,a.reference_1d))
