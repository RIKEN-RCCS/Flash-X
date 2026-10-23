#!/usr/bin/env python3
"""Compare one post-Hydro snapshot across rank-local CSV files; no HDF5 required."""
import argparse
import csv
import glob
import json
import math
import re
from pathlib import Path

FIELDS=['step','time','dt','level','i','j','k','x','y','z','dx','dy','dz',
        'density','velocity_x','velocity_y','velocity_z','pressure','internal_energy','total_energy']
QUANTITIES=FIELDS[13:]

def parse_real(value):
    """Read E/D notation and legacy Fortran fields with an omitted exponent letter."""
    value=value.strip().replace('D','E').replace('d','e')
    # ES formats without an explicit exponent width may write 1.234-146.
    match=re.fullmatch(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d{3,})',value)
    if match:
        value=match[1]+'E'+match[2]
    return float(value)

def load(patterns,digits=12):
    files=sorted({name for pattern in patterns for name in glob.glob(str(pattern))})
    if not files:
        raise ValueError('No dump files matched: '+str(patterns))
    cells={}; metadata=None
    for name in files:
        with open(name) as stream:
            reader=csv.DictReader(stream,skipinitialspace=True)
            if reader.fieldnames!=FIELDS:
                raise ValueError('Unexpected header: '+name)
            for record in reader:
                row={key:parse_real(value) for key,value in record.items()}
                if not all(math.isfinite(v) for v in row.values()):
                    raise ValueError('Nonfinite value: '+name)
                for key in ['step','level','i','j','k']:
                    if not row[key].is_integer():
                        raise ValueError('Noninteger index: '+name)
                    row[key]=int(row[key])
                if any(row[q]<=0 for q in ['dx','dy','dz','density']):
                    raise ValueError('Nonpositive cell volume or density: '+name)
                signature=(row['step'],row['time'],row['dt'])
                if metadata is None:
                    metadata=signature
                if signature!=metadata:
                    raise ValueError('Mixed steps/times/timesteps: select files for exactly one snapshot')
                key=(row['level'],*(round(row[q],digits) for q in ['x','y','z']))
                if key in cells:
                    raise ValueError('Duplicate/rounded-colliding cell key: '+str(key))
                cells[key]=row
    if not cells:
        raise ValueError('Snapshot has no owned cells')
    return cells,metadata,files

def totals(cells):
    rows=[cells[k] for k in sorted(cells)]
    masses=[r['density']*r['dx']*r['dy']*r['dz'] for r in rows]
    result={'mass':math.fsum(masses)}
    for key in ['velocity_x','velocity_y','velocity_z','total_energy']:
        name={'velocity_x':'momentum_x','velocity_y':'momentum_y','velocity_z':'momentum_z',
              'total_energy':'energy'}[key]
        result[name]=math.fsum(m*r[key] for m,r in zip(masses,rows))
    return result

def compare(reference,candidate,atol=1e-10,rtol=1e-8,time_atol=1e-12,digits=12):
    a,ma,fa=load(reference,digits); b,mb,fb=load(candidate,digits)
    if ma[0]!=mb[0] or abs(ma[1]-mb[1])>time_atol or abs(ma[2]-mb[2])>time_atol:
        raise ValueError(f'Snapshot step/time/dt differ: {ma} versus {mb}')
    if a.keys()!=b.keys():
        raise ValueError(f'Different mesh/cell keys: {len(a.keys()-b.keys())} reference-only, '
                         f'{len(b.keys()-a.keys())} candidate-only')
    maxima={q:0. for q in QUANTITIES}; bad={q:0 for q in QUANTITIES}
    for key in a:
        aa,bb=a[key],b[key]
        for q in ['dx','dy','dz']:
            if not math.isclose(aa[q],bb[q],rel_tol=1e-12,abs_tol=1e-14):
                raise ValueError('Different mesh widths at '+str(key))
        for q in QUANTITIES:
            difference=abs(aa[q]-bb[q]); maxima[q]=max(maxima[q],difference)
            if difference>atol+rtol*max(abs(aa[q]),abs(bb[q])):
                bad[q]+=1
    ta,tb=totals(a),totals(b)
    return {'passed':not any(bad.values()),'cells':len(a),'step':ma[0],'reference_time':ma[1],
            'candidate_time':mb[1],'reference_files':fa,'candidate_files':fb,
            'max_absolute_difference':maxima,'cells_outside_tolerance':bad,
            'reference_totals':ta,'candidate_totals':tb,
            'total_differences':{k:tb[k]-ta[k] for k in ta}}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',nargs='+',required=True,help='Filenames or quoted rank-file glob for ONE step')
    p.add_argument('--candidate',nargs='+',required=True)
    p.add_argument('--atol',type=float,default=1e-10)
    p.add_argument('--rtol',type=float,default=1e-8)
    p.add_argument('--time-atol',type=float,default=1e-12)
    p.add_argument('--coordinate-digits',type=int,default=12)
    p.add_argument('--report',type=Path)
    args=p.parse_args()
    if min(args.atol,args.rtol,args.time_atol)<0:
        p.error('Tolerances must be nonnegative')
    try:
        report=compare(args.reference,args.candidate,args.atol,args.rtol,args.time_atol,args.coordinate_digits)
    except ValueError as error:
        p.exit(1,str(error)+'\n')
    print(json.dumps(report,indent=2))
    if args.report:
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    if not report['passed']:
        p.exit(1,'Cell comparison failed\n')

if __name__=='__main__':
    main()
