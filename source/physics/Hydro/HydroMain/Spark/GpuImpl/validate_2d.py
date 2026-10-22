#!/usr/bin/env python3
"""Run 2D GPU validation on RIKYU. Does not build; requires an allocated GPU."""
import argparse
import os
from pathlib import Path
import subprocess
from compare_2d import compare, read

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,default=Path('build2d/periodic_batch_opt'))
    p.add_argument('--output',type=Path,default=Path('validation2d'))
    p.add_argument('--reference-1d',type=Path,help='64-cell SSPRK3 NewImpl result at t=0.02')
    a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy(); env['OMP_TARGET_OFFLOAD']='MANDATORY'; env.pop('GPU_ALLOW_HOST',None)
    exe=str(a.exe.resolve())
    def run(gx,gy,bx,by,test):
        path=a.output/f'{gx}x{gy}-{bx}x{by}-case{test}.csv'
        r=subprocess.run([exe,str(gx),str(gy),str(bx),str(by),'0.02',str(path),str(test)],
                         env=env,capture_output=True,text=True)
        path.with_suffix('.log').write_text(r.stdout+r.stderr)
        if r.returncode:
            raise RuntimeError('Solver failure; inspect '+str(path.with_suffix('.log')))
        rows=read(path)
        if len(rows)!=gx*gy or len({(r[0],r[1]) for r in rows})!=gx*gy:
            raise RuntimeError('Incomplete/duplicate 2D output: '+str(path))
        return path
    baseline={}
    for case in (1,2,3):
        baseline[case]=run(64,64,1,1,case)
        for bx,by in ((2,2),(4,1),(1,4)):
            result=run(64,64,bx,by,case)
            print(case,bx,by,compare(baseline[case],result),flush=True)
    print('x/y rotation:',compare(baseline[1],baseline[2],transpose=True))
    if a.reference_1d:
        for case,axis in ((1,'x'),(2,'y')):
            print('NewImpl planar reference',axis,compare(a.reference_1d,baseline[case],axis=axis))
    # Rectangular cells and asymmetric block grid exercise separate spacings and strides.
    rectangular=run(96,64,1,1,3)
    print('rectangular decomposition:',compare(rectangular,run(96,64,3,2,3)))
    print('2D validation passed; results in',a.output)

if __name__=='__main__':
    main()
