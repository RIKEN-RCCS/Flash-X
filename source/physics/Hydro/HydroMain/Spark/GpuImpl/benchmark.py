#!/usr/bin/env python3
"""Run optimized solver repetitions on the allocated RIKYU GPU; never build it."""
import argparse
import csv
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import time

FLOAT=r'([-+0-9.EeDd]+)'
def parse(text):
    def number(label):
        match=re.search(re.escape(label)+r'\s*'+FLOAT,text)
        if not match:
            raise ValueError('Missing timing: '+label)
        return float(match[1].replace('D','E').replace('d','e'))
    match=re.search(r'steps,time=\s*(\d+)',text)
    if not match or not re.search(r'initial_device=\s*0\b',text):
        raise ValueError('Missing step count or device execution confirmation')
    return int(match[1]),number('resident stepping seconds='),number('including transfer seconds=')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,default=Path('build/periodic_batch_opt'))
    p.add_argument('--cells',type=int,nargs='+',default=[4096,16384,65536])
    p.add_argument('--blocks',type=int,nargs='+',default=[1,16,64])
    p.add_argument('--steps',type=int,default=100)
    p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',type=Path,default=Path('benchmark-results'))
    a=p.parse_args()
    if a.steps<1 or a.repeats<1 or any(n<8 for n in a.cells) or any(b<1 for b in a.blocks):
        p.error('Positive steps/repeats/blocks and at least 8 cells are required')
    if any(n%b for n in a.cells for b in a.blocks):
        p.error('Every cell count must be divisible by every block count')
    exe=a.exe.resolve()
    if not exe.is_file():
        p.error('Build the optimized executable first: '+str(exe))
    a.output.mkdir(parents=True,exist_ok=True)
    # Keep each campaign separate; preserve prior results and its metadata.
    run=a.output/time.strftime('%Y%m%d-%H%M%S',time.gmtime())
    run.mkdir(exist_ok=False)
    env=os.environ.copy()
    env['OMP_TARGET_OFFLOAD']='MANDATORY'
    env.pop('GPU_ALLOW_HOST',None)
    metadata={'executable':str(exe),'host':platform.node(),'platform':platform.platform(),
              'arguments':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
              'environment':{k:v for k,v in env.items() if k.startswith(('OMP_','CUDA_','NVCOMPILER_'))}}
    for label,cmd in [('gpu',['nvidia-smi']),('compiler',['nvfortran','--version'])]:
        try:
            result=subprocess.run(cmd,capture_output=True,text=True)
            metadata[label]=result.stdout+result.stderr
        except OSError as e:
            metadata[label]=str(e)
    (run/'metadata.json').write_text(json.dumps(metadata,indent=2))
    rows=[]
    for n in a.cells:
        for blocks in a.blocks:
            for repeat in range(a.repeats):
                prefix=run/f'n{n}-b{blocks}-r{repeat+1}'
                # Large time cap: fixed step count governs the benchmark.
                cmd=[str(exe),str(n),str(blocks),'1.0',str(prefix.with_suffix('.csv')),str(a.steps)]
                result=subprocess.run(cmd,env=env,capture_output=True,text=True)
                text=result.stdout+result.stderr
                prefix.with_suffix('.log').write_text(text)
                if result.returncode:
                    raise SystemExit('Solver failed; inspect '+str(prefix.with_suffix('.log')))
                steps,resident,total=parse(text)
                if steps!=a.steps or resident<=0 or total<=0:
                    raise SystemExit('Unexpected timings/step count; regenerate and rebuild with step-limit support')
                rows.append(dict(cells=n,blocks=blocks,cells_per_block=n//blocks,repeat=repeat+1,
                                 steps=steps,resident_seconds=resident,transfer_inclusive_seconds=total,
                                 seconds_per_step=resident/steps,cell_updates_per_second=n*steps/resident))
                print(f'{n} cells, {blocks} blocks, repeat {repeat+1}: {resident/steps:.6g} s/step',flush=True)
    with (run/'timings.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print('Median resident seconds per step:')
    for n in a.cells:
        for blocks in a.blocks:
            group=[r['seconds_per_step'] for r in rows if r['cells']==n and r['blocks']==blocks]
            print(n,blocks,statistics.median(group),' min=',min(group),' max=',max(group))
    print('Results:',run)

if __name__=='__main__':
    main()
