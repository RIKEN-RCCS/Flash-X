#!/usr/bin/env python3
"""Generate a 1D batched OpenMP-target baseline; never compile or execute it."""
from pathlib import Path
import argparse
import importlib.util
import re
import sys
ROOT = Path(__file__).resolve().parent
REF = ROOT.parent / 'NewImpl'
spec = importlib.util.spec_from_file_location('gpu_macro_expand', REF/'macro_expand.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)

def generate(output):
    output.mkdir(parents=True, exist_ok=True)
    definitions = [REF/'hydro_layout.ini', REF/'hydro_helpers.ini',
                   REF/'update_helpers.ini', REF/'diagnostic_helpers.ini', REF/'driver_helpers.ini', ROOT/'gpu_layout.ini']
    exp = m.MacroExpander(m.load_macros(definitions))
    source = (ROOT/'batch_kernels.F90-mc').read_text()
    # Reuse the reference face algebra verbatim, excluding traversal/reduction/reporting.
    face = (REF/'hy_getFaceFlux.F90-mc').read_text()
    face = face.split('            @hy_assembleFaceStates', 1)[1]
    face = '@hy_assembleFaceStates' + face.split('\n          end do', 1)[0]
    face = face.replace('work%minus','minus').replace('work%plus','plus')
    face = face.replace('data,','state,').replace('maxSpeed=max(maxSpeed,speed)', '')
    face = face.replace('return', 'exit cell_work')
    source = source.replace('! INSERT_REFERENCE_FACE_BODY', face)
    source = m.lower_blocks(exp.expand_text(source))
    # Hoisted macro temporaries must be private to each collapsed iteration.
    def private(match):
        declarations = match.group(1)
        names=[]
        for line in declarations.splitlines():
            if '::' not in line or line.lstrip().startswith('!'):
                continue
            for item in m.split_arguments(line.split('::',1)[1]):
                names.append(re.match(r'\s*(\w+)', item).group(1))
        groups=[names[i:i+3] for i in range(0,len(names),3)]
        clauses=['!$omp& private(' + ','.join(group) + ') &' for group in groups]
        return declarations + '\n    error=0\n!$omp target teams distribute parallel do collapse(2) &\n' + \
            '\n'.join(clauses) + '\n!$omp& reduction(max:error)'
    source = re.sub(r'! PRIVATE_BEGIN\n(.*?)! PRIVATE_END', private, source, flags=re.S)
    # Shared small maps are initialized on the host before the face launch.
    start=source.index('  subroutine faces')
    end=source.index('  end subroutine',start)
    part=source[start:end].replace('    error=0', '    inputMap=[1,2,3,4,5,6,7]; outputMap=[1,2,3,4,5]\n    shockVariable=9; d=1; e=[1,0,0]\n    error=0')
    source=source[:start]+part+source[end:]
    resident={'flatten':'state,flat','reconstruct':'state,flat,minus,plus',
              'faces':'state,minus,plus,flux','update':'reference,state,flux,integrated',
              'acoustics':'state,soundWork','markers':'state,soundWork','timestep_rate':'state'}
    for name,arrays in resident.items():
        start=source.index('  subroutine '+name+'(')
        end=source.index('  end subroutine',start)
        part=source[start:end].replace('collapse(2)', 'collapse(2) map(alloc:'+arrays+')')
        source=source[:start]+part+source[end:]
    start=source.index('  subroutine timestep_rate(')
    part=source[start:].replace('    error=0', '    error=0; rate=0._rk').replace('reduction(max:error)', 'reduction(max:error,rate)')
    source=source[:start]+part
    (output/'batch_kernels.f90').write_text(m.wrap_fortran(source))
    # Only type declarations/constants from the reference module are needed.
    types=(REF/'hy_getFaceFlux.F90-mc').read_text().split('\ncontains',1)[0]
    types=types.replace('public :: cell_field, flux_field, hy_getFaceFlux', 'public :: cell_field, flux_field')
    types=m.MacroExpander(m.load_macros([REF/'hydro_layout.ini'])).expand_text(types)
    (output/'hydro_types.f90').write_text(types+'\nend module\n')
    for name in ['periodic_batch.F90','Makefile']:
        (output/name).write_text(m.wrap_fortran((ROOT/name).read_text()) if name.endswith('.F90') else (ROOT/name).read_text())
    print(output)

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'build')
    generate(parser.parse_args().output.resolve())
