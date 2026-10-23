#!/usr/bin/env python3
"""Generate a 1D/2D batched OpenMP-target baseline; never compile or execute it."""
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

# Use the reference expander for hygiene and control flow, but omit argument
# descriptor association in device code. Actuals in this baseline are pure
# references/expressions: unlike ASSOCIATE, substitution does not evaluate once.
_reference_scoped_definition = m.scoped_definition

def device_scoped_definition(macro, actuals, serial):
    text = _reference_scoped_definition(macro, actuals, serial)
    if not macro.args:
        return text
    start = text.index("associate ( &\n")
    # The reference emitter puts one selector per line, including the last ')'.
    lines = text[start:].splitlines(keepends=True)
    text = text[:start] + ''.join(lines[1+len(macro.args):])
    text = text.replace('end associate\n', '')
    mapping = {}
    for name, actual in zip(macro.args, actuals):
        if actual.strip().startswith('(/'):
            raise m.ExpansionError('Device helpers require named array actuals')
        identifier = re.fullmatch(r'[A-Za-z_]\w*(?:%[A-Za-z_]\w*)*', actual.strip())
        mapping[f'macro_scope_{serial}_arg_{name}'.lower()] = actual if identifier else '('+actual+')'
    return m.rename_fortran(text, mapping, f'macro_scope_{serial}')

m.scoped_definition = device_scoped_definition

def generate(output, dimension=1):
    if dimension not in (1,2):
        raise ValueError("dimension must be 1 or 2")
    output.mkdir(parents=True, exist_ok=True)
    definitions = [REF/'hydro_layout.ini', REF/'hydro_helpers.ini',
                   REF/'update_helpers.ini', REF/'diagnostic_helpers.ini', REF/'driver_helpers.ini', ROOT/('gpu_layout_2d.ini-template' if dimension==2 else 'gpu_layout.ini-template')]
    exp = m.MacroExpander(m.load_macros(definitions))
    source = (ROOT/('batch_kernels_2d.F90-mc' if dimension==2 else 'batch_kernels.F90-mc')).read_text()
    # Reuse the reference face algebra verbatim, excluding traversal/reduction/reporting.
    face = (REF/'hy_getFaceFlux.F90-mc').read_text()
    face = face.split('            @hy_assembleFaceStates', 1)[1]
    face = '@hy_assembleFaceStates' + face.split('\n          end do', 1)[0]
    face = face.replace('work%minus','minus').replace('work%plus','plus')
    face = face.replace('data,','state,').replace('maxSpeed=max(maxSpeed,speed)', '')
    face = face.replace('return', 'exit cell_work')
    source = source.replace('! INSERT_REFERENCE_FACE_BODY', face)
    source = m.lower_blocks(exp.expand_text(source))
    if re.search(r'\bassociate\s*\(', source, re.I):
        raise m.ExpansionError('Unexpected device ASSOCIATE remains')
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
        return declarations + '\n    error=0\n!$omp target teams distribute parallel do collapse('+str(dimension+1)+') &\n' + \
            '\n'.join(clauses) + '\n!$omp& reduction(max:error)'
    source = re.sub(r'! PRIVATE_BEGIN\n(.*?)! PRIVATE_END', private, source, flags=re.S)
    # Shared small maps are initialized on the host before the face launch.
    start=source.index('  subroutine faces')
    end=source.index('  end subroutine',start)
    init='    inputMap=[1,2,3,4,5,6,7]; outputMap=[1,2,3,4,5]\n'
    init+=('    shockVariable=9; e=0; e(d)=1\n' if dimension==2 else '    shockVariable=9; d=1; e=[1,0,0]\n')
    part=source[start:end].replace('    error=0',init+'    error=0')
    source=source[:start]+part+source[end:]
    for name in ['update','timestep_rate']:
        start=source.index('  subroutine '+name+'(')
        end=source.index('  end subroutine',start)
        part=source[start:end].replace('    error=0', ('    spacing=[dx,dy,1._rk]\n    error=0' if dimension==2 else '    spacing=[dx,1._rk,1._rk]\n    error=0'))
        source=source[:start]+part+source[end:]
    if dimension==2:
        start=source.index('  subroutine reconstruct(')
        end=source.index('  end subroutine',start)
        part=source[start:end].replace('    error=0','    e=0; e(d)=1\n    error=0')
        source=source[:start]+part+source[end:]
    resident={'flatten':'state,flat','reconstruct':'state,flat,minus,plus',
              'faces':'state,minus,plus,flux','update':'reference,state,flux,integrated',
              'acoustics':'state,soundWork','markers':'state,soundWork','timestep_rate':'state'}
    for name,arrays in resident.items():
        start=source.index('  subroutine '+name+'(')
        end=source.index('  end subroutine',start)
        part=source[start:end].replace('collapse('+str(dimension+1)+') &', 'collapse('+str(dimension+1)+') map(alloc:'+arrays+') &')
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
        (output/name).write_text(m.wrap_fortran((ROOT/('periodic_batch_2d.F90' if dimension==2 else name)).read_text()) if name.endswith('.F90') else (ROOT/name).read_text())
    print(output)

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--dimension',type=int,choices=[1,2],default=1)
    args=parser.parse_args()
    output=args.output or ROOT/('build2d' if args.dimension==2 else 'build')
    generate(output.resolve(),args.dimension)
