#!/usr/bin/env python3
"""Generate Flash-X AMR kernels using unchanged NewImpl arithmetic. No build/run."""
from pathlib import Path
import argparse
import importlib.util
import re
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('gpu_generator',ROOT.parent/'generate.py')
g=importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
m=g.m

def generate(output):
    output.mkdir(parents=True,exist_ok=True)
    defs=[g.REF/'hydro_layout.ini',g.REF/'hydro_helpers.ini',g.REF/'update_helpers.ini',ROOT/'amr_layout.ini-template']
    source=(ROOT/'amr_kernels.F90-template').read_text()
    face=(g.REF/'hy_getFaceFlux.F90-mc').read_text().split('            @hy_assembleFaceStates',1)[1]
    face='@hy_assembleFaceStates'+face.split('\n          end do',1)[0]
    face=face.replace('work%minus','minus').replace('work%plus','plus').replace('data,','state,')
    face=face.replace('maxSpeed=max(maxSpeed,speed)','').replace('return','exit cell_work')
    diagnostic_markers={
        '@hy_assembleFaceStates':'FACE_STATES',
        '@hy_estimateHydroWaveSpeeds':'WAVE_SPEEDS',
        '@hy_selectHLLEFlux':'HLLE',
        '@hy_computeHydroContact':'CONTACT',
        '@hy_computeHydroStarState(ql':'LEFT_STAR',
        '@hy_computeHydroStarState(qr':'RIGHT_STAR',
    }
    for token,label in diagnostic_markers.items():
        face=face.replace(token,'! DIAGNOSTIC_'+label+'\n            '+token)
    source=m.lower_blocks(m.MacroExpander(m.load_macros(defs)).expand_text(source.replace('! INSERT_REFERENCE_FACE_BODY',face)))
    def private(match):
        decl=match[1]; names=[]
        for line in decl.splitlines():
            if '::' in line and not line.lstrip().startswith('!'):
                names += [re.match(r'\s*(\w+)',s)[1] for s in m.split_arguments(line.split('::',1)[1])]
        clauses=['!$omp& private('+','.join(names[i:i+3])+') &' for i in range(0,len(names),3)]
        return decl+'\n    error=0\n!$omp target teams distribute parallel do collapse(4) &\n'+'\n'.join(clauses)+'\n!$omp& reduction(max:error)'
    source=re.sub(r'! PRIVATE_BEGIN\n(.*?)! PRIVATE_END',private,source,flags=re.S)
    for name,arrays in {'flatten':'state,flat','reconstruct':'state,flat,minus,plus','faces':'state,minus,plus,flux','update':'reference,state,flux'}.items():
        start=source.index('  subroutine '+name+'('); end=source.index('  end subroutine',start)
        part=source[start:end]
        init=''
        if name in ('reconstruct','faces'): init='    e=0; e(d)=1\n'
        if name=='faces': init+='    inputMap=primitiveMap; outputMap=conservedMap; shockVariable=9\n'
        if name=='update': init+='    spacing=[dx,dy,dz]\n'
        part=part.replace('    error=0',init+'    error=0').replace('collapse(4) &','collapse(4) map(alloc:'+arrays+') &')
        source=source[:start]+part+source[end:]
    # Replay the identical expanded face algebra serially on the host only after
    # a device failure. Successful execution pays no extra transfers or launches.
    start=source.index('  subroutine faces('); end=source.index('  end subroutine',start)+len('  end subroutine')
    diagnostic=source[start:end].replace('subroutine faces(', 'subroutine diagnose_faces(',1)
    diagnostic='\n'.join(line for line in diagnostic.splitlines() if not line.lstrip().startswith('!$omp'))
    diagnostic=diagnostic.replace('    logical::vacuum,inShock',
        '    logical::vacuum,inShock\n    character(32)::diagnosticPoint')
    diagnostic=diagnostic.replace('status=HY_OK\n        cell_work',
        "status=HY_OK; sl=0; sr=0; contact=0; pstar=0; ql=0; qr=0\n        diagnosticPoint='face states'\n        cell_work")
    for label in diagnostic_markers.values():
        diagnostic=diagnostic.replace('! DIAGNOSTIC_'+label, "diagnosticPoint='"+label+"'")
    diagnostic=diagnostic.replace('        error=max(error,status)', """        if(status/=HY_OK) then
          error=status
          print *, 'HY_GPU_AMR first failing face: batch slot,i,j,k,d=',block,i,j,k,d
          print *, 'HY_GPU_AMR failing calculation=',trim(diagnosticPoint),' status=',status
          print *, 'HY_GPU_AMR qL [rho,ux,uy,uz,p,gamma,rhoE]=',ql
          print *, 'HY_GPU_AMR qR [rho,ux,uy,uz,p,gamma,rhoE]=',qr
          print *, 'HY_GPU_AMR SL,SR,contact,Pstar=',sl,sr,contact,pstar
          print *, 'HY_GPU_AMR SR-SL,SL-contact,SR-contact=',sr-sl,sl-contact,sr-contact
          return
        end if""")
    source=source.replace('end module',diagnostic+'\nend module')
    assert not re.search(r'\bassociate\s*\(|macro_scope_\d+_arg_|@',source,re.I)
    (output/'hy_gpu_amr_kernels.F90').write_text(m.wrap_fortran(source))
    # Separate module name avoids colliding with existing NewImpl compatibility modules.
    types=(g.REF/'hy_getFaceFlux.F90-mc').read_text().split('\ncontains',1)[0]
    types=types.replace('newimpl_hydro','hy_gpu_amr_types').replace('rk=selected_real_kind(15,307)','rk=kind(1.0)')
    types=types.replace('public :: cell_field, flux_field, hy_getFaceFlux','public :: cell_field, flux_field')
    types=types.replace('flattening_eta=1.e-99_rk', 'flattening_eta=real(max(1.e-99_8,real(tiny(1._rk),8)),rk)')
    types=m.MacroExpander(m.load_macros([g.REF/'hydro_layout.ini'])).expand_text(types)
    (output/'hy_gpu_amr_types.F90').write_text(types+'\nend module\n')
    print(output)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT)
    generate(p.parse_args().output.resolve())
