"""Generation/setup-contract checks; never compile or execute a Fortran solver."""
from pathlib import Path
import sys
import tempfile
import unittest
import re
import shutil
import subprocess
import generate_amr as g

ROOT=Path(__file__).resolve().parent
SPARK=ROOT.parent.parent
REPO=next(p for p in ROOT.parents if (p/'bin/preProcess.py').is_file())
sys.path.insert(0,str(REPO/'bin'))
from preProcess import preProcess
from globals import SetupError

class AmrGenerationTests(unittest.TestCase):
    def test_device_sources_and_precision(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory); g.generate(output)
            text=(output/'hy_gpu_amr_kernels.F90').read_text()
            types=(output/'hy_gpu_amr_types.F90').read_text()
            self.assertNotRegex(text,r'(?i)\bassociate\s*\(|macro_scope_\d+_arg_|@')
            self.assertIn('state(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb)',text)
            self.assertIn('spacing=[dx,dy,dz]',text)
            self.assertIn('lower(2)-1,upper(2)+1',text)
            self.assertIn('lower(1)-e(1),upper(1)+e(1)',text)
            self.assertIn('coeff(3)*dt',text)
            self.assertIn('densityFloor',text)
            self.assertIn('rk=kind(1.0)',types)
            self.assertNotIn('flattening_eta=1.e-99_rk',types)
            self.assertTrue(all(len(line)<=132 for line in text.splitlines()))
            # Checked-in generated sources must match the reproducible generator output.
            for name in ['hy_gpu_amr_types.F90','hy_gpu_amr_kernels.F90']:
                self.assertEqual((output/name).read_text(),(ROOT/name).read_text())

    def test_2d_and_3d_preprocessing(self):
        cpp=shutil.which('clang') or shutil.which('cpp')
        if not cpp:
            self.skipTest('C preprocessor unavailable; no Fortran compilation is used')
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            for dim in (2,3):
                with self.subTest(dimension=dim):
                    (output/'Simulation.h').write_text(f'#define NDIM {dim}\n#define NSPECIES 0\n'
                        '#define NMASS_SCALARS 0\n#define NSTENCIL 4\n#define EOS_GAMMA_INCLUDED\n')
                    for name in ['constants.h','Spark.h','Eos.h']:
                        (output/name).write_text('')
                    texts={}
                    for name in ['hy_gpu_amr_kernels.F90','hy_gpu_amr_advance.F90','hy_gpu_amr_adapter.F90']:
                        result=subprocess.run([cpp,'-E','-P','-traditional-cpp','-x','c','-I',str(output),str(ROOT/name)],
                                              capture_output=True,text=True)
                        self.assertEqual(result.returncode,0,result.stderr)
                        text=result.stdout; texts[name]=text
                        # Balanced loops include generated named early-exit scopes.
                        stack=[]
                        for number,line in enumerate(text.splitlines(),1):
                            line=line.strip().lower()
                            if re.match(r'^(?:\w+:\s*)?do(?:\s|$)',line):
                                stack.append(number)
                            if re.match(r'^end\s*do\b',line):
                                self.assertTrue(stack,(name,number))
                                stack.pop()
                        self.assertEqual(stack,[],name)
                    kernel=texts['hy_gpu_amr_kernels.F90']
                    self.assertIn('do j=lower(2),upper(2)+e(2)',kernel)
                    self.assertIn('do k=lower(3),upper(3)+e(3)',kernel)
                    self.assertIn('do j=lower(2),upper(2)',kernel)
                    driver=texts['hy_gpu_amr_advance.F90']
                    self.assertIn(f'do d=1,{dim}',driver)
                    self.assertIn('do k=1,nz+ez',driver)
                    self.assertIn('collapse(4)',driver)
                    adapter=texts['hy_gpu_amr_adapter.F90']
                    if dim==3:
                        self.assertIn('do kz=1,nz+1',adapter)
                        self.assertIn('buffers(ix,jy,kz,v,3,q)',adapter)
                    else:
                        self.assertNotIn('buffers(ix,jy,kz,v,3,q)',adapter)
                        self.assertIn('spacing(3)=1.0',adapter)

    def test_setup_selection_and_rejections(self):
        names='useMP5 useTVD useFog useLim03 useGLM nontelescoping omp_offload oacc Milhoja threadWithinBlock GridIndexReorder'.split()
        flags={name:False for name in names}
        lines=preProcess(flags).processFile(str(ROOT/'Config'))
        self.assertIn('PPDEFINE HY_GPU_AMR','\n'.join(line for _,line in lines))
        for name in names:
            with self.subTest(option=name):
                with self.assertRaises(SetupError):
                    preProcess(dict(flags,**{name:True})).processFile(str(ROOT/'Config'))

    def test_host_reflux_remains_outside_gpu_branch(self):
        text=(SPARK/'Hydro.F90-mc').read_text()
        gpu=text.index('      call hy_gpu_amr_hydro(dt)')
        boundary=text.index('#endif',gpu)
        self.assertLess(boundary,text.index('call Grid_communicateFluxes(ALLDIR, UNSPEC_LEVEL)',gpu))
        self.assertIn('call hy_gpu_amr_check()',text[:gpu])
        driver=(ROOT/'hy_gpu_amr_advance.F90').read_text()
        self.assertNotIn('call Grid_',driver)
        self.assertIn('lower(1:NDIM)=1-(2-stage)*stencil',driver)
        self.assertIn('buffers(i,j,k,v,d,b)+wt*flux(i,j,k,v,d,b)',driver)
        self.assertNotIn('wt*dt',driver)
        self.assertIn('candidate>max(eosFloor,eintSwitch*ekin)',driver)
        self.assertIn('state(i,j,k,8,b)=internal+ekin',driver)
        adapter=(ROOT/'hy_gpu_amr_adapter.F90').read_text()
        self.assertIn('state(:,:,:,:,1:count)',adapter)
        self.assertIn('if(count>0) call flush_batch()',adapter)
        self.assertIn('call tile%releaseDataPtr(u,CENTER)',adapter)
        self.assertIn('call tiles(q)%getDataPtr(u,CENTER)',adapter)
        self.assertIn('call Eos_multiDim(MODE_DENS_EI',adapter)

if __name__=='__main__':
    unittest.main()
