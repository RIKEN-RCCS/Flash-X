"""Generation/setup-contract checks; never compile or execute a Fortran solver."""
from pathlib import Path
import sys
import tempfile
import unittest
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
            self.assertIn('state(1-ng:nx+ng,1-ng:ny+ng,9,nb)',text)
            self.assertIn('spacing=[dx,dy,1._rk]',text)
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
        self.assertIn('lower=1-(2-stage)*stencil',driver)
        self.assertIn('buffers(i,j,v,d,b)+wt*flux(i,j,v,d,b)',driver)
        self.assertNotIn('wt*dt',driver)
        self.assertIn('candidate>max(eosFloor,eintSwitch*ekin)',driver)
        self.assertIn('state(i,j,8,b)=internal+ekin',driver)
        adapter=(ROOT/'hy_gpu_amr_adapter.F90').read_text()
        self.assertIn('state(:,:,:,1:count)',adapter)
        self.assertIn('if(count>0) call flush_batch()',adapter)
        self.assertIn('call tile%releaseDataPtr(u,CENTER)',adapter)
        self.assertIn('call tiles(q)%getDataPtr(u,CENTER)',adapter)
        self.assertIn('call Eos_multiDim(MODE_DENS_EI',adapter)

if __name__=='__main__':
    unittest.main()
