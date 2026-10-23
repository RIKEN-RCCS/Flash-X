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
            diagnostic=text.split('  subroutine diagnose_faces(',1)[1]
            self.assertNotIn('!$omp',diagnostic)
            self.assertIn("diagnosticPoint='LEFT_STAR'",diagnostic)
            self.assertIn("diagnosticPoint='RIGHT_STAR'",diagnostic)
            self.assertIn('first failing face',diagnostic)
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
            for dim,cpu in ((2,False),(3,False),(2,True),(3,True)):
                with self.subTest(dimension=dim,cpu=cpu):
                    (output/'Simulation.h').write_text(f'#define NDIM {dim}\n#define NSPECIES 0\n'
                        '#define NMASS_SCALARS 0\n#define NSTENCIL 4\n#define EOS_GAMMA_INCLUDED\n')
                    if cpu:
                        with (output/'Simulation.h').open('a') as f:
                            f.write('#define HY_GPU_AMR_CPU_TEST\n')
                    for name in ['constants.h','Spark.h','Eos.h']:
                        (output/name).write_text('')
                    texts={}
                    for name in ['hy_gpu_amr_kernels.F90','hy_gpu_amr_advance.F90','hy_gpu_amr_adapter.F90','hy_gpu_amr_dump.F90','hy_gpu_amr_timing.F90']:
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
                    if cpu:
                        self.assertIn('CPU testing mode',adapter)
                        self.assertNotIn('target device unavailable',adapter)
                    else:
                        self.assertIn('target device unavailable',adapter)
                        self.assertNotIn('CPU testing mode',adapter)
                    if dim==3:
                        self.assertIn('do kz=1,nz+1',adapter)
                        self.assertIn('buffers(ix,jy,kz,v,3,q)',adapter)
                    else:
                        self.assertNotIn('buffers(ix,jy,kz,v,3,q)',adapter)
                        self.assertIn('spacing(3)=1.0',adapter)

    def test_timing_boundaries_and_build_dependencies(self):
        for name in ('hy_gpu_amr_adapter.F90','hy_gpu_amr_advance.F90'):
            text=(ROOT/name).read_text()
            stack=[]
            for line in text.splitlines():
                start=re.search(r'call profile_start\((P_\w+),',line)
                stop=re.search(r'call profile_stop\((P_\w+),',line)
                if start:
                    stack.append(start[1])
                if stop:
                    self.assertTrue(stack,(name,line))
                    self.assertEqual(stack.pop(),stop[1],(name,line))
            self.assertEqual(stack,[],name)
            # Instrumentation must not split a continued Fortran CALL.
            self.assertNotRegex(text,r'&\s*\n\s*call profile_')
            self.assertTrue(all(len(line)<=132 for line in text.splitlines()),name)
        make=(ROOT/'Makefile').read_text()
        self.assertIn('hy_gpu_amr_timing.o: Timers_interface.o',make)
        self.assertIn('hy_gpu_amr_kernels.o hy_gpu_amr_timing.o',make)
        adapter=(ROOT/'hy_gpu_amr_adapter.F90').read_text()
        self.assertIn('if(count==capacity) call flush_batch(1)',adapter)
        self.assertIn('if(count>0) call flush_batch(3)',adapter)
        self.assertIn('call flush_batch(2)',adapter)
        self.assertLess(adapter.index('call profile_stop(P_TOTAL'),adapter.index('call profile_write'))

    def test_native_macro_parser_does_not_discover_custom_layouts(self):
        # Flash-X loads *.ini from selected units, including AMR's GpuImpl parent.
        # Generated storage definitions must remain available only to our expander.
        from macroProcessor import macroProcessor
        native=macroProcessor()
        for directory in (ROOT.parent,ROOT):
            for path in directory.glob('*.ini'):
                native.loadDefs(str(path))
        for name in ('gpu_layout.ini','gpu_layout_2d.ini'):
            self.assertFalse((ROOT.parent/name).exists())
            self.assertTrue((ROOT.parent/(name+'-template')).is_file())
        self.assertTrue((ROOT/'amr_layout.ini-template').is_file())

    def test_setup_selection_and_rejections(self):
        names='useMP5 useTVD useFog useLim03 useGLM nontelescoping omp_offload oacc Milhoja threadWithinBlock GridIndexReorder'.split()
        flags={name:False for name in names}; flags['gpuAmrCpuTest']=False
        lines=preProcess(flags).processFile(str(ROOT/'Config'))
        self.assertIn('PPDEFINE HY_GPU_AMR','\n'.join(line for _,line in lines))
        for name in names:
            with self.subTest(option=name):
                with self.assertRaises(SetupError):
                    preProcess(dict(flags,**{name:True})).processFile(str(ROOT/'Config'))

    def test_standalone_makefile_preserves_flash_flags(self):
        make=shutil.which('make')
        if not make:
            self.skipTest('make unavailable')
        with tempfile.TemporaryDirectory() as directory:
            wrapper=Path(directory)/'Makefile'
            wrapper.write_text('FLASHBINARY = true\nFFLAGS = -c -O2\n'
                'include '+str(ROOT.parent/'Makefile')+'\n'
                '.PHONY: inspect\ninspect:\n\t@echo FLAGS=$(FFLAGS)\n')
            result=subprocess.run([make,'-f',str(wrapper),'inspect'],capture_output=True,text=True,check=True)
            self.assertIn('FLAGS=-c -O2',result.stdout)
            result=subprocess.run([make,'-f',str(wrapper),'-n','periodic_batch'],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            result=subprocess.run([make,'-f',str(ROOT.parent/'Makefile'),'-n','nvidia'],
                                  capture_output=True,text=True,check=True)
            self.assertIn('nvfortran -O3 -mp=gpu',result.stdout)

    def test_cpu_mode_library_selection(self):
        names='useMP5 useTVD useFog useLim03 useGLM nontelescoping omp_offload oacc Milhoja threadWithinBlock GridIndexReorder'.split()
        for cpu in (False,True):
            flags={name:False for name in names}; flags['gpuAmrCpuTest']=cpu
            lines=preProcess(flags).processFile(str(ROOT/'Config'))
            text='\n'.join(line.strip() for _,line in lines)
            self.assertIn('LIBRARY '+('OMP_HOST' if cpu else 'OMP_OL'),text)
            self.assertNotIn('LIBRARY '+('OMP_OL' if cpu else 'OMP_HOST'),text)
            self.assertEqual('PPDEFINE HY_GPU_AMR_CPU_TEST' in text,cpu)

    def test_host_reflux_remains_outside_gpu_branch(self):
        text=(SPARK/'Hydro.F90-mc').read_text()
        gpu=text.index('      call hy_gpu_amr_hydro(dt)')
        boundary=text.index('#endif',gpu)
        self.assertLess(boundary,text.index('call Grid_communicateFluxes(ALLDIR, UNSPEC_LEVEL)',gpu))
        self.assertIn('call hy_gpu_amr_check()',text[:gpu])
        self.assertGreater(text.index('call hy_gpu_amr_dump(timeEndAdv,dt)'),text.rindex('call hy_rk_correctFluxes'))
        driver=(ROOT/'hy_gpu_amr_advance.F90').read_text()
        self.assertNotIn('call Grid_',driver)
        self.assertIn('lower(1:NDIM)=1-(2-stage)*stencil',driver)
        self.assertIn('buffers(i,j,k,v,d,b)+wt*flux(i,j,k,v,d,b)',driver)
        self.assertNotIn('wt*dt',driver)
        self.assertIn('candidate>max(eosFloor,eintSwitch*ekin)',driver)
        self.assertIn("failedOperation='flattening'",driver)
        self.assertIn("'Riemann flux direction ',d",driver)
        self.assertIn('state(i,j,k,8,b)=internal+ekin',driver)
        adapter=(ROOT/'hy_gpu_amr_adapter.F90').read_text()
        self.assertIn('state(:,:,:,:,1:count)',adapter)
        self.assertIn('if(count>0) call flush_batch(3)',adapter)
        self.assertIn('call tile%releaseDataPtr(u,CENTER)',adapter)
        self.assertIn('call tiles(q)%getDataPtr(u,CENTER)',adapter)
        self.assertIn('call Eos_multiDim(MODE_DENS_EI',adapter)
        self.assertIn('call Driver_abort(trim(abortMessage))',adapter)
        self.assertIn('trim(failedOperation)',adapter)

if __name__=='__main__':
    unittest.main()
