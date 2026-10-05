#!/usr/bin/env python3
"""Macro contract tests and compiled numerical checks for three mesh layouts."""
import importlib.util
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('newimpl_expander', ROOT / 'macro_expand.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)

class Tests(unittest.TestCase):
    def test_arguments(self):
        macros = {
            'ref': m.Macro('ref', ('a', 'i'), '${a}(${i})', Path('test'), 1),
            'id': m.Macro('id', ('x',), '${x}', Path('test'), 2),
        }
        exp = m.MacroExpander(macros)
        self.assertEqual(exp.expand_text('@ref(i,i-1)@'), 'i(i-1)')
        self.assertEqual(exp.expand_text('@ref(field,max(i-1,0))@'), 'field(max(i-1,0))')
        self.assertEqual(exp.expand_text('@ref(field,@id(i+1)@)@'), 'field(i+1)')
        self.assertEqual(exp.expand_text('@id(@id(i+1)@)@'), 'i+1')
        for source in ['@unknown@', '@ref(a)@', '@ref(a,i)', '@ref(a,,i)@']:
            with self.assertRaises(m.ExpansionError):
                exp.expand_text(source)
        recursive = m.MacroExpander({'loop': m.Macro('loop', (), '@loop@', Path('test'), 1)})
        with self.assertRaises(m.ExpansionError):
            recursive.expand_text('@loop@')


    def test_scoped_helpers(self):
        import re
        macros = m.load_macros([ROOT / 'hydro_layout.ini', ROOT / 'hydro_helpers.ini'])
        exp = m.MacroExpander(macros)
        expanded = exp.expand_text('@hy_computeHydroContact(ql,qr,1,sl,sr,c,p,s)@\n'
                                   '@hy_computeHydroContact(ql,qr,1,sl,sr,c,p,s)@')
        names = re.findall(r'(macro_scope_\d+): block', expanded)
        self.assertEqual(len(names), 2)
        self.assertEqual(len(set(names)), 2)
        self.assertNotRegex(expanded, r'\breturn\b')
        self.assertIn('exit ' + names[0], expanded)
        source = (ROOT / 'hy_getFaceFlux.F90-mc').read_text()
        self.assertEqual(re.findall(r'^  subroutine (\w+)', source, re.M), ['hy_getFaceFlux'])
        self.assertNotRegex(source, r'\bcall hy_(?!getFaceFlux)')

    def test_block_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'scope.ini'
            path.write_text("simple\nargs=x,out,obj\nscope=block\ndeclarations=\n"
                            "  integer :: i\ndefinition=\n"
                            "  i=x\n  obj%i=i\n  out=i\n"
                            "  if (x<0) return\n  print *, 'return i' ! i and return\n")
            macros = m.load_macros([path])
            self.assertEqual(macros['simple'].scope, 'block')
            expanded = m.MacroExpander(macros).expand_text('@simple(next(),i,obj)@')
            self.assertEqual(expanded.count('next()'), 1)
            self.assertIn('macro_scope_1_arg_obj%i=macro_scope_1_local_i', expanded)
            self.assertIn("'return i' ! i and return", expanded)
            self.assertIn('exit macro_scope_1', expanded)
            for text in ['bad\nscope=unknown\ndefinition=x',
                         'bad\ndeclarations=integer :: i\ndefinition=x']:
                path.write_text(text)
                with self.assertRaises(m.DefinitionError):
                    m.load_macros([path])

    def test_hoisted_declarations(self):
        text = ("program p\n!$macro_imports\nimplicit none\n!$macro_declarations\n"
                "print *, 1\nmacro_scope_1: block\ninteger :: macro_scope_1_local_i\n"
                "macro_scope_1_local_i=2\nexit macro_scope_1\nend block macro_scope_1\nend program\n")
        lowered = m.lower_blocks(text)
        self.assertLess(lowered.index('integer ::'), lowered.index('print *, 1'))
        self.assertIn('macro_scope_1: do', lowered)
        self.assertNotIn(': block', lowered)
        with self.assertRaises(m.ExpansionError):
            m.lower_blocks('program p\nmacro_scope_1: block\nend block macro_scope_1\n')

    def test_flashx_wrapper(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is needed for wrapper tests')
        for ndim in [1, 2, 3]:
            for double in [False, True]:
                with self.subTest(ndim=ndim, double=double), tempfile.TemporaryDirectory() as tmp:
                    build = Path(tmp)
                    (build / 'Simulation.h').write_text(
                        f'#define NDIM {ndim}\n#define MDIM 3\n#define NSTENCIL 3\n'
                        '#define NSPECIES 0\n#define NMASS_SCALARS 0\n'
                        '#define DENS_VAR 8\n#define VELX_VAR 3\n#define VELY_VAR 10\n'
                        '#define VELZ_VAR 2\n#define PRES_VAR 6\n#define GAMC_VAR 9\n'
                        '#define EINT_VAR 4\n#define SHOK_VAR 11\n')
                    (build / 'constants.h').write_text('#define LOW 1\n#define HIGH 2\n')
                    shutil.copy(ROOT.parent / 'Spark.h', build / 'Spark.h')
                    (build / 'mock.f90').write_text(
                        'module Hydro_data\nimplicit none\nlogical :: hy_hybridRiemann=.true.\n'
                        'real :: hy_cvisc=0.2,hy_tiny=1.e-30,hy_smalldens=1.e-12,hy_smallpres=1.e-12\n'
                        'end module\nsubroutine Driver_abort(message)\ncharacter(*) :: message\n'
                        'print *,message\nstop 99\nend subroutine\n')
                    definitions = [ROOT / 'hydro_layout.ini', ROOT / 'hydro_helpers.ini']
                    m.expand_file(ROOT / 'hy_getFaceFlux.F90-mc', build / 'hydro.f90',
                                  m.MacroExpander(m.load_macros(definitions)))
                    definitions.append(ROOT / 'hydro_layout_flashx.ini')
                    m.expand_file(ROOT / 'hy_rk_getFaceFlux_wrapper.F90-mc', build / 'wrapper.F90',
                                  m.MacroExpander(m.load_macros(definitions)))
                    flags = ['-fdefault-real-8'] if double else []
                    subprocess.run([compiler, '-cpp', '-std=f2003', '-fcheck=all',
                                    '-ffpe-trap=invalid,zero,overflow', *flags, '-I', str(build),
                                    'hydro.f90', 'mock.f90', 'wrapper.F90',
                                    str(ROOT / 'test_wrapper.F90'), '-o', 'test'], cwd=build, check=True)
                    subprocess.run([str(build / 'test')], cwd=build, check=True)

    def test_update_wrapper(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is needed for wrapper tests')
        for ndim in [1, 2, 3]:
            for double in [False, True]:
                with self.subTest(ndim=ndim, double=double), tempfile.TemporaryDirectory() as tmp:
                    build = Path(tmp)
                    (build / 'Simulation.h').write_text(
                        f'#define NDIM {ndim}\n#define MDIM 3\n#define NSPECIES 0\n#define NMASS_SCALARS 0\n'
                        '#define DENS_VAR 8\n#define VELX_VAR 3\n#define VELY_VAR 10\n#define VELZ_VAR 2\n'
                        '#define PRES_VAR 6\n#define GAMC_VAR 9\n#define EINT_VAR 4\n#define ENER_VAR 7\n'
                        '#define GPOT_VAR 5\n#define GPOL_VAR 11\n#define GRAVITY\n')
                    (build / 'constants.h').write_text('#define LOW 1\n#define HIGH 2\n#define CARTESIAN 1\n')
                    shutil.copy(ROOT.parent / 'Spark.h', build / 'Spark.h')
                    (build / 'mock.f90').write_text(
                        'module Hydro_data\nimplicit none\ninteger :: hy_geometry=1\n'
                        'real :: hy_smallE=1.e-12,hy_smalldens=1.e-12\n'
                        'real :: hy_coeffArray(2,3)=reshape([0.,0.25,1.,0.75,0.,0.5],[2,3])\n'
                        'end module\nsubroutine Driver_abort(message)\ncharacter(*) :: message\n'
                        'print *, message\nstop 99\nend subroutine\n')
                    defs = [ROOT / 'hydro_layout.ini', ROOT / 'hydro_helpers.ini']
                    m.expand_file(ROOT / 'hy_getFaceFlux.F90-mc', build / 'hydro.f90',
                                  m.MacroExpander(m.load_macros(defs)))
                    defs += [ROOT / 'update_helpers.ini', ROOT / 'hydro_layout_flashx.ini']
                    if double:
                        defs.append(ROOT / 'update_flashx_nontelescoping.ini')
                    driver = (ROOT / 'test_update_wrapper.F90').read_text()
                    if double:
                        driver = driver.replace('tmpState(1:,loGC(1):,loGC(2):,loGC(3):)',
                                                'tmpState(1:,lo(1):,lo(2):,lo(3):)')
                        driver = driver.replace('MAXSTAGE,state,ref,gravity',
                                                'MAXSTAGE,state,ref(:,1:,1:,1:),gravity')
                    (build / 'test.F90').write_text(driver)
                    m.expand_file(ROOT / 'hy_rk_updateSoln_wrapper.F90-mc', build / 'wrapper.F90',
                                  m.MacroExpander(m.load_macros(defs)))
                    flags = ['-fdefault-real-8'] if double else []
                    subprocess.run([compiler, '-cpp', '-std=f2003', '-fcheck=all',
                                    '-ffpe-trap=invalid,zero,overflow', *flags, '-I', str(build),
                                    'hydro.f90', 'mock.f90', 'wrapper.F90',
                                    'test.F90', '-o', 'test'], cwd=build, check=True)
                    subprocess.run([str(build / 'test')], cwd=build, check=True)

    def test_update(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is needed for update tests')
        for layout in [None, 'hydro_layout_spatial_first.ini', 'hydro_layout_records.ini']:
            with self.subTest(layout=layout), tempfile.TemporaryDirectory() as tmp:
                build = Path(tmp)
                defs = [ROOT / 'hydro_layout.ini', ROOT / 'hydro_helpers.ini', ROOT / 'update_helpers.ini']
                if layout:
                    defs.append(ROOT / layout)
                exp = m.MacroExpander(m.load_macros(defs))
                for source, target in [('hy_getFaceFlux.F90-mc','hydro.f90'),
                                       ('hy_updateSolution.F90-mc','update.f90'),
                                       ('test_update.F90-mc','test.f90')]:
                    m.expand_file(ROOT / source, build / target, exp)
                subprocess.run([compiler, '-std=f2003', '-Wall', '-Wextra', '-Wno-compare-reals',
                                '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                                'hydro.f90', 'update.f90', 'test.f90', '-o', 'test'], cwd=build, check=True)
                subprocess.run([str(build / 'test')], cwd=build, check=True)

    def test_hydro(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is needed for numerical tests')
        for alternate in [None, 'hydro_layout_spatial_first.ini', 'hydro_layout_records.ini']:
            with self.subTest(layout=alternate), tempfile.TemporaryDirectory() as tmp:
                build = Path(tmp)
                definitions = [ROOT / 'hydro_layout.ini', ROOT / 'hydro_helpers.ini']
                if alternate:
                    definitions.append(ROOT / alternate)
                exp = m.MacroExpander(m.load_macros(definitions))
                for source, target in [('hy_getFaceFlux.F90-mc', 'hydro.f90'),
                                       ('test_hydro.F90-mc', 'test.f90')]:
                    m.expand_file(ROOT / source, build / target, exp)
                    generated = (build / target).read_text()
                    self.assertNotRegex(generated, r'(?im)^\s*(?:\w+:\s*)?(?:end\s+)?block\b')
                subprocess.run([compiler, '-std=f2003', '-Wall', '-Wextra', '-Wno-compare-reals',
                                '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                                '-fimplicit-none', 'hydro.f90', 'test.f90', '-o', 'test'],
                               cwd=build, check=True)
                subprocess.run([str(build / 'test')], cwd=build, check=True)

if __name__ == '__main__':
    unittest.main()
