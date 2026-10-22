"""Source-generation checks only: no Fortran compilation or device execution."""
import re
import tempfile
import unittest
from pathlib import Path
import generate as g
from benchmark import parse
from compare_2d import compare, FIELDS
import csv

class GenerationTests(unittest.TestCase):
    def test_direct_binding_preserves_scope_and_precedence(self):
        macro=g.m.Macro('sample',('a','factor','out'),
                        "out=factor*a\nif(a<0) return\nprint *, 'a' ! a\n",
                        Path('test'),1,scope='block',declarations='integer :: scratch')
        text=g.device_scoped_definition(macro,['q','coeff(3)*dt','result'],1)
        self.assertIn('result=(coeff(3)*dt)*q',text)
        self.assertIn('exit macro_scope_1',text)
        self.assertIn("'a' ! a",text)
        self.assertIn('macro_scope_1_local_scratch',text)
        self.assertNotIn('associate',text)

    def test_benchmark_output_parser(self):
        text=('initial_device= 0 target_devices= 1\n'
              'resident stepping seconds= 1.25D-01\n'
              'including transfer seconds= 1.30E-01\nsteps,time= 100 0.01\n')
        self.assertEqual(parse(text),(100,0.125,0.13))
        with self.assertRaises(ValueError):
            parse(text.replace('initial_device= 0','initial_device= 1'))
        with self.assertRaises(ValueError):
            parse('incomplete output')

    def test_rotated_2d_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            first=Path(directory)/'x.csv'; second=Path(directory)/'y.csv'
            def write(path,row):
                with path.open('w',newline='') as f:
                    writer=csv.writer(f); writer.writerow(FIELDS); writer.writerow(row)
            write(first,[0.25,0.75,1,0.2,0,1,2.5,2.52])
            write(second,[0.75,0.25,1,0,0.2,1,2.5,2.52])
            self.assertEqual(compare(first,second,transpose=True),0)
            with self.assertRaises(ValueError):
                compare(first,second)
            write(second,[0.75,0.25,-1,0,0.2,1,2.5,2.52])
            with self.assertRaises(ValueError):
                compare(first,second,transpose=True)

    def test_2d_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            g.generate(Path(directory),2)
            text=(Path(directory)/'batch_kernels.f90').read_text()
            driver=(Path(directory)/'periodic_batch.F90').read_text()
            self.assertIn('module gpu_batch_2d',text)
            self.assertIn('collapse(3)',text)
            self.assertIn('spacing=[dx,dy,1._rk]',text)
            self.assertIn('e=0; e(d)=1',text)
            self.assertIn('(-2:nx+3,-2:ny+3,9,nb)',text)
            self.assertNotRegex(text,r'macro_scope_\d+_arg_|(?i:associate)\s*\(')
            self.assertIn('do d=1,2',driver)
            self.assertIn('globaly/ny',driver)
            self.assertIn('flux(ii,ny+1,vv,2,bb)',driver)
            self.assertNotIn('@',text)
            self.assertTrue(all(len(line)<=132 for line in text.splitlines()))

    def test_complete_sources_have_no_argument_descriptors(self):
        with tempfile.TemporaryDirectory() as directory:
            g.generate(Path(directory))
            text=(Path(directory)/'batch_kernels.f90').read_text()
            self.assertNotRegex(text,r'(?i)\bassociate\s*\(')
            self.assertNotRegex(text,r'macro_scope_\d+_arg_')
            self.assertNotIn('@',text)
            self.assertNotIn('! PRIVATE_',text)
            self.assertIn('spacing=[dx,1._rk,1._rk]',text)
            self.assertIn('reduction(max:error,rate)',text)
            self.assertTrue(all(len(line)<=132 for line in text.splitlines()))

if __name__=='__main__':
    unittest.main()
