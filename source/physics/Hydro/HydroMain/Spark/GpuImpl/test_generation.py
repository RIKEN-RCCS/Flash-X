"""Source-generation checks only: no Fortran compilation or device execution."""
import re
import tempfile
import unittest
from pathlib import Path
import generate as g

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
