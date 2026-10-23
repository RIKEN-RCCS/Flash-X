"""Text dump comparison tests; synthetic data only, no Flash-X execution."""
import csv
from pathlib import Path
import tempfile
import unittest
from compare_dumps import FIELDS,compare,load

class DumpTests(unittest.TestCase):
    def write(self,path,rows):
        with path.open('w',newline='') as stream:
            writer=csv.writer(stream); writer.writerow(FIELDS); writer.writerows(rows)

    def row(self,x):
        return [10,0.1,0.01,2,1,1,1,x,0.5,0,0.5,1,1,1,0,0,0,1,2.5,2.5]

    def test_different_rank_partitions_and_indices(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); a=root/'a.csv'; b=root/'b1.csv'; c=root/'b2.csv'; empty=root/'b3.csv'
            left=self.row(0.25); right=self.row(0.75)
            self.write(a,[left,right]); right[4]=99
            self.write(b,[right]); self.write(c,[left]); self.write(empty,[])
            report=compare([a],[str(root/'b*.csv')])
            self.assertTrue(report['passed']); self.assertEqual(report['cells'],2)
            self.assertEqual(report['reference_totals']['mass'],1)
            self.assertEqual(report['reference_totals']['energy'],2.5)

    def test_mismatch_and_duplicate_detection(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); a=root/'a.csv'; b=root/'b.csv'
            self.write(a,[self.row(0.25)])
            row=self.row(0.25); row[FIELDS.index('pressure')]=2
            self.write(b,[row]); self.assertFalse(compare([a],[b])['passed'])
            row[1]=0.2; self.write(b,[row])
            with self.assertRaisesRegex(ValueError,'step/time/dt'):
                compare([a],[b])
            self.write(b,[self.row(0.25),self.row(0.25)])
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                load([b])

    def test_missing_cells_and_nonfinite_data(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); a=root/'a.csv'; b=root/'b.csv'
            self.write(a,[self.row(0.25),self.row(0.75)])
            self.write(b,[self.row(0.25)])
            with self.assertRaisesRegex(ValueError,'mesh/cell keys'):
                compare([a],[b])
            row=self.row(0.25); row[13]=float('nan'); self.write(b,[row])
            with self.assertRaisesRegex(ValueError,'Nonfinite'):
                load([b])

if __name__=='__main__':
    unittest.main()
