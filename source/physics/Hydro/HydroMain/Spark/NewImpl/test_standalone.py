#!/usr/bin/env python3
"""End-to-end periodic Sod checks: layout, block decomposition, RK, and conservation."""
import csv
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from build_standalone import build

class StandaloneTests(unittest.TestCase):
    def test_periodic_sod(self):
        if not shutil.which("gfortran"):
            self.skipTest("gfortran is required")
        baseline = {}
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for layout in ["variable-first", "spatial-first", "records"]:
                exe = build(root / layout, layout)
                for scheme in [2, 3]:
                    for blocks in [1, 4]:
                        with self.subTest(layout=layout, scheme=scheme, blocks=blocks):
                            output = root / f"{layout}-{scheme}-{blocks}.csv"
                            result = subprocess.run(
                                [str(exe), "128", str(blocks), "0.05", str(output), str(scheme)],
                                capture_output=True, text=True)
                            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                            with output.open() as stream:
                                rows = [[float(value) for value in row.values()]
                                        for row in csv.DictReader(stream)]
                            self.assertEqual(len(rows), 128)
                            self.assertTrue(all(math.isfinite(v) for row in rows for v in row))
                            self.assertTrue(all(row[1] > 0 and row[3] > 0 and row[4] > 0 for row in rows))
                            mass = sum(row[1] for row in rows) / 128
                            momentum = sum(row[1] * row[2] for row in rows) / 128
                            energy = sum(row[1] * row[5] for row in rows) / 128
                            self.assertAlmostEqual(mass, 0.5625, delta=1e-10)
                            self.assertAlmostEqual(momentum, 0, delta=1e-10)
                            self.assertAlmostEqual(energy, 1.375, delta=1e-10)
                            for row in rows:
                                self.assertAlmostEqual(row[3], 0.4 * row[1] * row[4], delta=1e-10)
                                self.assertAlmostEqual(row[5], row[4] + 0.5 * row[2]**2, delta=1e-10)
                            # Exact Sod left-star plateau, between rarefaction tail and contact.
                            star = [r for r in rows if 0.512 < r[0] < 0.530]
                            self.assertTrue(star)
                            for row in star:
                                self.assertAlmostEqual(row[1], 0.426319428, delta=0.045)
                                self.assertAlmostEqual(row[2], 0.927452620, delta=0.055)
                                self.assertAlmostEqual(row[3], 0.303130178, delta=0.030)
                            if scheme not in baseline:
                                baseline[scheme] = rows
                            else:
                                max_difference = max(abs(a-b) for row_a,row_b in zip(rows,baseline[scheme])
                                                     for a,b in zip(row_a,row_b))
                                self.assertLess(max_difference, 2e-11)
                bad = subprocess.run([str(exe), "127", "4", "0.01", str(root/"invalid.csv"), "2"],
                                     capture_output=True, text=True)
                self.assertNotEqual(bad.returncode, 0)

if __name__ == "__main__":
    unittest.main()
