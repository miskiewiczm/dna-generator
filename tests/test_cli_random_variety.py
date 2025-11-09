import unittest
import sys
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


class TestCLIRandomVariety(unittest.TestCase):
    def test_cli_random_mode_variety(self):
        cmd = [
            sys.executable, '-m', 'dna_generator',
            '--initial', 'CCTGTCATCACGCTAGTAAC',
            '--length', '120',
            '--sequences-only',
            '--mode', 'random',
            '--quiet',
            '--count', '10',
            '--profile', 'sequence_only',
        ]
        env = dict(**os.environ)
        env['PYTHONPATH'] = str(ROOT.parent)
        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT.parent), env=env)
        self.assertEqual(proc.returncode, 0, f"CLI failed: {proc.stderr}")
        lines = [ln.strip() for ln in proc.stdout.strip().splitlines() if ln.strip()]
        self.assertTrue(lines, "No output sequences produced")
        self.assertGreater(len(set(lines)), 1, "RANDOM mode CLI should output diverse sequences across --count")


if __name__ == '__main__':
    unittest.main()
