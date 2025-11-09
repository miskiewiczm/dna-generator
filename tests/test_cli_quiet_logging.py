import unittest
import sys
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


class TestCLIQuietLogging(unittest.TestCase):
    def test_quiet_suppresses_info_logs(self):
        cmd = [
            sys.executable, '-m', 'dna_generator',
            '--initial', 'CCTGTCATCACGCTAGTAAC',
            '--length', '120',
            '--sequences-only',
            '--quiet',
            '--seed', '3',
            '--profile', 'sequence_only',
        ]
        env = dict(**os.environ)
        env['PYTHONPATH'] = str(ROOT.parent)
        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT.parent), env=env)
        self.assertEqual(proc.returncode, 0, f"CLI failed: {proc.stderr}")
        out = proc.stdout.strip()
        # Expect only sequence characters and newlines (no log lines with timestamps or levels)
        self.assertTrue(out, "No sequence output produced")
        self.assertRegex(out, r'^[ATGC\n]+$', f"Unexpected characters in output (likely logs present): {out[:120]}")


if __name__ == '__main__':
    unittest.main()
