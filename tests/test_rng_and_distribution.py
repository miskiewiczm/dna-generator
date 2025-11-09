import unittest
import sys
import subprocess
import os
from pathlib import Path

# Ensure package import works when running tests from repo root
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent))

from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode


class TestRNGAndDiversity(unittest.TestCase):
    def setUp(self):
        self.initial = "CCTGTCATCACGCTAGTAAC"
        self.length = 200

    def test_different_seeds_yield_different_sequences_when_heuristics_disabled(self):
        cfg = GeneratorConfig(
            generation_mode=GenerationMode.DETERMINISTIC,
            validation_profile='sequence_only',  # avoid primer3 dependency
            enable_backtrack_heuristics=False,   # ensure RNG choices affect output
            enable_progress_logging=False,
            window_size=20,
            max_backtrack_attempts=50000,
        )
        gen = DNAGenerator(cfg)

        seeds = [1, 2, 3]
        sequences = []
        for s in seeds:
            res = gen.generate(self.initial, self.length, seed=s)
            self.assertTrue(res.success, f"Generation failed for seed={s}: {res.error_message}")
            sequences.append(res.sequence)

        # At least two different sequences across different seeds
        self.assertGreater(len(set(sequences)), 1, "Different seeds should produce different sequences when heuristics are disabled")

    def test_nucleotide_diversity_not_starved_in_A_and_G(self):
        cfg = GeneratorConfig(
            generation_mode=GenerationMode.DETERMINISTIC,
            validation_profile='sequence_only',  # sequence-only checks
            enable_progress_logging=False,
            window_size=20,
            max_backtrack_attempts=50000,
        )
        gen = DNAGenerator(cfg)

        seeds = [1, 2, 3, 4, 5]
        min_frac_A = 1.0
        min_frac_G = 1.0

        for s in seeds:
            res = gen.generate(self.initial, self.length, seed=s)
            self.assertTrue(res.success, f"Generation failed for seed={s}: {res.error_message}")
            seq = res.sequence
            frac_A = seq.count('A') / len(seq)
            frac_G = seq.count('G') / len(seq)
            min_frac_A = min(min_frac_A, frac_A)
            min_frac_G = min(min_frac_G, frac_G)

        # Sanity: do not starve A or G below 5% across these seeds
        self.assertGreaterEqual(min_frac_A, 0.05, f"Observed minimal A fraction too low: {min_frac_A:.3f}")
        self.assertGreaterEqual(min_frac_G, 0.05, f"Observed minimal G fraction too low: {min_frac_G:.3f}")

    def test_random_mode_multiple_sequences_not_identical(self):
        cfg = GeneratorConfig(
            generation_mode=GenerationMode.RANDOM,
            validation_profile='sequence_only',
            enable_progress_logging=False,
            window_size=20,
            max_backtrack_attempts=50000,
        )
        gen = DNAGenerator(cfg)
        results = gen.generate_multiple(self.initial, 120, count=10)
        seqs = [r.sequence for r in results if r.success]
        self.assertTrue(seqs, "No sequences generated in RANDOM mode")
        self.assertGreater(len(set(seqs)), 1, "RANDOM mode should produce diverse sequences across multiple generations")

    def test_cli_different_seeds_produce_different_sequences(self):
        # Smoke test via CLI to ensure seed propagates through CLI layer
        cmd_base = [
            sys.executable, "-m", "dna_generator",
            "--initial", self.initial,
            "--length", str(self.length),
            "--sequences-only",
            "--quiet",
            "--profile", "sequence_only",
            "--no-heuristics",
        ]

        env = dict(**os.environ)
        # Ensure the parent of the package is on sys.path in the subprocess
        env["PYTHONPATH"] = str(ROOT.parent)
        proc1 = subprocess.run(
            cmd_base + ["--seed", "1"],
            capture_output=True,
            text=True,
            cwd=str(ROOT.parent),
            env=env,
        )
        self.assertEqual(proc1.returncode, 0, f"CLI failed: {proc1.stderr}")
        out1 = proc1.stdout.strip()

        proc2 = subprocess.run(
            cmd_base + ["--seed", "2"],
            capture_output=True,
            text=True,
            cwd=str(ROOT.parent),
            env=env,
        )
        self.assertEqual(proc2.returncode, 0, f"CLI failed: {proc2.stderr}")
        out2 = proc2.stdout.strip()

        # Expect different sequences between different seeds
        self.assertNotEqual(out1, out2, "CLI produced identical sequences for different seeds with heuristics disabled")


if __name__ == '__main__':
    unittest.main()
