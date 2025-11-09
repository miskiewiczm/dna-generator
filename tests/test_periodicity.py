import unittest
import sys
from pathlib import Path

# Ensure package import works when running tests from repo root
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent))

from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode


class TestPeriodicity(unittest.TestCase):
    def test_no_strong_periodicity_at_window_size(self):
        """Long sequence should not exhibit strong repeating pattern at lag=window_size."""
        initial = "CCTGTCATCACGCTAGTAAC"
        length = 600
        window_size = 20

        cfg = GeneratorConfig(
            generation_mode=GenerationMode.DETERMINISTIC,
            validation_profile='sequence_only',  # avoid primer3
            window_size=window_size,
            enable_progress_logging=False,
            # heuristics enabled by default
        )
        gen = DNAGenerator(cfg)
        res = gen.generate(initial, length, seed=2)
        self.assertTrue(res.success, f"Generation failed: {res.error_message}")

        seq = res.sequence
        lag = window_size
        matches = 0
        total = len(seq) - lag
        for i in range(lag, len(seq)):
            if seq[i] == seq[i - lag]:
                matches += 1
        match_ratio = matches / total if total > 0 else 0.0

        # Expect no pathological periodicity at lag equal to window size
        # Threshold chosen to catch near-deterministic repetition but allow natural correlations
        self.assertLess(match_ratio, 0.7, f"Excessive periodicity at lag {lag}: ratio={match_ratio:.3f}")


if __name__ == '__main__':
    unittest.main()

