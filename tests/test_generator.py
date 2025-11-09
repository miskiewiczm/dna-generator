import unittest
import sys
from pathlib import Path

# Ensure package import works when running tests from repo root
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent))

from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode


class TestGeneratorDeterminism(unittest.TestCase):
    def test_deterministic_mode_reproducible(self):
        # Use minimal validation to test determinism without complexity
        cfg = GeneratorConfig(
            generation_mode=GenerationMode.DETERMINISTIC,
            seed=123,
            validation_profile='sequence_only',  # Basic sequence validation
            enable_progress_logging=False,
            max_backtrack_attempts=10000,
            window_size=12,
            max_homopolymer_length=6,  # More relaxed
            max_dinucleotide_repeats=4,  # More relaxed
        )

        gen = DNAGenerator(cfg)
        initial = 'ATGCTAGCTAGC'  # Longer initial sequence
        target = 50  # Shorter target
        res1 = gen.generate(initial, target)
        res2 = gen.generate(initial, target)
        self.assertTrue(res1.success, f"First generation failed: {res1.error_message}")
        self.assertTrue(res2.success, f"Second generation failed: {res2.error_message}")
        self.assertEqual(res1.sequence, res2.sequence, "Deterministic mode should produce identical sequences")


if __name__ == '__main__':
    unittest.main()

