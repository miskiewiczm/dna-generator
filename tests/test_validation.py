import unittest
from typing import Dict, Any

import sys
from pathlib import Path

# Ensure package import works when running tests from repo root
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent))

from dna_commons import DNAValidator
from dna_commons.validation import ValidationRules


class TestValidationRules(unittest.TestCase):
    def test_homopolymer_semantics_inclusive(self):
        rules = ValidationRules(
            gc_content=False,
            melting_temperature=False,
            hairpin_structures=False,
            homodimer_structures=False,
            homopolymer_runs=True,
            dinucleotide_repeats=False,
            three_prime_stability=False,
            max_homopolymer_length=4
        )
        v = DNAValidator(rules)
        # Run of exactly 4 is allowed
        self.assertTrue(v.validate_window('ATGCAAAATGCA'))
        # Run of 5 should be rejected
        self.assertFalse(v.validate_window('ATGCAAAAATGCA'))

    def test_dinucleotide_repeats_limit(self):
        rules = ValidationRules(
            gc_content=False,
            melting_temperature=False,
            hairpin_structures=False,
            homodimer_structures=False,
            homopolymer_runs=False,
            dinucleotide_repeats=True,
            three_prime_stability=False,
            max_dinucleotide_repeats=2
        )
        v = DNAValidator(rules)
        # AT repeated twice (ATAT) allowed
        self.assertTrue(v.validate_window('GGATATCC'))
        # AT repeated three times (ATATAT) rejected
        self.assertFalse(v.validate_window('GGATATATCC'))

    def test_three_prime_stability_limit(self):
        rules = ValidationRules(
            gc_content=False,
            melting_temperature=False,
            hairpin_structures=False,
            homodimer_structures=False,
            homopolymer_runs=False,
            dinucleotide_repeats=False,
            three_prime_stability=True
        )
        v = DNAValidator(rules)
        # Last 5 nt have 4 GC (GGGGC) -> should fail due to 3' stability
        self.assertFalse(v.validate_window('ATATAGGGGC'))
        # Last 5 nt have 3 GC -> should pass
        self.assertTrue(v.validate_window('ATATAGGGAT'))

    def test_gc_content_window_check_only(self):
        # Enable only GC check to avoid primer3 dependency and other checks
        rules = ValidationRules(
            gc_content=True,
            melting_temperature=False,
            hairpin_structures=False,
            homodimer_structures=False,
            homopolymer_runs=False,
            dinucleotide_repeats=False,
            three_prime_stability=False,
            min_gc=0.45,
            max_gc=0.55
        )
        v = DNAValidator(rules)
        # 50% GC should pass with default 0.45-0.55
        self.assertTrue(v.validate_window('ATGCCGTA'))
        # 100% GC should fail
        self.assertFalse(v.validate_window('GGGGCCCC'))

    def test_validate_sequence_without_primer3_checks(self):
        # Ensure validate_sequence works without termodynamic checks enabled
        rules = ValidationRules(
            gc_content=True,
            melting_temperature=False,
            hairpin_structures=False,
            homodimer_structures=False,
            homopolymer_runs=False,
            dinucleotide_repeats=False,
            three_prime_stability=False,
            min_gc=0.0,
            max_gc=1.0
        )
        v = DNAValidator(rules)
        # ATGC has 50% GC and should be valid under gc_content-only rules
        metrics = v.validate_sequence('ATGC')
        self.assertTrue(metrics.is_valid)


if __name__ == '__main__':
    unittest.main()

