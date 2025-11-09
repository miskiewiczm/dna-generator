"""
Backtracking engine for DNA sequence generation.

This module contains the core backtracking algorithm logic separated
from the main DNAGenerator class for better separation of concerns.
"""

from typing import Optional, List, Dict, Any, Tuple
import logging

from .config import GeneratorConfig, GenerationMode
from dna_commons import DNAValidator, DeterministicRandom

logger = logging.getLogger(__name__)


class BacktrackingEngine:
    """
    Engine implementing the backtracking algorithm for DNA generation.

    This class encapsulates the core algorithm logic, separated from
    the DNAGenerator API layer.
    """

    # Heuristic scoring constants
    HOMOPOLYMER_PENALTY = 1e6
    PRIME3_STABILITY_PENALTY = 5e5
    DINUCLEOTIDE_PENALTY = 2e5
    DIVERSITY_WEIGHT = 0.05
    RANDOM_MODE_MARGIN = 0.02
    SCORE_EPSILON = 1e-9

    def __init__(self, config: GeneratorConfig, validator: DNAValidator):
        """
        Initialize the backtracking engine.

        Args:
            config: Generation configuration
            validator: DNA sequence validator
        """
        self.config = config
        self.validator = validator

    def generate_sequence(self,
                         initial_sequence: str,
                         target_length: int,
                         random_gen: DeterministicRandom) -> Tuple[Optional[str], Dict[str, Any]]:
        """
        Generate a sequence using the backtracking algorithm.

        Args:
            initial_sequence: Initial sequence
            target_length: Target length
            random_gen: Random number generator

        Returns:
            Tuple (sequence, stats) or (None, stats) on failure
        """
        bases = ['A', 'T', 'G', 'C']
        seq_list = list(initial_sequence)

        # Initialize statistics
        stats = self._initialize_backtracking_stats(len(initial_sequence))

        # Initialize backtracking states
        states = []
        states.append(bases.copy())

        while len(seq_list) < target_length:
            stats['total_attempts'] += 1

            # Check stopping conditions
            if self._should_stop_backtracking(stats, states):
                break

            # Handle backtracking if needed
            if self._handle_backtrack_step(states, seq_list, len(initial_sequence), stats):
                continue

            # Choose and test candidate
            current_state = states[-1]
            candidate = self._choose_candidate_with_heuristics(
                seq_list,
                current_state,
                random_gen
            )
            current_state.remove(candidate)

            # Test sequence with candidate
            test_seq = seq_list + [candidate]
            window = self._get_analysis_window(test_seq)

            # Validate and accept if valid
            if self.validator.validate_window(window):
                self._accept_candidate(seq_list, candidate, states, stats, target_length)

        # Check if full sequence was generated
        if len(seq_list) == target_length:
            return "".join(seq_list), stats
        else:
            return None, stats

    def _initialize_backtracking_stats(self, initial_length: int) -> Dict[str, Any]:
        """Initialize statistics tracking for backtracking algorithm."""
        return {
            'backtrack_count': 0,
            'total_attempts': 0,
            'max_depth_reached': initial_length,
            'validation_failures': {
                'gc_content': 0,
                'melting_temp': 0,
                'homopolymers': 0,
                'dinucleotide_repeats': 0,
                'three_prime_stability': 0,
                'hairpin': 0,
                'homodimer': 0
            },
            'window_rollup': {
                'gc_min': None,
                'gc_max': None,
                'tm_min': None,
                'tm_max': None,
                'hairpin_tm_max': None,
                'homodimer_tm_max': None,
            }
        }

    def _should_stop_backtracking(self, stats: Dict[str, Any], states: List) -> bool:
        """Check if backtracking should stop due to limits or exhaustion."""
        # Check attempt limit
        if stats['total_attempts'] > self.config.max_backtrack_attempts:
            logger.warning(f"Exceeded backtracking attempt limit ({self.config.max_backtrack_attempts})")
            return True

        # Check if exhausted all options
        if not states:
            logger.warning("Exhausted all backtracking options")
            return True

        return False

    def _handle_backtrack_step(self, states: List, seq_list: List[str],
                              initial_length: int, stats: Dict[str, Any]) -> bool:
        """Handle one backtracking step. Returns True if should continue."""
        current_state = states[-1]

        if not current_state:
            # Backtrack - no more options at this position
            states.pop()
            if len(seq_list) > initial_length:
                seq_list.pop()
                stats['backtrack_count'] += 1
            return True  # Continue with next iteration

        return False  # Don't continue, proceed with candidate selection

    def _get_analysis_window(self, test_seq: List[str]) -> str:
        """Get analysis window from test sequence."""
        if len(test_seq) >= self.config.window_size:
            return "".join(test_seq[-self.config.window_size:])
        else:
            return "".join(test_seq)

    def _accept_candidate(self, seq_list: List[str], candidate: str,
                         states: List, stats: Dict[str, Any], target_length: int) -> None:
        """Accept a candidate nucleotide and update state."""
        bases = ['A', 'T', 'G', 'C']

        # Accept nucleotide
        seq_list.append(candidate)
        stats['max_depth_reached'] = max(stats['max_depth_reached'], len(seq_list))

        # Add state for next position (if needed)
        states.append(bases.copy())

        # Progress logging
        if self.config.enable_progress_logging and len(seq_list) % 50 == 0:
            progress = (len(seq_list) / target_length) * 100
            logger.info(f"Progress: {len(seq_list)}/{target_length} ({progress:.1f}%)")

    def _choose_candidate_with_heuristics(self,
                                          seq_list: List[str],
                                          options: List[str],
                                          random_gen: DeterministicRandom) -> str:
        """
        Choose a candidate using lightweight heuristics.

        Falls back to random choice when heuristics are disabled.
        """
        if not self.config.enable_backtrack_heuristics or len(options) <= 1:
            return random_gen.choice(options)

        # Only calculate target_gc if GC content validation is enabled
        if self.validator.rules.gc_content:
            target_gc = (self.config.min_gc + self.config.max_gc) / 2.0
        else:
            target_gc = None  # No GC targeting when validation disabled

        # Score each option
        scored = [(self._calculate_heuristic_score(base, seq_list, target_gc), base)
                 for base in options]
        scored.sort(key=lambda x: x[0])

        return self._select_candidate_from_scores(scored, random_gen)

    def _calculate_heuristic_score(self, base: str, seq_list: List[str],
                                  target_gc: float) -> float:
        """Calculate heuristic score for a candidate base."""
        test_seq = seq_list + [base]

        # Get analysis window
        window = self._get_analysis_window(test_seq)

        # Hard constraints (high penalties) - only for enabled validation rules
        if self.validator.rules.homopolymer_runs:
            has_homopolymer_violations, _ = self.validator._check_homopolymer_runs(window)
            if has_homopolymer_violations:
                return self.HOMOPOLYMER_PENALTY

        if self.validator.rules.three_prime_stability:
            if not self.validator._check_3_prime_stability(window):
                return self.PRIME3_STABILITY_PENALTY

        if self.validator.rules.dinucleotide_repeats:
            has_dinucleotide_violations, _ = self.validator._check_dinucleotide_repeats(window)
            if has_dinucleotide_violations:
                return self.DINUCLEOTIDE_PENALTY

        # Soft preference: move GC towards target in window (only if enabled)
        dist = 0.0
        if target_gc is not None:
            gc = self.validator._calculate_gc_content(window)
            dist = abs(gc - target_gc)

        # Diversity term
        if len(window) > 0:
            base_freq = window.count(base) / len(window)
        else:
            base_freq = 0.0

        # Novelty term to reduce periodic patterns
        recent_len = max(100, self.config.window_size)
        recent = ("".join(seq_list))[-recent_len:]
        novelty_penalty = 0.0

        for k in (3, 4, 5):
            if len(test_seq) >= k:
                kmer = ("".join(test_seq[-k:]))
                count = recent.count(kmer)
                novelty_penalty += 0.01 * count

        return dist + self.DIVERSITY_WEIGHT * base_freq + novelty_penalty

    def _select_candidate_from_scores(self, scored: List[tuple],
                                     random_gen: DeterministicRandom) -> str:
        """Select final candidate from scored options."""
        best_score = scored[0][0]

        # RANDOM mode: choose randomly among top candidates within margin
        if self.config.generation_mode == GenerationMode.RANDOM:
            eligible = [b for s, b in scored if s <= best_score + self.RANDOM_MODE_MARGIN]
            if len(eligible) > 1:
                return random_gen.choice(eligible)

        # Deterministic: if multiple exactly-best, choose randomly
        best_candidates = [b for s, b in scored if abs(s - best_score) < self.SCORE_EPSILON]
        if len(best_candidates) > 1:
            return random_gen.choice(best_candidates)
        return scored[0][1]