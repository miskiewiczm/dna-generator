"""
Main module for the DNA sequence generator.

Provides the DNAGenerator class that generates DNA sequences using a
backtracking algorithm with biochemical quality control.
"""

from typing import Optional, List, Dict, Any
from dataclasses import dataclass
import logging
import time
import secrets

from .config import GeneratorConfig, GenerationMode
from dna_commons import (
    DNAValidator,
    QualityMetrics,
    DeterministicRandom,
    generate_seed_from_string,
    SequenceAnalyzer,
    normalize_dna_sequence
)
from .exceptions import (
    GeneratorError,
    BacktrackingError,
    ValidationError,
    InputError
)
from .backtracking_engine import BacktrackingEngine

logger = logging.getLogger(__name__)


@dataclass
class GenerationResult:
    """
    Result of a DNA sequence generation operation.
    
    Attributes:
        success: Whether generation succeeded
        sequence: Generated DNA sequence (None on failure)
        initial_sequence: Initial sequence
        target_length: Target sequence length
        actual_length: Actual length of the generated sequence
        quality_metrics: Sequence quality metrics
        generation_stats: Generation process statistics
        error_message: Error message (None on success)
        generation_time: Generation time in seconds
    """
    success: bool
    sequence: Optional[str] = None
    initial_sequence: Optional[str] = None
    target_length: Optional[int] = None
    actual_length: Optional[int] = None
    quality_metrics: Optional[QualityMetrics] = None
    generation_stats: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    generation_time: Optional[float] = None
    
    def to_dict(self, raw: bool = False) -> Dict[str, Any]:
        """Convert result to a dictionary.
        
        Args:
            raw: If True, include numeric/raw quality metrics
        """
        return {
            'success': self.success,
            'sequence': self.sequence,
            'initial_sequence': self.initial_sequence,
            'target_length': self.target_length,
            'actual_length': self.actual_length,
            'quality_metrics': self.quality_metrics.to_dict(raw=raw) if self.quality_metrics else None,
            'generation_stats': self.generation_stats,
            'error_message': self.error_message,
            'generation_time': f"{self.generation_time:.2f}s" if self.generation_time else None
        }


class DNAGenerator:
    """
    Main class for generating DNA sequences.

    Provides high-level API for DNA sequence generation using
    the BacktrackingEngine for algorithm implementation.
    """

    def __init__(self, config: Optional[GeneratorConfig] = None):
        """
        Initialize the DNA generator.

        Args:
            config: Generator configuration (None = defaults)
        """
        self.config = config or GeneratorConfig()
        # Convert GeneratorConfig to ValidationRules and ThermodynamicParams for dna_commons
        from dna_commons import ValidationRules, ThermodynamicParams

        # Create validation rules from config - use validation_rules dict from config
        rules = ValidationRules(
            gc_content=self.config.validation_rules.get('gc_content', True),
            melting_temperature=self.config.validation_rules.get('melting_temperature', True),
            homopolymer_runs=self.config.validation_rules.get('homopolymer_runs', True),
            dinucleotide_repeats=self.config.validation_rules.get('dinucleotide_repeats', True),
            three_prime_stability=self.config.validation_rules.get('three_prime_stability', True),
            hairpin_structures=self.config.validation_rules.get('hairpin_structures', True),
            homodimer_structures=self.config.validation_rules.get('homodimer_structures', True),
            min_gc=self.config.min_gc,
            max_gc=self.config.max_gc,
            min_tm=self.config.min_tm,
            max_tm=self.config.max_tm,
            max_hairpin_tm=self.config.max_hairpin_tm,
            max_homodimer_tm=self.config.max_homodimer_tm,
            max_homopolymer_length=self.config.max_homopolymer_length,
            max_dinucleotide_repeats=self.config.max_dinucleotide_repeats
        )

        # Create thermodynamic parameters - load from default config
        thermoparams = ThermodynamicParams.load_default()

        self.validator = DNAValidator(rules, thermoparams)
        self.analyzer = SequenceAnalyzer()
        self.engine = BacktrackingEngine(self.config, self.validator)

        logger.info(f"Initialized DNAGenerator in {self.config.generation_mode.value} mode")
    
    def generate(self,
                 initial_sequence: str,
                 target_length: int,
                 seed: Optional[int] = None) -> GenerationResult:
        """
        Generate a DNA sequence of the requested length.
        
        Args:
            initial_sequence: Initial sequence
            target_length: Target sequence length
            seed: Optional seed (overrides config seed)
            
        Returns:
            GenerationResult with generation outcome
        """
        start_time = time.time()
        
        try:
            # Walidacja i normalizacja danych wejściowych
            initial_sequence = normalize_dna_sequence(initial_sequence)
            
            if target_length < len(initial_sequence):
                raise InputError(
                    f"Target length ({target_length}) must be >= initial sequence length ({len(initial_sequence)})",
                    parameter="target_length",
                    value=target_length
                )
            
            # Przygotowanie generatora losowego
            if self.config.generation_mode == GenerationMode.DETERMINISTIC:
                if seed is None:
                    seed = self.config.seed
                if seed is None:
                    # Generuj seed z parametrów
                    additional_data = f"{target_length}_{self.config.window_size}_{self.config.min_gc:.2f}_{self.config.max_gc:.2f}"
                    seed = generate_seed_from_string(initial_sequence, additional_data)
                logger.info(f"Using deterministic mode with seed={seed}")
            else:
                seed = None
                logger.info("Using random mode")
            
            deterministic = (self.config.generation_mode == GenerationMode.DETERMINISTIC)
            random_gen = DeterministicRandom(seed, deterministic)

            # Delegate generation to engine
            sequence, stats = self.engine.generate_sequence(
                initial_sequence,
                target_length,
                random_gen
            )
            
            if sequence is None:
                return GenerationResult(
                    success=False,
                    initial_sequence=initial_sequence,
                    target_length=target_length,
                    generation_stats=stats,
                    error_message="Failed to generate a sequence meeting window criteria",
                    generation_time=time.time() - start_time
                )
            
            # Walidacja końcowej sekwencji
            quality_metrics = self.validator.validate_sequence(sequence)
            
            # Analiza sekwencji
            analysis = self.analyzer.analyze_sequence(sequence)
            stats['sequence_analysis'] = analysis
            
            return GenerationResult(
                success=True,
                sequence=sequence,
                initial_sequence=initial_sequence,
                target_length=target_length,
                actual_length=len(sequence),
                quality_metrics=quality_metrics,
                generation_stats=stats,
                generation_time=time.time() - start_time
            )
            
        except Exception as e:
            logger.error(f"Error during generation: {e}")
            return GenerationResult(
                success=False,
                initial_sequence=initial_sequence,
                target_length=target_length,
                error_message=str(e),
                generation_time=time.time() - start_time
            )
    
    
    
    def generate_multiple(self,
                         initial_sequence: str,
                         target_length: int,
                         count: int = 5,
                         seed: Optional[int] = None) -> List[GenerationResult]:
        """
        Generate multiple sequences (useful for determinism testing).
        
        Args:
            initial_sequence: Initial sequence
            target_length: Target length
            count: Number of sequences to generate
            seed: Optional seed
            
        Returns:
            List of generation results
        """
        results = []
        
        logger.info(f"Generating {count} sequences in {self.config.generation_mode.value} mode")
        
        for i in range(count):
            logger.info(f"Generating sequence {i+1}/{count}")
            if self.config.generation_mode == GenerationMode.RANDOM:
                # Ensure different RNG seeds per sequence to promote diversity
                iter_seed = secrets.randbits(32)
            else:
                iter_seed = seed
            result = self.generate(initial_sequence, target_length, iter_seed)
            results.append(result)
        
        # Analiza determinizmu
        if self.config.generation_mode == GenerationMode.DETERMINISTIC:
            sequences = [r.sequence for r in results if r.success]
            if sequences:
                all_identical = all(seq == sequences[0] for seq in sequences)
                logger.info(f"All sequences identical: {all_identical}")
                if not all_identical:
                    unique_count = len(set(sequences))
                    logger.warning(f"Generated {unique_count} unique sequences in deterministic mode!")
        
        return results
