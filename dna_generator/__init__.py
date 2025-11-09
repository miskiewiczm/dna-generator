"""
Deterministic DNA sequence generator.

This package generates DNA sequences with biochemical quality control
using a backtracking algorithm. It supports deterministic (reproducible)
and purely random modes.

Main components:
    DNAGenerator: Main class to generate DNA sequences
    GeneratorConfig: Generation and validation configuration
    DNAValidator: Validation of DNA sequence quality
    GenerationMode: Generation mode enum (DETERMINISTIC/RANDOM)

Quick example:
    from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode
    
    # Deterministic mode
    config = GeneratorConfig(generation_mode=GenerationMode.DETERMINISTIC)
    generator = DNAGenerator(config)
    result = generator.generate("ATGCATGC", target_length=100)
    
    # Random mode
    config = GeneratorConfig(generation_mode=GenerationMode.RANDOM)
    generator = DNAGenerator(config)
    result = generator.generate("ATGCATGC", target_length=100)
    
    if result.success:
        print(f"DNA sequence: {result.sequence}")
        print(f"Quality: {result.quality_metrics}")
"""

from .config import GeneratorConfig, GenerationMode, DEFAULT_CONFIG
from .generator import DNAGenerator, GenerationResult
# Import from DNA Commons instead of local modules
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dna_commons import (
    DNAValidator,
    QualityMetrics,
    DeterministicRandom,
    SequenceAnalyzer,
    Primer3Adapter,
    ThermodynamicParams,
    PRIMER3_AVAILABLE,
    generate_seed_from_string
)

from .exceptions import (
    GeneratorError,
    ValidationError,
    BacktrackingError,
    ConfigurationError
)
from .backtracking_engine import BacktrackingEngine

__version__ = "1.0.0"
__author__ = "DNA Generator Team"

__all__ = [
    "DNAGenerator",
    "GeneratorConfig",
    "GenerationMode",
    "GenerationResult",
    "DNAValidator",
    "QualityMetrics",
    "DeterministicRandom",
    "SequenceAnalyzer",
    "GeneratorError",
    "ValidationError",
    "BacktrackingError",
    "ConfigurationError",
    "DEFAULT_CONFIG",
    "generate_seed_from_string",
    "Primer3Adapter",
    "ThermodynamicParams",
    "PRIMER3_AVAILABLE",
    "BacktrackingEngine"
]
